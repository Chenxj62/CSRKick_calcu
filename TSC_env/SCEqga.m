function [s_array, sigx_array, sigy_array, k1_values] = SCEqga(beamline, init_params, calc_k1_flag)
% SCEqga - 计算考虑空间电荷效应的束流包络演化
%
% 输入:
%   beamline     - 束线元件结构体数组
%   init_params  - 初始参数 [betax, betay, alphax, alphay, emitx, emity, sigz, Q, E0]
%   calc_k1_flag - 是否计算k1值 (默认0)
%
% 输出:
%   s_array      - 纵向位置数组
%   sigx_array   - 水平包络数组
%   sigy_array   - 垂直包络数组
%   k1_values    - 四极铁k1值数组(如果calc_k1_flag=1)

% 输入验证
if nargin < 3
    calc_k1_flag = false;
else
    calc_k1_flag = logical(calc_k1_flag);
end

% 初始化输出
k1_values = [];

% 解析初始参数（使用结构体更清晰）
params = struct('betax', init_params(1), 'betay', init_params(2), ...
                'alphax', init_params(3), 'alphay', init_params(4), ...
                'emitx', init_params(5), 'emity', init_params(6), ...
                'sigz', init_params(7), 'Q', init_params(8), ...
                'E0', init_params(9)+0.01e7);

% 物理常数定义为持久变量以减少重复赋值
persistent PHYS_CONST;
if isempty(PHYS_CONST)
    PHYS_CONST = struct('me_c2', 0.511e6, 'rc', 2.8179e-15, ...
                       'e', 1.602e-19, 'c', 2.998e8);
end

% 检查加速腔元素（优化的向量化操作）
element_types = cellfun(@(x) lower(x), {beamline.type}, 'UniformOutput', false);
has_cavity = any(strcmp(element_types, 'cavity'));

% 无加速腔时使用固定k1模式（避免重复代码）
if ~has_cavity
    [s_array, sigx_array, sigy_array] = processFixedK1Mode(beamline, params, PHYS_CONST);
    if calc_k1_flag
        quad_mask = strcmp(element_types, 'quadrupole');
        k1_values = [beamline(quad_mask).k1];
    end
    return;
end

% 有加速腔的完整计算
[s_array, sigx_array, sigy_array, k1_values] = ...
    processWithCavity(beamline, params, PHYS_CONST, calc_k1_flag);

end

%% 主处理函数 - 含加速腔
function [s_array, sigx_array, sigy_array, k1_values] = ...
    processWithCavity(beamline, params, PHYS_CONST, calc_k1_flag)

% 计算初始条件
gamma0 = params.E0 / PHYS_CONST.me_c2;
geo_emit = struct('x', params.emitx / gamma0, 'y', params.emity / gamma0);

% 初始包络向量 [σx, σx', σy, σy', E]
y0 = calculateInitialConditions(params, geo_emit);

% 预分配数组（动态增长以节省内存）
result = struct('s', [], 'sigx', [], 'sigy', []);
k1_values = [];

% 四极铁梯度缓存
quad_gradients = containers.Map();

% ODE求解器选项（平衡精度和速度）
ode_opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, ...
                  'MaxStep', 0.02, 'Stats', 'off');

% 状态变量
current_y = y0;
current_s = 0;

% 主循环：逐元件处理
for i = 1:length(beamline)
    elem = beamline(i);
    
    % 跳过零长度元件
    if elem.length <= 0
        continue;
    end
    
    % 计算元件区间
    s_span = [current_s, current_s + elem.length];
    
    % 根据元件类型处理
    [z_seg, y_seg, elem_k1] = processElement(elem, s_span, current_y, ...
        params, PHYS_CONST, quad_gradients, ode_opts, calc_k1_flag);
    
    % 存储结果（避免重复代码）
    result = appendResults(result, z_seg, y_seg);
    
    % 收集k1值
    if calc_k1_flag && ~isempty(elem_k1)
        k1_values(end+1) = elem_k1;
    end
    
    % 更新状态
    current_y = y_seg(end, :)';
    current_s = s_span(2);
end

% 输出结果
s_array = result.s;
sigx_array = result.sigx;
sigy_array = result.sigy;

end

%% 元件处理函数
function [z_seg, y_seg, elem_k1] = processElement(elem, s_span, y0, ...
    params, PHYS_CONST, quad_gradients, ode_opts, calc_k1_flag)

elem_k1 = [];
elem_type = lower(elem.type);

switch elem_type
    case {'drift', 'sextupole', 'marker'}
        % 无场区域
        [z_seg, y_seg] = ode45(@(s,y) envelopeODE_drift(s, y, params, PHYS_CONST), ...
                              s_span, y0, ode_opts);
        
    case 'quadrupole'
        % 四极铁处理（含梯度缓存）
        [z_seg, y_seg, elem_k1] = processQuadrupole(elem, s_span, y0, ...
            params, PHYS_CONST, quad_gradients, ode_opts, calc_k1_flag);
        
    case 'cavity'
        % 加速腔处理
        [z_seg, y_seg] = processCavity(elem, s_span, y0, params, PHYS_CONST, ode_opts);
        
    case {'dipole', 'bend', 'sbend', 'rbend'}
        error('SCEqga:UnsupportedElement', ...
              '空间电荷计算不支持色散元件: %s', elem.type);
        
    otherwise
        warning('SCEqga:UnknownElement', ...
                '未知元件类型 "%s"，作为漂移段处理', elem.type);
        [z_seg, y_seg] = ode45(@(s,y) envelopeODE_drift(s, y, params, PHYS_CONST), ...
                              s_span, y0, ode_opts);
end

end

%% 四极铁处理
function [z_seg, y_seg, avg_k1] = processQuadrupole(elem, s_span, y0, ...
    params, PHYS_CONST, quad_gradients, ode_opts, calc_k1_flag)

avg_k1 = [];

% 获取或计算梯度
if isKey(quad_gradients, elem.name)
    gradient = quad_gradients(elem.name);
else
    gamma_current = y0(5) / PHYS_CONST.me_c2;
    gradient = elem.k1 * gamma_current * PHYS_CONST.me_c2 / PHYS_CONST.c;
    quad_gradients(elem.name) = gradient;
end

% ODE求解
[z_seg, y_seg] = ode45(@(s,y) envelopeODE_quad(s, y, gradient, params, PHYS_CONST), ...
                      s_span, y0, ode_opts);

% 计算平均k1（如果需要）
if calc_k1_flag
    gamma_in = y0(5) / PHYS_CONST.me_c2;
    gamma_out = y_seg(end,5) / PHYS_CONST.me_c2;
    k1_in = gradient * PHYS_CONST.c / (gamma_in * PHYS_CONST.me_c2);
    k1_out = gradient * PHYS_CONST.c / (gamma_out * PHYS_CONST.me_c2);
    avg_k1 = (k1_in + k1_out) / 2;
end

end

%% 加速腔处理
function [z_seg, y_seg] = processCavity(elem, s_span, y0, params, PHYS_CONST, ode_opts)

% 计算腔长度（如果需要）
if ~isfield(elem, 'length') || elem.length == 0
    elem.length = elem.num_cells * 0.5 * PHYS_CONST.c / (elem.frequency * 1e6);
end

% 获取腔参数
A = getFieldValue(elem, 'gradient', 0) * 1e6;  % V/m
phi = getFieldValue(elem, 'phase', 0) * pi/180; % rad

% 入口边缘聚焦
y_current = y0;
if A * cos(phi) > 0
    dgamma_rel = A * cos(phi) / y_current(5);
    kick = -dgamma_rel / 2;
    y_current(2) = y_current(2) + kick * y_current(1);
    y_current(4) = y_current(4) + kick * y_current(3);
end

% 腔内传输
[z_seg, y_seg] = ode45(@(s,y) envelopeODE_cavity(s, y, A, phi, params, PHYS_CONST), ...
                      [s_span(1), s_span(2)], y_current, ode_opts);

% 出口边缘聚焦
y_exit = y_seg(end,:)';
if A * cos(phi) > 0
    dgamma_rel = A * cos(phi) / y_exit(5);
    kick = dgamma_rel / 2;
    y_exit(2) = y_exit(2) + kick * y_exit(1);
    y_exit(4) = y_exit(4) + kick * y_exit(3);
    y_seg(end,:) = y_exit';
end

end

%% 微分方程定义

% 漂移段
function dydt = envelopeODE_drift(~, y, params, PHYS_CONST)
dydt = envelopeODE_general(y, 0, 0, 0, params, PHYS_CONST);
end

% 四极铁
function dydt = envelopeODE_quad(~, y, gradient, params, PHYS_CONST)
gamma = y(5) / PHYS_CONST.me_c2;
k1 = gradient * PHYS_CONST.c / (gamma * PHYS_CONST.me_c2);
dydt = envelopeODE_general(y, k1, 0, 0, params, PHYS_CONST);
end

% 加速腔
function dydt = envelopeODE_cavity(~, y, A, phi, params, PHYS_CONST)
dydt = envelopeODE_general(y, 0, A, phi, params, PHYS_CONST);
end

% 通用包络方程
function dydt = envelopeODE_general(y, k1, A, phi, params, PHYS_CONST)

% 确保数值稳定性
sigx = max(y(1), 1e-10);
sigy = max(y(3), 1e-10);
E = max(y(5), PHYS_CONST.me_c2);
gamma = E / PHYS_CONST.me_c2;

% 几何发射度
geo_emitx = params.emitx / gamma;
geo_emity = params.emity / gamma;

% 发射度力
emit_force_x = geo_emitx^2 / sigx^3;
emit_force_y = geo_emity^2 / sigy^3;

% 空间电荷力
[sc_x, sc_y] = calculateSpaceChargeForces(sigx, sigy, gamma, params, PHYS_CONST);

% 加速效应
if A > 0
    dE_ds = A * cos(phi);
    dgamma_rel = dE_ds / E;
    if abs(cos(phi)) > 1e-6
        second_order = -dgamma_rel^2 / (7 * cos(phi)^2);
    else
        second_order = 0;
    end
else
    dE_ds = 0;
    dgamma_rel = 0;
    second_order = 0;
end

% 组装导数
dydt = [
    y(2);                                                          % σx'
    emit_force_x + sc_x - k1*sigx - dgamma_rel*y(2) + second_order*sigx;  % σx''
    y(4);                                                          % σy'
    emit_force_y + sc_y + k1*sigy - dgamma_rel*y(4) + second_order*sigy;  % σy''
    dE_ds                                                          % E'
];

end

%% 空间电荷计算（优化版）
function [sc_x, sc_y] = calculateSpaceChargeForces(sigx, sigy, gamma, params, PHYS_CONST)

if abs(params.Q) < 1e-20
    sc_x = 0;
    sc_y = 0;
    return;
end

% 空间电荷参数
N = params.Q / PHYS_CONST.e;
lambda = 1 / (params.sigz * gamma);
A0 = sqrt(pi) / (pi * gamma^2) * N * PHYS_CONST.rc * lambda;

if A0 <= 0
    sc_x = 0;
    sc_y = 0;
else
    sc_x = spaceChargeForce(sigx, sigy, A0, 1/sigx);
    sc_y = spaceChargeForce(sigy, sigx, A0, 1/sigy);
end

end

% 优化的空间电荷力计算
function force = spaceChargeForce(sig1, sig2, A0, sig1_inv)
% 使用更稳定的数值计算
k = sig2 * sig1_inv;
k = max(1e-6, min(k, 1e6));

% 避免数值溢出
k_power = k^(-1/5.5);
k_power = max(1e-6, min(k_power, 1e6));

f = 1 - 1/sqrt(1 + k_power);

% 防止指数下溢
exp_arg = -1/(8 * k^(1/5.5));
exp_arg = max(exp_arg, -50);
g = (1/4.25) / ((1+k) * (1 - exp(exp_arg)));

force = f * g * A0 * sig1_inv;
force = max(-1e6, min(force, 1e6));  % 限制范围
end

%% 辅助函数

% 计算初始条件
function y0 = calculateInitialConditions(params, geo_emit)
sig0x = sqrt(params.betax * geo_emit.x);
sigp0x = -params.alphax * sqrt(geo_emit.x / params.betax);
sig0y = sqrt(params.betay * geo_emit.y);
sigp0y = -params.alphay * sqrt(geo_emit.y / params.betay);
y0 = [sig0x; sigp0x; sig0y; sigp0y; params.E0];
end

% 结果追加
function result = appendResults(result, z_seg, y_seg)
if isempty(result.s)
    result.s = z_seg;
    result.sigx = y_seg(:,1);
    result.sigy = y_seg(:,3);
else
    result.s = [result.s; z_seg(2:end)];
    result.sigx = [result.sigx; y_seg(2:end,1)];
    result.sigy = [result.sigy; y_seg(2:end,3)];
end
end

% 获取结构体字段值
function val = getFieldValue(s, field, default)
if isfield(s, field)
    val = s.(field);
else
    val = default;
end
end

% 固定K1模式处理（简化版）
function [s_array, sigx_array, sigy_array] = processFixedK1Mode(beamline, params, PHYS_CONST)

% 类似processWithCavity但使用固定k1
gamma0 = params.E0 / PHYS_CONST.me_c2;
geo_emit = struct('x', params.emitx / gamma0, 'y', params.emity / gamma0);
y0 = calculateInitialConditions(params, geo_emit);

result = struct('s', [], 'sigx', [], 'sigy', []);
ode_opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'MaxStep', 0.02);

current_y = y0;
current_s = 0;

for i = 1:length(beamline)
    elem = beamline(i);
    if elem.length <= 0, continue; end
    
    s_span = [current_s, current_s + elem.length];
    
    switch lower(elem.type)
        case 'quadrupole'
            k1 = elem.k1;
            [z_seg, y_seg] = ode45(@(s,y) envelopeODE_fixedK1(s, y, k1, params, PHYS_CONST), ...
                                  s_span, current_y, ode_opts);
        otherwise
            [z_seg, y_seg] = ode45(@(s,y) envelopeODE_drift(s, y, params, PHYS_CONST), ...
                                  s_span, current_y, ode_opts);
    end
    
    result = appendResults(result, z_seg, y_seg);
    current_y = y_seg(end,:)';
    current_s = s_span(2);
end

s_array = result.s;
sigx_array = result.sigx;
sigy_array = result.sigy;

end

% 固定K1的ODE
function dydt = envelopeODE_fixedK1(~, y, k1, params, PHYS_CONST)
dydt = envelopeODE_general(y, k1, 0, 0, params, PHYS_CONST);
end