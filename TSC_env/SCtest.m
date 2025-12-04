function [s_array, sigx_array, sigy_array, k1_values] = SCEnv(beamline, init_params, E0, calc_k1_flag, beamfile1, num_slice1)
% SCtest - 统一的空间电荷包络计算入口函数
% 当输入参数小于5时：使用单一包络计算
% 当输入参数大于等于5时：使用基于切片的多包络计算

% 检查输入参数个数，决定计算模式
if nargin < 5
    % 单一包络模式
    [s_array, sigx_array, sigy_array, k1_values] = SCenvelopeEq(beamline, init_params, E0, calc_k1_flag);
    return;
else
    % 基于切片的多包络计算模式
    % 检查输入参数
    if nargin < 6
        num_slice1=20;
    end
    
    beamfile = beamfile1;
    num_slice = num_slice1;
    
    % 读取切片粒子分布
    slice_counts = read_z_distribution(beamfile, num_slice);
    total_particles = slice_counts(1);  % 总粒子数
    slice_particle_counts = slice_counts(2:end-1);  % 各切片粒子数
    
    % 解析初始参数
    betax = init_params(1);
    betay = init_params(2);
    alphax = init_params(3);
    alphay = init_params(4);
    emitx = init_params(5);
    emity = init_params(6);
    sigz = init_params(7);
    Q_total = init_params(8);  % 总电荷量
    
    % 物理常数
    me_c2 = 0.511e6;
    rc = 2.8179e-15;
    e = 1.602e-19;
    c = 2.998e8;
    
    % 计算总的实际粒子数
    N_total_real = abs(Q_total) / e;  % 根据总电荷量计算的实际粒子数
    
    % 初始gamma和几何发射度
    gamma0 = E0 / me_c2;
    geo_emitx0 = emitx / gamma0;
    geo_emity0 = emity / gamma0;
    
    z_val = slice_counts(end)/num_slice;  % 每个切片的长度
    
    % 初始包络参数
    sig0x = sqrt(betax * geo_emitx0);
    sigp0x = -alphax * sqrt(geo_emitx0 / betax);
    sig0y = sqrt(betay * geo_emity0);
    sigp0y = -alphay * sqrt(geo_emity0 / betay);
    
    % 初始化k1值存储
    if calc_k1_flag
        k1_values = [];
    end
    
    % 四极铁梯度存储
    quad_gradients = containers.Map();
    
    % 初始条件
    current_y = [sig0x, sigp0x, sig0y, sigp0y, E0];
    
    % 优化的ODE选项
    options = odeset('RelTol', 5e-6, 'AbsTol', 5e-8, 'MaxStep', 0.03, 'Stats', 'off');
    
    % 存储所有切片的结果
    all_slice_results = cell(num_slice, 1);
    s_common = [];  % 公共的s坐标
    
    % 逐切片计算
    for slice_idx = 1:num_slice
        % 计算当前切片的电荷量和实际粒子数
        Q_slice = Q_total * slice_particle_counts(slice_idx) / total_particles;
        N_slice_real = abs(Q_slice) / e;  % 根据切片电荷量计算的实际粒子数
        
        % 重置到初始位置（每个切片都从头开始传输）
        current_y_slice = current_y;
        current_s_slice = 0;
        
        % 为当前切片创建临时存储
        s_array_slice = [];
        sigx_array_slice = [];
        sigy_array_slice = [];
        
        % 逐元素求解当前切片
        for i = 1:length(beamline)
            element = beamline(i);
            
            % 跳过长度为0的元素
            if element.length == 0
                continue;
            end
            
            s_start = current_s_slice;
            s_end = current_s_slice + element.length;
            
            % 根据元素类型设置微分方程
            switch lower(element.type)
                case 'drift'
                    [z_seg, y_seg] = ode45(@(z,y) envelopeODE_sliced(z, y, 0, 0, 0, N_slice_real, z_val, rc, sigz, emitx, emity, me_c2), ...
                                           [s_start, s_end], current_y_slice, options);
                    
                case 'quadrupole'
                    k1_design = element.k1;
                    
                    if quad_gradients.isKey(element.name)
                        gradient = quad_gradients(element.name);
                    else
                        gamma_current = current_y_slice(5) / me_c2;
                        gradient = k1_design * gamma_current * me_c2 * 1.602e-19 / (e * c);
                        quad_gradients(element.name) = gradient;
                    end
                    
                    if calc_k1_flag && slice_idx == 1  % 只在第一个切片时计算k1值
                        gamma_entrance = current_y_slice(5) / me_c2;
                        k1_entrance = (e * gradient * c) / (gamma_entrance * me_c2 * 1.602e-19);
                        
                        [z_temp, y_temp] = ode45(@(z,y) envelopeODE_quad_sliced(z, y, gradient, N_slice_real, z_val, rc, sigz, emitx, emity, me_c2, e, c), ...
                                               [s_start, s_end], current_y_slice, options);
                        
                        gamma_exit = y_temp(end,5) / me_c2;
                        k1_exit = (e * gradient * c) / (gamma_exit * me_c2 * 1.602e-19);
                        k1_average = (k1_entrance + k1_exit) / 2;
                        k1_values = [k1_values, k1_average];
                        
                        z_seg = z_temp;
                        y_seg = y_temp;
                    else
                        [z_seg, y_seg] = ode45(@(z,y) envelopeODE_quad_sliced(z, y, gradient, N_slice_real, z_val, rc, sigz, emitx, emity, me_c2, e, c), ...
                                               [s_start, s_end], current_y_slice, options);
                    end
                    
                case 'sextupole'
                    [z_seg, y_seg] = ode45(@(z,y) envelopeODE_sliced(z, y, 0, 0, 0, N_slice_real, z_val, rc, sigz, emitx, emity, me_c2), ...
                                           [s_start, s_end], current_y_slice, options);
                    
                case 'cavity'
                    if isfield(element, 'num_cells') && isfield(element, 'frequency')
                        element.length = element.num_cells * 0.5 * 3e8 / (element.frequency * 1e6);
                    end
                    
                    A = 0;
                    if isfield(element, 'gradient')
                        A = element.gradient * 1e6;
                    end

                    phi = 0;
                    if isfield(element, 'phase')
                        phi = element.phase * pi/180;
                    end
                    
                    energy_gain = A * element.length * cos(phi);
                    
                    % 入口edge focusing kick
                    if energy_gain > 0
                        dgamma_over_gamma_entrance = A*cos(phi) / current_y_slice(5);
                        M21_entrance = -dgamma_over_gamma_entrance / 2;
                        
                        current_y_slice(2) = current_y_slice(2) + M21_entrance * current_y_slice(1);
                        current_y_slice(4) = current_y_slice(4) + M21_entrance * current_y_slice(3);
                    end
                    
                    [z_seg, y_seg] = ode45(@(z,y) envelopeODE_sliced(z, y, 0, A, phi, N_slice_real, z_val, rc, sigz, emitx, emity, me_c2), ...
                                           [s_start, s_end], current_y_slice, options);
                    
                    current_y_slice = y_seg(end,:)';
                    
                    % 出口edge focusing kick
                    if energy_gain > 0
                        dgamma_over_gamma_exit = A*cos(phi) / current_y_slice(5);
                        M21_exit = dgamma_over_gamma_exit / 2;
                        
                        current_y_slice(2) = current_y_slice(2) + M21_exit * current_y_slice(1);
                        current_y_slice(4) = current_y_slice(4) + M21_exit * current_y_slice(3);
                    end
                    
                case {'dipole', 'bend', 'sbend', 'rbend'}
                    error('SC计算不支持色散元件: %s', element.type);
                    
                case 'marker'
                    continue;
                    
                otherwise
                    [z_seg, y_seg] = ode45(@(z,y) envelopeODE_sliced(z, y, 0, 0, 0, N_slice_real, z_val, rc, sigz, emitx, emity, me_c2), ...
                                           [s_start, s_end], current_y_slice, options);
            end
            
            % 存储当前切片的结果
            if strcmpi(element.type, 'cavity')
                if isempty(s_array_slice)
                    s_array_slice = z_seg;
                    sigx_array_slice = y_seg(:,1);
                    sigy_array_slice = y_seg(:,3);
                else
                    s_array_slice = [s_array_slice; z_seg(2:end)];
                    sigx_array_slice = [sigx_array_slice; y_seg(2:end,1)];
                    sigy_array_slice = [sigy_array_slice; y_seg(2:end,3)];
                end
            else
                if isempty(s_array_slice)
                    s_array_slice = z_seg;
                    sigx_array_slice = y_seg(:,1);
                    sigy_array_slice = y_seg(:,3);
                else
                    s_array_slice = [s_array_slice; z_seg(2:end)];
                    sigx_array_slice = [sigx_array_slice; y_seg(2:end,1)];
                    sigy_array_slice = [sigy_array_slice; y_seg(2:end,3)];
                end
                current_y_slice = y_seg(end,:)';
            end
            
            current_s_slice = s_end;
        end
        
        % 存储当前切片的完整结果
        all_slice_results{slice_idx} = struct('s', s_array_slice, ...
                                              'sigx', sigx_array_slice, ...
                                              'sigy', sigy_array_slice, ...
                                              'particle_count', slice_particle_counts(slice_idx));
        
        % 设置公共的s坐标（使用第一个切片的s坐标）
        if slice_idx == 1
            s_common = s_array_slice;
        end
    end
    
    % 减少RMS求解点数以加速计算
    num_points_original = length(s_common);
    downsample_factor = 3;  % 降采样因子
    num_points_reduced = ceil(num_points_original / downsample_factor);
    
    % 创建降采样的索引
    sample_indices = round(linspace(1, num_points_original, num_points_reduced));
    s_reduced = s_common(sample_indices);
    
    % 计算加权RMS结果
    sigx_weighted_rms = zeros(num_points_reduced, 1);
    sigy_weighted_rms = zeros(num_points_reduced, 1);
    
    for i = 1:num_points_reduced
        weighted_sum_x = 0;
        weighted_sum_y = 0;
        original_idx = sample_indices(i);
        
        for slice_idx = 1:num_slice
            slice_result = all_slice_results{slice_idx};
            
            if length(slice_result.s) == num_points_original
                sigx_val = slice_result.sigx(original_idx);
                sigy_val = slice_result.sigy(original_idx);
            else
                sigx_val = interp1(slice_result.s, slice_result.sigx, s_reduced(i), 'linear', 'extrap');
                sigy_val = interp1(slice_result.s, slice_result.sigy, s_reduced(i), 'linear', 'extrap');
            end
            
            weight = slice_result.particle_count;
            weighted_sum_x = weighted_sum_x + sigx_val^2 * weight;
            weighted_sum_y = weighted_sum_y + sigy_val^2 * weight;
        end
        
        sigx_weighted_rms(i) = sqrt(weighted_sum_x / total_particles);
        sigy_weighted_rms(i) = sqrt(weighted_sum_y / total_particles);
    end
    
    % 输出结果
    s_array = s_reduced;
    sigx_array = sigx_weighted_rms;
    sigy_array = sigy_weighted_rms;
    
    if calc_k1_flag && ~exist('k1_values', 'var')
        k1_values = [];
    end
end

end

% 单一包络计算函数
function [s_array, sigx_array, sigy_array, k1_values] = SCenvelopeEq(beamline, init_params, E0, calc_k1_flag)
% SCenvelopeEq - 使用新格式beamline计算包络演化（优化版）

% 检查输入参数
if nargin < 4
    calc_k1_flag = 0;
end

% 解析初始参数
betax = init_params(1);
betay = init_params(2);
alphax = init_params(3);
alphay = init_params(4);
emitx = init_params(5);
emity = init_params(6);
sigz = init_params(7);
Q = init_params(8);

% 物理常数
me_c2 = 0.511e6;
rc = 2.8179e-15;
e = 1.602e-19;
c = 2.998e8;

% 检查是否有加速腔元素（向量化操作）
element_types = {beamline.type};
has_cavity = any(strcmpi(element_types, 'cavity'));

% 如果没有加速腔，使用固定k1模式
if ~has_cavity
    [s_array, sigx_array, sigy_array] = SCenvelopeEq_fixedK1(beamline, init_params, E0);
    if calc_k1_flag
        % 直接输出固定k1值
        quad_indices = strcmpi(element_types, 'quadrupole');
        k1_values = [beamline(quad_indices).k1];
    else
        k1_values = [];
    end
    return;
end

% 初始gamma和几何发射度
gamma0 = E0 / me_c2;
geo_emitx0 = emitx / gamma0;
geo_emity0 = emity / gamma0;

% 初始包络参数
sig0x = sqrt(betax * geo_emitx0);
sigp0x = -alphax * sqrt(geo_emitx0 / betax);
sig0y = sqrt(betay * geo_emity0);
sigp0y = -alphay * sqrt(geo_emity0 / betay);

% 预分配输出数组
s_array = [];
sigx_array = [];
sigy_array = [];

% 初始化k1值存储（如果需要）
if calc_k1_flag
    k1_values = [];
end

% 四极铁梯度存储
quad_gradients = containers.Map();

% 初始条件
current_y = [sig0x, sigp0x, sig0y, sigp0y, E0];
current_s = 0;

% 优化的ODE选项
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'MaxStep', 0.02, 'Stats', 'off');

% 逐元素求解
for i = 1:length(beamline)
    element = beamline(i);
    
    if element.length == 0
        continue;
    end
    
    s_start = current_s;
    s_end = current_s + element.length;
    
    % 根据元素类型设置微分方程
    switch lower(element.type)
        case 'drift'
            [z_seg, y_seg] = ode45(@(z,y) envelopeODE(z, y, 0, 0, 0, Q, rc, sigz, emitx, emity, me_c2), ...
                                   [s_start, s_end], current_y, options);
            
        case 'quadrupole'
            k1_design = element.k1;
            
            if quad_gradients.isKey(element.name)
                gradient = quad_gradients(element.name);
            else
                gamma_current = current_y(5) / me_c2;
                gradient = k1_design * gamma_current * me_c2 * 1.602e-19 / (e * c);
                quad_gradients(element.name) = gradient;
            end
            
            if calc_k1_flag
                gamma_entrance = current_y(5) / me_c2;
                k1_entrance = (e * gradient * c) / (gamma_entrance * me_c2 * 1.602e-19);
                
                [z_temp, y_temp] = ode45(@(z,y) envelopeODE_quad(z, y, gradient, Q, rc, sigz, emitx, emity, me_c2, e, c), ...
                                       [s_start, s_end], current_y, options);
                
                gamma_exit = y_temp(end,5) / me_c2;
                k1_exit = (e * gradient * c) / (gamma_exit * me_c2 * 1.602e-19);
                k1_average = (k1_entrance + k1_exit) / 2;
                k1_values = [k1_values, k1_average];
                
                z_seg = z_temp;
                y_seg = y_temp;
            else
                [z_seg, y_seg] = ode45(@(z,y) envelopeODE_quad(z, y, gradient, Q, rc, sigz, emitx, emity, me_c2, e, c), ...
                                       [s_start, s_end], current_y, options);
            end
            
        case 'sextupole'
            [z_seg, y_seg] = ode45(@(z,y) envelopeODE(z, y, 0, 0, 0, Q, rc, sigz, emitx, emity, me_c2), ...
                                   [s_start, s_end], current_y, options);
            
        case 'cavity'
            element.length = element.num_cells * 0.5 * 3e8 / (element.frequency * 1e6);
            
            if isfield(element, 'gradient')
                A = element.gradient * 1e6;
            else
                A = 0;
            end

            if isfield(element, 'phase')
                phi = element.phase * pi/180;
            else
                phi = 0;
            end
            
            energy_gain = A * element.length * cos(phi);
            
            % 入口edge focusing kick
            if energy_gain > 0
                dgamma_over_gamma_entrance = A*cos(phi) / current_y(5);
                M21_entrance = -dgamma_over_gamma_entrance / 2;
                
                current_y(2) = current_y(2) + M21_entrance * current_y(1);
                current_y(4) = current_y(4) + M21_entrance * current_y(3);
            end
            
            [z_seg, y_seg] = ode45(@(z,y) envelopeODE(z, y, 0, A, phi, Q, rc, sigz, emitx, emity, me_c2), ...
                                   [s_start, s_end], current_y, options);
            
            current_y = y_seg(end,:)';
            
            % 出口edge focusing kick
            if energy_gain > 0
                dgamma_over_gamma_exit = A*cos(phi) / current_y(5);
                M21_exit = dgamma_over_gamma_exit / 2;
                
                current_y(2) = current_y(2) + M21_exit * current_y(1);
                current_y(4) = current_y(4) + M21_exit * current_y(3);
            end
            
        case {'dipole', 'bend', 'sbend', 'rbend'}
            error('SC计算不支持色散元件: %s', element.type);
            
        case 'marker'
            continue;
            
        otherwise
            [z_seg, y_seg] = ode45(@(z,y) envelopeODE(z, y, 0, 0, 0, Q, rc, sigz, emitx, emity, me_c2), ...
                                   [s_start, s_end], current_y, options);
    end
    
    % 存储结果
    if strcmpi(element.type, 'cavity')
        if isempty(s_array)
            s_array = z_seg;
            sigx_array = y_seg(:,1);
            sigy_array = y_seg(:,3);
        else
            s_array = [s_array; z_seg(2:end)];
            sigx_array = [sigx_array; y_seg(2:end,1)];
            sigy_array = [sigy_array; y_seg(2:end,3)];
        end
    else
        if isempty(s_array)
            s_array = z_seg;
            sigx_array = y_seg(:,1);
            sigy_array = y_seg(:,3);
        else
            s_array = [s_array; z_seg(2:end)];
            sigx_array = [sigx_array; y_seg(2:end,1)];
            sigy_array = [sigy_array; y_seg(2:end,3)];
        end
        current_y = y_seg(end,:)';
    end
    
    current_s = s_end;
end

if calc_k1_flag && ~exist('k1_values', 'var')
    k1_values = [];
end

end

% 其余辅助函数保持不变...
function [s_array, sigx_array, sigy_array] = SCenvelopeEq_fixedK1(beamline, init_params, E0)
% 没有加速腔时使用固定k1的快速计算

% 解析初始参数
betax = init_params(1);
betay = init_params(2);
alphax = init_params(3);
alphay = init_params(4);
emitx = init_params(5);
emity = init_params(6);
sigz = init_params(7);
Q = init_params(8);

% 物理常数
me_c2 = 0.511e6;
rc = 2.8179e-15;

% 初始gamma和几何发射度
gamma0 = E0 / me_c2;
geo_emitx0 = emitx / gamma0;
geo_emity0 = emity / gamma0;

% 初始包络参数
sig0x = sqrt(betax * geo_emitx0);
sigp0x = -alphax * sqrt(geo_emitx0 / betax);
sig0y = sqrt(betay * geo_emity0);
sigp0y = -alphay * sqrt(geo_emity0 / betay);

% 初始化输出数组
s_array = [];
sigx_array = [];
sigy_array = [];

% 初始条件
current_y = [sig0x, sigp0x, sig0y, sigp0y, E0];
current_s = 0;

% 优化的ODE选项
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'MaxStep', 0.02, 'Stats', 'off');

% 逐元素求解
for i = 1:length(beamline)
    element = beamline(i);
    
    if element.length == 0
        continue;
    end
    
    s_start = current_s;
    s_end = current_s + element.length;
    
    switch lower(element.type)
        case 'drift'
            [z_seg, y_seg] = ode45(@(z,y) envelopeODE(z, y, 0, 0, 0, Q, rc, sigz, emitx, emity, me_c2), ...
                                   [s_start, s_end], current_y, options);
            
        case 'quadrupole'
            k1 = element.k1;
            [z_seg, y_seg] = ode45(@(z,y) envelopeODE(z, y, k1, 0, 0, Q, rc, sigz, emitx, emity, me_c2), ...
                                   [s_start, s_end], current_y, options);
            
        case 'sextupole'
            [z_seg, y_seg] = ode45(@(z,y) envelopeODE(z, y, 0, 0, 0, Q, rc, sigz, emitx, emity, me_c2), ...
                                   [s_start, s_end], current_y, options);
            
        case {'dipole', 'bend', 'sbend', 'rbend'}
            error('SC计算不支持色散元件: %s', element.type);
            
        case 'marker'
            continue;
            
        otherwise
            [z_seg, y_seg] = ode45(@(z,y) envelopeODE(z, y, 0, 0, 0, Q, rc, sigz, emitx, emity, me_c2), ...
                                   [s_start, s_end], current_y, options);
    end
    
    % 存储结果
    if isempty(s_array)
        s_array = z_seg;
        sigx_array = y_seg(:,1);
        sigy_array = y_seg(:,3);
    else
        s_array = [s_array; z_seg(2:end)];
        sigx_array = [sigx_array; y_seg(2:end,1)];
        sigy_array = [sigy_array; y_seg(2:end,3)];
    end
    
    current_y = y_seg(end,:)';
    current_s = s_end;
end

end

% 所有ODE函数保持之前的定义...
function dydt = envelopeODE_sliced(s, y, k1, A, phi, N_slice_real, z_val, rc, sigz, norm_emitx, norm_emity, me_c2)
% 切片模式的ODE函数
sigx = max(y(1), 1e-10);
sigy = max(y(3), 1e-10);
E = max(y(5), me_c2);
gamma = E / me_c2;

geo_emitx = norm_emitx / gamma;
geo_emity = norm_emity / gamma;

emittance_x = geo_emitx^2 / (sigx^3);
emittance_y = geo_emity^2 / (sigy^3);

if N_slice_real > 0 && z_val > 0
    e = 1.602e-19;
    N_density = N_slice_real / z_val;
    lambda = 1 / (sigz * gamma);
    A0 = 1/pi*sqrt(pi)/gamma/gamma*N_density * rc * lambda/(1/(2*sqrt(pi)*sigz));
    
    if A0 > 0
        sc_x = spaceChargeXFast(sigx, sigy, A0, 1/sigx);
        sc_y = spaceChargeYFast(sigx, sigy, A0, 1/sigy);
    else
        sc_x = 0; sc_y = 0;
    end
else
    sc_x = 0;
    sc_y = 0;
end

if nargin >= 5 && A > 0
    dE_ds = A * cos(phi);
    dgamma_over_gamma = dE_ds / E;
    if abs(cos(phi)) > 1e-6
        second_order_coeff = -1 / (8* cos(phi)^2) * dgamma_over_gamma^2;
    else
        second_order_coeff = 0;
    end
else
    dE_ds = 0;
    dgamma_over_gamma = 0;
    second_order_coeff = 0;
end

dydt = [
    y(2);
    emittance_x + sc_x - k1 * sigx - dgamma_over_gamma * y(2) + second_order_coeff * sigx;
    y(4);
    emittance_y + sc_y + k1 * sigy - dgamma_over_gamma * y(4) + second_order_coeff * sigy;
    dE_ds
];
end

function dydt = envelopeODE_quad_sliced(s, y, gradient, N_slice_real, z_val, rc, sigz, norm_emitx, norm_emity, me_c2, e, c)
% 切片模式的四极铁ODE函数
sigx = max(y(1), 1e-10);
sigy = max(y(3), 1e-10);
E = max(y(5), me_c2);
gamma = E / me_c2;

k1 = (e * gradient * c) / (gamma * me_c2 * 1.602e-19);

geo_emitx = norm_emitx / gamma;
geo_emity = norm_emity / gamma;

emittance_x = geo_emitx^2 / (sigx^3);
emittance_y = geo_emity^2 / (sigy^3);

if N_slice_real > 0 && z_val > 0
    N_density = N_slice_real / z_val;
    lambda = 1 / (sigz * gamma);
    A0 = 1/pi*sqrt(pi)/gamma/gamma*N_density * rc * lambda/(1/(2*sqrt(pi)*sigz));
    
    if A0 > 0
        sc_x = spaceChargeXFast(sigx, sigy, A0, 1/sigx);
        sc_y = spaceChargeYFast(sigx, sigy, A0, 1/sigy);
    else
        sc_x = 0; sc_y = 0;
    end
else
    sc_x = 0;
    sc_y = 0;
end

dydt = [
    y(2);
    emittance_x + sc_x - k1 * sigx;
    y(4);
    emittance_y + sc_y + k1 * sigy;
    0
];
end

function dydt = envelopeODE_quad(s, y, gradient, Q, rc, sigz, norm_emitx, norm_emity, me_c2, e, c)
% 单一包络模式的四极铁ODE函数
sigx = max(y(1), 1e-10);
sigy = max(y(3), 1e-10);
E = max(y(5), me_c2);
gamma = E / me_c2;

k1 = (e * gradient * c) / (gamma * me_c2 * 1.602e-19);

geo_emitx = norm_emitx / gamma;
geo_emity = norm_emity / gamma;

emittance_x = geo_emitx^2 / (sigx^3);
emittance_y = geo_emity^2 / (sigy^3);

if abs(Q) > 1e-20
    N = Q / e;
    lambda = 1 / (sigz * gamma);
    A0 = 1/pi*sqrt(pi)/gamma/gamma*N * rc * lambda;
    
    if A0 > 0
        sc_x = spaceChargeXFast(sigx, sigy, A0, 1/sigx);
        sc_y = spaceChargeYFast(sigx, sigy, A0, 1/sigy);
    else
        sc_x = 0; sc_y = 0;
    end
else
    sc_x = 0;
    sc_y = 0;
end

dydt = [
    y(2);
    emittance_x + sc_x - k1 * sigx;
    y(4);
    emittance_y + sc_y + k1 * sigy;
    0
];
end

function dydt = envelopeODE(s, y, k1, A, phi, Q, rc, sigz, norm_emitx, norm_emity, me_c2)
% 单一包络模式的通用ODE函数
sigx = max(y(1), 1e-10);
sigy = max(y(3), 1e-10);
E = max(y(5), me_c2);
gamma = E / me_c2;

geo_emitx = norm_emitx / gamma;
geo_emity = norm_emity / gamma;

emittance_x = geo_emitx^2 / (sigx^3);
emittance_y = geo_emity^2 / (sigy^3);

if abs(Q) > 1e-20
    e = 1.602e-19;
    N = Q / e;
    lambda = 1 / (sigz * gamma);
    A0 = 1/pi*sqrt(pi)/gamma/gamma*N * rc * lambda;
    
    if A0 > 0
        sc_x = spaceChargeXFast(sigx, sigy, A0, 1/sigx);
        sc_y = spaceChargeYFast(sigx, sigy, A0, 1/sigy);
    else
        sc_x = 0; sc_y = 0;
    end
else
    sc_x = 0;
    sc_y = 0;
end

if nargin >= 5 && A > 0
    dE_ds = A * cos(phi);
    dgamma_over_gamma = dE_ds / E;
    if abs(cos(phi)) > 1e-6
        second_order_coeff = -1 / (8* cos(phi)^2) * dgamma_over_gamma^2;
    else
        second_order_coeff = 0;
    end
else
    dE_ds = 0;
    dgamma_over_gamma = 0;
    second_order_coeff = 0;
end

dydt = [
    y(2);
    emittance_x + sc_x - k1 * sigx - dgamma_over_gamma * y(2) + second_order_coeff * sigx;
    y(4);
    emittance_y + sc_y + k1 * sigy - dgamma_over_gamma * y(4) + second_order_coeff * sigy;
    dE_ds
];
end

function val = spaceChargeXFast(sigx, sigy, A0, sigx_inv)
    k = max(min(sigy * sigx_inv, 1e6), 1e-6);
    k_power = max(min(k^(-1/5.5), 1e6), 1e-6);
    f = 1 - 1/sqrt(1 + k_power);
    exp_arg = max(-1/(8*k^(1/5.5)), -50);
    g = (1/4.25) / (1+k) / (1 - max(exp(exp_arg), 1e-15));
    val = max(min(f * g * A0 * sigx_inv, 1e6), -1e6);
end

function val = spaceChargeYFast(sigx, sigy, A0, sigy_inv)
    k = max(min(sigx * sigy_inv, 1e6), 1e-6);
    k_power = max(min(k^(-1/5.5), 1e6), 1e-6);
    f = 1 - 1/sqrt(1 + k_power);
    exp_arg = max(-1/(8*k^(1/5.5)), -50);
    g = (1/4.25) / (1+k) / (1 - max(exp(exp_arg), 1e-15));
    val = max(min(f * g * A0 * sigy_inv, 1e6), -1e6);
end