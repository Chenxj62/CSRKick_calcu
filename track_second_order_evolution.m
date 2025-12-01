function [s_positions, T126_evolution, T116_evolution, T216_evolution, T226_evolution, ...
          T346_evolution, T336_evolution, T436_evolution, T446_evolution] = ...
    track_second_order_evolution(beamline, delta)
    % 计算二阶传输矩阵元素沿束线的演化
    %
    % 输入参数:
    %   beamline - 元件数组 [ELE.D1, ELE.Q1, ...] 
    %   delta - 计算时使用的偏差值 (默认 0.001)
    %
    % 输出参数:
    %   s_positions - 沿束线的位置数组
    %   
    %   X方向的二阶项:
    %   T126_evolution - T126(x对px和delta的二阶项)沿s的演化
    %   T116_evolution - T116(x对x和delta的二阶项)沿s的演化
    %   T216_evolution - T216(px对x和delta的二阶项)沿s的演化
    %   T226_evolution - T226(px对px和delta的二阶项)沿s的演化
    %
    %   Y方向的二阶项:
    %   T346_evolution - T346(y对py和delta的二阶项)沿s的演化
    %   T336_evolution - T336(y对y和delta的二阶项)沿s的演化
    %   T436_evolution - T436(py对y和delta的二阶项)沿s的演化
    %   T446_evolution - T446(py对py和delta的二阶项)沿s的演化
    
    % 默认参数
    if nargin < 2
        delta = 0.001;
    end
    
    % 初始化存储数组
    num_elements = length(beamline);
    s_positions = zeros(num_elements + 1, 1);
    
    % X方向的二阶项
    T126_evolution = zeros(num_elements + 1, 1);
    T116_evolution = zeros(num_elements + 1, 1);
    T216_evolution = zeros(num_elements + 1, 1);
    T226_evolution = zeros(num_elements + 1, 1);
    
    % Y方向的二阶项
    T346_evolution = zeros(num_elements + 1, 1);
    T336_evolution = zeros(num_elements + 1, 1);
    T436_evolution = zeros(num_elements + 1, 1);
    T446_evolution = zeros(num_elements + 1, 1);
    
    % 初始位置和初始值
    current_s = 0;
    s_positions(1) = current_s;
    
    % 初始值都为0
    T126_evolution(1) = 0;
    T116_evolution(1) = 0;
    T216_evolution(1) = 0;
    T226_evolution(1) = 0;
    
    T346_evolution(1) = 0;
    T336_evolution(1) = 0;
    T436_evolution(1) = 0;
    T446_evolution(1) = 0;
    
  
    
    for i = 1:num_elements
        % 计算当前位置的二阶项
        current_beamline = beamline(1:i);
        
        % X方向的二阶项
        T126_evolution(i+1) = calculate_Tijk(current_beamline, 1, 2, 6, delta);
        T116_evolution(i+1) = calculate_Tijk(current_beamline, 1, 1, 6, delta);
        T216_evolution(i+1) = calculate_Tijk(current_beamline, 2, 1, 6, delta);
        T226_evolution(i+1) = calculate_Tijk(current_beamline, 2, 2, 6, delta);
        
        % Y方向的二阶项
        T346_evolution(i+1) = calculate_Tijk(current_beamline, 3, 4, 6, delta);
        T336_evolution(i+1) = calculate_Tijk(current_beamline, 3, 3, 6, delta);
        T436_evolution(i+1) = calculate_Tijk(current_beamline, 4, 3, 6, delta);
        T446_evolution(i+1) = calculate_Tijk(current_beamline, 4, 4, 6, delta);
        
        % 更新s位置
        element = beamline(i);
        if isfield(element, 'length')
            current_s = current_s + element.length;
        end
        s_positions(i+1) = current_s;
        
        % 显示进度
        
    end
    
    
end

function Tijk = calculate_Tijk(beamline, iindex, jindex, kindex, delta)
    % 计算二阶传输矩阵元素 T_ijk
    % 
    % 参数:
    %   iindex - 输出坐标索引 (1=x, 2=px, 3=y, 4=py, 5=z, 6=delta)
    %   jindex - 第一个输入坐标索引
    %   kindex - 第二个输入坐标索引
    %   delta - 偏差值
    %
    % 返回:
    %   Tijk - 二阶传输矩阵元素值
    
    if jindex == kindex
        % 对角项: T_ijj
        % 使用公式: T_ijj = [f(δ) + f(-δ) - 2f(0)] / (2δ²)
        p1 = zeros(6, 1);
        p2 = zeros(6, 1);
        p0 = zeros(6, 1);
        
        p1(jindex) = delta;
        p2(jindex) = -delta;
        
        % 追踪粒子
        p1 = track_particle(p1, beamline);
        p2 = track_particle(p2, beamline);
        p0 = track_particle(p0, beamline);
        
        % 计算二阶项
        % 使用更精确的公式,包含零点
        Tijk = (p1(iindex) + p2(iindex) - 2*p0(iindex)) / (2 * delta^2);
        
    else
        % 非对角项: T_ijk (j≠k)
        % 使用混合偏导数公式
        p1 = zeros(6, 1);
        p2 = zeros(6, 1);
        p3 = zeros(6, 1);
        p4 = zeros(6, 1);
        
        % 配置四个测试粒子
        p1(jindex) = delta;
        p1(kindex) = delta;
        
        p2(jindex) = delta;
        p2(kindex) = -delta;
        
        p3(jindex) = -delta;
        p3(kindex) = delta;
        
        p4(jindex) = -delta;
        p4(kindex) = -delta;
        
        % 追踪粒子
        p1 = track_particle(p1, beamline);
        p2 = track_particle(p2, beamline);
        p3 = track_particle(p3, beamline);
        p4 = track_particle(p4, beamline);
        
        % 计算混合偏导数: ∂²f/∂j∂k
        % T_ijk = [f(δ,δ) - f(δ,-δ) - f(-δ,δ) + f(-δ,-δ)] / (4δ²)
        Tijk = (p1(iindex) - p2(iindex) - p3(iindex) + p4(iindex)) / (4 * delta^2);
    end
end

function p_final = track_particle(p_initial, beamline)
    % 追踪粒子穿过束线
    
    p = p_initial;
    
    for i = 1:length(beamline)
        element = beamline(i);
        
        if strcmpi(element.type, 'drift')
            p = trackDrift(p, element.length);
            
        elseif strcmpi(element.type, 'marker')
            % 标记点不影响粒子轨迹
            continue;
            
        elseif strcmpi(element.type, 'quadrupole')
            if abs(element.k1) < 1e-4
                p = trackDrift(p, element.length);
            else
                p = trackQuadrupole(p, element.length, element.k1);
            end
            
        elseif strcmpi(element.type, 'sextupole')
            p = trackSextupole(p, element.length, element.k2);
            
        elseif strcmpi(element.type, 'cavity')
            % 如果束线中有腔体,作为漂移处理(中心粒子)
            if isfield(element, 'length')
                p = trackDrift(p, element.length);
            end
            
        else
            % 未知元件类型,作为漂移处理
            if isfield(element, 'length') && element.length > 0
                p = trackDrift(p, element.length);
            end
        end
    end
    
    p_final = p;
end

function p_final = trackDrift(p_initial, L)
    % 漂移段跟踪函数
    
    x1 = p_initial(1);
    px1 = p_initial(2);
    y1 = p_initial(3);
    py1 = p_initial(4);
    z1 = p_initial(5);
    delta = p_initial(6);
    
    % 计算纵向动量
    pl = sqrt(1 - (px1^2 + py1^2)/(1 + delta)^2);
    
    % 计算漂移后的坐标
    x2 = x1 + L * px1 / ((1 + delta) * pl);
    px2 = px1;
    
    y2 = y1 + L * py1 / ((1 + delta) * pl);
    py2 = py1;
    
    z2 = z1 + (1.0 - 1.0/pl) * L;
    delta2 = delta;
    
    p_final = [x2; px2; y2; py2; z2; delta2];
end

function p_final = trackQuadrupole(p_initial, L, k1)
    % 四极铁跟踪函数
    
    x1 = p_initial(1);
    px1 = p_initial(2);
    y1 = p_initial(3);
    py1 = p_initial(4);
    z1 = p_initial(5);
    delta = p_initial(6);
    
    % 计算ω值
    wx = sqrt(abs(k1) / (1 + delta));
    wy = wx;
    
    % 处理x方向的参数
    if k1 > 0
        cx = cos(wx * L);
        sx = sin(wx * L) / wx;
        tx = -1;
    else
        cx = cosh(wx * L);
        sx = sinh(wx * L) / wx;
        tx = 1;
    end
    
    % 处理y方向的参数
    if k1 > 0
        cy = cosh(wy * L);
        sy = sinh(wy * L) / wy;
        ty = 1;
    else
        cy = cos(wy * L);
        sy = sin(wy * L) / wy;
        ty = -1;
    end
    
    % 计算传输后的坐标
    x2 = cx * x1 + sx * px1 / (1 + delta);
    px2 = tx * wx^2 * (1 + delta) * sx * x1 + cx * px1;
    
    y2 = cy * y1 + sy * py1 / (1 + delta);
    py2 = ty * wy^2 * (1 + delta) * sy * y1 + cy * py1;
    
    % 计算m参数
    m511 = tx * wx^2 * (L - cx * sx) / 4;
    m512 = -tx * wx^2 * sx^2 / (2 * (1 + delta));
    m522 = -1 * (L + cx * sx) / (4 * (1 + delta)^2);
    
    m533 = ty * wy^2 * (L - cy * sy) / 4;
    m534 = -ty * wy^2 * sy^2 / (2 * (1 + delta));
    m544 = -1 * (L + cy * sy) / (4 * (1 + delta)^2);
    
    % 计算纵向位置
    z2 = z1 + m511 * x1^2 + m512 * x1 * px1 + m522 * px1^2 + ...
         m533 * y1^2 + m534 * y1 * py1 + m544 * py1^2;
    
    delta2 = delta;
    
    p_final = [x2; px2; y2; py2; z2; delta2];
end

function p_final = trackSextupole(p_initial, L, k2)
    % 六极铁跟踪函数 - 使用drift-kick-drift模型
    
    x = p_initial(1);
    px = p_initial(2);
    y = p_initial(3);
    py = p_initial(4);
    z = p_initial(5);
    delta = p_initial(6);
    
    L_half = L / 2;
    p_final = p_initial;
    
    % 第一个漂移(L/2)
    [x_new, y_new, z_new] = bmadDrift(x, px, y, py, z, delta, L_half);
    p_final(1) = x_new;
    p_final(3) = y_new;
    p_final(5) = z_new;
    
    % 六极铁踢度
    x = p_final(1);
    px = p_final(2);
    y = p_final(3);
    py = p_final(4);
    
    px_new = px - 0.5 * k2 * (x^2 - y^2) * L;
    py_new = py + k2 * x * y * L;
    
    p_final(2) = px_new;
    p_final(4) = py_new;
    
    % 第二个漂移(L/2)
    x = p_final(1);
    px = p_final(2);
    y = p_final(3);
    py = p_final(4);
    z = p_final(5);
    
    [x_new, y_new, z_new] = bmadDrift(x, px, y, py, z, delta, L_half);
    p_final(1) = x_new;
    p_final(3) = y_new;
    p_final(5) = z_new;
end

function [x_out, y_out, z_out] = bmadDrift(x, px, y, py, z, delta, L)
    % BMAD漂移段计算
    
    pl = sqrt(1 - (px^2 + py^2)/(1 + delta)^2);
    
    x_out = x + L * px / ((1 + delta) * pl);
    y_out = y + L * py / ((1 + delta) * pl);
    z_out = z + (1.0 - 1.0/pl) * L;
end