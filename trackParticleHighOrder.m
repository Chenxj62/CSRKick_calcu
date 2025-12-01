function p_final = trackParticleHighOrder(p_initial, beamline, order)
    % 使用高阶传输矩阵追踪粒子（仅支持四极铁和漂移段）
    %
    % 输入参数:
    %   p_initial - 初始粒子坐标 [x; px; y; py; z; delta]
    %   beamline  - 元件数组
    %   order     - 传输矩阵阶数 (1=仅一阶, 2=包含二阶, 默认为2)
    %
    % 输出参数:
    %   p_final   - 最终粒子坐标 [x; px; y; py; z; delta]
    
    % 默认使用二阶
    if nargin < 3
        order = 2;
    end
    
    % 初始化当前粒子状态
    p_current = p_initial;
    
    % 逐元件追踪
    for i = 1:length(beamline)
        element = beamline(i);
        
        % 根据元件类型选择追踪方法
        if strcmpi(element.type, 'drift')
            p_current = trackDriftHighOrder(p_current, element.length, order);
            
        elseif strcmpi(element.type, 'marker')
            % 标记点不影响粒子轨迹
            continue;
            
        elseif strcmpi(element.type, 'quadrupole')
            if abs(element.k1) < 1e-4
                % 弱四极铁当作漂移
                p_current = trackDriftHighOrder(p_current, element.length, order);
            else
                p_current = trackQuadrupoleHighOrder(p_current, element.length, element.k1, order);
            end
            
        else
            warning('元件类型 %s 不支持，仅支持drift和quadrupole，跳过该元件', element.type);
        end
    end
    
    % 返回最终粒子坐标
    p_final = p_current;
end

%% ========== 漂移段高阶追踪 ==========

function p_final = trackDriftHighOrder(p_initial, L, order)
    % 漂移段高阶追踪 - 基于BMAD公式24.60-24.61
    
    x = p_initial(1);
    px = p_initial(2);
    y = p_initial(3);
    py = p_initial(4);
    z = p_initial(5);
    delta = p_initial(6);
    
    % 计算纵向动量
    pl = sqrt(1 - (px^2 + py^2)/(1 + delta)^2);
    
    % 一阶传输矩阵元素
    if order >= 1
        % x方向
        x_new = x + L * px / ((1 + delta) * pl);
        px_new = px;
        
        % y方向
        y_new = y + L * py / ((1 + delta) * pl);
        py_new = py;
        
        % z方向
        z_new = z + (1.0 - 1.0/pl) * L;
        delta_new = delta;
    end
    
    % 二阶修正
    if order >= 2
        % 计算二阶项 T_ijk
        % 基于公式: pl ≈ 1 - (px^2 + py^2)/(2(1+delta)^2) + ...
        
        % x的二阶修正
        % T_166: x对delta的二阶依赖
        T_166 = L * px / ((1 + delta)^3 * pl^3) * (px^2 + py^2) / 2;
        
        % T_126: x对px和delta的耦合
        T_126 = -L * px / ((1 + delta)^2 * pl);
        
        % T_122: x对px^2的依赖
        T_122 = L * px^3 / ((1 + delta)^3 * pl^3);
        
        % T_144: x对py^2的依赖
        T_144 = L * px * py^2 / ((1 + delta)^3 * pl^3);
        
        x_new = x_new + T_166 * delta^2 + T_126 * px * delta + ...
                T_122 * px^2 + T_144 * py^2;
        
        % y的二阶修正（类似x）
        T_366 = L * py / ((1 + delta)^3 * pl^3) * (px^2 + py^2) / 2;
        T_346 = -L * py / ((1 + delta)^2 * pl);
        T_344 = L * py^3 / ((1 + delta)^3 * pl^3);
        T_322 = L * py * px^2 / ((1 + delta)^3 * pl^3);
        
        y_new = y_new + T_366 * delta^2 + T_346 * py * delta + ...
                T_344 * py^2 + T_322 * px^2;
        
        % z的二阶修正
        % T_566: z对delta的二阶依赖
        T_566 = -L / (2 * pl^3 * (1 + delta)^2);
        
        % T_522: z对px^2的依赖
        T_522 = L * px^2 / (2 * (1 + delta)^2 * pl^3);
        
        % T_544: z对py^2的依赖
        T_544 = L * py^2 / (2 * (1 + delta)^2 * pl^3);
        
        z_new = z_new + T_566 * delta^2 + T_522 * px^2 + T_544 * py^2;
    end
    
    p_final = [x_new; px_new; y_new; py_new; z_new; delta_new];
end

%% ========== 四极铁高阶追踪 ==========

function p_final = trackQuadrupoleHighOrder(p_initial, L, k1, order)
    % 四极铁高阶追踪 - 基于BMAD公式24.81-24.84
    
    x = p_initial(1);
    px = p_initial(2);
    y = p_initial(3);
    py = p_initial(4);
    z = p_initial(5);
    delta = p_initial(6);
    
    % 计算omega参数
    wx = sqrt(abs(k1) / (1 + delta));
    wy = wx;
    
    % x方向参数
    if k1 > 0
        cx = cos(wx * L);
        sx = sin(wx * L) / wx;
        tx = -1;
    else
        cx = cosh(wx * L);
        sx = sinh(wx * L) / wx;
        tx = 1;
    end
    
    % y方向参数
    if k1 > 0
        cy = cosh(wy * L);
        sy = sinh(wy * L) / wy;
        ty = 1;
    else
        cy = cos(wy * L);
        sy = sin(wy * L) / wy;
        ty = -1;
    end
    
    % 一阶传输（公式24.81）
    if order >= 1
        x_new = cx * x + sx * px / (1 + delta);
        px_new = tx * wx^2 * (1 + delta) * sx * x + cx * px;
        
        y_new = cy * y + sy * py / (1 + delta);
        py_new = ty * wy^2 * (1 + delta) * sy * y + cy * py;
        
        % z的一阶项
        m511 = tx * wx^2 * (L - cx * sx) / 4;
        m512 = -tx * wx^2 * sx^2 / (2 * (1 + delta));
        m522 = -1 * (L + cx * sx) / (4 * (1 + delta)^2);
        
        m533 = ty * wy^2 * (L - cy * sy) / 4;
        m534 = -ty * wy^2 * sy^2 / (2 * (1 + delta));
        m544 = -1 * (L + cy * sy) / (4 * (1 + delta)^2);
        
        z_new = z + m511 * x^2 + m512 * x * px + m522 * px^2 + ...
                m533 * y^2 + m534 * y * py + m544 * py^2;
        
        delta_new = delta;
    end
    
    % 二阶修正
    if order >= 2
        % x方向的二阶项
        % T_166: x对delta的二阶依赖
        T_166 = -sx * px * delta / (1 + delta)^2;
        
        % T_126: x对px和delta的耦合
        T_126 = -sx / (1 + delta);
        
        x_new = x_new + T_166 * delta^2 + T_126 * px * delta;
        
        % px方向的二阶项
        % T_266: px对delta的二阶依赖
        T_266 = tx * wx^2 * sx * x * delta;
        
        % T_216: px对x和delta的耦合  
        T_216 = tx * k1 * sx / (1 + delta);
        
        px_new = px_new + T_266 * delta^2 + T_216 * x * delta;
        
        % y方向的二阶项（类似x）
        T_366 = -sy * py * delta / (1 + delta)^2;
        T_346 = -sy / (1 + delta);
        
        y_new = y_new + T_366 * delta^2 + T_346 * py * delta;
        
        % py方向的二阶项
        T_466 = ty * wy^2 * sy * y * delta;
        T_436 = ty * k1 * sy / (1 + delta);
        
        py_new = py_new + T_466 * delta^2 + T_436 * y * delta;
        
        % z的二阶项
        % T_566: z对delta的二阶依赖
        T_566 = m522 * px^2 / (1 + delta)^2 + m544 * py^2 / (1 + delta)^2;
        
        % T_516: z对x和delta的耦合
        T_516 = 2 * m511 * x / (1 + delta);
        
        % T_536: z对y和delta的耦合
        T_536 = 2 * m533 * y / (1 + delta);
        
        z_new = z_new + T_566 * delta^2 + T_516 * x * delta + T_536 * y * delta;
    end
    
    p_final = [x_new; px_new; y_new; py_new; z_new; delta_new];
end