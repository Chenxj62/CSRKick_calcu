function [ Tijk] = nonlinear_calculator(beamline,iindex,jindex, kindex,delta)
    % 计算加速器束线的非线性色散项T166, T566, T266
    %
    % 输入参数:
    %   beamline - 元件数组 [ELE.D1, ELE.Q1, ...] 
    %   delta - 计算色散时使用的能散值
    %
    % 输出参数:
    %   T166 - x相对于delta的二阶传输矩阵项
    %   T566 - z相对于delta的二阶传输矩阵项
    %   T266 - px相对于delta的二阶传输矩阵项
    
    % 默认参数
    if nargin < 2
        delta = 0.001;  % 默认能散值
    end

    
    
    % 创建三个粒子：中心粒子和两个能散粒子
    p1 = [0; 0; 0; 0; 0; 0];        % 中心粒子 (x, px, y, py, z, delta)
    p2 = [0; 0; 0; 0; 0;0];    % 正能散粒子
    p3 = [0; 0; 0; 0; 0;0];   % 负能散粒子
    p4 = [0; 0; 0; 0; 0;0]; 



    % 显示开始计算
if jindex==kindex
    p1(jindex)=delta;
    p2(jindex)=-delta;
    % 追踪粒子穿过每个元件
    for i = 1:length(beamline)
        element = beamline(i);
        
        % 根据元件类型跟踪粒子
        if strcmpi(element.type, 'drift')
        
            p1 = trackDrift(p1, element.length);
            p2 = trackDrift(p2, element.length);
           
        elseif strcmpi(element.type, 'edge')
        
            p1 = trackEdge(p1, element.h,element.angle);
            p2 = trackEdge(p2, element.h,element.angle);
         
        elseif strcmpi(element.type, 'marker')
            % 标记点不影响粒子轨迹，直接跳过
            continue;
            
        elseif strcmpi(element.type, 'quadrupole')

            if abs(element.k1)<1e-4

            p1 = trackDrift(p1, element.length);
            p2 = trackDrift(p2, element.length);

            else

            p1 = trackQuadrupole(p1, element.length, element.k1);
            p2 = trackQuadrupole(p2, element.length, element.k1);
            end
            
        elseif strcmpi(element.type, 'sextupole')

            p1 = trackSextupole(p1, element.length, element.k2);
            p2 = trackSextupole(p2, element.length, element.k2);
            
            
        elseif strcmpi(element.type, 'dipole') 
          % 调试 - 打印弯铁元件的所有字段

fields = fieldnames(element);


% 从h字段提取曲率和半径
if isfield(element, 'h') && abs(element.h) > 1e-10
    g = element.h;         % 曲率 (m^-1)
    rho = 1/g;             % 半径 (m)
 
else
    % 如果没有h字段或h接近零，使用默认值
    rho = 5.0;             % 默认半径 5米
    g = 1/rho;             % 计算曲率
    warning('找不到有效的曲率字段h，使用默认半径%.2f米', rho);
end

% 检查角度字段
if isfield(element, 'angle') && element.angle ~= 0
    angle = element.angle;
else
    % 如果没有角度字段，从曲率和长度计算
    angle = g * element.length;
end

% 确保我们有合理的弧长
L =rho*angle;


% 获取k1值（如果有）
k1 = 0;
if isfield(element, 'k1')
    k1 = element.k1;
end


            
            % 使用BMAD方法跟踪粒子
            p1 = trackParticleK1Zero(p1, g, L, k1);
            p2 = trackParticleK1Zero(p2, g, L, k1);
           
        else
            warning('未知元件类型 %s，跳过该元件', element.type);
        end
    end
    

    % 计算二阶色散项
Tijk=(p1(iindex)+p2(iindex))/2/delta^2;
    
else
   
    p1(jindex)=delta; p1(kindex)=delta;
    p2(jindex)=0; p2(kindex)=delta;
    p3(jindex)=delta; p3(kindex)=0;
   
    % 追踪粒子穿过每个元件
    for i = 1:length(beamline)
        element = beamline(i);
        
        % 根据元件类型跟踪粒子
        if strcmpi(element.type, 'drift')
        
            p1 = trackDrift(p1, element.length);
            p2 = trackDrift(p2, element.length);
            p3 = trackDrift(p3, element.length);
         
        elseif strcmpi(element.type, 'marker')
            % 标记点不影响粒子轨迹，直接跳过
            continue;
        elseif strcmpi(element.type, 'edge')
                     
         
        
            p1 = trackEdge(p1, element.h,element.angle);
            p2 = trackEdge(p2, element.h,element.angle);
            p3 = trackEdge(p3, element.h,element.angle);
            
        elseif strcmpi(element.type, 'quadrupole')

            if abs(element.k1)<1e-4

            p1 = trackDrift(p1, element.length);
            p2 = trackDrift(p2, element.length);
            p3 = trackDrift(p3, element.length);

            else

            p1 = trackQuadrupole(p1, element.length, element.k1);
            p2 = trackQuadrupole(p2, element.length, element.k1);
            p3 = trackQuadrupole(p3, element.length, element.k1);
            end


           
            
        elseif strcmpi(element.type, 'sextupole')

            p1 = trackSextupole(p1, element.length, element.k2);
            p2 = trackSextupole(p2, element.length, element.k2);
            p3 = trackSextupole(p3, element.length, element.k2);
            
        elseif strcmpi(element.type, 'dipole') 
          % 调试 - 打印弯铁元件的所有字段

fields = fieldnames(element);


% 从h字段提取曲率和半径
if isfield(element, 'h') && abs(element.h) > 1e-10
    g = element.h;         % 曲率 (m^-1)
    rho = 1/g;             % 半径 (m)
 
else
    % 如果没有h字段或h接近零，使用默认值
    rho = 5.0;             % 默认半径 5米
    g = 1/rho;             % 计算曲率
    warning('找不到有效的曲率字段h，使用默认半径%.2f米', rho);
end

% 检查角度字段
if isfield(element, 'angle') && element.angle ~= 0
    angle = element.angle;
else
    % 如果没有角度字段，从曲率和长度计算
    angle = g * element.length;
end

% 确保我们有合理的弧长
L =rho*angle;


% 获取k1值（如果有）
k1 = 0;
if isfield(element, 'k1')
    k1 = element.k1;
end


            
            % 使用BMAD方法跟踪粒子
            p1 = trackParticleK1Zero(p1, g, L, k1);
            p2 = trackParticleK1Zero(p2, g, L, k1);
            p3 = trackParticleK1Zero(p3, g, L, k1);
           
        else
            warning('未知元件类型 %s，跳过该元件', element.type);
        end
    end
    

    % 计算二阶色散项
T1D1=(p1(iindex)-p2(iindex))/delta;
t2D1=(p3(iindex)-p4(iindex))/delta;

Tijk=(T1D1-t2D1)/delta;

end
end

function p_final = trackDrift(p_initial, L)
    % 漂移段跟踪函数 - 基于BMAD手册公式24.60-24.61
    
    % 提取初始相空间坐标
    x1 = p_initial(1);
    px1 = p_initial(2);
    y1 = p_initial(3);
    py1 = p_initial(4);
    z1 = p_initial(5);
    delta = p_initial(6);
    
    % 计算纵向动量 pl (公式24.61)
    pl = sqrt(1 - (px1^2 + py1^2)/(1 + delta)^2);
    
    % 计算标准化速度 β (假设参考粒子 βref = 1)
    beta = 1;  % 相对论β = v/c
    beta_ref = 1.0;                          % 参考粒子β = 1（相对论极限）
    
    % 使用公式24.60计算漂移后的坐标
    x2 = x1 + L * px1 / ((1 + delta) * pl);
    px2 = px1;
    
    y2 = y1 + L * py1 / ((1 + delta) * pl);
    py2 = py1;
    
    % 计算纵向位置
    z2 = z1 + (beta/beta_ref - 1/pl) * L;
    
    % 能散保持不变
    delta2 = delta;
    
    % 返回最终坐标
    p_final = [x2; px2; y2; py2; z2; delta2];
end

function p_final = trackEdge(p_initial, h,angle)
    % 漂移段跟踪函数 - 基于BMAD手册公式24.60-24.61
    
    
    % 提取初始相空间坐标
    x1 = p_initial(1);
    px1 = p_initial(2);
    y1 = p_initial(3);
    py1 = p_initial(4);
    z1 = p_initial(5);
    delta = p_initial(6);
    

    
    % 使用公式24.60计算漂移后的坐标
    x2 = x1 ;
    px2 = x1*h*tan(angle)+px1;

    y2 = y1 ;
    py2 = -y1*h*tan(angle)+py1;
    

    % 计算纵向位置
    z2 = z1 ;
    
    % 能散保持不变
    delta2 = delta;
    
    % 返回最终坐标
    p_final = [x2; px2; y2; py2; z2; delta2];
end

function p_final = trackQuadrupole(p_initial, L, k1)
    % 四极铁跟踪函数 - 基于BMAD手册公式24.81-24.84
    
    % 提取初始相空间坐标
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
    
    % 计算传输后的坐标 (公式24.81)
    x2 = cx * x1 + sx * px1 / (1 + delta);
    px2 = tx * wx^2 * (1 + delta) * sx * x1 + cx * px1;
    
    y2 = cy * y1 + sy * py1 / (1 + delta);
    py2 = ty * wy^2 * (1 + delta) * sy * y1 + cy * py1;
    
    % 计算m参数 (公式24.84)
    m511 = tx * wx^2 * (L - cx * sx) / 4;
    m512 = -tx * wx^2 * sx^2 / (2 * (1 + delta));
    m522 = -1 * (L + cx * sx) / (4 * (1 + delta)^2);
    
    m533 = ty * wy^2 * (L - cy * sy) / 4;
    m534 = -ty * wy^2 * sy^2 / (2 * (1 + delta));
    m544 = -1 * (L + cy * sy) / (4 * (1 + delta)^2);
    
    % 计算纵向位置
    z2 = z1 + m511 * x1^2 + m512 * x1 * px1 + m522 * px1^2 + ...
         m533 * y1^2 + m534 * y1 * py1 + m544 * py1^2;
    
    % 能散保持不变
    delta2 = delta;
    
    % 返回最终相空间坐标
    p_final = [x2; px2; y2; py2; z2; delta2];
end

function p_final = trackSextupole(p_initial, L, k2)
    % 六极铁跟踪函数 - 使用drift-kick-drift模型
    
    % 提取初始相空间坐标
    x = p_initial(1);
    px = p_initial(2);
    y = p_initial(3);
    py = p_initial(4);
    z = p_initial(5);
    delta = p_initial(6);
    
    % 半长度
    L_half = L / 2;
    
    % 初始化输出
    p_final = p_initial;
    
    % 第一个漂移(L/2)
    [x_new, y_new, z_new] = bmadDrift(x, px, y, py, z, delta, L_half);
    
    p_final(1) = x_new;
    p_final(3) = y_new;
    p_final(5) = z_new;
    
    % 六极铁踢度(完整强度)
    x = p_final(1);
    px = p_final(2);
    y = p_final(3);
    py = p_final(4);
    
    % 六极铁踢度方程(薄透镜近似)，使用全长
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
    % BMAD漂移段计算方法
    % beta_ref 和 beta 都设为1
    
    % 计算纵向动量
    pl = sqrt(1 - (px^2 + py^2)/(1 + delta)^2);
    
    % 计算位置变化
    x_out = x + L * px / ((1 + delta) * pl);
    y_out = y + L * py / ((1 + delta) * pl);
    
    % 计算纵向位置变化 (beta=1, beta_ref=1)
    z_out = z + (1.0 - 1.0/pl) * L;
end

function p_final = trackParticleK1Zero(p_initial, g, L, k1)
    % 使用公式24.47实现BMAD手册中的弯铁跟踪方法(k1=0)
    k1=0;
    % 提取初始相空间坐标
    x = p_initial(1);    % 初始x位置
    px = p_initial(2);   % 初始px动量
    y = p_initial(3);    % 初始y位置
    py = p_initial(4);   % 初始py动量
    z = p_initial(5);    % 初始z位置
    delta = p_initial(6); % 能散
    
    % 计算实际曲率gtot - 正确处理能散效应
    gtot = g;  % 实际曲率与能散相关，能量越高曲率越小
    
    % 计算中间参数
    [m5, m51, m52, m511, m512, m522, m533, m534, m544, xc] = calculateMParameters(g, gtot, L, delta);
    
    % 计算kx
    kx = k1 + g*gtot;  % 设计曲率
    
    % 计算omega和相关参数
    if kx > 0
        wx = sqrt(abs(kx) / (1 + delta));
        cx = cos(wx * L);
        sx = sin(wx * L) / wx;
        tx = -1;
    else
        wx = sqrt(abs(kx) / (1 + delta));
        cx = cosh(wx * L);
        sx = sinh(wx * L) / wx;
        tx = 1;
    end
    
    % 对于k1=0，处理ky
    ky = k1;  % k1=0
    
        wy = 0;
        cy = 1;
        sy = L;
        ty = 1;
    
    % 计算最终坐标（公式24.47）
    x_final = cx * (x - xc) + sx * px / (1 + delta) + xc;
    px_final = tx * wx^2 * (1 + delta) * sx * (x - xc) + cx * px;
    
    y_final = cy * y + sy * py / (1 + delta);
    py_final = ty * wy^2 * (1 + delta) * sy * y + cy * py;
    
    % 计算z坐标
    z_final = z + m5 + m51 * (x - xc) + m52 * px + m511 * (x - xc)^2 + ...
              m512 * (x - xc) * px + m522 * px^2 + m533 * y^2 + m534 * y * py + m544 * py^2;
    
    % 能散保持不变
    delta_final = delta;
    
    % 返回最终坐标
    p_final = [x_final; px_final; y_final; py_final; z_final; delta_final];
end

function [m5, m51, m52, m511, m512, m522, m533, m534, m544, xc] = calculateMParameters(g, gtot, L, delta)
    % 计算公式24.50中的m参数
    % g是设计曲率，gtot是实际曲率
    
    % 处理k1=0的特殊情况
    gtot = g;  % 实际曲率与能散相关，能量越高曲率越小
    k1 = 0;
    kx = k1 + g*gtot;
    
    % 对于kx接近零的情况进行特殊处理
    if abs(kx) < 1e-10
        % 当kx接近零时，使用泰勒展开近似
        xc = (g * (1 + delta) - gtot) / 1e-10;  % 避免除以零
        
        % 漂移段近似
        wx = 1e-5;
        cx = 1;
        sx = L;
        tx = 1;
    else
        xc = (g * (1 + delta) - gtot) / kx;
        
        if kx > 0
            wx = sqrt(abs(kx) / (1 + delta));
            cx = cos(wx * L);
            sx = sin(wx * L) / wx;
            tx = -1;
        else
            wx = sqrt(abs(kx) / (1 + delta));
            cx = cosh(wx * L);
            sx = sinh(wx * L) / wx;
            tx = 1;
        end
    end
    
    % 对于k1=0，处理ky
    ky = k1;  % k1=0
    
        wy = 0;
        cy = 1;
        sy = L;
        ty = 1;
    
    % 计算m参数（根据公式24.50）
    m5 = -g * xc * L;
    m51 = -g * sx;
    
    m52 = tx * g * (1 - cx) / ((1 + delta) * wx^2);
    
    m511 = tx * wx^2 * (L - cx * sx) / 4;
    
    m512 = -tx * wx^2 * sx^2 / (2 * (1 + delta));
    
    m522 = -1 * (L + cx * sx) / (4 * (1 + delta)^2);
    
    m533 = ty * wy^2 * (L - cy * sy) / 4;
    
    m534 = -ty * wy^2 * sy^2 / (2 * (1 + delta));
    
    m544 = -1 * (L + cy * sy) / (4 * (1 + delta)^2);
end