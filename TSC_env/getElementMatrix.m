function R = getElementMatrix(element, gamma)
    % 计算元件的传输矩阵
    % 输入:
    %   element - 元件结构体
    %   gamma - 相对论因子 (可选，仅cavity需要)
    % 输出:
    %   R - 6x6传输矩阵
    
    R = eye(6);
    L = element.length;
    
    % 物理常数
    m0 = 0.511;  % 电子静止质量 (MeV)
    
    switch lower(element.type)
        case 'edge'
            h = element.h;
            angle = element.angle;
            R(2,1) = h*tan(angle);
            R(4,3) = -R(2,1);

        case 'drift'
            R(1,2) = L;
            R(3,4) = L;
            
        case 'dipole'
            h = element.h;
            if abs(h) < 1e-10
                R(1,2) = L;
                R(3,4) = L;
            else
                theta = L * h;
                cos_t = cos(theta);
                sin_t = sin(theta);
                R(1,1:2) = [cos_t, sin_t/h];
                R(3,4) = L;
                R(1,6) = (1 - cos_t)/h;
                R(2,1:2) = [-h*sin_t, cos_t];
                R(2,6) = sin_t;
                R(5,1) = -sin_t;
                R(5,2) = -(1 - cos_t)/h;
                R(5,6) = -(theta - sin_t)/h;
            end
            
        case {'quadrupole', 'quad'}
            k_actual = element.k1;
            if abs(k_actual) < 1e-10
                R(1,2) = L;
                R(3,4) = L;
            else
                if k_actual > 0
                    % 水平聚焦，垂直散焦
                    kl = sqrt(k_actual) * L;
                    R(1,1) = cos(kl);
                    R(1,2) = sin(kl)/sqrt(k_actual);
                    R(2,1) = -sqrt(k_actual)*sin(kl);
                    R(2,2) = cos(kl);
                    
                    R(3,3) = cosh(kl);
                    R(3,4) = sinh(kl)/sqrt(k_actual);
                    R(4,3) = sqrt(k_actual)*sinh(kl);
                    R(4,4) = cosh(kl);
                elseif k_actual < 0
                    % 水平散焦，垂直聚焦
                    kl = sqrt(-k_actual) * L;
                    R(1,1) = cosh(kl);
                    R(1,2) = sinh(kl)/sqrt(-k_actual);
                    R(2,1) = sqrt(-k_actual)*sinh(kl);
                    R(2,2) = cosh(kl);
                    
                    R(3,3) = cos(kl);
                    R(3,4) = sin(kl)/sqrt(-k_actual);
                    R(4,3) = -sqrt(-k_actual)*sin(kl);
                    R(4,4) = cos(kl);
                end
            end
            
        case {'sextupole', 'sext'}
            % 简化处理，只考虑长度效应
            R(1,2) = L;
            R(3,4) = L;
            
        case 'marker'
            % Marker不改变矩阵

        case {'cavity', 'rfcavity', 'rf'}
            % 检查是否提供了gamma参数
            if nargin < 2 || isempty(gamma)
                error('计算cavity传输矩阵需要提供gamma（相对论因子）参数');
            end
            
            % 获取加速腔参数
            if isfield(element, 'gradient')
                gradient = element.gradient;  % 梯度 (MV/m)
            else
                gradient = 0;
            end
           matlab
% 计算 krf (frequency 单位: MHz)
krf = 2 * pi * element.frequency / 299.792458;

% 计算 R65

            if isfield(element, 'phase')
                phase_deg = element.phase;    % 相位 (度)
            else
                phase_deg = 0;
            end
            
            % 转换相位为弧度
            Delta_phi = phase_deg * pi/180;
            
            % 计算能量增益
            if isfield(element, 'voltage')
                % 如果直接给出电压
                energy_gain = element.voltage * cos(Delta_phi);
            else
                % 从梯度和长度计算
                energy_gain = gradient * L * cos(Delta_phi);
            end
            
            % 计算入口和出口的相对论因子
            gamma_i = gamma;
            energy_i = gamma_i * m0;
            energy_f = energy_i + energy_gain;
            gamma_f = energy_f / m0;
            
            % 防止能量变为负值或零
            if gamma_f <= 0
                warning('Cavity: 出口能量非正值，gamma_f=%f，设置为入口值', gamma_f);
                gamma_f = gamma_i;
            end
            
            % 计算能量增长梯度 (单位长度的gamma增量)
            if L > 0
                gamma_prime = (gamma_f - gamma_i) / L;
            else
                gamma_prime = 0;
            end
            
            % 计算alpha参数 - 相位推进
            if abs(cos(Delta_phi)) > 1e-10 && abs(gamma_f/gamma_i - 1) > 1e-10
                alpha = (1/(2*sqrt(2)*cos(Delta_phi))) * log(gamma_f/gamma_i);
            else
                alpha = 0;
            end
            
            % eta参数 (聚焦强度因子)
            eta = 1;
            
            % 构建2x2横向传输矩阵
            if abs(alpha) > 1e-10 && abs(gamma_prime) > 1e-10
                % 使用完整公式
                M11 = cos(alpha) - sqrt(2/eta) * cos(Delta_phi) * sin(alpha);
                M12 = sqrt(8/eta) * (gamma_i/gamma_prime) * cos(Delta_phi) * sin(alpha);
                M21 = -(gamma_prime/gamma_i) * (cos(Delta_phi)/sqrt(2*eta) + sqrt(eta/8)/cos(Delta_phi)) * sin(alpha);
                M22 = (gamma_i/gamma_f) * (cos(alpha) + sqrt(2/eta) * cos(Delta_phi) * sin(alpha));
            else
                % 退化为绝热阻尼近似
                r = sqrt(gamma_i / gamma_f);
                M11 = r;
                M12 = 0;
                M21 = 0;
                M22 = 1/r;
            end
            
            % 构建6x6传输矩阵
            % 能量增益因子
            dgam = gamma_f / gamma_i;
            
            % 横向传输 (包含能量修正的平方根)
            % 注意：这里不应用sqrt(dgam)，因为完整公式中M22已包含gamma_i/gamma_f
            R(1:2, 1:2) =sqrt(1)*[M11, M12; (gamma_i / gamma_f)*M21, M22];
            R(3:4, 3:4) =sqrt(1)*[M11, M12; (gamma_i / gamma_f)*M21, M22];  % x和y方向相同
            
            % 纵向传输
            R(5,5) = 1;
            R(5,6) = 0;  % 简化处理
            R(6,5) = (1 - gamma_i / gamma_f) * krf * tan(Delta_phi);
            
            % 能量变化信息存储在R(6,6)
            % 注意：这里用倒数表示，调用者需要用1/R(6,6)获取能量增益比
            R(6,6) = gamma_i / gamma_f;
            
        case 'tgu' 
            gam = element.gamma;  % 能量
            alp = element.alp;    % 磁场梯度
            alpc = element.alpc;  % 矫正线圈梯度
            K0 = element.k0;
            ku = element.ku;

            if abs(alp - alpc) < 1e-10
                % 当 alpc ≈ alp 时，TGU退化为漂移段
                R(1,2) = L;
                R(3,4) = L;
            else
                % 计算辅助参数
                eta_c = 1 / (alp - alpc);
                eta_alpha = (2 + K0^2) / (alp * K0^2);
                
                % 计算kx和ky
                kx = sqrt((alp * K0^2) / (2 * gam^2) * (alp - alpc));
                ky = sqrt(K0^2 / (2 * gam^2) * (ku^2 + alp^2 + alp * alpc));
                
                % 检查kx是否为实数
                if (alp - alpc) < 0
                    warning('TGU: alp - alpc < 0, kx becomes imaginary. Check parameters.');
                end
                
                % 计算三角函数值
                cos_kxz = cos(kx * L);
                sin_kxz = sin(kx * L);
                cos_kyz = cos(ky * L);
                sin_kyz = sin(ky * L);
                
                % 构建传输矩阵R的各个元素
                % 第一行
                R(1,1) = cos_kxz;
                R(1,2) = (1/kx) * sin_kxz;
                R(1,6) = eta_c * (1 - cos_kxz);
                
                % 第二行
                R(2,1) = -kx * sin_kxz;
                R(2,2) = cos_kxz;
                R(2,6) = eta_c * kx * sin_kxz;
                
                % 第三行
                R(3,3) = cos_kyz;
                R(3,4) = (1/ky) * sin_kyz;
                
                % 第四行
                R(4,3) = -ky * sin_kyz;
                R(4,4) = cos_kyz;
                
                % 第五行
                R(5,1) = -eta_c * kx * sin_kxz;
                R(5,2) = -eta_c * (1 - cos_kxz);
                R(5,6) = kx^2 * eta_c * eta_alpha * L - eta_c^2 * kx * (kx * L - sin_kxz);
                
                % 第六行 (R(6,6) = 1，其他元素保持默认值0)
            end

        otherwise
            % 默认当作漂移段处理
            R(1,2) = L;
            R(3,4) = L;
    end
end