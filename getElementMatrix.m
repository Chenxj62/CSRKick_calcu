function R = getElementMatrix(element)
    % 计算元件的传输矩阵
    R = eye(6);
    L = element.length;
    switch lower(element.type)
        case 'edge'
         h = element.h;
         angle=element.angle;
     
            R(2,1)=h*tan(angle);
            R(4,3)=-R(2,1);

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

                 fprintf('待补充');
            % % 使用与原函数相同的方法处理加速腔
            % 
            % % 获取参数
            % if isfield(element, 'gradient')
            %     gradient = element.gradient;  % 梯度 (MV/m)
            % else
            %     gradient = 0;
            % end
            % 
            % if isfield(element, 'phase')
            %     phase_deg = element.phase;    % 相位 (度)
            % else
            %     phase_deg = 0;
            % end
            % 
            % % 转换相位为弧度
            % Delta_phi = phase_deg * pi/180;
            % 
            % % 计算能量增益 (对于单个cell或整个腔)
            % if isfield(element, 'voltage')
            %     energy_gain = element.voltage * cos(Delta_phi);
            % else
            %     energy_gain = gradient * L * cos(Delta_phi);
            % end
            % 
            % % 计算相对论因子
            % gamma_i = energy / m0;
            % gamma_f = (energy + energy_gain) / m0;
            % 
            % % 计算能量增长梯度
            % gamma_prime = gradient * cos(Delta_phi) / m0; % 单位长度gamma增量
            % 
            % % 计算alpha参数 - 当前cell中的能量增益
            % alpha = (1/(2*sqrt(2)*cos(Delta_phi))) * log(gamma_f/gamma_i);
            % 
            % % 计算eta值 (η = 1)
            % eta = 1;
            % 
            % % 构建更精确的传输矩阵
            % M11 = cos(alpha) - (2/eta)^0.5 * cos(Delta_phi) * sin(alpha);
            % M12 = (8/eta)^0.5 * (gamma_i/gamma_prime) * cos(Delta_phi) * sin(alpha);
            % M21 = -(gamma_prime/gamma_i) * (cos(Delta_phi)/sqrt(2*eta) + (eta/8)^0.5/cos(Delta_phi)) * sin(alpha);
            % M22 = (gamma_i/gamma_f) * (cos(alpha) + (2/eta)^0.5 * cos(Delta_phi) * sin(alpha));
            % 
            % % 构建2x2传输矩阵并扩展到6x6
            % M_2x2 = [M11, M12; M21, M22];
            % 
            %   dgam= gamma_f / gamma_i;
            % 
            % % 同样的矩阵用于x和y平面
            % R(1:2, 1:2) = sqrt(dgam)*M_2x2;
            % R(3:4, 3:4) = sqrt(dgam)*M_2x2;
            % 
            % % 纵向传输部分
            % R(5,6) = 1;
            % 
            % % dgam存储在额外信息中，以便调用者应用能量修正
            % R(6,6) = gamma_i / gamma_f;  % 能量增益比例，用于dgam修正
        case 'tgu' 
            gam=element.gamma; %能量
            alp=element.alp;  %磁场梯度
            alpc=element.alpc; %矫正线圈梯度
            K0=element.k0;
            ku=element.ku;

         if abs(alp - alpc) < 1e-10
        % 当 alpc ≈ alp 时，TGU退化为漂移段
        R(1,2) = L;
        R(3,4) = L;
    else
        % 检查是否有其他必要参数


        
        % 计算辅助参数
        eta_c = 1 / (alp - alpc);
        eta_alpha = (2 + K0^2) / (alp * K0^2);
        
        % 计算kx和ky
        kx = sqrt((alp * K0^2) / (2 * gam^2) * (alp - alpc));
        ky = sqrt(K0^2 / (2 * gam^2) * (ku^2 + alp^2 + alp * alpc));
        
        % 检查kx是否为实数（当alp - alpc > 0时）
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
        R(5,6) =  kx^2 * eta_c * eta_alpha * L - eta_c^2 * kx * (kx * L - sin_kxz);
        
        % 第六行 (R(6,6) = 1，其他元素保持默认值0)
        % R(6,6) 已经在 R = eye(6) 中设置为1
         end

        otherwise
            % 默认当作漂移段处理
            R(1,2) = L;
            R(3,4) = L;
    end

    
end