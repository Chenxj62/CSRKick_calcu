function [final_beta, evolution, periodic_twiss] = calcBetaEvolution(beamline, initial_energy, initial_beta, plot_flag,calc_periodic)
% calcBetaEvolution - 计算束线Beta函数演化和周期解
%
% 输入参数:
%   beamline - 元件结构体数组
%   initial_energy - 初始能量(MeV)
%   initial_beta - 初始Beta参数 [betax, betay, alphax, alphay]
%   calc_periodic - 是否计算周期解 (可选，默认false)
%   plot_flag - 是否绘图 (可选，默认true)
%
% 输出参数:
%   final_beta - 最终Beta参数 [betax, betay, alphax, alphay]
%   evolution - 束线Beta演化记录结构体
%   periodic_twiss - 周期Twiss参数 [betax, betay, alphax, alphay]

% 默认参数
if nargin < 4
    plot_flag = false;
end
if nargin < 5
   calc_periodic = true;
end

% 物理常数
m0 = 0.511;    % 电子静止质量 (MeV)

% 初始化演化记录 - 增加容量以容纳更多采样点
num_elements = length(beamline);
num_slices = 10;  % 每个元件的切片数
max_points = num_elements * (num_slices + 2) + 100; % 额外空间用于cavity等复杂元件

evolution = struct();
evolution.s = zeros(max_points, 1);
evolution.energy = zeros(max_points, 1);
evolution.beta_x = zeros(max_points, 1);
evolution.beta_y = zeros(max_points, 1);
evolution.alpha_x = zeros(max_points, 1);
evolution.alpha_y = zeros(max_points, 1);
evolution.element = cell(max_points, 1);

% 初始值
evolution.s(1) = 0;
evolution.energy(1) = initial_energy;
evolution.beta_x(1) = initial_beta(1);
evolution.beta_y(1) = initial_beta(2);
evolution.alpha_x(1) = initial_beta(3);
evolution.alpha_y(1) = initial_beta(4);
evolution.element{1} = 'Start';

% 当前状态
current_beta = initial_beta;
current_energy = initial_energy;
current_s = 0;
point_index = 2;

% 初始化总传输矩阵
M_total_x = eye(2);
M_total_y = eye(2);

% 直接处理beamline结构体数组中的每个元件
for i = 1:length(beamline)
    element = beamline(i);
    
    % 确保元素有必需的字段
    if ~isfield(element, 'type')
        warning('元素 #%d 没有type字段，将被跳过', i);
        continue;
    end
    element_type = element.type;
    
    if ~isfield(element, 'length')
        warning('元素 %s 没有length字段，将设为0', element_type);
        element.length = 0;
    end
    element_length = element.length;
    
    if ~isfield(element, 'name')
        element.name = sprintf('Unnamed_%s_%d', element_type, i);
    end
    element_name = element.name;
    
    % 提取当前Beta参数
    beta_x = current_beta(1);
    beta_y = current_beta(2);
    alpha_x = current_beta(3);
    alpha_y = current_beta(4);
    gamma_x = (1 + alpha_x^2) / beta_x;
    gamma_y = (1 + alpha_y^2) / beta_y;
    
    % 记录元件起点
    evolution.s(point_index) = current_s;
    evolution.energy(point_index) = current_energy;
    evolution.beta_x(point_index) = beta_x;
    evolution.beta_y(point_index) = beta_y;
    evolution.alpha_x(point_index) = alpha_x;
    evolution.alpha_y(point_index) = alpha_y;
    evolution.element{point_index} = [element_name '_start'];
    point_index = point_index + 1;
    
    % 根据元件类型处理
    switch element_type
        case 'cavity'
            % 处理加速腔 - 保持原有实现
            gradient = element.gradient;      % 梯度 (MV/m)
            n_cells = element.num_cells;      % cell数
            frequency = element.frequency;    % 频率 (MHz)
            phase_deg = element.phase;        % 相位 (度)
            
            % 计算单个cell长度 (半个波长)
            wavelength = 299.792458 / frequency; % 波长 (m), c=299,792,458 m/s
            cell_length = wavelength / 2;     % 半个波长
            
            % 计算每个cell的能量增益
            Delta_phi = phase_deg * pi/180;   % 相位(弧度)
            energy_gain_per_cell = gradient * cell_length * cos(Delta_phi); % MeV/cell
            
            % 记录初始值
            cell_s = 0;
            cell_energy = current_energy;
            cell_beta_x = beta_x;
            cell_beta_y = beta_y;
            cell_alpha_x = alpha_x;
            cell_alpha_y = alpha_y;
            cell_gamma_x = gamma_x;
            cell_gamma_y = gamma_y;
            
            % 初始化腔体总传输矩阵
            M_cavity_x = eye(2);
            M_cavity_y = eye(2);
            
            % 常数
            rest_energy = m0;  % 电子静质能 (MeV)
            
            % 遍历每个cell
            for j = 1:n_cells
                % 计算相对论因子
                gamma_i = cell_energy / rest_energy;
                gamma_f = (cell_energy + energy_gain_per_cell) / rest_energy;
                
                % 计算能量增长梯度
                gamma_prime = gradient * cos(Delta_phi) / rest_energy; % 单位长度gamma增量
                
                % 计算alpha参数 - 当前cell中的能量增益
                alpha = (1/(2*sqrt(2)*cos(Delta_phi))) * log(gamma_f/gamma_i);
                
                % 计算eta值 (η = 1)
                eta = 1;
                
                % 构建更精确的传输矩阵 (单个cell的矩阵)
                M11 = cos(alpha) - (2/eta)^0.5 * cos(Delta_phi) * sin(alpha);
                M12 = (8/eta)^0.5 * (gamma_i/gamma_prime) * cos(Delta_phi) * sin(alpha);
                M21 = -(gamma_prime/gamma_i) * (cos(Delta_phi)/sqrt(2*eta) + (eta/8)^0.5/cos(Delta_phi)) * sin(alpha);
                M22 = (gamma_i/gamma_f) * (cos(alpha) + (2/eta)^0.5 * cos(Delta_phi) * sin(alpha));
                
                M = [M11, M12; M21, M22];
                
                % 更新能量
                cell_energy_old = cell_energy;
                cell_energy = cell_energy + energy_gain_per_cell;
                dgam = cell_energy / cell_energy_old;
                
                % 构建Twiss矩阵
                T_x = [cell_beta_x, -cell_alpha_x; -cell_alpha_x, cell_gamma_x];
                T_y = [cell_beta_y, -cell_alpha_y; -cell_alpha_y, cell_gamma_y];
                
                % 应用传输矩阵
                T_x = M * T_x * M';
                T_y = M * T_y * M';
                
                % 应用能量修正
                T_x = dgam * T_x;
                T_y = dgam * T_y;
                
                % 提取新的Beta参数
                cell_beta_x = T_x(1,1);
                cell_alpha_x = -T_x(1,2);
                cell_gamma_x = T_x(2,2);
                
                cell_beta_y = T_y(1,1);
                cell_alpha_y = -T_y(1,2);
                cell_gamma_y = T_y(2,2);
                
                % 更新腔体总传输矩阵
                M_cavity_x = M * M_cavity_x;
                M_cavity_y = M * M_cavity_y;
                
                % 记录每个cell的结果
                cell_s = cell_s + cell_length;
                evolution.s(point_index) = current_s + cell_s;
                evolution.energy(point_index) = cell_energy;
                evolution.beta_x(point_index) = cell_beta_x;
                evolution.beta_y(point_index) = cell_beta_y;
                evolution.alpha_x(point_index) = cell_alpha_x;
                evolution.alpha_y(point_index) = cell_alpha_y;
                evolution.element{point_index} = [element_name '_cell' num2str(j)];
                point_index = point_index + 1;
            end
            
            % 更新总传输矩阵
            M_total_x = M_cavity_x * M_total_x;
            M_total_y = M_cavity_y * M_total_y;
            
            % 更新当前状态
            current_energy = cell_energy;
            current_beta = [cell_beta_x, cell_beta_y, cell_alpha_x, cell_alpha_y];
            current_s = current_s + cell_s;
            
        otherwise
            % 对于除cavity外的所有元件，使用getElementMatrix函数
            
            % 如果是零长度元件，直接记录
            if element_length <= 0
                evolution.s(point_index) = current_s;
                evolution.energy(point_index) = current_energy;
                evolution.beta_x(point_index) = beta_x;
                evolution.beta_y(point_index) = beta_y;
                evolution.alpha_x(point_index) = alpha_x;
                evolution.alpha_y(point_index) = alpha_y;
                evolution.element{point_index} = element_name;
                point_index = point_index + 1;
            else
                % 使用切片法处理有长度的元件
                slice_length = element_length / num_slices;
                
                % 当前切片的状态
                slice_beta_x = beta_x;
                slice_beta_y = beta_y;
                slice_alpha_x = alpha_x;
                slice_alpha_y = alpha_y;
                slice_gamma_x = gamma_x;
                slice_gamma_y = gamma_y;
                slice_s = 0;
                slice_energy = current_energy;
                
                % 初始化元件总传输矩阵
                M_element_x = eye(2);
                M_element_y = eye(2);
                
                % 处理每个切片
                for slice = 1:num_slices
                    % 创建切片元件
                    slice_element = element;
                    slice_element.length = slice_length;
                    
                    % 使用getElementMatrix获取传输矩阵
                    R = getElementMatrix(slice_element);
                    
                    % 提取2x2传输矩阵
                    M_x = R(1:2, 1:2);
                    M_y = R(3:4, 3:4);
                    
                    % 处理能量变化（如果R(6,6) != 1，说明有能量变化）
                    if abs(R(6,6) - 1) > 1e-10
                        energy_ratio = 1 / R(6,6);  % 能量增益比例
                        slice_energy = slice_energy * energy_ratio;
                    end
                    
                    % 更新元件传输矩阵
                    M_element_x = M_x * M_element_x;
                    M_element_y = M_y * M_element_y;
                    
                    % 构建Twiss矩阵
                    T_x = [slice_beta_x, -slice_alpha_x; -slice_alpha_x, slice_gamma_x];
                    T_y = [slice_beta_y, -slice_alpha_y; -slice_alpha_y, slice_gamma_y];
                    
                    % 使用传输矩阵更新Twiss矩阵
                    T_x_new = M_x * T_x * M_x';
                    T_y_new = M_y * T_y * M_y';
                    
                    % 应用能量修正（如果有能量变化）
                    if abs(R(6,6) - 1) > 1e-10
                        dgam = energy_ratio;
                        T_x_new = dgam * T_x_new;
                        T_y_new = dgam * T_y_new;
                    end
                    
                    % 提取新的Beta参数
                    slice_beta_x = T_x_new(1,1);
                    slice_alpha_x = -T_x_new(1,2);
                    slice_gamma_x = T_x_new(2,2);
                    
                    slice_beta_y = T_y_new(1,1);
                    slice_alpha_y = -T_y_new(1,2);
                    slice_gamma_y = T_y_new(2,2);
                    
                    % 更新切片位置
                    slice_s = slice_s + slice_length;
                    
                    % 记录除了最后一个切片以外的中间切片
                    if slice < num_slices
                        evolution.s(point_index) = current_s + slice_s;
                        evolution.energy(point_index) = slice_energy;
                        evolution.beta_x(point_index) = slice_beta_x;
                        evolution.beta_y(point_index) = slice_beta_y;
                        evolution.alpha_x(point_index) = slice_alpha_x;
                        evolution.alpha_y(point_index) = slice_alpha_y;
                        evolution.element{point_index} = [element_name '_slice' num2str(slice)];
                        point_index = point_index + 1;
                    end
                end
                
                % 更新总传输矩阵
                M_total_x = M_element_x * M_total_x;
                M_total_y = M_element_y * M_total_y;
                
                % 记录元件的终点
                evolution.s(point_index) = current_s + element_length;
                evolution.energy(point_index) = slice_energy;
                evolution.beta_x(point_index) = slice_beta_x;
                evolution.beta_y(point_index) = slice_beta_y;
                evolution.alpha_x(point_index) = slice_alpha_x;
                evolution.alpha_y(point_index) = slice_alpha_y;
                evolution.element{point_index} = [element_name '_end'];
                point_index = point_index + 1;
                
                % 更新当前状态
                current_beta = [slice_beta_x, slice_beta_y, slice_alpha_x, slice_alpha_y];
                current_energy = slice_energy;
                current_s = current_s + element_length;
            end
    end
end

% 裁剪数组到实际大小
evolution.s = evolution.s(1:point_index-1);
evolution.energy = evolution.energy(1:point_index-1);
evolution.beta_x = evolution.beta_x(1:point_index-1);
evolution.beta_y = evolution.beta_y(1:point_index-1);
evolution.alpha_x = evolution.alpha_x(1:point_index-1);
evolution.alpha_y = evolution.alpha_y(1:point_index-1);
evolution.element = evolution.element(1:point_index-1);

% 返回最终Beta参数
final_beta = current_beta;

% 计算周期解或对称解
periodic_twiss = [];
if calc_periodic
    % 使用直接数值求解方法计算周期/对称Twiss参数
    [beta_x_per, beta_y_per, alpha_x_per, alpha_y_per] = calculatePeriodicTwiss(M_total_x, M_total_y);
    
    % 检查是否至少有一个方向有解
    has_x_solution = ~isempty(beta_x_per);
    has_y_solution = ~isempty(beta_y_per);
    
    if has_x_solution || has_y_solution
        % 构建Twiss参数数组
        if has_x_solution && has_y_solution
            % 两个方向都有解
            periodic_twiss = [beta_x_per, beta_y_per, alpha_x_per, alpha_y_per];
        else
            % 只有一个方向有解，将另一个方向设为NaN以保持数组大小一致
            if has_x_solution
                periodic_twiss = [beta_x_per, NaN, alpha_x_per, NaN];
            else
                periodic_twiss = [NaN, beta_y_per, NaN, alpha_y_per];
            end
        end
        
        % 输出结果信息
        trace_x = trace(M_total_x);
        trace_y = trace(M_total_y);
        
        if has_x_solution && has_y_solution
            fprintf('水平和垂直方向均有解:\n');
        elseif has_x_solution
            fprintf('仅水平方向有解, 垂直方向无解\n');
        elseif has_y_solution
            fprintf('仅垂直方向有解, 水平方向无解\n');
        else
            fprintf('水平和垂直方向均无解\n');
        end
        
        if has_x_solution
            if abs(trace_x) < 2
                solution_type = '周期解';
            else
                solution_type = '对称解';
            end
            fprintf('X: beta=%f, alpha=%f (%s, trace=%f)\n', beta_x_per, alpha_x_per, solution_type, trace_x);
        end
        if has_y_solution
            if abs(trace_y) < 2
                solution_type = '周期解';
            else
                solution_type = '对称解';
            end
            fprintf('Y: beta=%f, alpha=%f (%s, trace=%f)\n', beta_y_per, alpha_y_per, solution_type, trace_y);
        end
    else
        warning('无法计算有效的Twiss参数');
    end
end

% 绘图 - 只保留Beta函数绘图
if plot_flag
    % 设置自定义颜色
    deep_blue = [0, 0, 0.7];   % 深蓝色
    deep_red = [0.7, 0, 0];    % 深红色
    
    % 创建图形窗口
    figure('Position', [100, 100, 700, 400]);
    
    % 绘制beta_x和beta_y演化
    plot(evolution.s, evolution.beta_x, 'Color', deep_blue, 'LineWidth', 1.5);
    hold on;
    plot(evolution.s, evolution.beta_y, 'Color', deep_red, 'LineWidth', 1.5);
    
    % 设置图表格式
    xlabel('位置 (m)', 'FontSize', 10);
    ylabel('Beta 函数 (m)', 'FontSize', 10);
    title('Beta 函数沿束线的演化', 'FontSize', 11);
    legend('\beta_x', '\beta_y', 'Location', 'best', 'FontSize', 9);
    grid on;
    
    % 设置坐标轴范围，起始于零
    xlim([0, max(evolution.s)]);
    ylim([0, max([max(evolution.beta_x), max(evolution.beta_y)]) * 1.05]);
    
    % 设置背景为白色
    set(gcf, 'Color', 'w');
end

end

function [beta_x, beta_y, alpha_x, alpha_y] = calculatePeriodicTwiss(M_x, M_y)
    % 使用数值方法计算对称Twiss参数
    % 输入: M_x, M_y - 水平和垂直方向的传输矩阵
    % 输出: 对称Twiss参数
    
    % 计算水平方向参数
    [beta_x, alpha_x] = calculateSymmetricTwissFsolve(M_x, 'x');
    
    % 计算垂直方向参数
    [beta_y, alpha_y] = calculateSymmetricTwissFsolve(M_y, 'y');
end

function [beta, alpha] = calculateSymmetricTwissFsolve(M, direction)
    % 使用fsolve数值方法求解对称Twiss参数
    % 输入: 
    %   M - 传输矩阵
    %   direction - 方向 ('x'或'y')，用于信息输出
    % 输出: 
    %   beta, alpha - 对称Twiss参数
    
    % 默认值
    beta = [];
    alpha = [];
    
    % 检查轨迹进行系统稳定性分析
    trace_val = trace(M);
    is_stable = abs(trace_val) < 2;
    
    if is_stable
        % 稳定系统：计算周期解
        try
            % 提取矩阵元素
            c = M(1,1);
            s = M(1,2);
            cp = M(2,1);
            sp = M(2,2);
            
            % 计算通用参数
            mu = acos(0.5 * trace_val);  % 相位推进
            sin_mu = sin(mu);
            
            % 只有当sin(mu)不接近零时才有稳定解
            if abs(sin_mu) < 1e-10
                return;
            end
            
            % 计算beta
            beta_val = abs(s / sin_mu);
            
            % 计算alpha
            alpha_val = (c - sp) / (2 * sin_mu);
            
            % 验证解
            gamma_val = (1 + alpha_val^2) / beta_val;
            T = [beta_val, -alpha_val; -alpha_val, gamma_val];
            T_transformed = M * T * M';
            
            beta_transformed = T_transformed(1,1);
            alpha_transformed = -T_transformed(1,2);
            
            % 检查是否满足周期条件
            if abs(beta_transformed - beta_val) < 1e-6 && abs(alpha_transformed - alpha_val) < 1e-6
                beta = beta_val;
                alpha = alpha_val;
            end
        catch
            % 捕获计算错误
        end
    else
        % 不稳定系统：尝试计算对称解
        try
            % 检查是否有fsolve函数
            if ~exist('fsolve', 'file')
                warning('%s方向：需要Optimization Toolbox中的fsolve函数来求解对称解', upper(direction));
                return;
            end
            
            % 初始猜测值
            % 提取矩阵元素用于更好的初始猜测
            c = M(1,1);
            s = M(1,2);
            cp = M(2,1);
            sp = M(2,2);
            
            % 使用系统特性为初始猜测提供更好的值
            if abs(s) > 1e-10
                alpha0 = 0;  % 对称系统可能alpha接近于0
                beta0 = abs(s); 
            else
                % 备用猜测
                alpha0 = 0;
                beta0 = 1;
            end
            
            % 定义方程组：对于对称解，beta1=beta2, alpha1=-alpha2
            % 使用传输矩阵关系：M * T * M' = T'，其中T=[beta,-alpha;-alpha,gamma]，T'=[beta,alpha;alpha,gamma]
            % 其中gamma = (1+alpha^2)/beta
            % 定义方程组
        symTwissEq = @(params) symmetricTwissEquations(params, M);
        
        % 设置fsolve选项 - 增加最大迭代次数和函数评估次数
        options = optimoptions('fsolve', ...
            'Display', 'off', ...                         % 不显示迭代过程
            'Algorithm', 'levenberg-marquardt', ...       % 算法选择
            'MaxIterations', 1000000, ...                 % 增加最大迭代次数
            'MaxFunctionEvaluations', 1000000, ...        % 增加最大函数评估次数
            'FunctionTolerance', 1e-8, ...                % 提高函数容差精度
            'OptimalityTolerance', 1e-8, ...              % 提高最优性容差
            'StepTolerance', 1e-10);                      % 提高步长容差
        
        % 使用fsolve求解
        [result, ~, exitflag] = fsolve(symTwissEq, [beta0, alpha0], options);
            
            % 检查是否有解
            if exitflag > 0
                beta_val = result(1);
                alpha_val = result(2);
                
                % 确保beta为正值
                if beta_val <= 0
                    return;
                end
                
                % 验证解
                gamma_val = (1 + alpha_val^2) / beta_val;
                T = [beta_val, -alpha_val; -alpha_val, gamma_val];
                T_transformed = M * T * M';
                
                beta_transformed = T_transformed(1,1);
                alpha_transformed = -T_transformed(1,2);
                
                % 对称解条件：beta_transformed ≈ beta, alpha_transformed ≈ -alpha
                if abs(beta_transformed - beta_val) < 1e-4 && abs(alpha_transformed + alpha_val) < 1e-4
                    beta = beta_val;
                    alpha = alpha_val;
                end
            end
        catch e
            warning('%s方向：计算对称解时出错', upper(direction));
        end
        
        % 如果fsolve方法失败，尝试使用解析方法
        if isempty(beta) && abs(M(1,2)) > 1e-10
            try
                % 使用二次方程：M12*alpha^2 + (M11-M22)*alpha - M21 = 0
                A = M(1,2);
                B = M(1,1) - M(2,2);
                C = -M(2,1);
                
                discriminant = B^2 - 4*A*C;
                
                if discriminant >= 0
                    % 计算两个解
                    alpha1 = (-B + sqrt(discriminant))/(2*A);
                    alpha2 = (-B - sqrt(discriminant))/(2*A);
                    
                    % 计算对应的beta值
                    beta_val = abs(M(1,2) / (alpha2 - alpha1));
                    
                    % 验证解
                    if beta_val > 0
                        gamma_val = (1 + alpha1^2) / beta_val;
                        
                        % 验证关系式
                        test1 = M(1,1) * beta_val - M(1,2) * alpha1;
                        test2 = beta_val / gamma_val;
                        
                        if abs(test1 - test2) < 1e-4
                            beta = beta_val;
                            alpha = alpha1;
                        end
                    end
                end
            catch
                % 解析方法也失败
            end
        end
    end
end

function F = symmetricTwissEquations(params, M)
    % 定义对称Twiss参数方程组
    % 输入:
    %   params - [beta, alpha]
    %   M - 传输矩阵
    % 输出:
    %   F - 残差向量
    
    beta = params(1);
    alpha = params(2);
    
    % 检查参数有效性
    if beta <= 0
        % 如果beta无效，返回大的残差
        F = [1e6; 1e6];
        return;
    end
    
    % 计算gamma
    gamma = (1 + alpha^2) / beta;
    
    % 构建Twiss矩阵
    T = [beta, -alpha; -alpha, gamma];
    
    % 应用传输矩阵变换
    T_transformed = M * T * M';
    
    % 提取变换后的参数
    beta_transformed = T_transformed(1,1);
    alpha_transformed = -T_transformed(1,2);
    
    % 对称解条件：beta_transformed = beta, alpha_transformed = -alpha
    F = [beta_transformed - beta; 
         alpha_transformed + alpha];
end