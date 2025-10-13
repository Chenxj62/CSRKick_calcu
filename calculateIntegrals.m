function [bend_integrals, segment_integrals, min_R16] = calculateIntegrals(beamline, i1, i2, h_value, calculate_min_R16)
    % 计算特定积分，从起点(s_i)到终点(s_f)，使用解析和数值结合的方法计算积分
    %
    % 输入参数:
    %   beamline: 束线结构体数组
    %   i1: 弯铁内的i值
    %   i2: 段间的i值 
    %   h_value: h值，用于计算C(s)=1/(1-h*R56(s))
    %   calculate_min_R16: 是否计算R16最小值（1表示计算，其他值表示不计算）
    %
    % 输出参数:
    %   bend_integrals: 弯铁内的积分结果 [R51积分, R52积分]
    %   segment_integrals: 段间的积分结果 [R51积分, R52积分]
    %   min_R16: R16的最小值（如果calculate_min_R16=1），否则为0
    
    % 处理默认参数
    if nargin < 5
        calculate_min_R16 = 0;
    end
    
    % 初始化返回值
    min_R16 = 0;
    
    % 初始化积分结果
    intB_R51 = 0; intB_R52 = 0;  % 弯铁积分
    intD_R51 = 0; intD_R52 = 0;  % 段间积分
    
    % 初始化传输矩阵和位置
    R_cumulative = eye(6);
    cumulative_s = 0;
    
    % 用于存储R16值（如果需要计算最小值）
    if calculate_min_R16 == 1
        all_R16_values = [];
    end
    
    % 遍历每个元件
    for elem_idx = 1:length(beamline)
        element = beamline(elem_idx);
        
        % 判断当前元件是否为弯铁
        is_bend = isBendElement(element);
        
        % 根据元件类型选择积分方法
        if is_bend
            % 弯铁：使用数值积分
            [dR51, dR52, R16_vals] = integrateBendElement(element, R_cumulative, h_value, i1, calculate_min_R16);
            intB_R51 = intB_R51 + dR51;
            intB_R52 = intB_R52 + dR52;
        else
            % 非弯铁：使用解析积分
            [dR51, dR52, R16_vals] = integrateDriftElement(element, R_cumulative, h_value, i2, calculate_min_R16);
            intD_R51 = intD_R51 + dR51;
            intD_R52 = intD_R52 + dR52;
        end
        
        % 收集R16值
        if calculate_min_R16 == 1 && ~isempty(R16_vals)
            all_R16_values = [all_R16_values; R16_vals];
        end
        
        % 更新累积传输矩阵
        R_element = getElementMatrix(element);
        R_cumulative = R_element * R_cumulative;
        cumulative_s = cumulative_s + element.length;
    end
    
    % 设置返回值
    bend_integrals = [intB_R51, intB_R52];
    segment_integrals = [intD_R51, intD_R52];
    
    % 计算R16最小值
    if calculate_min_R16 == 1 && ~isempty(all_R16_values)
        min_R16 = min(abs(all_R16_values));
    end
end

function is_bend = isBendElement(element)
    % 判断元件是否为弯铁
    % 可以根据element的类型字段或弯转角度来判断
    
    if isfield(element, 'type')
        % 如果有类型字段
        bend_types = {'bend', 'dipole', 'sbend', 'rbend', 'sector', 'rectangular'};
        is_bend = any(strcmpi(element.type, bend_types));
    elseif isfield(element, 'angle') && element.angle ~= 0
        % 如果有非零弯转角度
        is_bend = true;
    elseif isfield(element, 'rho') && isfinite(element.rho) && element.rho ~= 0
        % 如果有有限非零弯转半径
        is_bend = true;
    else
        % 默认为非弯铁
        is_bend = false;
    end
end

function [intR51, intR52, R16_values] = integrateBendElement(element, R_initial, h_value, power, calculate_R16)
    % 弯铁元件的数值积分
    
    % 细分弯铁元件
    n_steps = 200;  % 积分步数
    ds = element.length / n_steps;
    
    intR51 = 0;
    intR52 = 0;
    R16_values = [];
    
    % 创建部分元件用于积分
    part_element = element;
    R_current = R_initial;
    
    for i = 1:n_steps
        % 设置当前步长
        part_element.length = ds;
        
        % 计算步长中点的传输矩阵（更准确的积分）
        part_element.length = ds/2;
        R_mid = getElementMatrix(part_element) * R_current;
        
        % 在中点计算被积函数值
        R51_mid = R_mid(5, 1);
        R52_mid = R_mid(5, 2);
        R56_mid = R_mid(5, 6);
        C_mid = 1 / (1 - h_value * R56_mid);
        
        % 累积积分
        intR51 = intR51 + R51_mid * C_mid^power * ds;
        intR52 = intR52 + R52_mid * C_mid^power * ds;
        
        % 收集R16值
        if calculate_R16 == 1
            R16_values = [R16_values; R_mid(1, 6)];
        end
        
        % 更新到下一步
        part_element.length = ds;
        R_current = getElementMatrix(part_element) * R_current;
    end
end

function [intR51, intR52, R16_values] = integrateDriftElement(element, R_initial, h_value, power, calculate_R16)
    % 非弯铁元件的解析或简化积分
    
    % 对于漂移空间等线性元件，传输矩阵参数是线性变化的
    % 可以使用解析积分或简化的数值积分
    
    R16_values = [];
    
    % 计算初始和最终传输矩阵
    R_start = R_initial;
    R_element = getElementMatrix(element);
    R_end = R_element * R_initial;
    
    % 收集R16值
    if calculate_R16 == 1
        % 对于非弯铁，可以只取几个采样点
        n_samples = 10;
        part_element = element;
        R_current = R_initial;
        for i = 1:n_samples
            part_element.length = element.length * i / n_samples;
            R_temp = getElementMatrix(part_element) * R_initial;
            R16_values = [R16_values; R_temp(1, 6)];
        end
    end
    
    if strcmp(element.type, 'drift') || element.length == 0
        % 漂移空间：传输矩阵线性变化，可以解析积分
        % 使用中点法则作为近似
        L = element.length;
        part_element = element;
        part_element.length = L/2;
        R_mid = getElementMatrix(part_element) * R_initial;
        
        R51_mid = R_mid(5, 1);
        R52_mid = R_mid(5, 2);
        R56_mid = R_mid(5, 6);
        C_mid = 1 / (1 - h_value * R56_mid);
        
        intR51 = R51_mid * C_mid^power * L;
        intR52 = R52_mid * C_mid^power * L;
    else
        % 其他非弯铁元件：使用数值积分但步数较少
        n_steps = 50;
        ds = element.length / n_steps;
        
        intR51 = 0;
        intR52 = 0;
        
        part_element = element;
        R_current = R_initial;
        
        for i = 1:n_steps
            part_element.length = ds/2;
            R_mid = getElementMatrix(part_element) * R_current;
            
            R51_mid = R_mid(5, 1);
            R52_mid = R_mid(5, 2);
            R56_mid = R_mid(5, 6);
            C_mid = 1 / (1 - h_value * R56_mid);
            
            intR51 = intR51 + R51_mid * C_mid^power * ds;
            intR52 = intR52 + R52_mid * C_mid^power * ds;
            
            part_element.length = ds;
            R_current = getElementMatrix(part_element) * R_current;
        end
    end
end