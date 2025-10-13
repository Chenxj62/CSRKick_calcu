function [bend_integrals, segment_integrals, min_R16] = CalcCSRkick(beamline, h_value, calculate_min_R16, sigz)
    % 计算特定积分，从起点(s_i)到终点(s_f)，使用解析和数值结合的方法计算积分
    %
    % 输入参数:
    %   beamline: 束线结构体数组
    %   h_value: h值，用于计算C(s)=1/(1-h*R56(s))
    %   calculate_min_R16: 是否计算R16最小值（1表示计算，其他值表示不计算）
    %   sigz: 束团长度
    %
    % 输出参数:
    %   bend_integrals: 弯铁内的积分结果 [R51积分, R52积分]
    %   segment_integrals: 段间的积分结果 [R51积分, R52积分]
    %   min_R16: R16的最小值（如果calculate_min_R16=1），否则为0
    
    % 处理默认参数
    if nargin < 3
        calculate_min_R16 = 0;
    end
    if nargin < 4
        sigz = 1e-3;  % 默认束团长度
    end
    
    % 提取所有弯铁的曲率半径
    bend_rhos = extractAllBendRadii(beamline);
    
    % 初始化返回值
    min_R16 = 0;
    
    % 初始化积分结果
    intB_R51 = 0; intB_R52 = 0;  % 弯铁积分
    sum_segment_R51 = 0; sum_segment_R52 = 0;  % 段间积分求和（最后取绝对值）
    
    % 初始化传输矩阵和位置
    R_cumulative = eye(6);
    cumulative_s = 0;
    
    % 用于存储R16值（如果需要计算最小值）
    if calculate_min_R16 == 1
        all_R16_values = [];
    end
    
    % 用于跟踪当前段间状态
    in_segment = false;
    segment_start_R = [];
    segment_length = 0;
    current_bend_index = 0;  % 当前弯铁索引
    
    % 遍历每个元件
    for elem_idx = 1:length(beamline)
        element = beamline(elem_idx);
        
        % 判断当前元件是否为弯铁
        is_bend = isBendElement(element);
        
        if is_bend
            % 弯铁元件
            current_bend_index = current_bend_index + 1;
            
            % 如果之前在段间，先处理完这个段间
            if in_segment && segment_length > 0
                % 第i个段间使用第i+1个弯铁的半径（即当前弯铁）
                if current_bend_index <= length(bend_rhos)
                    rho_for_segment = bend_rhos(current_bend_index);
                    [dR51, dR52] = calculateSegmentContribution(segment_start_R, segment_length, h_value, sigz, rho_for_segment);
                    sum_segment_R51 = sum_segment_R51 + dR51;
                    sum_segment_R52 = sum_segment_R52 + dR52;
                end
            end
            
            % 处理弯铁积分
            [dR51, dR52, R16_vals] = integrateBendElement(element, R_cumulative, h_value, calculate_min_R16);
            intB_R51 = intB_R51 + dR51;
            intB_R52 = intB_R52 + dR52;
            
            % 收集R16值
            if calculate_min_R16 == 1 && ~isempty(R16_vals)
                all_R16_values = [all_R16_values; R16_vals];
            end
            
            % 重置段间状态
            in_segment = false;
            segment_length = 0;
            
        else
            % 非弯铁元件
            
            % 收集R16值（只处理有长度的元件）
            if calculate_min_R16 == 1 && element.length > 0
                n_samples = 10;
                part_element = element;
                for i = 1:n_samples
                    part_element.length = element.length * i / n_samples;
                    R_temp = getElementMatrix(part_element) * R_cumulative;
                    all_R16_values = [all_R16_values; R_temp(1, 6)];
                end
            end
            
            % 只有当后面还有弯铁时才处理段间
            if current_bend_index < length(bend_rhos)
                % 开始新的段间或继续当前段间
                if ~in_segment
                    % 开始新的段间
                    in_segment = true;
                    segment_start_R = R_cumulative;
                    segment_length = element.length;
                else
                    % 继续当前段间
                    segment_length = segment_length + element.length;
                end
            end
        end
        
        % 更新累积传输矩阵
        R_element = getElementMatrix(element);
        R_cumulative = R_element * R_cumulative;
        cumulative_s = cumulative_s + element.length;
    end
    
    % 最后一个段间不处理（后面没有弯铁）
    
    % 设置返回值
    bend_integrals = [intB_R51, intB_R52];
    % 对段间求和结果取绝对值
    segment_integrals = [sum_segment_R51, sum_segment_R52];
    
    % 计算R16最小值
    if calculate_min_R16 == 1 && ~isempty(all_R16_values)
        min_R16 = min(all_R16_values);
    end
end

function bend_rhos = extractAllBendRadii(beamline)
    % 提取所有弯铁的曲率半径
    bend_rhos = [];
    
    for i = 1:length(beamline)
        element = beamline(i);
        if isBendElement(element)
            if isfield(element, 'h') && element.h ~= 0
                rho = 1 / element.h;  % 保持原始符号
                bend_rhos = [bend_rhos; rho];
            else
                warning('Bend element without valid h field, using default rho = 1.0');
                bend_rhos = [bend_rhos; 1.0];
            end
        end
    end
end

function [dR51, dR52] = calculateSegmentContribution(segment_start_R, segment_length, h_value, sigz, rho_typical)
    % 计算单个段间的贡献（不取绝对值，留给上层函数求和后再取绝对值）
    
    % 段间起点的传输矩阵元素
    R51 = segment_start_R(5, 1);
    R52 = segment_start_R(5, 2);
    R56 = segment_start_R(5, 6);
    
    % 压缩因子C在段间不变，按起点计算
    C = 1 / (1 - h_value * R56);
    
    % 计算phi值
    % phi0 = (6*sigz/abs(rho_typ))^(1/3)  使用绝对值
    phi0 = (6 * sigz / abs(rho_typical))^(1/3);
    % phi_i = phi0 / C_i^(1/3)
    phi_i = phi0 / (C^(1/3));
    
    % 计算CSR修正项，使用下游弯铁的半径绝对值
    log_term = log(2 * segment_length / (abs(rho_typical) * phi_i) + 1);
    
    % 计算贡献（不取绝对值）
    dR51 = R51 * C * log_term;
    dR52 = R52 * C * log_term;
end

function is_bend = isBendElement(element)
    % 判断元件是否为弯铁
    % 可以根据element的类型字段或弯转角度来判断
    
    if isfield(element, 'type')
        % 如果有类型字段
        bend_types = {'bend', 'dipole', 'sbend', 'rbend', 'sector', 'rectangular'};
        is_bend = any(strcmpi(element.type, bend_types));
    elseif isfield(element, 'h') && abs(element.h) > 1e-10
        % 如果有非零曲率
        is_bend = true;
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

function [intR51, intR52, R16_values] = integrateBendElement(element, R_initial, h_value, calculate_R16)
    % 弯铁元件的数值积分
    
    % 获取弯铁的曲率h值
    if isfield(element, 'h')
        h_bend = element.h;
    else
        h_bend = 0;
    end
    
    % 细分弯铁元件
    n_steps = 500;  % 积分步数
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
        
        % 添加CSR修正因子：|h|^(2/3)  使用绝对值
        if abs(h_bend) > 1e-10
            csr_factor = abs(h_bend)^(2/3);
        else
            csr_factor = 1;
        end
        
        % 累积积分，包含CSR修正
        intR51 = intR51 + R51_mid * C_mid * csr_factor * ds;
        intR52 = intR52 + R52_mid * C_mid * csr_factor * ds;
        
        % 收集R16值
        if calculate_R16 == 1
            R16_values = [R16_values; R_mid(1, 6)];
        end
        
        % 更新到下一步
        part_element.length = ds;
        R_current = getElementMatrix(part_element) * R_current;
    end
end