function [min_Rij, max_Rij, s_min, s_max] = findRijmax(beamline, i, j)
    % 查找束线中传输矩阵元素Rij的最小值和最大值，以及它们在束线中的位置
    %
    % 输入参数:
    %   beamline: 束线结构体数组
    %   i: 传输矩阵的行索引 (1-6)
    %   j: 传输矩阵的列索引 (1-6)
    %
    % 输出参数:
    %   min_Rij: Rij的最小值
    %   max_Rij: Rij的最大值
    %   s_min: 最小值所在的位置 (m)
    %   s_max: 最大值所在的位置 (m)
    
    % 检查输入参数是否有效
    if i < 1 || i > 6 || j < 1 || j > 6
        error('索引 i 和 j 必须在 1 到 6 之间');
    end
    
    % 计算总长度
    total_length = 0;
    for k = 1:length(beamline)
        total_length = total_length + beamline(k).length;
    end
    
    % 设置扫描步数
    scan_steps = 1000;  % 增加扫描点数确保不遗漏
    
    % 初始化存储位置和Rij值的数组
    s_positions = [];
    Rij_values = [];
    
    % 初始化累积传输矩阵和长度
    R_cumulative = eye(6);
    cumulative_length = 0;
    
    % 添加起始点 (s=0)
    s_positions(end+1) = 0;
    Rij_values(end+1) = R_cumulative(i, j);
    
    % 处理每个元件
    for k = 1:length(beamline)
        element = beamline(k);
        element_length = element.length;
        
        % 添加元件入口位置（如果不是第一个元件）
        if k > 1
            s_positions(end+1) = cumulative_length;
            Rij_values(end+1) = R_cumulative(i, j);
        end
        
        if element_length > 0
            % 确定当前元件的扫描步数
            num_element_steps = max(10, round(element_length / total_length * scan_steps));
            
            % 在元件内部进行扫描（不包括入口和出口，它们单独处理）
            for scan_step = 1:(num_element_steps-1)
                s_local = scan_step * element_length / num_element_steps;
                
                % 计算当前位置的传输矩阵
                R_element_partial = getPartialElementMatrix(element, s_local);
                R_scan = R_element_partial * R_cumulative;
                
                % 记录位置和Rij值
                s_positions(end+1) = cumulative_length + s_local;
                Rij_values(end+1) = R_scan(i, j);
            end
        end
        
        % 更新累积传输矩阵和长度
        R_element = getElementMatrix(element);
        R_cumulative = R_element * R_cumulative;
        cumulative_length = cumulative_length + element_length;
        
        % 添加元件出口位置
        s_positions(end+1) = cumulative_length;
        Rij_values(end+1) = R_cumulative(i, j);
    end
    
    % 转换为列向量
    s_positions = s_positions(:);
    Rij_values = Rij_values(:);
    
    % 移除重复位置（保留数值较大的那个）
    [unique_s, unique_idx] = unique(s_positions);
    if length(unique_s) < length(s_positions)
        % 有重复位置，需要处理
        s_clean = [];
        Rij_clean = [];
        
        for idx = 1:length(unique_s)
            s_val = unique_s(idx);
            matching_indices = find(abs(s_positions - s_val) < 1e-12);
            
            if length(matching_indices) == 1
                s_clean(end+1) = s_positions(matching_indices);
                Rij_clean(end+1) = Rij_values(matching_indices);
            else
                % 有多个相同位置，取Rij绝对值最大的
                [~, max_abs_idx] = max(abs(Rij_values(matching_indices)));
                best_idx = matching_indices(max_abs_idx);
                s_clean(end+1) = s_positions(best_idx);
                Rij_clean(end+1) = Rij_values(best_idx);
            end
        end
        
        s_positions = s_clean(:);
        Rij_values = Rij_clean(:);
    else
        s_positions = s_positions(unique_idx);
        Rij_values = Rij_values(unique_idx);
    end
    
    % 找出最小值和最大值
    [min_Rij, min_idx] = min(Rij_values);
    [max_Rij, max_idx] = max(Rij_values);
    
    % 记录最小值和最大值的位置
    s_min = s_positions(min_idx);
    s_max = s_positions(max_idx);
 
end

%
% 计算部分元件的传输矩阵（对于长度s的子元件）
function R = getPartialElementMatrix(element, s)
    % 检查s是否超过元件长度
    if s > element.length
        s = element.length;
    end
    
    if s <= 0
        R = eye(6);
        return;
    end
    
    % 创建一个与原元件类型相同但长度为s的新元件
    partial_element = element;
    partial_element.length = s;
    
    % 如果是偏转磁铁，需要调整角度和其他参数
    if strcmp(element.type, 'dipole') && element.length > 0
        % 按比例调整角度
        partial_element.angle = s * element.angle / element.length;
        % h保持不变，因为它是曲率（1/rho）
    elseif strcmp(element.type, 'quadrupole') && element.length > 0
        % 四极铁的k1保持不变，因为它是梯度
    elseif strcmp(element.type, 'sextupole') && element.length > 0
        % 六极铁的k2保持不变
    end
    
    % 计算这个部分元件的传输矩阵
    R = getElementMatrix(partial_element);
end