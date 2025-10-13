function [final_value, evolution, bend_averages] = calcsliene(beamline, initial_energy, initial_beta, epsilon_x, h0, sigmaz,sigmae, plot_flag)
% calcDispersionFunction - 计算色散函数随束线的演化

% 默认参数
if nargin < 7
    plot_flag = true;
end

% 计算sigma55
sigma55 = sigmaz^2;

sigma66=sigmae^2;
epsilon_x = epsilon_x/(initial_energy/0.511);

% 初始Beta参数（用于色散函数计算）
beta_x_initial = initial_beta(1);
alpha_x_initial = initial_beta(2);
gamma_x_initial = (1 + alpha_x_initial^2) / beta_x_initial;

% 初始化弯铁平均值记录
bend_averages = struct();
bend_averages.names = {};
bend_averages.average_values = [];
bend_averages.start_positions = [];
bend_averages.end_positions = [];
bend_averages.lengths = [];

if plot_flag
    % 需要绘图时，计算完整演化
    num_elements = length(beamline);
    num_slices = 10;
    max_points = num_elements * (num_slices + 2) + 100;
    
    evolution = struct();
    evolution.s = zeros(max_points, 1);
    evolution.energy = zeros(max_points, 1);
    evolution.beta_x = zeros(max_points, 1);
    evolution.alpha_x = zeros(max_points, 1);
    evolution.R51 = zeros(max_points, 1);
    evolution.R52 = zeros(max_points, 1);
    evolution.R56 = zeros(max_points, 1);
    evolution.R16 = zeros(max_points, 1);
    evolution.dispersion_function = zeros(max_points, 1);
    evolution.numerator = zeros(max_points, 1);
    evolution.denominator = zeros(max_points, 1);
    evolution.element = cell(max_points, 1);
    
    % 初始值
    evolution.s(1) = 0;
    evolution.energy(1) = initial_energy;
    evolution.beta_x(1) = beta_x_initial;
    evolution.alpha_x(1) = alpha_x_initial;
    evolution.R51(1) = 0;
    evolution.R52(1) = 0;
    evolution.R56(1) = 0;
    evolution.R16(1) = 0;
    evolution.element{1} = 'Start';
    
    % 计算初始函数值
    [disp_func_initial, num_initial, den_initial] = calculateDispersionValue(beta_x_initial, alpha_x_initial, gamma_x_initial, ...
        0, 0, 0, epsilon_x, h0, sigma55, sigma66);
    evolution.dispersion_function(1) = disp_func_initial;
    evolution.numerator(1) = num_initial;
    evolution.denominator(1) = den_initial;
    
    % 当前状态
    current_beta_x = beta_x_initial;
    current_alpha_x = alpha_x_initial;
    current_energy = initial_energy;
    current_s = 0;
    point_index = 2;
    
    % 初始化累积传输矩阵
    R_total = eye(6);
    
    % 先识别弯铁组
    bend_groups = identifyBendGroups(beamline);
    
    % 为每个弯铁组存储所有取值
    bend_group_values = cell(length(bend_groups.start_indices), 1);
    for k = 1:length(bend_groups.start_indices)
        bend_group_values{k} = [];
    end
    
    % 遍历束线元件
    for i = 1:length(beamline)
        element = beamline(i);
        
        if ~isfield(element, 'type')
            warning('元素 #%d 没有type字段，将被跳过', i);
            continue;
        end
        element_type = element.type;
        
        if ~isfield(element, 'length')
            element.length = 0;
        end
        element_length = element.length;
        
        if ~isfield(element, 'name')
            element.name = sprintf('Unnamed_%s_%d', element_type, i);
        end
        element_name = element.name;
        
        % 跳过cavity
        if strcmp(element_type, 'cavity')
            warning('跳过cavity元件: %s', element_name);
            continue;
        end
        
        % 检查是否为弯铁组的起始
        group_info = findBendGroupInfo(i, bend_groups);
        
        if group_info.is_start
            bend_start_s = current_s;
            bend_start_point = point_index;
        end
        
        % 记录元件起点
        evolution.s(point_index) = current_s;
        evolution.energy(point_index) = current_energy;
        evolution.beta_x(point_index) = current_beta_x;
        evolution.alpha_x(point_index) = current_alpha_x;
        evolution.R51(point_index) = R_total(5,1);
        evolution.R52(point_index) = R_total(5,2);
        evolution.R56(point_index) = R_total(5,6);
        evolution.R16(point_index) = R_total(1,6);
        evolution.element{point_index} = [element_name '_start'];
        
        % 计算色散函数
        [disp_func, num_val, den_val] = calculateDispersionValue(beta_x_initial, alpha_x_initial, gamma_x_initial, ...
            R_total(5,1), R_total(5,2), R_total(5,6), epsilon_x, h0, sigma55, sigma66);
        evolution.dispersion_function(point_index) = disp_func;
        evolution.numerator(point_index) = num_val;
        evolution.denominator(point_index) = den_val;
        
        % 如果在弯铁组内，存储取值
        if group_info.in_group
            sqrt_value = sqrt(disp_func);
            bend_group_index = find(group_info.group_index == 1:length(bend_groups.start_indices));
            bend_group_values{bend_group_index}(end+1) = sqrt_value;
        end
        
        point_index = point_index + 1;
        
        if element_length <= 0
            % 零长度元件
            R_element = getElementMatrix(element);
            R_total = R_element * R_total;
            
            evolution.s(point_index) = current_s;
            evolution.energy(point_index) = current_energy;
            evolution.beta_x(point_index) = current_beta_x;
            evolution.alpha_x(point_index) = current_alpha_x;
            evolution.R51(point_index) = R_total(5,1);
            evolution.R52(point_index) = R_total(5,2);
            evolution.R56(point_index) = R_total(5,6);
            evolution.R16(point_index) = R_total(1,6);
            evolution.element{point_index} = element_name;
            
            [disp_func, num_val, den_val] = calculateDispersionValue(beta_x_initial, alpha_x_initial, gamma_x_initial, ...
                R_total(5,1), R_total(5,2), R_total(5,6), epsilon_x, h0, sigma55, sigma66);
            evolution.dispersion_function(point_index) = disp_func;
            evolution.numerator(point_index) = num_val;
            evolution.denominator(point_index) = den_val;
            
            % 如果在弯铁组内，存储取值
            if group_info.in_group
                sqrt_value = sqrt(disp_func);
                bend_group_index = find(group_info.group_index == 1:length(bend_groups.start_indices));
                bend_group_values{bend_group_index}(end+1) = sqrt_value;
            end
            
            point_index = point_index + 1;
        else
            % 有长度的元件
            slice_length = element_length / num_slices;
            slice_beta_x = current_beta_x;
            slice_alpha_x = current_alpha_x;
            slice_energy = current_energy;
            R_element = eye(6);
            
            for slice = 1:num_slices
                slice_element = element;
                slice_element.length = slice_length;
                
                R_slice = getElementMatrix(slice_element);
                R_element = R_slice * R_element;
                
                % 更新Beta参数
                M_x = R_slice(1:2, 1:2);
                T_x = [slice_beta_x, -slice_alpha_x; -slice_alpha_x, (1+slice_alpha_x^2)/slice_beta_x];
                T_x_new = M_x * T_x * M_x';
                slice_beta_x = T_x_new(1,1);
                slice_alpha_x = -T_x_new(1,2);
                
                if abs(R_slice(6,6) - 1) > 1e-10
                    slice_energy = slice_energy / R_slice(6,6);
                end
                
                % 记录中间切片
                if slice < num_slices
                    evolution.s(point_index) = current_s + slice * slice_length;
                    evolution.energy(point_index) = slice_energy;
                    evolution.beta_x(point_index) = slice_beta_x;
                    evolution.alpha_x(point_index) = slice_alpha_x;
                    
                    R_temp = R_element * R_total;
                    evolution.R51(point_index) = R_temp(5,1);
                    evolution.R52(point_index) = R_temp(5,2);
                    evolution.R56(point_index) = R_temp(5,6);
                    evolution.R16(point_index) = R_temp(1,6);
                    evolution.element{point_index} = [element_name '_slice' num2str(slice)];
                    
                    [disp_func, num_val, den_val] = calculateDispersionValue(beta_x_initial, alpha_x_initial, gamma_x_initial, ...
                        R_temp(5,1), R_temp(5,2), R_temp(5,6), epsilon_x, h0, sigma55, sigma66);
                    evolution.dispersion_function(point_index) = disp_func;
                    evolution.numerator(point_index) = num_val;
                    evolution.denominator(point_index) = den_val;
                    
                    % 如果在弯铁组内，存储取值
                    if group_info.in_group
                        sqrt_value = sqrt(disp_func);
                        bend_group_index = find(group_info.group_index == 1:length(bend_groups.start_indices));
                        bend_group_values{bend_group_index}(end+1) = sqrt_value;
                    end
                    
                    point_index = point_index + 1;
                end
            end
            
            % 记录元件终点
            evolution.s(point_index) = current_s + element_length;
            evolution.energy(point_index) = slice_energy;
            evolution.beta_x(point_index) = slice_beta_x;
            evolution.alpha_x(point_index) = slice_alpha_x;
            
            R_total = R_element * R_total;
            evolution.R51(point_index) = R_total(5,1);
            evolution.R52(point_index) = R_total(5,2);
            evolution.R56(point_index) = R_total(5,6);
            evolution.R16(point_index) = R_total(1,6);
            evolution.element{point_index} = [element_name '_end'];
            
            [disp_func, num_val, den_val] = calculateDispersionValue(beta_x_initial, alpha_x_initial, gamma_x_initial, ...
                R_total(5,1), R_total(5,2), R_total(5,6), epsilon_x, h0, sigma55, sigma66);
            evolution.dispersion_function(point_index) = disp_func;
            evolution.numerator(point_index) = num_val;
            evolution.denominator(point_index) = den_val;
            
            % 如果在弯铁组内，存储取值
            if group_info.in_group
                sqrt_value = sqrt(disp_func);
                bend_group_index = find(group_info.group_index == 1:length(bend_groups.start_indices));
                bend_group_values{bend_group_index}(end+1) = sqrt_value;
            end
            
            point_index = point_index + 1;
            
            % 更新当前状态
            current_beta_x = slice_beta_x;
            current_alpha_x = slice_alpha_x;
            current_energy = slice_energy;
            current_s = current_s + element_length;
        end
        
        % 检查是否为弯铁组的结束
        if group_info.is_end
            % 计算弯铁组简单平均值
            bend_group_index = find(group_info.group_index == 1:length(bend_groups.start_indices));
            if ~isempty(bend_group_values{bend_group_index})
                average_value = mean(bend_group_values{bend_group_index});
                
                bend_averages.names{end+1} = group_info.name;
                bend_averages.average_values(end+1) = average_value;
                bend_averages.start_positions(end+1) = bend_start_s;
                bend_averages.end_positions(end+1) = current_s;
                bend_averages.lengths(end+1) = current_s - bend_start_s;
            end
        end
    end
    
    % 裁剪数组
    evolution.s = evolution.s(1:point_index-1);
    evolution.energy = evolution.energy(1:point_index-1);
    evolution.beta_x = evolution.beta_x(1:point_index-1);
    evolution.alpha_x = evolution.alpha_x(1:point_index-1);
    evolution.R51 = evolution.R51(1:point_index-1);
    evolution.R52 = evolution.R52(1:point_index-1);
    evolution.R56 = evolution.R56(1:point_index-1);
    evolution.R16 = evolution.R16(1:point_index-1);
    evolution.dispersion_function = evolution.dispersion_function(1:point_index-1);
    evolution.numerator = evolution.numerator(1:point_index-1);
    evolution.denominator = evolution.denominator(1:point_index-1);
    evolution.element = evolution.element(1:point_index-1);
    
    % 绘图
    set(0, 'DefaultAxesFontName', 'Cambria');
    set(0, 'DefaultTextFontName', 'Cambria');
    
    figure('Position', [100, 100, 700, 500]);
    
    % 子图1
    subplot(2,1,1);
    yyaxis left
    plot(evolution.s, evolution.beta_x,'-', 'LineWidth', 2.5,'Color', [0.8,0,0]);
    ylabel('\beta_x (m)', 'FontName', 'Cambria', 'FontSize', 15, 'Color',  [0.8,0,0]);
    ylim([0, max(evolution.beta_x) * 1.1]);
    set(gca, 'YColor',  [0.8,0,0]);
    
    yyaxis right
    plot(evolution.s, evolution.R16, '-', 'LineWidth', 2.5,'Color',[0.1,0,0.9]);
    ylabel('\eta_x (m)', 'FontName', 'Cambria', 'FontSize', 14, 'Color', [0.1,0,0.9]);
    set(gca, 'YColor',  [0.1,0,0.9]);
    
    xlabel('s (m)', 'FontName', 'Cambria', 'FontSize', 14);
    grid off;
    set(gca, 'LineWidth', 2, 'FontName', 'Cambria', 'FontSize', 14);
    
    % 子图2
    subplot(2,1,2);
    hold on;
    y_data = 100*sqrt(evolution.dispersion_function);
    y_lim = [min(y_data) * 0.9, max(y_data) * 1.5];
    
    % 绘制弯铁区域背景
    for i = 1:length(bend_averages.names)
        fill([bend_averages.start_positions(i), bend_averages.end_positions(i), ...
              bend_averages.end_positions(i), bend_averages.start_positions(i)], ...
             [y_lim(1), y_lim(1), y_lim(2), y_lim(2)], ...
             [0.7, 0.7, 0.7], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
    end
    
    % 绘制主曲线
    plot(evolution.s, y_data, 'k-', 'LineWidth', 2.5);
    
    % 绘制平均值线和标注
    for i = 1:length(bend_averages.names)
        plot([bend_averages.start_positions(i), bend_averages.end_positions(i)], ...
             [1e2*bend_averages.average_values(i), 1e2*bend_averages.average_values(i)], ...
             'k-', 'LineWidth', 2);
        
        text_x = (bend_averages.start_positions(i) + bend_averages.end_positions(i)) / 2;
        text_y = y_lim(2) - (y_lim(2) - y_lim(1)) * 0.15;
        text(text_x, text_y, sprintf('%.1e', 100*bend_averages.average_values(i)), ...
             'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', 'k', ...
              'FontName', 'Cambria');
    end
    
    hold off;
    ylim(y_lim);
    xlim([-2,max(evolution.s)+2])
    xlabel('s (m)', 'FontName', 'Cambria', 'FontSize', 14);
    ylabel('\sigma_e[%]', 'FontName', 'Cambria', 'FontSize', 14);
    grid off;
    box on;
    set(gca, 'LineWidth', 2, 'FontName', 'Cambria', 'FontSize', 14);
    
else
    % 简化的非绘图版本
    R_total = eye(6);
    current_s = 0;
    
    bend_groups = identifyBendGroups(beamline);
    
    % 为每个弯铁组存储所有取值
    bend_group_values = cell(length(bend_groups.start_indices), 1);
    for k = 1:length(bend_groups.start_indices)
        bend_group_values{k} = [];
    end
    
    for i = 1:length(beamline)
        element = beamline(i);
        
        if ~isfield(element, 'type') || strcmp(element.type, 'cavity')
            continue;
        end
        
        if ~isfield(element, 'length')
            element.length = 0;
        end
        
        if ~isfield(element, 'name')
            element.name = sprintf('Unnamed_%s_%d', element.type, i);
        end
        
        group_info = findBendGroupInfo(i, bend_groups);
        
        if group_info.is_start
            bend_start_s = current_s;
        end
        
        if element.length <= 0
            R_element = getElementMatrix(element);
            R_total = R_element * R_total;
            
            % 如果在弯铁组内，存储取值
            if group_info.in_group
                [disp_func, ~, ~] = calculateDispersionValue(beta_x_initial, alpha_x_initial, gamma_x_initial, ...
                    R_total(5,1), R_total(5,2), R_total(5,6), epsilon_x, h0, sigma55, sigma66);
                sqrt_value = sqrt(disp_func);
                bend_group_index = find(group_info.group_index == 1:length(bend_groups.start_indices));
                bend_group_values{bend_group_index}(end+1) = sqrt_value;
            end
        else
            num_slices = 10;
            slice_length = element.length / num_slices;
            R_element = eye(6);
            
            for slice = 1:num_slices
                slice_element = element;
                slice_element.length = slice_length;
                R_slice = getElementMatrix(slice_element);
                R_element = R_slice * R_element;
                
                if group_info.in_group
                    R_temp = R_element * R_total;
                    [disp_func, ~, ~] = calculateDispersionValue(beta_x_initial, alpha_x_initial, gamma_x_initial, ...
                        R_temp(5,1), R_temp(5,2), R_temp(5,6), epsilon_x, h0, sigma55, sigma66);
                    sqrt_value = sqrt(disp_func);
                    bend_group_index = find(group_info.group_index == 1:length(bend_groups.start_indices));
                    bend_group_values{bend_group_index}(end+1) = sqrt_value;
                end
            end
            
            R_total = R_element * R_total;
            current_s = current_s + element.length;
        end
        
        if group_info.is_end
            % 计算弯铁组简单平均值
            bend_group_index = find(group_info.group_index == 1:length(bend_groups.start_indices));
            if ~isempty(bend_group_values{bend_group_index})
                average_value = mean(bend_group_values{bend_group_index});
                
                bend_averages.names{end+1} = group_info.name;
                bend_averages.average_values(end+1) = average_value;
                bend_averages.start_positions(end+1) = bend_start_s;
                bend_averages.end_positions(end+1) = current_s;
                bend_averages.lengths(end+1) = current_s - bend_start_s;
            end
        end
    end
    
    evolution = struct();
    evolution.R51 = R_total(5,1);
    evolution.R52 = R_total(5,2);
    evolution.R56 = R_total(5,6);
end

% 计算最终值
[final_disp_func, ~, ~] = calculateDispersionValue(beta_x_initial, alpha_x_initial, gamma_x_initial, ...
    R_total(5,1), R_total(5,2), R_total(5,6), epsilon_x, h0, sigma55,sigma66);
final_value = sqrt(final_disp_func);

end

function bend_groups = identifyBendGroups(beamline)
% 识别弯铁组
bend_groups = struct();
bend_groups.start_indices = [];
bend_groups.end_indices = [];
bend_groups.names = {};

current_group_start = [];
group_names = {};

for i = 1:length(beamline)
    element = beamline(i);
    
    if ~isfield(element, 'type')
        continue;
    end
    
    if ~isfield(element, 'length')
        element.length = 0;
    end
    
    if ~isfield(element, 'name')
        element.name = sprintf('Unnamed_%s_%d', element.type, i);
    end
    
    is_bend = ismember(lower(element.type), {'bend', 'sbend', 'dipole'}) && element.length > 0;
    
    if is_bend
        if isempty(current_group_start)
            % 开始新组
            current_group_start = i;
            group_names = {element.name};
        else
            % 添加到当前组
            group_names{end+1} = element.name;
        end
    else
        if ~isempty(current_group_start)
            if element.length == 0
                % 零长度元件，继续当前组
            else
                % 结束当前组
                bend_groups.start_indices(end+1) = current_group_start;
                bend_groups.end_indices(end+1) = i-1;
                bend_groups.names{end+1} = strjoin(group_names, '+');
                current_group_start = [];
                group_names = {};
            end
        end
    end
end

% 处理最后一个组
if ~isempty(current_group_start)
    bend_groups.start_indices(end+1) = current_group_start;
    bend_groups.end_indices(end+1) = length(beamline);
    bend_groups.names{end+1} = strjoin(group_names, '+');
end

end

function group_info = findBendGroupInfo(element_index, bend_groups)
% 查找元素所属的弯铁组信息
group_info = struct();
group_info.in_group = false;
group_info.is_start = false;
group_info.is_end = false;
group_info.name = '';
group_info.group_index = 0;

if isempty(bend_groups.start_indices)
    return;
end

for i = 1:length(bend_groups.start_indices)
    if element_index >= bend_groups.start_indices(i) && element_index <= bend_groups.end_indices(i)
        group_info.in_group = true;
        group_info.is_start = (element_index == bend_groups.start_indices(i));
        group_info.is_end = (element_index == bend_groups.end_indices(i));
        group_info.name = bend_groups.names{i};
        group_info.group_index = i;
        break;
    end
end

end

function [dispersion_value, numerator, denominator] = calculateDispersionValue(beta_x, alpha_x, gamma_x, R51, R52, R56, epsilon_x, h0, sigma55, sigma66)
% 计算 s11, s12, s22
s11 = beta_x * epsilon_x;
s12 = -alpha_x * epsilon_x;  
s22 = gamma_x * epsilon_x;
s55 = sigma55;
s66 = sigma66;

% 计算公共项
common_term = R51^2 * s11 + 2 * R51 * R52 * s12 + R52^2 * s22;

% 计算分子
numerator = h0^2 * common_term * s55 + (common_term + s55) * s66;

% 计算分母  
denominator = common_term + (1 + h0 * R56)^2 * s55 + R56^2 * s66;

if abs(denominator) < 1e-15
    dispersion_value = 0;
    warning('分母接近零，色散函数值设为0');
else
    dispersion_value = numerator / denominator;
end
end