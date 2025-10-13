function floorplan = plotFloorPlan(beamline, options)
% plotFloorPlan - Plot beamline floor plan with track deflection
%
% Input parameters:
%   beamline - Array of element structures
%   options - Optional parameter structure:
%     .figure_size - Figure size [width, height] (default [800, 600])
%     .rect_width_factor - Rectangle width factor (default 0.75) - relative to element length
%     .rect_length - Rectangle length (default 0.2) - absolute length in meters
%     .x_range - X-axis range [min, max] (default calculated from data)
%     .y_range - Y-axis range [min, max] (default calculated from data)
%     .axis_line_width - Axis line width (default 1.5)
%     .font_name - Font name (default 'Cambria Math')
%     .font_size - Axis and tick font size (default 10)
%     .title_text - Custom title text (default 'Beamline Floor Plan')
%     .title_font_size - Title font size (default 14)
%     .rotation_angle - Overall rotation angle in degrees (default 0)

% Set default parameters
if nargin < 2
    options = struct();
end

if ~isfield(options, 'figure_size')
    options.figure_size = [600, 400];
end
if ~isfield(options, 'line_width')
    options.line_width = 2;
end
if ~isfield(options, 'axis_line_width')
    options.axis_line_width = 1.5;
end
if ~isfield(options, 'font_name')
    options.font_name = 'Cambria Math';
end
if ~isfield(options, 'font_size')
    options.font_size = 15;
end
if ~isfield(options, 'title_text')
    options.title_text = 'Beamline Floor Plan';  % 默认标题文字
end
if ~isfield(options, 'title_font_size')
    options.title_font_size = 18;
end
if ~isfield(options, 'rect_width_factor')
    options.rect_width_factor = 1;  % 元件宽度为长度的0.75倍（沿束线方向）
end
if ~isfield(options, 'rect_length')
    options.rect_length = .6;  % 元件长度为绝对值0.2米（垂直于束线方向）
end
if ~isfield(options, 'focus_elements')
    options.focus_elements = {'quadrupole', 'dipole', 'sextupole', 'cavity', 'TGU'};
end
if ~isfield(options, 'ang')
    options.ang = 0;  % 默认不旋转
end

% 将旋转角度转换为弧度
rotation_rad = options.ang * pi / 180;



% Initialize beamline layout data
max_points = length(beamline) * 20;
floorplan = struct();
floorplan.x = zeros(max_points, 1);
floorplan.y = zeros(max_points, 1);
floorplan.s = zeros(max_points, 1);
floorplan.theta = zeros(max_points, 1);
floorplan.elements = struct('name', {}, 'type', {}, 'pos_start', [], 'pos_end', [], 'x_center', [], 'y_center', [], 's_center', [], 'theta_center', []);

% Initialize current position and direction
current_x = 0;
current_y = 0;
current_s = 0;
current_theta = 0;
point_index = 1;

% Record initial point
floorplan.x(point_index) = current_x;
floorplan.y(point_index) = current_y;
floorplan.s(point_index) = current_s;
floorplan.theta(point_index) = current_theta;
point_index = point_index + 1;

% Track element types and positions
element_index = 1;

% Process each element
for i = 1:length(beamline)
    element = beamline(i);
    
    % Ensure element has required fields
    if ~isfield(element, 'type')
        warning('Element #%d has no type field, will be skipped', i);
        continue;
    end
    
    if ~isfield(element, 'length')
        warning('Element %s has no length field, will be set to 0', element.type);
        element.length = 0;
    end
    
    if ~isfield(element, 'name')
        element.name = sprintf('Unnamed_%s_%d', element.type, i);
    end
    
    element_type = element.type;
    element_length = element.length;
    element_name = element.name;
    
    % Record element start position
    element_start_x = current_x;
    element_start_y = current_y;
    element_start_s = current_s;
    element_start_theta = current_theta;
    
    % Process based on element type
    switch element_type
        case 'dipole'
            % Handle dipole magnets
            angle = 0;
            
            % Get angle from h or angle field
            if isfield(element, 'angle')
                angle = element.angle;
            elseif isfield(element, 'h') && element.h ~= 0
                % Calculate angle from curvature
                rho = 1 / element.h;
                angle = element_length / rho;
            else
                warning('Dipole %s has no angle information, will be treated as straight section', element_name);
            end
            
            if abs(angle) < 1e-10 || element_length <= 0
                % If angle is close to zero or length is zero, process as straight section
                new_x = current_x + element_length * cos(current_theta);
                new_y = current_y + element_length * sin(current_theta);
                
                % Record element endpoint
                floorplan.x(point_index) = new_x;
                floorplan.y(point_index) = new_y;
                floorplan.s(point_index) = current_s + element_length;
                floorplan.theta(point_index) = current_theta;
                point_index = point_index + 1;
            else
                % Calculate dipole parameters
                radius = element_length / angle;
                
                % Calculate circle center
                circle_center_x = current_x - radius * sin(current_theta);
                circle_center_y = current_y + radius * cos(current_theta);
                
                % Calculate start and end angles
                start_angle = current_theta - pi/2;
                end_angle = start_angle + angle;
                
                % Sample multiple points on the dipole to draw the arc
                num_samples = 20;
                angles = linspace(start_angle, end_angle, num_samples);
                
                for j = 1:num_samples
                    % Calculate current sample point position
                    sample_x = circle_center_x + radius * cos(angles(j));
                    sample_y = circle_center_y + radius * sin(angles(j));
                    sample_s = current_s + (j-1) * element_length / (num_samples-1);
                    sample_theta = angles(j) + pi/2;
                    
                    % Record sample point
                    floorplan.x(point_index) = sample_x;
                    floorplan.y(point_index) = sample_y;
                    floorplan.s(point_index) = sample_s;
                    floorplan.theta(point_index) = sample_theta;
                    point_index = point_index + 1;
                end
                
                % Update current direction to dipole exit direction
                current_theta = current_theta + angle;
            end
            
        otherwise
            % Handle all other elements (as straight sections)
            new_x = current_x + element_length * cos(current_theta);
            new_y = current_y + element_length * sin(current_theta);
            
            % Record element endpoint
            floorplan.x(point_index) = new_x;
            floorplan.y(point_index) = new_y;
            floorplan.s(point_index) = current_s + element_length;
            floorplan.theta(point_index) = current_theta;
            point_index = point_index + 1;
    end
    
    % Update current position
    current_x = floorplan.x(point_index-1);
    current_y = floorplan.y(point_index-1);
    current_s = current_s + element_length;
    
    % Record element information (for plotting)
    if ismember(element_type, options.focus_elements)
        % Calculate element center position
        element_end_x = current_x;
        element_end_y = current_y;
        element_end_s = current_s;
        element_end_theta = current_theta;
        
        % For straight elements, center position is midpoint of start and end
        element_center_x = (element_start_x + element_end_x) / 2;
        element_center_y = (element_start_y + element_end_y) / 2;
        element_center_s = (element_start_s + element_end_s) / 2;
        element_center_theta = element_start_theta;
        
        % If dipole, calculate more accurate center position
        if strcmp(element_type, 'dipole')
            angle = 0;
            if isfield(element, 'angle')
                angle = element.angle;
            elseif isfield(element, 'h') && element.h ~= 0
                rho = 1 / element.h;
                angle = element_length / rho;
            end
            
            if abs(angle) > 1e-10
                % For dipoles, center position is at the midpoint of the arc
                radius = element_length / angle;
                circle_center_x = element_start_x - radius * sin(element_start_theta);
                circle_center_y = element_start_y + radius * cos(element_start_theta);
                
                % Calculate midpoint angle
                mid_angle = element_start_theta - pi/2 + angle/2;
                
                % Calculate midpoint coordinates
                element_center_x = circle_center_x + radius * cos(mid_angle);
                element_center_y = circle_center_y + radius * sin(mid_angle);
                
                % Dipole center point tangent direction
                element_center_theta = mid_angle + pi/2;
            end
        end
        
        floorplan.elements(element_index).name = element_name;
        floorplan.elements(element_index).type = element_type;
        floorplan.elements(element_index).pos_start = [element_start_x, element_start_y];
        floorplan.elements(element_index).pos_end = [element_end_x, element_end_y];
        floorplan.elements(element_index).x_center = element_center_x;
        floorplan.elements(element_index).y_center = element_center_y;
        floorplan.elements(element_index).s_center = element_center_s;
        floorplan.elements(element_index).theta_center = element_center_theta;
        floorplan.elements(element_index).length = element_length;
        
        % 为cavity和TGU添加额外的信息
        if strcmp(element_type, 'cavity')
            if isfield(element, 'num_cells')
                floorplan.elements(element_index).num_cells = element.num_cells;
            else
                floorplan.elements(element_index).num_cells = 1; % 默认值
            end
        end
        
        element_index = element_index + 1;
    end
end

% Trim arrays to actual size
floorplan.x = floorplan.x(1:point_index-1);
floorplan.y = floorplan.y(1:point_index-1);
floorplan.s = floorplan.s(1:point_index-1);
floorplan.theta = floorplan.theta(1:point_index-1);

% 应用整体旋转
if abs(rotation_rad) > 1e-10
    % 旋转束线轨迹
    x_rot = floorplan.x * cos(rotation_rad) - floorplan.y * sin(rotation_rad);
    y_rot = floorplan.x * sin(rotation_rad) + floorplan.y * cos(rotation_rad);
    floorplan.x = x_rot;
    floorplan.y = y_rot;
    
    % 旋转束线角度
    floorplan.theta = floorplan.theta + rotation_rad;
    
    % 旋转元件位置
    for i = 1:length(floorplan.elements)
        % 旋转起始位置
        x_start_rot = floorplan.elements(i).pos_start(1) * cos(rotation_rad) - floorplan.elements(i).pos_start(2) * sin(rotation_rad);
        y_start_rot = floorplan.elements(i).pos_start(1) * sin(rotation_rad) + floorplan.elements(i).pos_start(2) * cos(rotation_rad);
        floorplan.elements(i).pos_start = [x_start_rot, y_start_rot];
        
        % 旋转结束位置
        x_end_rot = floorplan.elements(i).pos_end(1) * cos(rotation_rad) - floorplan.elements(i).pos_end(2) * sin(rotation_rad);
        y_end_rot = floorplan.elements(i).pos_end(1) * sin(rotation_rad) + floorplan.elements(i).pos_end(2) * cos(rotation_rad);
        floorplan.elements(i).pos_end = [x_end_rot, y_end_rot];
        
        % 旋转中心位置
        x_center_rot = floorplan.elements(i).x_center * cos(rotation_rad) - floorplan.elements(i).y_center * sin(rotation_rad);
        y_center_rot = floorplan.elements(i).x_center * sin(rotation_rad) + floorplan.elements(i).y_center * cos(rotation_rad);
        floorplan.elements(i).x_center = x_center_rot;
        floorplan.elements(i).y_center = y_center_rot;
        
        % 旋转元件角度
        floorplan.elements(i).theta_center = floorplan.elements(i).theta_center + rotation_rad;
    end
end

% Calculate layout boundaries (after rotation)
floorplan.x_min = min(floorplan.x);
floorplan.x_max = max(floorplan.x);
floorplan.y_min = min(floorplan.y);
floorplan.y_max = max(floorplan.y);
floorplan.total_length = max(floorplan.s);

% Plot beamline layout
fig = figure('Position', [100, 100, options.figure_size]);
set(fig, 'Color', 'white');

% Plot beamline track
plot(floorplan.x, floorplan.y, 'k-', 'LineWidth', options.line_width);
hold on;

% Prepare for legend
legend_handles = [];
legend_texts = {};
element_types_shown = {};

% Draw special element markers
for i = 1:length(floorplan.elements)
    element = floorplan.elements(i);
    x_center = element.x_center;
    y_center = element.y_center;
    element_type = element.type;
    element_length = element.length;
    element_theta = element.theta_center;
    
    % Choose different colors and shapes based on element type
    switch element_type
        case 'quadrupole'
            fill_color = [.6, 0, 0];
            element_display_name = 'Quadrupole';
            
            % Calculate rectangle parameters
            rect_width = element_length * options.rect_width_factor;
            rect_length = options.rect_length;
            rect_width = max(rect_width, 0.05);
            
            % Draw rectangle
            h = drawRectangle(x_center, y_center, rect_width, rect_length, element_theta, fill_color);
            
        case 'dipole'
            fill_color = [0.1, 0, 0.8];
            element_display_name = 'Dipole';
            
            % Calculate rectangle parameters
            rect_width = element_length * options.rect_width_factor*0.5;
            rect_length = options.rect_length;
            rect_width = max(rect_width, 0.05);
            
            % Draw rectangle
            h = drawRectangle(x_center, y_center, rect_width, rect_length, element_theta, fill_color);
            
        case 'sextupole'
            fill_color = [0.2, 0, .5];
            element_display_name = 'Sextupole';
            
            % Calculate rectangle parameters
            rect_width = element_length * options.rect_width_factor;
            rect_length = options.rect_length;
            rect_width = max(rect_width, 0.05);
            
            % Draw rectangle
            h = drawRectangle(x_center, y_center, rect_width, rect_length, element_theta, fill_color);
            
        case 'cavity'
            fill_color = [0.5, 0, 0.5];  % 紫色
            element_display_name = 'Cavity';
            
            % 获取cell数量
            if isfield(element, 'num_cells')
                num_cells = element.num_cells;
            else
                num_cells = 1;
            end
            
            % 计算每个cell的位置
            cell_length = element_length / num_cells;
            start_pos_x = element.pos_start(1);
            start_pos_y = element.pos_start(2);
            
            % 为图例创建一个组合对象
            h_group = hggroup;
            
            % 绘制每个cell为椭圆
            for j = 1:num_cells
                % 计算当前cell的中心位置
                cell_offset = (j - 0.5) * cell_length;
                cell_center_x = start_pos_x + cell_offset * cos(element_theta);
                cell_center_y = start_pos_y + cell_offset * sin(element_theta);
                
                % 椭圆参数（沿束线方向稍长，垂直方向稍短）
                ellipse_width = cell_length ;  % 沿束线方向
                ellipse_height = options.rect_length * 1.2;  % 垂直束线方向
                
                % 绘制椭圆
                h_cell = drawEllipse(cell_center_x, cell_center_y, ellipse_width, ellipse_height, element_theta, fill_color);
                set(h_cell, 'Parent', h_group);
            end
            
            h = h_group;
            
        case 'TGU'
            element_display_name = 'Undulator';
            
            % TGU用9个长方形拼接，红蓝交替
            num_segments = 9;
            segment_length = element_length / num_segments;
            segment_width = options.rect_length ;  % 垂直束线方向的宽度
            
            start_pos_x = element.pos_start(1);
            start_pos_y = element.pos_start(2);
            
            % 红蓝交替的颜色
            colors = {[0.8, 0, 0], [0, 0,0]};  % 红色和黑色
            
            % 为图例创建一个组合对象
            h_group = hggroup;
            
            for j = 1:num_segments
                % 计算当前段的中心位置
                segment_offset = (j - 0.5) * segment_length;
                segment_center_x = start_pos_x + segment_offset * cos(element_theta);
                segment_center_y = start_pos_y + segment_offset * sin(element_theta);
                
                % 选择颜色（红黑交替）
                color_index = mod(j-1, 2) + 1;
                segment_color = colors{color_index};
                
                % 绘制长方形段
                h_segment = drawRectangle(segment_center_x, segment_center_y, segment_length*0.9, segment_width, element_theta, segment_color);
                set(h_segment, 'Parent', h_group);
            end
            
            h = h_group;
            
        otherwise
            fill_color = [0.5, 0.5, 0.5];
            element_display_name = element_type;
            
            % Calculate rectangle parameters
            rect_width = element_length * options.rect_width_factor;
            rect_length = options.rect_length;
            rect_width = max(rect_width, 0.05);
            
            % Draw rectangle
            h = drawRectangle(x_center, y_center, rect_width, rect_length, element_theta, fill_color);
    end
    
    % Collect unique element types for legend
    if ~ismember(element_type, element_types_shown)
        element_types_shown{end+1} = element_type;
        legend_handles(end+1) = h;
        legend_texts{end+1} = element_display_name;
    end
end

% Add legend
if ~isempty(legend_handles)
    lgd = legend(legend_handles, legend_texts, 'Location', 'northeast');
    % 图例字体大小与坐标轴刻度字体大小一致
    set(lgd, 'FontName', options.font_name, 'FontSize', options.font_size);
end

% Set figure properties - 使用自定义标题
title(options.title_text, 'FontSize', options.title_font_size, 'FontName', options.font_name);
xlabel('X (m)', 'FontSize', options.font_size, 'FontName', options.font_name);
ylabel('Y (m)', 'FontSize', options.font_size, 'FontName', options.font_name);

% Set X and Y axis range
if isfield(options, 'x_range')
    xlim(options.x_range);
else
    % Default: add a small margin to the data range
    xlim([floorplan.x_min - 2.0, floorplan.x_max + 2.0]);
end

if isfield(options, 'y_range')
    ylim(options.y_range);
else
    % Default: add a small margin to the data range
    ylim([floorplan.y_min - 2.0, floorplan.y_max + 2.0]);
end

% Set axis equal scale and appearance
axis equal;
grid off;
box on;

% Set axis line width and font
ax = gca;
set(ax, 'FontName', options.font_name, 'FontSize', options.font_size, 'LineWidth', options.axis_line_width);

end

% 辅助函数：绘制长方形
function h = drawRectangle(x_center, y_center, rect_width, rect_length, theta, fill_color)
    % Calculate rectangle vertices
    perp_dir = theta + pi/2;
    
    p1_x = x_center + (rect_length/2) * cos(perp_dir) + (rect_width/2) * cos(theta);
    p1_y = y_center + (rect_length/2) * sin(perp_dir) + (rect_width/2) * sin(theta);
    
    p2_x = x_center + (rect_length/2) * cos(perp_dir) - (rect_width/2) * cos(theta);
    p2_y = y_center + (rect_length/2) * sin(perp_dir) - (rect_width/2) * sin(theta);
    
    p3_x = x_center - (rect_length/2) * cos(perp_dir) - (rect_width/2) * cos(theta);
    p3_y = y_center - (rect_length/2) * sin(perp_dir) - (rect_width/2) * sin(theta);
    
    p4_x = x_center - (rect_length/2) * cos(perp_dir) + (rect_width/2) * cos(theta);
    p4_y = y_center - (rect_length/2) * sin(perp_dir) + (rect_width/2) * sin(theta);
    
    % Draw rectangle
    h = patch([p1_x, p2_x, p3_x, p4_x], [p1_y, p2_y, p3_y, p4_y], fill_color, 'EdgeColor', 'none');
end

% 辅助函数：绘制椭圆
function h = drawEllipse(x_center, y_center, width, height, theta, fill_color)
    % 创建椭圆参数方程的参数
    t = linspace(0, 2*pi, 50);
    
    % 椭圆在局部坐标系中的坐标
    x_local = (width/2) * cos(t);
    y_local = (height/2) * sin(t);
    
    % 旋转到全局坐标系
    x_global = x_center + x_local * cos(theta) - y_local * sin(theta);
    y_global = y_center + x_local * sin(theta) + y_local * cos(theta);
    
    % 绘制椭圆
    h = patch(x_global, y_global, fill_color, 'EdgeColor', 'none');
end