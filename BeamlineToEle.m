function BeamlineToEle(beamline)
% convertBeamlineToElements - 将beamline转换为元素定义并在命令行显示
% 
% Input:
%   beamline - 从readBmad得到的beamline结构数组
%
% 该函数会生成类似你原来代码的元素定义，并在命令行显示

% 计算各类元素数量
clc;
num_quads = sum(strcmp({beamline.type}, 'quadrupole'));
num_drifts = sum(strcmp({beamline.type}, 'drift'));
num_sexts = sum(strcmp({beamline.type}, 'sextupole'));
num_dipoles = sum(strcmp({beamline.type}, 'dipole') & [beamline.length] > 0);

% 重新分配索引，避免重叠
quad_start = 1;
drift_start = num_quads + 1;
sext_start = num_quads + num_drifts + 1;
rho_start = num_quads + num_drifts + num_sexts + 1;

total_vars = rho_start + num_dipoles - 1;

fprintf('adjust_vars = zeros(1, %d);\n\n', total_vars);

fprintf('%% 从beamline转换的元素定义\n');
fprintf('%% 使用adjust_vars变量进行参数调整（无重叠索引）\n\n');

% 分析beamline，收集基础参数
drift_lengths = [];
quad_k1s = [];
quad_lengths = [];
sext_k2s = [];
sext_lengths = [];
dipole_angles = [];
dipole_rhos = [];

% 检查每个dipole是否有对应的edge元素
dipole_has_e1 = [];
dipole_has_e2 = [];
dipole_e1_angles = [];
dipole_e2_angles = [];

% 元素计数器
drift_idx = 1;
quad_idx = 1;
sext_idx = 1;
dipole_idx = 1;

% 第一遍：收集所有基础参数
for i = 1:length(beamline)
    element = beamline(i);
    
    switch lower(element.type)
        case 'drift'
            drift_lengths(end+1) = element.length;
            
        case 'quadrupole'
            quad_k1s(end+1) = element.k1;
            quad_lengths(end+1) = element.length;
            
        case 'sextupole'
            sext_k2s(end+1) = element.k2;
            sext_lengths(end+1) = element.length;
            
        case 'dipole'
            if element.length > 0  % 只记录非thin dipole
                dipole_angles(end+1) = element.angle;
                if element.h ~= 0
                    dipole_rhos(end+1) = 1/element.h;  % rho = 1/h
                else
                    dipole_rhos(end+1) = element.length/element.angle;
                end
                
                % 检查前后是否有edge元素（基于位置判断）
                has_e1 = false;
                e1_angle = 0;
                if i > 1 && strcmp(beamline(i-1).type, 'edge')
                    has_e1 = true;
                    e1_angle = beamline(i-1).angle;
                end
                
                has_e2 = false;
                e2_angle = 0;
                if i < length(beamline) && strcmp(beamline(i+1).type, 'edge')
                    has_e2 = true;
                    e2_angle = beamline(i+1).angle;
                end
                
                dipole_has_e1(end+1) = has_e1;
                dipole_has_e2(end+1) = has_e2;
                dipole_e1_angles(end+1) = e1_angle;
                dipole_e2_angles(end+1) = e2_angle;
                
                dipole_idx = dipole_idx + 1;
            end
    end
end

% 显示基础参数定义
fprintf('%% 基础参数\n');
if ~isempty(dipole_angles)
    for i = 1:length(dipole_angles)
        fprintf('B%d_angle = %.6f;  %% rad\n', i, dipole_angles(i));
        fprintf('rho%d_base = %.5f;  %% B%d基础弯转半径\n', i, dipole_rhos(i), i);
    end
end

% 显示实际长度
if ~isempty(quad_lengths)
    fprintf('\n%% 四极铁长度\n');
    for i = 1:length(quad_lengths)
        fprintf('Lq%d = %.6f;  %% Q%d长度\n', i, quad_lengths(i), i);
    end
end

if ~isempty(sext_lengths)
    fprintf('\n%% 六极铁长度\n');
    for i = 1:length(sext_lengths)
        fprintf('Ls%d = %.6f;  %% S%d长度\n', i, sext_lengths(i), i);
    end
end

fprintf('\ncons = [];  %% 约束初始化\n\n');

% 重置计数器
drift_idx = 1;
quad_idx = 1;
sext_idx = 1;
dipole_idx = 1;

fprintf('%% 使用基础值+adjust_vars创建元件（无重叠索引）\n');

% 第二遍：生成createElement调用
for i = 1:length(beamline)
    element = beamline(i);
    
    % 跳过edge元素，它们会在dipole处理中生成
    if strcmp(element.type, 'edge')
        continue;
    end
    
    switch lower(element.type)
        case 'drift'
            adjust_idx = drift_start + drift_idx - 1;
            fprintf('ELE.D%d = createElement(''drift'', %.6f + adjust_vars(%d), 0, ''D%d'');\n', ...
                   drift_idx, drift_lengths(drift_idx), adjust_idx, drift_idx);
            drift_idx = drift_idx + 1;
            
        case 'quadrupole'
            adjust_idx = quad_start + quad_idx - 1;
            fprintf('ELE.Q%d = createElement(''quadrupole'', Lq%d, %.5f + adjust_vars(%d), ''Q%d'');\n', ...
                   quad_idx, quad_idx, quad_k1s(quad_idx), adjust_idx, quad_idx);
            quad_idx = quad_idx + 1;
            
        case 'sextupole'
            adjust_idx = sext_start + sext_idx - 1;
            fprintf('ELE.S%d = createElement(''sextupole'', Ls%d, %.5f + adjust_vars(%d), ''S%d'');\n', ...
                   sext_idx, sext_idx, sext_k2s(sext_idx), adjust_idx, sext_idx);
            sext_idx = sext_idx + 1;
            
        case 'dipole'
            if element.length > 0
                rho_adjust_idx = rho_start + dipole_idx - 1;
                fprintf('rho%d = %.5f + adjust_vars(%d);  %% B%d调整后的弯转半径\n', ...
                       dipole_idx, dipole_rhos(dipole_idx), rho_adjust_idx, dipole_idx);
                
                % 如果有E1，创建E1 edge元素作为变量
                if dipole_has_e1(dipole_idx)
                    fprintf('E1_%d = createElement(''edge'', 0, 1/rho%d, %.6f, ''E1_%d'');\n', ...
                           dipole_idx, dipole_idx, dipole_e1_angles(dipole_idx), dipole_idx);
                end
                
                % 创建dipole元素
                fprintf('ELE.B%d = createElement(''dipole'', B%d_angle*rho%d, 1/rho%d, ''B%d'');\n', ...
                       dipole_idx, dipole_idx, dipole_idx, dipole_idx, dipole_idx);
                
                % 如果有E2，创建E2 edge元素作为变量
                if dipole_has_e2(dipole_idx)
                    fprintf('E2_%d = createElement(''edge'', 0, 1/rho%d, %.6f, ''E2_%d'');\n', ...
                           dipole_idx, dipole_idx, dipole_e2_angles(dipole_idx), dipole_idx);
                end
                
                dipole_idx = dipole_idx + 1;
            end
    end
end

fprintf('\n%% 定义完整束线\n');
fprintf('beamline_elements = [');

% 生成beamline数组
beamline_elements = {};
drift_counter = 0;
quad_counter = 0;
sext_counter = 0;
dipole_counter = 0;

for i = 1:length(beamline)
    element = beamline(i);
    
    % 跳过edge元素，它们已经包含在dipole处理中
    if strcmpi(element.type, 'edge')
        continue;
    end
    
    % 根据元素类型添加到beamline
    switch lower(element.type)
        case 'drift'
            drift_counter = drift_counter + 1;
            beamline_elements{end+1} = sprintf('ELE.D%d', drift_counter);
            
        case 'quadrupole'
            quad_counter = quad_counter + 1;
            beamline_elements{end+1} = sprintf('ELE.Q%d', quad_counter);
            
        case 'sextupole'
            sext_counter = sext_counter + 1;
            beamline_elements{end+1} = sprintf('ELE.S%d', sext_counter);
            
        case 'dipole'
            if element.length > 0
                dipole_counter = dipole_counter + 1;
                
                % 添加E1（如果存在）
                if dipole_has_e1(dipole_counter)
                    beamline_elements{end+1} = sprintf('E1_%d', dipole_counter);
                end
                
                % 添加dipole本身
                beamline_elements{end+1} = sprintf('ELE.B%d', dipole_counter);
                
                % 添加E2（如果存在）
                if dipole_has_e2(dipole_counter)
                    beamline_elements{end+1} = sprintf('E2_%d', dipole_counter);
                end
            end
    end
end

% 输出beamline定义
for i = 1:length(beamline_elements)
    if i < length(beamline_elements)
        fprintf('%s, ', beamline_elements{i});
        if mod(i, 5) == 0  % 每5个元素换行
            fprintf('...\n                     ');
        end
    else
        fprintf('%s];\n\n', beamline_elements{i});
    end
end

% 显示变量分配说明
fprintf('%% adjust_vars变量的分配说明（无重叠索引）：\n');
fprintf('%% adjust_vars(%d-%d): 四极铁K1强度的调整量\n', quad_start, quad_start+num_quads-1);
fprintf('%% adjust_vars(%d-%d): 漂移段长度的调整量\n', drift_start, drift_start+num_drifts-1);
if num_sexts > 0
    fprintf('%% adjust_vars(%d-%d): 六极铁K2强度的调整量\n', sext_start, sext_start+num_sexts-1);
end
fprintf('%% adjust_vars(%d-%d): 弯转半径的调整量 (每个dipole独立)\n', rho_start, rho_start+num_dipoles-1);

% 生成约束定义
fprintf('\n%% 变量约束定义\n');
fprintf('D_lengths = [');
for i = 1:length(drift_lengths)
    if i < length(drift_lengths)
        fprintf('%.6f, ', drift_lengths(i));
    else
        fprintf('%.6f];\n', drift_lengths(i));
    end
end

% 生成弯转半径约束
if ~isempty(dipole_rhos)
    fprintf('rho_bases = [');
    for i = 1:length(dipole_rhos)
        if i < length(dipole_rhos)
            fprintf('%.5f, ', dipole_rhos(i));
        else
            fprintf('%.5f];\n', dipole_rhos(i));
        end
    end
    
    fprintf('\noptions.lb = [-10*ones(1,%d), -D_lengths, 0*ones(1,%d), -0.5*rho_bases];\n', num_quads, num_sexts);
    fprintf('options.ub = [10*ones(1,%d), D_lengths, 0*ones(1,%d), 2*rho_bases];\n', num_quads, num_sexts);
else
    fprintf('\noptions.lb = [-10*ones(1,%d), -D_lengths, 0*ones(1,%d)];\n', num_quads, num_sexts);
    fprintf('options.ub = [10*ones(1,%d), D_lengths, 0*ones(1,%d)];\n', num_quads, num_sexts);
end

end