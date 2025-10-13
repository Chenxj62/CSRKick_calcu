function bmadCode = writeBmad(beamline)
% writeBmad - Convert beamline elements to Bmad lattice format
% 处理重复元素，为重复出现的元素添加序号后缀
%
% Input parameters:
%   beamline - Array of element structures with fields:
%     .type - Element type (e.g., 'quadrupole', 'dipole', 'drift', 'sextupole')
%     .length - Element length [m]
%     .k1 - Quadrupole strength [1/m^2] (for quadrupoles)
%     .angle - Bend angle [rad] (for dipoles)
%     .k2 - Sextupole strength [1/m^3] (for sextupoles)
%     .name - Element name
%     .h - Curvature [1/m] (for dipoles)
%     .frequency - RF frequency [Hz] (for cavities)
%     .phase - RF phase [rad] (for cavities)
%
% Output:
%   bmadCode - String containing the Bmad lattice definition
%   (结果也会直接输出到命令窗口)

% Initialize output string

clc;
bmadCode = '';

% 创建一个用于跟踪元素名称出现次数的映射
elementNameCounts = containers.Map('KeyType', 'char', 'ValueType', 'int32');
% 创建一个用于存储修改后元素名称的数组
modifiedNames = cell(size(beamline));

% 第一遍：为每个元素分配一个唯一的名称
for i = 1:length(beamline)
    element = beamline(i);
    
    % 跳过没有类型的元素
    if ~isfield(element, 'type')
        warning('Element #%d has no type field, will be skipped', i);
        modifiedNames{i} = '';
        continue;
    end
    
    % 获取元素名称
    if ~isfield(element, 'name') || isempty(element.name)
        elementType = lower(element.type);
        element.name = sprintf('%s%d', elementType(1:min(3,length(elementType))), i);
    end
    
    baseName = element.name;
    
    % 检查这个名称之前是否出现过
    if isKey(elementNameCounts, baseName)
        % 更新计数并创建新名称
        count = elementNameCounts(baseName) + 1;
        elementNameCounts(baseName) = count;
        uniqueName = sprintf('%s_%d', baseName, count);
    else
        % 首次出现，记录并使用原始名称
        elementNameCounts(baseName) = 1;
        uniqueName = baseName;
    end
    
    % 存储修改后的名称
    modifiedNames{i} = uniqueName;
end

% 辅助函数：检查指定位置的元素是否为edge
function [isEdge, edgeAngle] = checkEdgeElement(beamline, index)
    isEdge = false;
    edgeAngle = 0;
    
    if index >= 1 && index <= length(beamline)
        element = beamline(index);
        if isfield(element, 'type') && strcmpi(element.type, 'edge')
            isEdge = true;
            if isfield(element, 'angle')
                edgeAngle = element.angle;
            end
        end
    end
end

% 第二遍：生成每个唯一元素的定义
for i = 1:length(beamline)
    element = beamline(i);
    
    % 跳过没有类型或已处理过的元素
    if ~isfield(element, 'type') || isempty(modifiedNames{i})
        continue;
    end
    
    elementType = lower(element.type);
    uniqueName = modifiedNames{i};
    
    % 跳过edge元素，它们会被合并到相邻的dipole中
    if strcmpi(elementType, 'edge')
        continue;
    end
    
    % 检查必要的长度字段
    if ~isfield(element, 'length') && ~strcmp(elementType, 'marker')
        warning('Element %s has no length field, will be set to 0', uniqueName);
        element.length = 0;
    end
    
    % 根据元素类型格式化定义
    switch elementType
        case {'quadrupole', 'quad'}
            % 检查k1是否提供
            if ~isfield(element, 'k1')
                warning('Quadrupole %s has no k1 field, will be set to 0', uniqueName);
                element.k1 = 0;
            end
            
            bmadLine = sprintf('%s: QUAD, L=%g, K1=%g', ...
                uniqueName, element.length, element.k1);
            
        case {'dipole', 'sbend', 'rbend', 'bend'}
            % 标准化为sbend用于Bmad
            bmadType = 'SBEND';
            dipoleLength = element.length;

                % 从曲率计算角度，使用dipole本身的长度
             angle = dipoleLength * element.h;

            % 检查前后是否有edge元素
            [hasE1, E1] = checkEdgeElement(beamline, i-1);
            [hasE2, E2] = checkEdgeElement(beamline, i+1);
            
            % 如果没有edge元素，设为0
            if ~hasE1
                E1 = 0;
            end
            if ~hasE2
                E2 = 0;
            end
            
            % 保存dipole本身的长度，不要被edge影响
            
            
            % 获取弯转角度，使用dipole本身的参数
       
            
            % 优先使用angle字段

            
            % 生成Bmad代码，使用dipole本身的长度
            if abs(angle) > 1e-10
                % 有弯转角度的情况
                if E1 ~= 0 || E2 ~= 0
                    bmadLine = sprintf('%s: %s, L=%g, ANGLE=%g, E1=%g, E2=%g', ...
                        uniqueName, bmadType, dipoleLength, angle, E1, E2);
                else
                    bmadLine = sprintf('%s: %s, L=%g, ANGLE=%g', ...
                        uniqueName, bmadType, dipoleLength, angle);
                end
            else
                % 没有弯转的情况
                bmadLine = sprintf('%s: %s, L=%g, ANGLE=0', ...
                    uniqueName, bmadType, dipoleLength);
            end
            
        case {'drift', 'drft'}
            bmadLine = sprintf('%s: DRIFT, L=%g', ...
                uniqueName, element.length);
            
        case {'sextupole', 'sext'}
            % 检查k2是否提供
            if ~isfield(element, 'k2')
                warning('Sextupole %s has no k2 field, will be set to 0', uniqueName);
                element.k2 = 0;
            end
            
            bmadLine = sprintf('%s: SEXT, L=%g, K2=%g', ...
                uniqueName, element.length, element.k2);
            
        case 'marker'
            bmadLine = sprintf('%s: MARKER', uniqueName);
            
        case {'cavity', 'rfcavity', 'rf'}
            % RF腔参数
            voltage = 0;
            freq = 0;
            phase = 0;
            
            if isfield(element, 'voltage')
                voltage = element.voltage;
            end
            if isfield(element, 'frequency')
                freq = element.frequency;
            end
            if isfield(element, 'phase')
                phase = element.phase;
            end
            
            bmadLine = sprintf('%s: RFCAVITY, L=%g, VOLTAGE=%g, FREQ=%g, PHASE=%g', ...
                uniqueName, element.length, voltage, freq, phase);
            
        otherwise
            % 通用元素，带长度
            bmadLine = sprintf('%s: DRIFT, L=%g  ! Originally %s type', ...
                uniqueName, element.length, elementType);
    end
    
    % 将元素添加到bmad代码，后面跟一个空行
    bmadCode = [bmadCode, bmadLine, sprintf('\n\n')];
    
    % 直接输出当前元素定义
    fprintf('%s\n\n', bmadLine);
end

% 收集所有有效的元素名称，保持束线顺序（排除edge元素）
validElementNames = {};
for i = 1:length(beamline)
    if ~isempty(modifiedNames{i}) && ~strcmpi(beamline(i).type, 'edge')
        validElementNames{end+1} = modifiedNames{i};
    end
end

% 如果有有效元素，添加完整束线定义
if ~isempty(validElementNames)
    % 手动构建元件名称列表，带有换行以提高可读性
    lineElements = '';
    charsPerLine = 80;
    currentLineLength = 0;
    
    for i = 1:length(validElementNames)
        if i > 1
            if currentLineLength + length(validElementNames{i}) + 2 > charsPerLine
                % 添加换行以提高可读性
                lineElements = [lineElements, ',&' sprintf('\n    ')];
                currentLineLength = 4; % 换行后的4个空格
            else
                lineElements = [lineElements, ', '];
                currentLineLength = currentLineLength + 2;
            end
        end
        
        lineElements = [lineElements, validElementNames{i}];
        currentLineLength = currentLineLength + length(validElementNames{i});
    end
    
    beamlineDef = sprintf('full_line: LINE = (%s)', lineElements);
    bmadCode = [bmadCode, beamlineDef, sprintf('\n\n')];
    
    % 直接输出完整束线定义
    fprintf('%s\n\n', beamlineDef);
    
    % 添加USE语句
    bmadCode = [bmadCode, 'USE, full_line'];
    fprintf('USE, full_line\n');
end

% 为避免显示ans='...'，使函数不返回任何值
if nargout == 0
    clear bmadCode
end
end