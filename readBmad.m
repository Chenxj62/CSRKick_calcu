function [beamline] = readBmad(bmadFile)
% readBmad - Convert Bmad lattice format to beamline elements
% 从Bmad格式文件转换为beamline结构数组
%
% Input parameters:
%   bmadFile - Filename of .bmad file
%
% Output:
%   beamline - Array of element structures

% 读取文件内容
if ~exist(bmadFile, 'file')
    error('File not found: %s', bmadFile);
end

fid = fopen(bmadFile, 'r');
if fid == -1
    error('Cannot open file: %s', bmadFile);
end

% 读取所有行
allLines = {};
while ~feof(fid)
    line = fgetl(fid);
    if ischar(line)
        allLines{end+1} = line;
    end
end
fclose(fid);

% 预处理：合并续行（以&结尾的行或以逗号结尾但下一行缩进的行）
processedLines = {};
i = 1;
while i <= length(allLines)
    currentLine = allLines{i};
    
    % 检查是否有续行符&或者行末有逗号且下一行是元素列表的一部分
    while i < length(allLines)
        % 检查续行符&
        if contains(currentLine, '&')
            currentLine = strrep(currentLine, '&', ' ');
            i = i + 1;
            if i <= length(allLines)
                nextLine = strtrim(allLines{i});
                currentLine = [currentLine, ' ', nextLine];
            end
        % 检查行末逗号且下一行可能是续行
        elseif endsWith(strtrim(currentLine), ',') && i < length(allLines)
            nextLine = strtrim(allLines{i+1});
            % 如果下一行不是新的定义（不包含冒号）且不是空行或注释
            if ~isempty(nextLine) && nextLine(1) ~= '!' && nextLine(1) ~= '%' && ~contains(nextLine, ':')
                i = i + 1;
                currentLine = [currentLine, ' ', nextLine];
            else
                break;
            end
        else
            break;
        end
    end
    
    processedLines{end+1} = currentLine;
    i = i + 1;
end

% 初始化
elementDefs = containers.Map('KeyType', 'char', 'ValueType', 'any');
lineDefs = containers.Map('KeyType', 'char', 'ValueType', 'any');
useLine = '';  % 用于存储USE命令指定的line

% 第一遍：解析所有定义
for i = 1:length(processedLines)
    line = strtrim(processedLines{i});
    
    % 跳过空行、注释行
    if isempty(line) || line(1) == '!' || line(1) == '%'
        continue;
    end
    
    % 检查USE命令
    if contains(upper(line), 'USE,')
        useLine = parseUseCommand(line);
        continue;
    end
    
    % 跳过参数行、BEGINNING行、CALL行
    if contains(upper(line), 'PARAMETER[') || ...
       contains(upper(line), 'BEGINNING[') || ...
       contains(upper(line), 'CALL,')
        continue;
    end
    
    % 检查是否包含冒号（元素或line定义）
    if ~contains(line, ':')
        continue;
    end
    
    % 分离名称和定义
    colonPos = strfind(line, ':');
    elementName = strtrim(line(1:colonPos(1)-1));
    definition = strtrim(line(colonPos(1)+1:end));
    
    % 检查是否为LINE定义 - 更宽松的匹配
    if contains(upper(definition), 'LINE') || (contains(definition, '(') && contains(definition, ')'))
        % 解析LINE定义
        lineElements = parseLineDefinition(definition);
        lineDefs(elementName) = lineElements;
        fprintf('Parsed LINE: %s with %d elements\n', elementName, length(lineElements));
    else
        % 解析普通元素定义
        element = parseElementDefinition(elementName, definition);
        if ~isempty(element)
            elementDefs(element.name) = element;
        end
    end
end

% 确定要使用的beamline - 只使用USE指定的line
finalBeamlineOrder = {};
if ~isempty(useLine)
    if isKey(lineDefs, useLine)
        finalBeamlineOrder = lineDefs(useLine);
        fprintf('Using LINE: %s\n', useLine);
    else
        error('USE命令指定的line "%s" 未找到！请检查文件中是否定义了该line。', useLine);
    end
else
    error('未找到USE命令！Bmad文件必须包含USE命令来指定要使用的beamline。');
end

% 展开beamline
beamline = expandBeamline(finalBeamlineOrder, elementDefs, lineDefs);

fprintf('Successfully parsed %d elements into beamline\n', length(beamline));

% 显示beamline摘要
if ~isempty(beamline)
    fprintf('\nBeamline summary:\n');
    elementTypes = {};
    for i = 1:length(beamline)
        elementTypes{end+1} = beamline(i).type;
    end
    [uniqueTypes, ~, idx] = unique(elementTypes);
    for i = 1:length(uniqueTypes)
        count = sum(idx == i);
        fprintf('  %s: %d elements\n', uniqueTypes{i}, count);
    end
    
    % 检查重复元素名称
    elementNames = {beamline.name};
    [uniqueNames, ~, nameIdx] = unique(elementNames);
    duplicateCount = 0;
    for i = 1:length(uniqueNames)
        count = sum(nameIdx == i);
        if count > 1
            duplicateCount = duplicateCount + count;
        end
    end
    if duplicateCount > 0
        fprintf('  Total elements with instance numbers: %d\n', duplicateCount);
    end
end

end

function useLine = parseUseCommand(line)
% 解析USE命令
useLine = '';

% 查找USE,
usePos = strfind(upper(line), 'USE,');
if isempty(usePos)
    return;
end

% 提取USE后面的内容
content = strtrim(line(usePos(1)+4:end));

% 移除注释
commentPos = strfind(content, '!');
if ~isempty(commentPos)
    content = content(1:commentPos(1)-1);
end
content = strtrim(content);

% 可能的格式：
% USE, linename
% USE, LINE=linename
if contains(upper(content), 'LINE=')
    equalPos = strfind(upper(content), 'LINE=');
    useLine = strtrim(content(equalPos(1)+5:end));
else
    useLine = content;
end

% 移除可能的分号
useLine = strrep(useLine, ';', '');
useLine = strtrim(useLine);

end

function lineElements = parseLineDefinition(definition)
% 解析LINE定义
lineElements = {};

% 找到等号
equalPos = strfind(definition, '=');
if isempty(equalPos)
    return;
end

content = definition(equalPos(1)+1:end);

% 移除LINE关键字（如果存在）
content = regexprep(content, '^\s*LINE\s*', '', 'ignorecase');

% 找到括号
openParen = strfind(content, '(');
closeParen = strfind(content, ')');

if isempty(openParen) || isempty(closeParen)
    % 如果没有括号，可能是简单格式：LINE = element1, element2, ...
    elementList = content;
else
    % 提取括号内的内容
    elementList = content(openParen(1)+1:closeParen(end)-1);
end

% 清理空格和换行
elementList = strrep(elementList, newline, ' ');
elementList = regexprep(elementList, '\s+', ' ');
elementList = strtrim(elementList);

% 移除可能的分号
elementList = strrep(elementList, ';', '');

% 按逗号分割
elements = strsplit(elementList, ',');

% 处理每个元素（包括负号表示反向）
for i = 1:length(elements)
    elementName = strtrim(elements{i});
    
    % 移除可能的括号
    elementName = strrep(elementName, '(', '');
    elementName = strrep(elementName, ')', '');
    elementName = strtrim(elementName);
    
    if ~isempty(elementName)
        lineElements{end+1} = elementName;
    end
end

% 调试输出
fprintf('DEBUG: Parsed %d elements from LINE definition\n', length(lineElements));

end

function elements = processDipoleEdges(element, instanceNum)
% 处理dipole元素的边缘角，创建edge元素
% instanceNum: 实例编号，用于区分同名元素的不同使用

elements = [];

% 创建元素的副本，避免修改原始定义
elementCopy = copyElement(element);

% 为元素添加实例编号（如果需要）
if instanceNum > 1
    elementCopy.name = sprintf('%s#%d', elementCopy.name, instanceNum);
end

if strcmpi(elementCopy.type, 'dipole')
    % 计算1/rho（曲率）
    if elementCopy.length > 0
        rho_inv = elementCopy.angle / elementCopy.length;  % 1/rho = angle/length
    elseif elementCopy.length == 0
        rho_inv = 0;
    end
    
    % 如果有E1，创建入口edge元素
    if elementCopy.e1 ~= 0
        edgeElement1 = createElement('edge', 0, rho_inv, elementCopy.e1, [elementCopy.name, '_e1']);
        elements = [elements, edgeElement1];
    end
    
    % 添加dipole本身（清除E1, E2）
    dipoleElement = elementCopy;
    dipoleElement.e1 = 0;
    dipoleElement.e2 = 0;
    dipoleElement = standardizeElement(dipoleElement);
    elements = [elements, dipoleElement];
    
    % 如果有E2，创建出口edge元素
    if elementCopy.e2 ~= 0
        edgeElement2 = createElement('edge', 0, rho_inv, elementCopy.e2, [elementCopy.name, '_e2']);
        elements = [elements, edgeElement2];
    end
else
    % 非dipole元素直接添加
    elementCopy = standardizeElement(elementCopy);
    elements = [elements, elementCopy];
end

end

function elementCopy = copyElement(element)
% 创建元素的深度副本
elementCopy = struct();
fields = fieldnames(element);
for i = 1:length(fields)
    elementCopy.(fields{i}) = element.(fields{i});
end
end

function element = createElement(type, length, h, angle, name)
% 创建标准化的元素结构
element = struct();
element.name = name;
element.type = type;
element.length = length;
element.k1 = 0;
element.k2 = 0;
element.angle = angle;
element.h = h;
element.e1 = 0;
element.e2 = 0;
element.voltage = 0;
element.frequency = 0;
element.phase = 0;
end

function element = standardizeElement(element)
% 标准化元素结构，确保所有元素都有相同的字段
% 定义所有可能的字段
allFields = {'name', 'type', 'length', 'k1', 'k2', 'angle', 'h', 'e1', 'e2', 'voltage', 'frequency', 'phase'};

% 创建标准化的元素结构
standardElement = struct();

for i = 1:length(allFields)
    fieldName = allFields{i};
    if isfield(element, fieldName)
        standardElement.(fieldName) = element.(fieldName);
    else
        % 为缺失的字段设置默认值
        switch fieldName
            case {'name', 'type'}
                standardElement.(fieldName) = '';
            otherwise
                standardElement.(fieldName) = 0;  % 数值字段默认为0
        end
    end
end

element = standardElement;
end

function element = parseElementDefinition(elementName, definition)
% 解析单个元素定义
element = [];

% 移除注释
commentPos = strfind(definition, '!');
if ~isempty(commentPos)
    definition = definition(1:commentPos(1)-1);
end
definition = strtrim(definition);

% 按逗号分割参数，但要注意处理空格
% 首先统一格式化：移除多余空格
definition = regexprep(definition, '\s*=\s*', '=');
definition = regexprep(definition, '\s*,\s*', ',');

parts = strsplit(definition, ',');
if isempty(parts)
    return;
end

% 第一部分是元素类型
bmadType = upper(strtrim(parts{1}));

% 初始化元素结构
element.name = elementName;
element.type = convertBmadType(bmadType);
element.length = 0;
element.angle = 0;
element.rho = 0;
element.h = 0;  % 曲率 (1/rho)

% 解析参数
for i = 2:length(parts)
    param = strtrim(parts{i});
    
    % 查找等号
    equalPos = strfind(param, '=');
    if isempty(equalPos)
        continue;
    end
    
    paramName = upper(strtrim(param(1:equalPos(1)-1)));
    paramValue = strtrim(param(equalPos(1)+1:end));
    
    % 转换参数值
    numValue = str2double(paramValue);
    if isnan(numValue)
        continue;
    end
    
    % 根据参数名称设置字段
    switch paramName
        case 'L'
            element.length = numValue;
        case 'K1'
            element.k1 = numValue;
        case 'K2'
            element.k2 = numValue;
        case 'ANGLE'
            element.angle = numValue;
        case {'RHO', 'RADIUS'}
            element.rho = numValue;
            if numValue ~= 0
                element.h = 1 / numValue;  % h = 1/rho
            end
        case 'G'
            element.h = numValue;  % G就是曲率 1/rho
            if numValue ~= 0
                element.rho = 1 / numValue;  % rho = 1/G
            end
        case 'E1'
            element.e1 = numValue;
        case 'E2'
            element.e2 = numValue;
        case 'VOLTAGE'
            element.voltage = numValue;
        case 'FREQ'
            element.frequency = numValue;
        case 'PHASE'
            element.phase = numValue;
    end
end

% 对于dipole元素，需要处理三种可能的定义方式
if strcmpi(element.type, 'dipole')
    % 情况1: RHO + ANGLE (最常见)
    if element.rho ~= 0 && element.angle ~= 0
        element.length = abs(element.angle * element.rho);
        element.h = 1 / element.rho;
        
    % 情况2: L + RHO  
    elseif element.length ~= 0 && element.rho ~= 0
        element.angle = element.length / element.rho;
        element.h = 1 / element.rho;
        
    % 情况3: L + ANGLE
    elseif element.length ~= 0 && element.angle ~= 0
        element.rho = element.length / element.angle;
        element.h = 1 / element.rho;
        
    % 情况4: G(曲率) + ANGLE
    elseif element.h ~= 0 && element.angle ~= 0
        element.rho = 1 / element.h;
        element.length = abs(element.angle * element.rho);
        
    % 情况5: L + G(曲率)
    elseif element.length ~= 0 && element.h ~= 0
        element.rho = 1 / element.h;
        element.angle = element.length / element.rho;
        

        
    % 情况7: 只有G，没有其他参数（thin lens情况）
    elseif element.h ~= 0 && element.length == 0 && element.angle == 0
        element.rho = 1 / element.h;
        % thin lens情况，长度为0
        
    % 错误情况：参数不足
    else
        if element.length == 0 && element.angle == 0 && element.rho == 0 && element.h == 0
            warning('Dipole %s: insufficient parameters for definition', elementName);
        end
    end
    
    % 如果只有L和G，计算角度
    if element.length ~= 0 && element.h ~= 0 && element.angle == 0
        element.angle = element.h * element.length;
        element.rho = element.length / element.angle;
    end
    
    % 确保所有参数的一致性检查
    if element.rho ~= 0 && element.h ~= 0
        % 检查rho和h是否一致
        if abs(element.h - 1/element.rho) > 1e-10
            warning('Dipole %s: inconsistent RHO and G values', elementName);
        end
    end
    
    if element.length ~= 0 && element.angle ~= 0 && element.rho ~= 0
        % 检查L, ANGLE, RHO是否一致
        expectedLength = abs(element.angle * element.rho);
        if abs(element.length - expectedLength) > 1e-10
            warning('Dipole %s: inconsistent L, ANGLE, RHO values', elementName);
        end
    end
end

% 标准化元素
element = standardizeElement(element);

end

function beamline = expandBeamline(beamlineOrder, elementDefs, lineDefs)
% 递归展开beamline，处理嵌套的LINE定义和反向元素
beamline = [];

% 用于跟踪每个元素名称的使用次数
elementUsageCount = containers.Map('KeyType', 'char', 'ValueType', 'int32');

for i = 1:length(beamlineOrder)
    elementName = beamlineOrder{i};
    
    % 检查是否为反向元素（以-开头）
    isReversed = false;
    if length(elementName) > 1 && elementName(1) == '-'
        isReversed = true;
        elementName = elementName(2:end);
    end
    
    % 检查是否为LINE定义
    if isKey(lineDefs, elementName)
        % 递归展开LINE
        subElements = lineDefs(elementName);
        if isReversed
            % 反向：只是颠倒顺序，不改变元素本身的属性
            subElements = subElements(end:-1:1);
            for j = 1:length(subElements)
                if length(subElements{j}) > 1 && subElements{j}(1) == '-'
                    subElements{j} = subElements{j}(2:end); % 双重负号变正
                else
                    subElements{j} = ['-', subElements{j}];
                end
            end
        end
        
        % 递归展开
        subBeamline = expandBeamline(subElements, elementDefs, lineDefs);
        beamline = [beamline, subBeamline];
        
    elseif isKey(elementDefs, elementName)
        % 普通元素
        
        % 更新使用计数
        if isKey(elementUsageCount, elementName)
            elementUsageCount(elementName) = elementUsageCount(elementName) + 1;
        else
            elementUsageCount(elementName) = 1;
        end
        
        instanceNum = elementUsageCount(elementName);
        
        % 获取元素定义的副本
        element = elementDefs(elementName);
        
        if isReversed
            % 创建反向版本的名称
            if instanceNum > 1
                element.name = sprintf('%s#%d_rev', elementName, instanceNum);
            else
                element.name = [elementName, '_rev'];
            end
        else
            % 如果是重复使用，添加实例编号
            if instanceNum > 1
                element.name = sprintf('%s#%d', elementName, instanceNum);
            end
        end
        
        % 处理dipole的edge元素
        elementWithEdges = processDipoleEdges(element, 1); % 这里传1因为已经在上面处理了实例编号
        
        % 如果是反向，需要颠倒edge的顺序
        if isReversed && length(elementWithEdges) > 1
            elementWithEdges = elementWithEdges(end:-1:1);
        end
        
        beamline = [beamline, elementWithEdges];
        
    else
        error('Element %s not found in definitions', elementName);
    end
end

end

function elementType = convertBmadType(bmadType)
% 将Bmad元素类型转换为标准类型
switch bmadType
    case 'QUAD'
        elementType = 'quadrupole';
    case 'QUADRUPOLE'
        elementType = 'quadrupole';
    case {'SBEND', 'RBEND'}
        elementType = 'dipole';
    case 'DRIFT'
        elementType = 'drift';
    case 'SEXT'
        elementType = 'sextupole';
    case 'SEXTUPOLE'
        elementType = 'sextupole';
    case 'MARKER'
        elementType = 'marker';
    case 'RFCAVITY'
        elementType = 'cavity';
    case 'EDGE'
        elementType = 'edge';
    case 'KICKER'
        elementType = 'kicker';
    case 'MONITOR'
        elementType = 'monitor';
    case 'MULTIPOLE'
        elementType = 'multipole';
    otherwise
        elementType = lower(bmadType);
end
end