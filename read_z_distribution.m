function slice_data = read_z_distribution(infile, num_slices)
    % 如果没有指定切片数量，默认使用20
    if nargin < 2
        num_slices = 20;
    end
    
    % 读取粒子总数
    [np] = textread(infile,'%f%*s%*s%*s',1,'headerlines',3);
    
    % 读取所有粒子坐标 (x, px, y, py, z, pz)
    [x, px, y, py, z, pz] = textread(infile,'%n%n%n%n%n%n%*f%*f%*f%*f%*f%*f%*f',np,'headerlines',9);
    
    % 计算能量偏差 delta = pz (假设pz已经是相对动量偏差)
    delta = pz;
    
    % 计算z的范围并创建等间距切片
    z_min = min(z);
    z_max = max(z);
    z_length = z_max - z_min;  % 计算z的总长度
    z_edges = linspace(z_min, z_max, num_slices + 1);
    
    % 初始化结果矩阵 (num_slices+2) × 17
    % 列定义：
    % 1: 粒子数（第一行为总粒子数，最后一行为z_length）
    % 2: betax
    % 3: alphax
    % 4: gammax
    % 5: emitx_geo（几何发射度）
    % 6: betay
    % 7: alphay
    % 8: gammay
    % 9: emity_geo（几何发射度）
    % 10: x_centroid（切片x中心）
    % 11: px_centroid（切片px中心）
    % 12: y_centroid（切片y中心）
    % 13: py_centroid（切片py中心）
    % 14: z_centroid（切片z中心）
    % 15: delta_centroid（切片能量偏差中心）
    % 16: <x*pz>（二阶关联项）
    % 17: <px*pz>（二阶关联项）
    slice_data = zeros(num_slices + 2, 17);
    
    % 第一行：总粒子数和全局中心坐标及全局二阶关联项
    slice_data(1, 1) = np;  % 总粒子数
    slice_data(1, 10) = mean(x);      % 全局x中心
    slice_data(1, 11) = mean(px);     % 全局px中心
    slice_data(1, 12) = mean(y);      % 全局y中心
    slice_data(1, 13) = mean(py);     % 全局py中心
    slice_data(1, 14) = mean(z);      % 全局z中心
    slice_data(1, 15) = mean(delta);  % 全局delta中心
    slice_data(1, 16) = mean(x .* delta);   % 全局<x*pz>
    slice_data(1, 17) = mean(px .* delta);  % 全局<px*pz>
    
    % 最后一行：z_length
    slice_data(end, 1) = z_length;  % z_length放在第一列最后一行
    
    % 计算每个切片的光学参数
    for i = 1:num_slices
        % 找到当前切片的粒子
        slice_idx = (z >= z_edges(i) & z < z_edges(i+1));
        
        % 如果最后一个切片，包含右边界
        if i == num_slices
            slice_idx = (z >= z_edges(i) & z <= z_edges(i+1));
        end
        
        slice_particles = sum(slice_idx);
        slice_data(i+1, 1) = slice_particles;  % 切片粒子数
        
        % 如果粒子数太少，给出警告但继续计算
        if slice_particles < 20
            warning('切片 %d 粒子数较少 (%d < 20)，统计精度可能不足。', i, slice_particles);
        end
        
        if slice_particles > 0  % 只要有粒子就尝试计算
            % 提取当前切片的坐标
            x_slice = x(slice_idx);
            px_slice = px(slice_idx);
            y_slice = y(slice_idx);
            py_slice = py(slice_idx);
            z_slice = z(slice_idx);
            delta_slice = delta(slice_idx);
            
            % 计算并存储切片中心位置
            x_mean = mean(x_slice);
            px_mean = mean(px_slice);
            y_mean = mean(y_slice);
            py_mean = mean(py_slice);
            z_mean = mean(z_slice);
            delta_mean = mean(delta_slice);
            
            slice_data(i+1, 10) = x_mean;     % x中心
            slice_data(i+1, 11) = px_mean;    % px中心
            slice_data(i+1, 12) = y_mean;     % y中心
            slice_data(i+1, 13) = py_mean;    % py中心
            slice_data(i+1, 14) = z_mean;     % z中心
            slice_data(i+1, 15) = delta_mean; % delta中心
            
            % 计算二阶关联项 <x*pz> 和 <px*pz>
            xpz_correlation = mean(x_slice .* delta_slice);
            pxpz_correlation = mean(px_slice .* delta_slice);
            
            slice_data(i+1, 16) = xpz_correlation;  % <x*pz>
            slice_data(i+1, 17) = pxpz_correlation; % <px*pz>
            
            % 计算相对于切片中心的坐标（用于发射度计算）
            x_rel = x_slice - x_mean;
            px_rel = px_slice - px_mean;
            y_rel = y_slice - y_mean;
            py_rel = py_slice - py_mean;
            
            % 计算二阶矩（相对于切片中心）
            % X方向
            sigma11_x = mean(x_rel.^2);      % <x²>
            sigma12_x = mean(x_rel .* px_rel); % <x·px>
            sigma22_x = mean(px_rel.^2);     % <px²>
            
            % Y方向
            sigma11_y = mean(y_rel.^2);      % <y²>
            sigma12_y = mean(y_rel .* py_rel); % <y·py>
            sigma22_y = mean(py_rel.^2);     % <py²>
            
            % 计算几何发射度（添加数值保护）
            emit_x_geo = sqrt(max(0, sigma11_x * sigma22_x - sigma12_x^2));
            emit_y_geo = sqrt(max(0, sigma11_y * sigma22_y - sigma12_y^2));
            
            % 避免数值问题
            if emit_x_geo < 1e-15
                emit_x_geo = 1e-15;
            end
            if emit_y_geo < 1e-15
                emit_y_geo = 1e-15;
            end
            
            % 计算Twiss参数（相对于切片中心）
            % β = <x²>/ε, α = -<x·px>/ε, γ = <px²>/ε
            % X方向
            beta_x = sigma11_x / emit_x_geo;
            alpha_x = -sigma12_x / emit_x_geo;
            gamma_x = sigma22_x / emit_x_geo;
            
            % Y方向
            beta_y = sigma11_y / emit_y_geo;
            alpha_y = -sigma12_y / emit_y_geo;
            gamma_y = sigma22_y / emit_y_geo;
            
            % 存储结果
            slice_data(i+1, 2) = beta_x;
            slice_data(i+1, 3) = alpha_x;
            slice_data(i+1, 4) = gamma_x;
            slice_data(i+1, 5) = emit_x_geo;
            slice_data(i+1, 6) = beta_y;
            slice_data(i+1, 7) = alpha_y;
            slice_data(i+1, 8) = gamma_y;
            slice_data(i+1, 9) = emit_y_geo;
            
        else
            % 如果切片中没有粒子，设为默认值
            warning('切片 %d 无粒子，设为默认值。', i);
            slice_data(i+1, 2:9) = [1.0, 0.0, 1.0, 1e-9, 1.0, 0.0, 1.0, 1e-9];
            slice_data(i+1, 10:17) = [0, 0, 0, 0, 0, 0, 0, 0];  % 中心位置和关联项设为0
        end
    end
    
    % 生成输出文件名
    [filepath, name, ext] = fileparts(infile);
    output_filename = fullfile(filepath, [name, '_slice_params.txt']);
    
    % 写入txt文件
    fid = fopen(output_filename, 'w');
    if fid == -1
        error('无法创建输出文件: %s', output_filename);
    end
    
    % 写入文件头
    fprintf(fid, '%% 切片光学参数分析结果\n');
    fprintf(fid, '%% 源文件: %s\n', infile);
    fprintf(fid, '%% 切片数: %d\n', num_slices);
    fprintf(fid, '%% 总粒子数: %d\n', np);
    fprintf(fid, '%% z方向长度: %.6f m\n', z_length);
    fprintf(fid, '%% 生成时间: %s\n', datestr(now));
    fprintf(fid, '%% \n');
    fprintf(fid, '%% 矩阵格式说明:\n');
    fprintf(fid, '%% 第1行: 总粒子数, 0~8列=0, 全局x_cent, px_cent, y_cent, py_cent, z_cent, delta_cent, <x*pz>, <px*pz>\n');
    fprintf(fid, '%% 第2~%d行: 切片粒子数, betax, alphax, gammax, emitx_geo, betay, alphay, gammay, emity_geo, x_cent, px_cent, y_cent, py_cent, z_cent, delta_cent, <x*pz>, <px*pz>\n', num_slices+1);
    fprintf(fid, '%% 最后一行: z_length, 0, ..., 0\n');
    fprintf(fid, '%% \n');
    fprintf(fid, '%% 列定义:\n');
    fprintf(fid, '%% 1: 粒子数(第1行=总数,最后行=z_length)  2: betax(m)      3: alphax       4: gammax(1/m)\n');
    fprintf(fid, '%% 5: emitx_geo(m*rad)                    6: betay(m)      7: alphay       8: gammay(1/m)\n');
    fprintf(fid, '%% 9: emity_geo(m*rad)                    10: x_cent(m)    11: px_cent     12: y_cent(m)\n');
    fprintf(fid, '%% 13: py_cent                            14: z_cent(m)    15: delta_cent  16: <x*pz>\n');
    fprintf(fid, '%% 17: <px*pz>\n');
    fprintf(fid, '%% \n');
    fprintf(fid, '%% 注意: 发射度和Twiss参数是相对于各切片中心计算的\n');
    fprintf(fid, '%% 注意: <x*pz>和<px*pz>是绝对值（未减去中心）\n');
    fprintf(fid, '%% \n');
    
    % 写入数据
    for i = 1:size(slice_data, 1)
        fprintf(fid, '%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n', slice_data(i, :));
    end
    
    fclose(fid);
    
    fprintf('读取粒子文件完成:\n');
    fprintf('总粒子数: %d\n', np);
    fprintf('z方向长度: %.6f m\n', z_length);
    fprintf('切片数: %d\n', num_slices);
    fprintf('全局 <x*pz>: %.6e\n', slice_data(1, 16));
    fprintf('全局 <px*pz>: %.6e\n', slice_data(1, 17));
    fprintf('结果已保存到: %s\n', output_filename);
end