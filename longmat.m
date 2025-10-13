function [energy_spread, rms_z, skewness, kurtosis, min_derivative] = longmat(beam_file, EA, E0, L, n_slices, phi01, phi02, r1, r2, T566_1, T566_2, U5666_1, U5666_2, plot_flag, write_flag)
% BEAM_TRANSPORT_ANALYSIS 束流传输分析函数

    if nargin < 14
        plot_flag = 0;  % 默认不绘图
    end
    if nargin < 15
        write_flag = 0;  % 默认不写文件
    end
    
    freq = 1300;
    EA1 = EA;
    EA2 = EA;
    phi1 = phi01*pi/180;
    phi2 = phi02*pi/180;
    
    %% 读取束流数据
    [np] = textread(beam_file, '%f%*s%*s%*s', 1, 'headerlines', 3);
    [nc] = textread(beam_file, '%f%*s%*s%*s', 1, 'headerlines', 6);
    [xs, xp, y, yp, z, p] = textread(beam_file, '%n%n%n%n%n%n%*f%*f%*f%*f%*f%*f%*f', np, 'headerlines', 9);
    
    %% 如果需要写文件，用所有粒子进行计算
    if write_flag == 1
        % 生成输出文件名
        [filepath, name, ext] = fileparts(beam_file);
        outputfile = fullfile(filepath, [name '_after_ac2' ext]);
        
        % 用所有粒子计算AC2出口状态并写入文件
        write_all_particles_after_ac2(xs, xp, y, yp, z, p, T566_1, T566_2, U5666_1, U5666_2, phi1, phi2, r1, r2, EA1, EA2, E0, L, freq, outputfile, nc);
    end
    
    %% 创建纵向切片和权重宏粒子（用于统计分析）
    z_min = min(z);
    z_max = max(z);
    z_edges = linspace(z_min, z_max, n_slices+1);
    z_centers = (z_edges(1:end-1) + z_edges(2:end)) / 2;
    
    macro_z = [];
    macro_p = [];
    macro_weight = [];
    
    for i = 1:n_slices
        in_slice = (z >= z_edges(i)) & (z < z_edges(i+1));
        
        if sum(in_slice) > 0
            macro_z = [macro_z; z_centers(i)];
            macro_p = [macro_p; mean(p(in_slice))];
            macro_weight = [macro_weight; sum(in_slice)];
        end
    end
    
    %% 用切片计算束流特性
    [energy_spread, rms_z, skewness, kurtosis, min_derivative, final_z, final_p] = calculate_beam_properties(macro_z, macro_p, macro_weight, T566_1, T566_2, U5666_1, U5666_2, phi1, phi2, r1, r2, EA1, EA2, E0, L, freq);
    
    %% 绘图 - 使用和写文件相同的宏粒子追踪方法
    if plot_flag == 1
        plot_beam_distribution_macro_tracking(xs, xp, y, yp, z, p, T566_1, T566_2, U5666_1, U5666_2, phi1, phi2, r1, r2, EA1, EA2, E0, L, freq, nc, energy_spread, rms_z, skewness, kurtosis, min_derivative);
    end
end

function write_all_particles_after_ac2(xs, xp, y, yp, z, p, T566_1, T566_2, U5666_1, U5666_2, phi1, phi2, r1, r2, EA1, EA2, E0, L, freq, outputfile, nc)
    % 用所有粒子进行追踪，计算AC2出口状态并写入文件
    
    krf_track = 2*pi/(3*10^8/freq/10^6);
    mean_p = mean(p);
    
    % 向量化计算，避免循环
    gamma = E0/0.511;
    z_current = z + L/gamma^2.*p;
    
    % 第一次加速 - 向量化
    mer2 = (p+1)*E0;
    linac2 = mer2 + EA1*cos(krf_track*z_current + phi1);
    mean_ene = E0*(1+mean_p) + EA1*cos(phi1);
    
    % 第一次加速后的相对能量偏差
    phase_6 = (linac2 - mean_ene)/mean_ene;
    
    % 第一级压缩
    z_after_bc1 = z_current + r1*phase_6 + T566_1*phase_6.^2 + U5666_1*phase_6.^3;
    
    % 第二次加速 - 向量化
    mer2_2 = (phase_6+1)*mean_ene;
    linac2_2 = mer2_2 + EA2*cos(krf_track*z_after_bc1 + phi2);
    mean_ene_2 = mean_ene + EA2*cos(phi2);
    
    % 第二次加速后的相对能量偏差
    linac22_2 = (linac2_2 - mean_ene_2)/mean_ene_2;
    
    % 创建输出数组
    phase_after_ac2 = [xs, xp, y, yp, z_after_bc1, linac22_2];
    
    % 归一化横向动量
    phase_after_ac2(:,2) = phase_after_ac2(:,2) ./ mean_ene_2 .* E0;
    phase_after_ac2(:,4) = phase_after_ac2(:,4) ./ mean_ene_2 .* E0;
    
    % 写入文件
    write_particles_to_file(phase_after_ac2, outputfile, nc);
end

function [energy_spread, rms_z, skewness, kurtosis, min_derivative, final_z, final_p] = calculate_beam_properties(macro_z, macro_p, macro_weight, T566_1, T566_2, U5666_1, U5666_2, phi1, phi2, r1, r2, EA1, EA2, E0, L, freq)
    
    % 先计算稳定性（最小导数）
    min_derivative = calculate_min_derivative(macro_z, macro_p, T566_1, T566_2, U5666_1, U5666_2, phi1, phi2, r1, r2, EA1, EA2, E0, L, freq);
    
    % 计算最终位置和能量
    try
        [final_z, final_p] = transport_function_full(macro_z, macro_p, T566_1, T566_2, U5666_1, U5666_2, phi1, phi2, r1, r2, EA1, EA2, E0, L, freq);
        
        % 检查结果有效性
        if any(~isfinite(final_z)) || any(~isfinite(final_p))
            error('传输结果包含无效值');
        end
        
    catch
        % 如果传输计算失败，用默认值
        energy_spread = 1e6;  % 很大的能散表示发散
        rms_z = 1e6;         % 很大的束长表示发散
        skewness = 0;
        kurtosis = 3;
        final_z = macro_z;   % 返回原始位置
        final_p = macro_p;   % 返回原始动量
        return;
    end
    
    % 计算加权统计量
    total_weight = sum(macro_weight);
    
    % 能散计算
    weighted_mean_p = sum(final_p .* macro_weight) / total_weight;
    energy_spread = sqrt(sum(macro_weight .* (final_p - weighted_mean_p).^2) / total_weight);
    
    % 位置统计量
    weighted_mean_z = sum(final_z .* macro_weight) / total_weight;
    weighted_var_z = sum(macro_weight .* (final_z - weighted_mean_z).^2) / total_weight;
    rms_z = sqrt(weighted_var_z);
    
    % 偏度和峰度
    if rms_z > 0
        moment3 = sum(macro_weight .* (final_z - weighted_mean_z).^3) / total_weight;
        skewness = moment3 / (rms_z^3);
        
        moment4 = sum(macro_weight .* (final_z - weighted_mean_z).^4) / total_weight;
        kurtosis = moment4 / (rms_z^4);
    else
        skewness = 0;
        kurtosis = 3;
    end
    
    % 确保所有输出都是有限的
    if ~isfinite(energy_spread)
        energy_spread = 1e6;
    end
    if ~isfinite(rms_z)
        rms_z = 1e6;
    end
    if ~isfinite(skewness)
        skewness = 0;
    end
    if ~isfinite(kurtosis)
        kurtosis = 3;
    end
end

function [final_z, final_p] = transport_function_full(z_in, p_in, T566_1, T566_2, U5666_1, U5666_2, phi1, phi2, r1, r2, EA1, EA2, E0, L, freq)
    
    krf_track = 2*pi/(3*10^8/freq/10^6);
    mean_p = mean(p_in);
    
    final_z = zeros(size(z_in));
    final_p = zeros(size(p_in));
    
    for i = 1:length(z_in)
        gamma = E0/0.511;
        z_current = z_in(i) + L/gamma^2*p_in(i);
        
        % 第一次加速
        mer2 = (p_in(i)+1)*E0;
        linac2 = mer2 + EA1*cos(krf_track*z_current + phi1);
        mean_ene = E0*(1+mean_p) + EA1*cos(phi1);
        
        % 第一次加速后的相对能量偏差
        phase_6 = (linac2 - mean_ene)/mean_ene;
        
        % 第一级压缩（包含U5666项）
        z_after_bc1 = z_current + r1*phase_6 + T566_1*phase_6^2 + U5666_1*phase_6^3;
        
        % 第二次加速
        mer2_2 = (phase_6+1)*mean_ene;
        linac2_2 = mer2_2 + EA2*cos(krf_track*z_after_bc1 + phi2);
        mean_ene_2 = mean_ene + EA2*cos(phi2);
        
        % 第二次加速后的相对能量偏差
        linac22_2 = (linac2_2 - mean_ene_2)/mean_ene_2;
        
        % 第二级压缩（包含U5666项）
        final_z(i) = z_after_bc1 + r2*linac22_2 + T566_2*linac22_2^2 + U5666_2*linac22_2^3;
        final_p(i) = linac22_2;
    end
end

function write_particles_to_file(phase, outputfile, nc)
    
    fid = fopen(outputfile, 'w');
    if fid == -1
        error('无法创建输出文件: %s', outputfile);
    end
    
    % 写入文件头
    a = {'!','ASCII::3';
         '0' '   ! ix_ele'; 
         '1' '   ! n_bunch';
         size(phase,1) '   ! n_particle';
         [] 'BEGIN_BUNCH';
         [] 'Electron';
         nc '   ! bunch_charge_tot';
         -0 '        ! z_center'; 
         0 '        ! t_center'};
    
    for i = 1:3
        fprintf(fid, '%s', a{i,1});
        fprintf(fid, '%s\n', a{i,2});
    end
    
    for i = 4:9
        fprintf(fid, '%d         ', a{i,1});
        fprintf(fid, '%s\n', a{i,2});
    end
    
    % 写入粒子数据
    for j = 1:size(phase,1)
        fprintf(fid, '%.15f                 ', phase(j,:));
        fprintf(fid, '\n');
    end
    
    fprintf(fid, 'END_BUNCH \n');
    fclose(fid);
    
    fprintf('已将 %d 个粒子写入文件: %s\n', size(phase,1), outputfile);
end

function min_derivative = calculate_min_derivative(macro_z, macro_p, T566_1, T566_2, U5666_1, U5666_2, phi1, phi2, r1, r2, EA1, EA2, E0, L, freq)
    
    min_derivative = Inf;
    dz = 1e-9;
    
    % 检查所有宏粒子点的导数
    for k = 1:length(macro_z)
        try
            % 原始点
            [final_z_orig, ~] = transport_function_full(macro_z(k), macro_p(k), T566_1, T566_2, U5666_1, U5666_2, phi1, phi2, r1, r2, EA1, EA2, E0, L, freq);
            
            % 扰动点
            z_pert = macro_z(k) + dz;
            [final_z_pert, ~] = transport_function_full(z_pert, macro_p(k), T566_1, T566_2, U5666_1, U5666_2, phi1, phi2, r1, r2, EA1, EA2, E0, L, freq);
            
            % 数值导数
            derivative = (final_z_pert - final_z_orig) / dz;
            
            % 更新最小导数
            if isfinite(derivative)
                min_derivative = min(min_derivative, derivative);
            end
            
        catch
            % 如果计算失败，认为导数为负无穷（极不稳定）
            min_derivative = -Inf;
            return;
        end
    end
    
    % 如果所有计算都失败
    if min_derivative == Inf
        min_derivative = -1e6;
    end
end

function plot_beam_distribution_macro_tracking(xs, xp, y, yp, z, p, T566_1, T566_2, U5666_1, U5666_2, phi1, phi2, r1, r2, EA1, EA2, E0, L, freq, nc, energy_spread, rms_z, skewness, kurtosis, min_derivative)
    % 使用宏粒子追踪进行绘图 - 与写文件完全相同的物理计算
    
    fprintf('开始宏粒子追踪计算用于绘图，粒子数: %d\n', length(z));
    tic;
    
    krf_track = 2*pi/(3*10^8/freq/10^6);
    mean_p = mean(p);
    
    % 向量化计算以提高速度
    gamma = E0/0.511;
    z_current = z + L/gamma^2.*p;
    
    % 第一次加速 - 向量化
    mer2 = (p+1)*E0;
    linac2 = mer2 + EA1*cos(krf_track*z_current + phi1);
    mean_ene = E0*(1+mean_p) + EA1*cos(phi1);
    
    % 第一次加速后的相对能量偏差
    phase_6 = (linac2 - mean_ene)/mean_ene;
    
    % 第一级压缩
    z_after_bc1 = z_current + r1*phase_6 + T566_1*phase_6.^2 + U5666_1*phase_6.^3;
    
    % 第二次加速 - 向量化
    mer2_2 = (phase_6+1)*mean_ene;
    linac2_2 = mer2_2 + EA2*cos(krf_track*z_after_bc1 + phi2);
    mean_ene_2 = mean_ene + EA2*cos(phi2);
    
    % 第二次加速后的相对能量偏差
    linac22_2 = (linac2_2 - mean_ene_2)/mean_ene_2;
    
    % 第二级压缩 - 最终位置
    final_z_all = z_after_bc1 + r2*linac22_2 + T566_2*linac22_2.^2 + U5666_2*linac22_2.^3;
    final_p_all = linac22_2;
    
    fprintf('宏粒子追踪计算完成，耗时: %.2f秒\n', toc);
    
    % 开始绘图
    figure('Color', 'white', 'Position', [100, 100, 800, 600]);
    
    % 设置Cambridge学术风格
    set(0, 'DefaultAxesFontName', 'Cambria');
    set(0, 'DefaultTextFontName', 'Cambria');
    set(0, 'DefaultAxesFontSize', 16);
    set(0, 'DefaultTextFontSize', 16);
    
    % 使用所有宏粒子的最终位置创建直方图
    z_final_mm = final_z_all * 1e3;  % 转换为mm
    
    % 中心化
    z_final_mm = z_final_mm - mean(z_final_mm);
    
    % 创建直方图 - 自适应bin数量
    n_particles = length(z_final_mm);
    n_bins = min(100, max(20, round(sqrt(n_particles))));  % 自适应bin数量
    
    % 计算直方图
    [counts, edges] = histcounts(z_final_mm, n_bins);
    bin_centers = (edges(1:end-1) + edges(2:end)) / 2;
    bin_width = edges(2) - edges(1);
    
    % 计算流强密度
    % 流强 = 每个bin的粒子数 / 总粒子数 * 总电荷 / (bin宽度/光速)
    current_density = counts / n_particles * nc / (bin_width * 1e-3 / 3e8);  % A
    
    % 绘制流强截面图
    plot(bin_centers, current_density, 'k-', 'LineWidth', 2.5);
    hold on;
    
    % 填充区域
    fill([bin_centers, fliplr(bin_centers)], [current_density, zeros(size(current_density))], ...
         [0.7 0.7 0.7], 'EdgeColor', 'k', 'LineWidth', 1.5, 'FaceAlpha', 0.3);
    
    % 设置坐标轴
    xlabel('Longitudinal Position z (mm)', 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Cambria');
    ylabel('Current (A)', 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Cambria');
    
    % 设置图形属性
    grid off;
    box on;
    set(gca, 'LineWidth', 2, 'FontSize', 18, 'FontName', 'Cambria');
    
    % 计算统计量
    peak_current = max(current_density);
    
    % 添加统计信息文本框
    stability_status = '';
    if min_derivative > 0
        stability_status = 'Stable';
    else
        stability_status = 'Unstable';
    end
    
    info_text = sprintf(['Peak Current: %.2f A\n' ...
                        'Energy Spread: %.3f %%\n' ...
                        'RMS Bunch Length: %.2f μm\n' ...
                        'Skewness: %.3f\n' ...
                        'Kurtosis: %.3f\n' ...
                        'Status: %s'], ...
                        peak_current, 100*energy_spread, rms_z*1e6, ...
                        skewness, kurtosis,  stability_status);
    
    % 添加文本框
    text(0.02, 0.98, info_text, 'Units', 'normalized', ...
         'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
         'BackgroundColor', 'white', 'EdgeColor', 'black', 'LineWidth', 1.5, ...
         'FontSize', 10, 'FontName', 'Cambria', 'FontWeight', 'normal');
    
    % 设置坐标轴范围
    axis tight;
    xlim_range = xlim;
    ylim_range = ylim;
    xlim([xlim_range(1) - 0.02*diff(xlim_range), xlim_range(2) + 0.02*diff(xlim_range)]);
    ylim([0, ylim_range(2) * 1.05]);
    
    % 设置刻度
    set(gca, 'TickDir', 'in', 'TickLength', [0.02, 0.02]);
    
    % 标记峰值点
    if ~isempty(current_density) && max(current_density) > 0
        [max_current, max_idx] = max(current_density);
        peak_position = bin_centers(max_idx);
        
        plot(peak_position, max_current, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
        
        % 添加峰值标注
        text(peak_position, max_current + 0.05*max_current, ...
             sprintf('Peak: %.2f mm\n%.2f A', peak_position, max_current), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
             'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Cambria');
    end
    
    
    
    hold off;
    
    % 恢复默认字体设置
    set(0, 'DefaultAxesFontName', 'remove');
    set(0, 'DefaultTextFontName', 'remove');
    set(0, 'DefaultAxesFontSize', 'remove');
    set(0, 'DefaultTextFontSize', 'remove');
    
    fprintf('绘图完成\n');
end