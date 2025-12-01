function [results] = SCEnvwe2(beamline, init_params, beamfile1)
% SCEnv - 统一的空间电荷包络计算入口函数（优化版本，整合输出）
%
% 输出：
%   results - 包含所有结果的结构体：
%     .s_array - 纵向位置
%     .sigx_array - x方向包络
%     .sigy_array - y方向包络
%     .sigpx_array - x方向包络导数
%     .sigpy_array - y方向包络导数
%     .x_centroid - x方向中心轨迹位置
%     .y_centroid - y方向中心轨迹位置
%     .px_centroid - x方向中心轨迹角度
%     .py_centroid - y方向中心轨迹角度
%     .gamma_array - 相对论因子γ
%     .k1_values - 四极铁归一化梯度（如适用）
%     .emitx_proj_array - x方向投影发射度（如适用）
%     .emity_proj_array - y方向投影发射度（如适用）

betax = init_params(1);
betay = init_params(2);
alphax = init_params(3);
alphay = init_params(4);
emitx = init_params(5);
emity = init_params(6);
sigz1 = init_params(7);

sigz = 100;
Q_total = init_params(8);
E0 = init_params(9) + 0.0e7;
calc_k1_flag = 1;

me_c2 = 0.511e6;  % 提前定义，两个分支都需要

if nargin < 3
    % 单一包络模式
    init_params_single = [betax, betay, alphax, alphay, emitx, emity, sigz1, Q_total, E0];
    [s_array, sigx_array, sigy_array, k1_values] = SCEqga(beamline, init_params_single, 1);
    
    % ========== 单包络模式：计算能量演化 ==========
    gamma_array = zeros(size(s_array));
    current_E = E0;
    gamma_array(1) = current_E / me_c2;
    
    current_s = 0;
    point_idx = 1;
    
    for i = 1:length(beamline)
        element = beamline(i);
        
        if strcmpi(element.type, 'cavity')
            % 计算cavity的能量增益
            if isfield(element, 'num_cells') && isfield(element, 'frequency')
                element.length = element.num_cells * 0.5 * 3e8 / (element.frequency * 1e6);
            end
            
            A = 0;
            if isfield(element, 'gradient')
                A = element.gradient * 1e6;
            end
            
            phi = 0;
            if isfield(element, 'phase')
                phi = element.phase * pi/180;
            end
            
            energy_gain = A * element.length * cos(phi);
            current_E = current_E + energy_gain;
        end
        
        % 更新s范围内的所有点的gamma值
        s_end = current_s + element.length;
        while point_idx <= length(s_array) && s_array(point_idx) <= s_end
            gamma_array(point_idx) = current_E / me_c2;
            point_idx = point_idx + 1;
        end
        
        current_s = s_end;
    end
    
    % 填充剩余点
    gamma_array(point_idx:end) = current_E / me_c2;
    
    % 整合输出
    results = struct();
    results.s_array = s_array;
    results.sigx_array = sigx_array;
    results.sigy_array = sigy_array;
    results.sigpx_array = [];
    results.sigpy_array = [];
    results.x_centroid = [];
    results.y_centroid = [];
    results.px_centroid = [];
    results.py_centroid = [];
    results.gamma_array = gamma_array;
    results.k1_values = k1_values;
    results.emitx_proj_array = [];
    results.emity_proj_array = [];
    return;
else
    % 基于切片的多包络计算模式
    beamfile = beamfile1;
    
    try
        slice_data_matrix = load(beamfile);
    catch
        error('无法读取文件: %s', beamfile);
    end
    
    % 解析切片数据
    total_particles = slice_data_matrix(1, 1);
    z_length = 1.2*slice_data_matrix(end, 1);
    slice_info = slice_data_matrix(2:end-1, :);
    
    num_slice = size(slice_info, 1);
    slice_particle_counts = slice_info(:, 1);
    betax_slices = slice_info(:, 2);
    alphax_slices = slice_info(:, 3);
    emitx_geo_slices = slice_info(:, 5);
    betay_slices = slice_info(:, 6);
    alphay_slices = slice_info(:, 7);
    emity_geo_slices = slice_info(:, 9);
    
    % 读取切片中心位置
    x_centroid_slices = slice_info(:, 10);
    px_centroid_slices = slice_info(:, 11);
    y_centroid_slices = slice_info(:, 12);
    py_centroid_slices = slice_info(:, 13);
    
    % 物理常数（预计算）
    rc = 2.8179e-15;
    e = 1.602e-19;
    c = 2.998e8;
    
    gamma0 = E0 / me_c2;
    z_val = z_length / num_slice;
    
    if calc_k1_flag
        k1_values = [];
    end
    
    quad_gradients = containers.Map();
    options = odeset('RelTol', 5e-6, 'AbsTol', 5e-8, 'MaxStep', 0.03, 'Stats', 'off');
    
    % **预分配所有切片的存储空间**
    total_length = sum([beamline.length]);
    estimated_points = ceil(total_length / 0.03) + length(beamline) * 10;
    
    % 使用结构体数组而非cell数组，提高访问速度
    all_slice_results(num_slice) = struct('s', [], 'sigx', [], 'sigpx', [], ...
                                          'sigy', [], 'sigpy', [], 'E', [], ...
                                          'norm_emitx', 0, 'norm_emity', 0, ...
                                          'particle_count', 0, 'valid', false);
    
    all_centroid_results(num_slice) = struct('s', [], 'x', [], 'px', [], ...
                                              'y', [], 'py', [], ...
                                              'particle_count', 0, 'valid', false);
    
    s_common = [];
    
    % **预处理：筛选有效切片**
    valid_slices = slice_particle_counts > 0 & ...
                   emitx_geo_slices > 1e-15 & ...
                   emity_geo_slices > 1e-15 & ...
                   betax_slices > 1e-10 & ...
                   betay_slices > 1e-10;
    
    valid_slice_indices = find(valid_slices);
    
    % 逐切片计算（仅处理有效切片）
    for idx = 1:length(valid_slice_indices)
        slice_idx = valid_slice_indices(idx);
        
        betax_slice = betax_slices(slice_idx);
        betay_slice = betay_slices(slice_idx);
        alphax_slice = alphax_slices(slice_idx);
        alphay_slice = alphay_slices(slice_idx);
        emitx_geo_slice = emitx_geo_slices(slice_idx);
        emity_geo_slice = emity_geo_slices(slice_idx);
        slice_particles = slice_particle_counts(slice_idx);
        
        x0_slice = x_centroid_slices(slice_idx);
        px0_slice = px_centroid_slices(slice_idx);
        y0_slice = y_centroid_slices(slice_idx);
        py0_slice = py_centroid_slices(slice_idx);
        
        norm_emitx_slice = emitx_geo_slice * gamma0;
        norm_emity_slice = emity_geo_slice * gamma0;
        
        % 初始包络参数
        sig0x_slice = sqrt(betax_slice * emitx_geo_slice);
        sigp0x_slice = -alphax_slice * sqrt(emitx_geo_slice / betax_slice);
        sig0y_slice = sqrt(betay_slice * emity_geo_slice);
        sigp0y_slice = -alphay_slice * sqrt(emity_geo_slice / betay_slice);
        
        % 安全性检查（合并到一个条件）
        if ~all(isfinite([sig0x_slice, sigp0x_slice, sig0y_slice, sigp0y_slice])) || ...
           sig0x_slice <= 0 || sig0y_slice <= 0
            continue;
        end
        
        Q_slice = Q_total * slice_particles / total_particles;
        N_slice_real = abs(Q_slice) / e;
        
        current_y_slice = [sig0x_slice, sigp0x_slice, sig0y_slice, sigp0y_slice, E0];
        current_s_slice = 0;
        current_centroid_slice = [x0_slice, px0_slice, y0_slice, py0_slice];
        
        % 预分配数组
        s_array_slice = zeros(estimated_points, 1);
        sigx_array_slice = zeros(estimated_points, 1);
        sigpx_array_slice = zeros(estimated_points, 1);
        sigy_array_slice = zeros(estimated_points, 1);
        sigpy_array_slice = zeros(estimated_points, 1);
        E_array_slice = zeros(estimated_points, 1);
        
        x_cent_array_slice = zeros(estimated_points, 1);
        px_cent_array_slice = zeros(estimated_points, 1);
        y_cent_array_slice = zeros(estimated_points, 1);
        py_cent_array_slice = zeros(estimated_points, 1);
        
        current_idx = 0;
        ode_failed = false;
        
        % **预计算元素类型（避免重复strcmp）**
        num_elements = length(beamline);
        element_types = zeros(num_elements, 1); % 1=drift, 2=quad, 3=cavity, 4=其他
        for i = 1:num_elements
            switch lower(beamline(i).type)
                case 'drift'
                    element_types(i) = 1;
                case 'quadrupole'
                    element_types(i) = 2;
                case 'cavity'
                    element_types(i) = 3;
                case {'sextupole', 'marker'}
                    element_types(i) = 0; % skip
                otherwise
                    element_types(i) = 4;
            end
        end
        
        % 逐元素求解
        for i = 1:num_elements
            element = beamline(i);
            
            if element.length == 0 || element_types(i) == 0
                continue;
            end
            
            s_start = current_s_slice;
            s_end = current_s_slice + element.length;
            
            % 合并有效性检查
            if any(~isfinite(current_y_slice)) || current_y_slice(1) <= 0 || current_y_slice(3) <= 0
                ode_failed = true;
                break;
            end
            
            % **根据预计算的类型分支（减少strcmp调用）**
            switch element_types(i)
                case 1 % drift
                    [z_seg, y_seg] = ode45(@(z,y) envelopeODE_sliced(z, y, 0, 0, 0, N_slice_real, z_val, rc, sigz, norm_emitx_slice, norm_emity_slice, me_c2), ...
                                           [s_start, s_end], current_y_slice, options);
                    
                case 2 % quadrupole
                    if quad_gradients.isKey(element.name)
                        gradient = quad_gradients(element.name);
                    else
                        gamma_current = current_y_slice(5) / me_c2;
                        gradient = element.k1 * gamma_current * me_c2 * 1.602e-19 / (e * c);
                        quad_gradients(element.name) = gradient;
                    end
                    
                    if calc_k1_flag && slice_idx == valid_slice_indices(1)
                        gamma_entrance = current_y_slice(5) / me_c2;
                        k1_entrance = (e * gradient * c) / (gamma_entrance * me_c2 * 1.602e-19);
                        
                        [z_seg, y_seg] = ode45(@(z,y) envelopeODE_quad_sliced(z, y, gradient, N_slice_real, z_val, rc, sigz, norm_emitx_slice, norm_emity_slice, me_c2, e, c), ...
                                               [s_start, s_end], current_y_slice, options);
                        
                        gamma_exit = y_seg(end,5) / me_c2;
                        k1_exit = (e * gradient * c) / (gamma_exit * me_c2 * 1.602e-19);
                        k1_values = [k1_values, (k1_entrance + k1_exit) / 2];
                    else
                        [z_seg, y_seg] = ode45(@(z,y) envelopeODE_quad_sliced(z, y, gradient, N_slice_real, z_val, rc, sigz, norm_emitx_slice, norm_emity_slice, me_c2, e, c), ...
                                               [s_start, s_end], current_y_slice, options);
                    end
                    
                case 3 % cavity
                    if isfield(element, 'num_cells') && isfield(element, 'frequency')
                        element.length = element.num_cells * 0.5 * 3e8 / (element.frequency * 1e6);
                    end
                    
                    A = 0;
                    if isfield(element, 'gradient')
                        A = element.gradient * 1e6;
                    end

                    phi = 0;
                    if isfield(element, 'phase')
                        phi = element.phase * pi/180;
                    end
                    
                    energy_gain = A * element.length * cos(phi);
                    
                    % ========== Entrance kick（包络 + 中心轨迹）==========
                    if energy_gain > 0
                        dgamma_over_gamma_entrance = A*cos(phi) / current_y_slice(5);
                        M21_entrance = -dgamma_over_gamma_entrance / 2;
                        
                        % 包络 entrance kick
                        current_y_slice(2) = current_y_slice(2) + M21_entrance * current_y_slice(1);
                        current_y_slice(4) = current_y_slice(4) + M21_entrance * current_y_slice(3);
                        
                        % 中心轨迹 entrance kick
                        current_centroid_slice(2) = current_centroid_slice(2) + M21_entrance * current_centroid_slice(1);
                        current_centroid_slice(4) = current_centroid_slice(4) + M21_entrance * current_centroid_slice(3);
                    end
                    
                    % ========== ODE 求解 ==========
                    [z_seg, y_seg] = ode45(@(z,y) envelopeODE_sliced(z, y, 0, A, phi, N_slice_real, z_val, rc, sigz, norm_emitx_slice, norm_emity_slice, me_c2), ...
                                           [s_start, s_end], current_y_slice, options);
                    
                    if isempty(z_seg) || any(any(~isfinite(y_seg))) || y_seg(end,1) <= 0 || y_seg(end,3) <= 0
                        ode_failed = true;
                        break;
                    end
                    
                    % ========== 中心运动计算（cavity内部）==========
                    n_points = length(z_seg);
                    cent_seg = zeros(n_points, 4);
                    cent_seg(1, :) = current_centroid_slice;  % 使用已经施加entrance kick的状态
                    
                    sigx_vec = y_seg(:, 1);
                    sigy_vec = y_seg(:, 3);
                    E_vec = y_seg(:, 5);
                    
                    for j = 2:n_points
                        ds_cent = z_seg(j) - z_seg(j-1);
                        dcdt = centroidODE_cavity_KV(cent_seg(j-1,:)', sigx_vec(j-1), sigy_vec(j-1), E_vec(j-1), A, phi, N_slice_real, z_val, me_c2);
                        cent_seg(j, :) = cent_seg(j-1, :) + dcdt' * ds_cent;
                    end
                    
                    % ========== Exit kick（包络 + 中心轨迹）==========
                    current_y_slice = y_seg(end,:)';
                    
                    if energy_gain > 0
                        dgamma_over_gamma_exit = A*cos(phi) / current_y_slice(5);
                        M21_exit = dgamma_over_gamma_exit / 2;
                        
                        % 包络 exit kick
                        current_y_slice(2) = current_y_slice(2) + M21_exit * current_y_slice(1);
                        current_y_slice(4) = current_y_slice(4) + M21_exit * current_y_slice(3);
                        
                        % 中心轨迹 exit kick
                        cent_seg(end, 2) = cent_seg(end, 2) + M21_exit * cent_seg(end, 1);
                        cent_seg(end, 4) = cent_seg(end, 4) + M21_exit * cent_seg(end, 3);
                    end
                    
                    current_centroid_slice = cent_seg(end, :)';
                    
                case 4 % 其他元素作为drift处理
                    [z_seg, y_seg] = ode45(@(z,y) envelopeODE_sliced(z, y, 0, 0, 0, N_slice_real, z_val, rc, sigz, norm_emitx_slice, norm_emity_slice, me_c2), ...
                                           [s_start, s_end], current_y_slice, options);
            end
            
            if isempty(z_seg) || any(any(~isfinite(y_seg))) || y_seg(end,1) <= 0 || y_seg(end,3) <= 0
                ode_failed = true;
                break;
            end
            
            % **中心运动计算（非cavity元素）**
            if element_types(i) ~= 3
                n_points = length(z_seg);
                cent_seg = zeros(n_points, 4);
                cent_seg(1, :) = current_centroid_slice;
                
                % 预提取包络参数（避免重复索引）
                sigx_vec = y_seg(:, 1);
                sigy_vec = y_seg(:, 3);
                E_vec = y_seg(:, 5);
                
                for j = 2:n_points
                    ds_cent = z_seg(j) - z_seg(j-1);
                    
                    % 根据元素类型计算导数
                    switch element_types(i)
                        case {1, 4} % drift
                            dcdt = centroidODE_drift_KV(cent_seg(j-1,:)', sigx_vec(j-1), sigy_vec(j-1), E_vec(j-1), N_slice_real, z_val, me_c2);
                        case 2 % quad
                            k1_cur = (e * gradient * c) / ((E_vec(j-1)/me_c2) * me_c2 * 1.602e-19);
                            dcdt = centroidODE_quad_KV(cent_seg(j-1,:)', sigx_vec(j-1), sigy_vec(j-1), k1_cur, E_vec(j-1), N_slice_real, z_val, me_c2);
                    end
                    
                    cent_seg(j, :) = cent_seg(j-1, :) + dcdt' * ds_cent;
                end
                
                current_centroid_slice = cent_seg(end, :)';
            end
            
            % **存储数据（优化索引操作）**
            n_new = length(z_seg);
            if current_idx == 0
                idx_range = 1:n_new;
            else
                idx_range = current_idx+1:current_idx+n_new-1;
                z_seg = z_seg(2:end);
                y_seg = y_seg(2:end, :);
                cent_seg = cent_seg(2:end, :);
            end
            
            s_array_slice(idx_range) = z_seg;
            sigx_array_slice(idx_range) = y_seg(:,1);
            sigpx_array_slice(idx_range) = y_seg(:,2);
            sigy_array_slice(idx_range) = y_seg(:,3);
            sigpy_array_slice(idx_range) = y_seg(:,4);
            E_array_slice(idx_range) = y_seg(:,5);
            
            x_cent_array_slice(idx_range) = cent_seg(:,1);
            px_cent_array_slice(idx_range) = cent_seg(:,2);
            y_cent_array_slice(idx_range) = cent_seg(:,3);
            py_cent_array_slice(idx_range) = cent_seg(:,4);
            
            current_idx = idx_range(end);
            
            if element_types(i) ~= 3 % 非cavity元素
                current_y_slice = y_seg(end,:)';
            end
            
            current_s_slice = s_end;
        end
        
        if ode_failed
            continue;
        end
        
        % 裁剪到实际长度
        s_array_slice = s_array_slice(1:current_idx);
        sigx_array_slice = sigx_array_slice(1:current_idx);
        sigpx_array_slice = sigpx_array_slice(1:current_idx);
        sigy_array_slice = sigy_array_slice(1:current_idx);
        sigpy_array_slice = sigpy_array_slice(1:current_idx);
        E_array_slice = E_array_slice(1:current_idx);
        
        x_cent_array_slice = x_cent_array_slice(1:current_idx);
        px_cent_array_slice = px_cent_array_slice(1:current_idx);
        y_cent_array_slice = y_cent_array_slice(1:current_idx);
        py_cent_array_slice = py_cent_array_slice(1:current_idx);
        
        % 存储结果
        all_slice_results(slice_idx).s = s_array_slice;
        all_slice_results(slice_idx).sigx = sigx_array_slice;
        all_slice_results(slice_idx).sigpx = sigpx_array_slice;
        all_slice_results(slice_idx).sigy = sigy_array_slice;
        all_slice_results(slice_idx).sigpy = sigpy_array_slice;
        all_slice_results(slice_idx).E = E_array_slice;
        all_slice_results(slice_idx).norm_emitx = norm_emitx_slice;
        all_slice_results(slice_idx).norm_emity = norm_emity_slice;
        all_slice_results(slice_idx).particle_count = slice_particles;
        all_slice_results(slice_idx).valid = true;
        
        all_centroid_results(slice_idx).s = s_array_slice;
        all_centroid_results(slice_idx).x = x_cent_array_slice;
        all_centroid_results(slice_idx).px = px_cent_array_slice;
        all_centroid_results(slice_idx).y = y_cent_array_slice;
        all_centroid_results(slice_idx).py = py_cent_array_slice;
        all_centroid_results(slice_idx).particle_count = slice_particles;
        all_centroid_results(slice_idx).valid = true;
        
        if isempty(s_common)
            s_common = s_array_slice;
        end
    end
    
    % **降采样（向量化）**
    num_points_original = length(s_common);
    downsample_factor = 1;
    num_points_reduced = ceil(num_points_original / downsample_factor);
    sample_indices = round(linspace(1, num_points_original, num_points_reduced));
    s_reduced = s_common(sample_indices);
    
    % **预分配最终结果数组**
    sigx_weighted_rms = zeros(num_points_reduced, 1);
    sigy_weighted_rms = zeros(num_points_reduced, 1);
    sigpx_weighted_rms = zeros(num_points_reduced, 1);
    sigpy_weighted_rms = zeros(num_points_reduced, 1);
    emitx_proj_array = zeros(num_points_reduced, 1);
    emity_proj_array = zeros(num_points_reduced, 1);
    x_centroid_weighted = zeros(num_points_reduced, 1);
    y_centroid_weighted = zeros(num_points_reduced, 1);
    px_centroid_weighted_array = zeros(num_points_reduced, 1);
    py_centroid_weighted_array = zeros(num_points_reduced, 1);
    gamma_weighted_array = zeros(num_points_reduced, 1);  % 新增gamma数组
    
    N_total = sum(slice_particle_counts);
    
    % **预插值所有切片数据（如果需要）**
    for slice_idx = valid_slice_indices'
        if all_slice_results(slice_idx).valid && ...
           length(all_slice_results(slice_idx).s) ~= num_points_original
            
            % 插值到统一网格
            all_slice_results(slice_idx).sigx_interp = interp1(all_slice_results(slice_idx).s, ...
                all_slice_results(slice_idx).sigx, s_common, 'linear', 'extrap');
            all_slice_results(slice_idx).sigpx_interp = interp1(all_slice_results(slice_idx).s, ...
                all_slice_results(slice_idx).sigpx, s_common, 'linear', 'extrap');
            all_slice_results(slice_idx).sigy_interp = interp1(all_slice_results(slice_idx).s, ...
                all_slice_results(slice_idx).sigy, s_common, 'linear', 'extrap');
            all_slice_results(slice_idx).sigpy_interp = interp1(all_slice_results(slice_idx).s, ...
                all_slice_results(slice_idx).sigpy, s_common, 'linear', 'extrap');
            all_slice_results(slice_idx).E_interp = interp1(all_slice_results(slice_idx).s, ...
                all_slice_results(slice_idx).E, s_common, 'linear', 'extrap');
            
            all_centroid_results(slice_idx).x_interp = interp1(all_centroid_results(slice_idx).s, ...
                all_centroid_results(slice_idx).x, s_common, 'linear', 'extrap');
            all_centroid_results(slice_idx).px_interp = interp1(all_centroid_results(slice_idx).s, ...
                all_centroid_results(slice_idx).px, s_common, 'linear', 'extrap');
            all_centroid_results(slice_idx).y_interp = interp1(all_centroid_results(slice_idx).s, ...
                all_centroid_results(slice_idx).y, s_common, 'linear', 'extrap');
            all_centroid_results(slice_idx).py_interp = interp1(all_centroid_results(slice_idx).s, ...
                all_centroid_results(slice_idx).py, s_common, 'linear', 'extrap');
        else
            % 无需插值，直接使用原始数据
            all_slice_results(slice_idx).sigx_interp = all_slice_results(slice_idx).sigx;
            all_slice_results(slice_idx).sigpx_interp = all_slice_results(slice_idx).sigpx;
            all_slice_results(slice_idx).sigy_interp = all_slice_results(slice_idx).sigy;
            all_slice_results(slice_idx).sigpy_interp = all_slice_results(slice_idx).sigpy;
            all_slice_results(slice_idx).E_interp = all_slice_results(slice_idx).E;
            
            all_centroid_results(slice_idx).x_interp = all_centroid_results(slice_idx).x;
            all_centroid_results(slice_idx).px_interp = all_centroid_results(slice_idx).px;
            all_centroid_results(slice_idx).y_interp = all_centroid_results(slice_idx).y;
            all_centroid_results(slice_idx).py_interp = all_centroid_results(slice_idx).py;
        end
    end
    
    % **向量化加权计算**
    for i = 1:num_points_reduced
        original_idx = sample_indices(i);
        
        % 收集所有有效切片的数据（向量化）
        valid_mask = [all_centroid_results.valid];
        valid_indices = find(valid_mask);
        n_valid = length(valid_indices);
        
        if n_valid == 0
            continue;
        end
        
        % 预分配临时数组
        weights = zeros(n_valid, 1);
        x_vals = zeros(n_valid, 1);
        y_vals = zeros(n_valid, 1);
        px_vals = zeros(n_valid, 1);
        py_vals = zeros(n_valid, 1);
        E_vals = zeros(n_valid, 1);  % 新增能量数组
        
        for j = 1:n_valid
            idx = valid_indices(j);
            weights(j) = all_centroid_results(idx).particle_count;
            x_vals(j) = all_centroid_results(idx).x_interp(original_idx);
            y_vals(j) = all_centroid_results(idx).y_interp(original_idx);
            px_vals(j) = all_centroid_results(idx).px_interp(original_idx);
            py_vals(j) = all_centroid_results(idx).py_interp(original_idx);
            E_vals(j) = all_slice_results(idx).E_interp(original_idx);  % 获取能量
        end
        
        % 向量化计算加权中心（包括角度中心和能量）
        total_weight = sum(weights);
        if total_weight > 0
            x_centroid_weighted(i) = sum(x_vals .* weights) / total_weight;
            y_centroid_weighted(i) = sum(y_vals .* weights) / total_weight;
            px_centroid_weighted = sum(px_vals .* weights) / total_weight;
            py_centroid_weighted = sum(py_vals .* weights) / total_weight;
            E_weighted = sum(E_vals .* weights) / total_weight;  % 加权平均能量
            
            % 计算加权平均gamma
            gamma_weighted_array(i) = E_weighted / me_c2;
            
            % 保存到数组
            px_centroid_weighted_array(i) = px_centroid_weighted;
            py_centroid_weighted_array(i) = py_centroid_weighted;
        end
        
        % 第二次循环：计算RMS和投影发射度
        weighted_sum_x = 0;
        weighted_sum_y = 0;
        weighted_sum_px = 0;
        weighted_sum_py = 0;
        
        x2_sum = 0;
        xxp_sum = 0;
        xp2_sum = 0;
        y2_sum = 0;
        yyp_sum = 0;
        yp2_sum = 0;
        
        for j = 1:n_valid
            idx = valid_indices(j);
            
            sigx_val = all_slice_results(idx).sigx_interp(original_idx);
            sigpx_val = all_slice_results(idx).sigpx_interp(original_idx);
            sigy_val = all_slice_results(idx).sigy_interp(original_idx);
            sigpy_val = all_slice_results(idx).sigpy_interp(original_idx);
            E_val = all_slice_results(idx).E_interp(original_idx);
            
            if ~isfinite(sigx_val) || ~isfinite(sigy_val) || sigx_val <= 0 || sigy_val <= 0
                continue;
            end
            
            weight = weights(j);
            
            x_cent_val = all_centroid_results(idx).x_interp(original_idx);
            px_cent_val = all_centroid_results(idx).px_interp(original_idx);
            y_cent_val = all_centroid_results(idx).y_interp(original_idx);
            py_cent_val = all_centroid_results(idx).py_interp(original_idx);
            
            % ========== 修正：计算相对于总体加权中心的偏移 ==========
            dx = x_cent_val - x_centroid_weighted(i);
            dy = y_cent_val - y_centroid_weighted(i);
            dpx = px_cent_val - px_centroid_weighted;
            dpy = py_cent_val - py_centroid_weighted;
            
            sigx_total_sq = sigx_val^2;
            sigy_total_sq = sigy_val^2;
            
            weighted_sum_x = weighted_sum_x + sigx_total_sq * weight;
            weighted_sum_y = weighted_sum_y + sigy_total_sq * weight;
            weighted_sum_px = weighted_sum_px + sigpx_val^2 * weight;
            weighted_sum_py = weighted_sum_py + sigpy_val^2 * weight;
            
            gamma_val = E_val / me_c2;
            emitx_geo = all_slice_results(idx).norm_emitx / gamma_val;
            emity_geo = all_slice_results(idx).norm_emity / gamma_val;
            
            betax_slice = sigx_val^2 / emitx_geo;
            alphax_slice = -sigx_val * sigpx_val / emitx_geo;
            gammax_slice = (1 + alphax_slice^2) / betax_slice;
            
            betay_slice = sigy_val^2 / emity_geo;
            alphay_slice = -sigy_val * sigpy_val / emity_geo;
            gammay_slice = (1 + alphay_slice^2) / betay_slice;
            
            % 投影发射度计算使用相对偏移（平行轴定理）
            x2_slice = betax_slice * emitx_geo + dx^2;
            xxp_slice = -alphax_slice * emitx_geo + dx * dpx;
            xp2_slice = gammax_slice * emitx_geo + dpx^2;
            
            y2_slice = betay_slice * emity_geo + dy^2;
            yyp_slice = -alphay_slice * emity_geo + dy * dpy;
            yp2_slice = gammay_slice * emity_geo + dpy^2;
            
            x2_sum = x2_sum + weight * x2_slice;
            xxp_sum = xxp_sum + weight * xxp_slice;
            xp2_sum = xp2_sum + weight * xp2_slice;
            
            y2_sum = y2_sum + weight * y2_slice;
            yyp_sum = yyp_sum + weight * yyp_slice;
            yp2_sum = yp2_sum + weight * yp2_slice;
        end
        
        if total_weight > 0
            sigx_weighted_rms(i) = sqrt(weighted_sum_x / total_weight);
            sigy_weighted_rms(i) = sqrt(weighted_sum_y / total_weight);
            sigpx_weighted_rms(i) = sqrt(weighted_sum_px / total_weight);
            sigpy_weighted_rms(i) = sqrt(weighted_sum_py / total_weight);
            
            x2_proj = x2_sum / N_total;
            xxp_proj = xxp_sum / N_total;
            xp2_proj = xp2_sum / N_total;
            
            y2_proj = y2_sum / N_total;
            yyp_proj = yyp_sum / N_total;
            yp2_proj = yp2_sum / N_total;
            
            emitx_proj_array(i) = sqrt(max(x2_proj * xp2_proj - xxp_proj^2, 0));
            emity_proj_array(i) = sqrt(max(y2_proj * yp2_proj - yyp_proj^2, 0));
        else
            sigx_weighted_rms(i) = NaN;
            sigy_weighted_rms(i) = NaN;
            sigpx_weighted_rms(i) = NaN;
            sigpy_weighted_rms(i) = NaN;
            emitx_proj_array(i) = NaN;
            emity_proj_array(i) = NaN;
        end
    end
    
    % **整合所有输出到一个结构体**
    results = struct();
    results.s_array = s_reduced;
    results.sigx_array = sigx_weighted_rms;
    results.sigy_array = sigy_weighted_rms;
    results.sigpx_array = sigpx_weighted_rms;
    results.sigpy_array = sigpy_weighted_rms;
    results.x_centroid = x_centroid_weighted;
    results.y_centroid = y_centroid_weighted;
    results.px_centroid = px_centroid_weighted_array;
    results.py_centroid = py_centroid_weighted_array;
    results.gamma_array = gamma_weighted_array;  % 新增gamma输出
    results.emitx_proj_array = emitx_proj_array;
    results.emity_proj_array = emity_proj_array;
    
    if calc_k1_flag && ~exist('k1_values', 'var')
        k1_values = [];
    end
    results.k1_values = k1_values;
end

end

% ========== ODE函数 ==========
function dydt = envelopeODE_sliced(s, y, k1, A, phi, N_slice_real, z_val, rc, sigz, norm_emitx, norm_emity, me_c2)
sigx = max(y(1), 1e-12);
sigy = max(y(3), 1e-12);
E = max(y(5), me_c2);
gamma = E / me_c2;

if ~isfinite(sigx) || ~isfinite(sigy) || sigx <= 0 || sigy <= 0
    dydt = [0; 0; 0; 0; 0];
    return;
end

geo_emitx = norm_emitx / gamma;
geo_emity = norm_emity / gamma;

emittance_x = geo_emitx^2 / (sigx^3);
emittance_y = geo_emity^2 / (sigy^3);

if N_slice_real > 0 && z_val > 0
    N_density = N_slice_real / z_val;
    lambda = 1 / (sigz * gamma);
    A0 = 1/pi*sqrt(pi)/gamma/gamma*N_density * rc * lambda/(1/(2*sqrt(pi)*sigz));
    
    if A0 > 0
        sc_x = spaceChargeXFast(sigx, sigy, A0, 1/sigx);
        sc_y = spaceChargeYFast(sigx, sigy, A0, 1/sigy);
    else
        sc_x = 0; sc_y = 0;
    end
else
    sc_x = 0;
    sc_y = 0;
end

if nargin >= 5 && A > 0
    dE_ds = A * cos(phi);
    dgamma_over_gamma = dE_ds / E;
    if abs(cos(phi)) > 1e-6
        second_order_coeff = -1 / (8* cos(phi)^2) * dgamma_over_gamma^2;
    else
        second_order_coeff = 0;
    end
else
    dE_ds = 0;
    dgamma_over_gamma = 0;
    second_order_coeff = 0;
end

dsigx_ds = y(2);
dsigpx_ds = emittance_x + sc_x - k1 * sigx - dgamma_over_gamma * y(2) + second_order_coeff * sigx;
dsigy_ds = y(4);
dsigpy_ds = emittance_y + sc_y + k1 * sigy - dgamma_over_gamma * y(4) + second_order_coeff * sigy;

dsigpx_ds = max(min(dsigpx_ds, 1e10), -1e10);
dsigpy_ds = max(min(dsigpy_ds, 1e10), -1e10);

dydt = [dsigx_ds; dsigpx_ds; dsigy_ds; dsigpy_ds; dE_ds];
end

function dydt = envelopeODE_quad_sliced(s, y, gradient, N_slice_real, z_val, rc, sigz, norm_emitx, norm_emity, me_c2, e, c)
sigx = max(y(1), 1e-12);
sigy = max(y(3), 1e-12);
E = max(y(5), me_c2);
gamma = E / me_c2;

k1 = (e * gradient * c) / (gamma * me_c2 * 1.602e-19);

dydt = envelopeODE_sliced(s, y, k1, 0, 0, N_slice_real, z_val, rc, sigz, norm_emitx, norm_emity, me_c2);
dydt(5) = 0;
end

% ========== 修正后的中心运动ODE函数 ==========
function dcdt = centroidODE_drift_KV(cent, sigx, sigy, E, N_slice_real, z_val, me_c2)
    x = cent(1);
    xp = cent(2);
    y = cent(3);
    yp = cent(4);
    
    gamma = E / me_c2;
    beta_rel = sqrt(1 - 1/gamma^2);
    
    if N_slice_real > 0 && sigx > 0 && sigy > 0 && z_val > 0
        r0 = 2.8179e-15;
        % 线电荷密度
        lambda = N_slice_real / z_val;
        % KV分布空间电荷力系数
        K_coeff = 2 * lambda * r0 / (gamma^3 * beta_rel^2);
        
        % x和y方向分别计算（椭圆束）
        Ksc_x = K_coeff / (sigx * (sigx + sigy));
        Ksc_y = K_coeff / (sigy * (sigx + sigy));
    else
        Ksc_x = 0;
        Ksc_y = 0;
    end
    
    dcdt = [xp; Ksc_x * x; yp; Ksc_y * y];
end

function dcdt = centroidODE_quad_KV(cent, sigx, sigy, k1, E, N_slice_real, z_val, me_c2)
    x = cent(1);
    xp = cent(2);
    y = cent(3);
    yp = cent(4);
    
    % 计算空间电荷力
    gamma = E / me_c2;
    beta_rel = sqrt(1 - 1/gamma^2);
    
    if N_slice_real > 0 && sigx > 0 && sigy > 0 && z_val > 0
        r0 = 2.8179e-15;
        lambda = N_slice_real / z_val;
        K_coeff = 2 * lambda * r0 / (gamma^3 * beta_rel^2);
        
        Ksc_x = K_coeff / (sigx * (sigx + sigy));
        Ksc_y = K_coeff / (sigy * (sigx + sigy));
    else
        Ksc_x = 0;
        Ksc_y = 0;
    end
    
    % 四极铁力 + 空间电荷力
    dcdt = [xp; 
            -k1 * x + Ksc_x * x; 
            yp; 
            k1 * y + Ksc_y * y];
end

function dcdt = centroidODE_cavity_KV(cent, sigx, sigy, E, A, phi, N_slice_real, z_val, me_c2)
    x = cent(1);
    xp = cent(2);
    y = cent(3);
    yp = cent(4);
    
    gamma = E / me_c2;
    beta_rel = sqrt(1 - 1/gamma^2);
    
    if N_slice_real > 0 && sigx > 0 && sigy > 0 && z_val > 0
        r0 = 2.8179e-15;
        lambda = N_slice_real / z_val;
        K_coeff = 2 * lambda * r0 / (gamma^3 * beta_rel^2);
        
        Ksc_x = K_coeff / (sigx * (sigx + sigy));
        Ksc_y = K_coeff / (sigy * (sigx + sigy));
    else
        Ksc_x = 0;
        Ksc_y = 0;
    end
    
    dgamma_over_gamma = A * cos(phi) / E;
    
    % 二阶色散系数
    if abs(cos(phi)) > 1e-6
        second_order_coeff = -dgamma_over_gamma^2 / (8 * cos(phi)^2);
    else
        second_order_coeff = 0;
    end
    
    dcdt = [xp; 
            -dgamma_over_gamma * xp + second_order_coeff * x + Ksc_x * x;
            yp;
            -dgamma_over_gamma * yp + second_order_coeff * y + Ksc_y * y];
end

function val = spaceChargeXFast(sigx, sigy, A0, sigx_inv)
    k = min(max(sigy * sigx_inv, 1e-6), 1e6);
    k_power = min(max(k^(-1/5.5), 1e-6), 1e6);
    f = 1 - 1/sqrt(1 + k_power);
    exp_arg = max(-1/(8*k^(1/5.5)), -50);
    g = (1/4.25) / (1+k) / (1 - max(exp(exp_arg), 1e-15));
    val = min(max(f * g * A0 * sigx_inv, -1e6), 1e6);
end

function val = spaceChargeYFast(sigx, sigy, A0, sigy_inv)
    k = min(max(sigx * sigy_inv, 1e-6), 1e6);
    k_power = min(max(k^(-1/5.5), 1e-6), 1e6);
    f = 1 - 1/sqrt(1 + k_power);
    exp_arg = max(-1/(8*k^(1/5.5)), -50);
    g = (1/4.25) / (1+k) / (1 - max(exp(exp_arg), 1e-15));
    val = min(max(f * g * A0 * sigy_inv, -1e6), 1e6);
end