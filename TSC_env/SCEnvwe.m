function [results] = SCEnvwe(beamline, init_params, beamfile1)
% SCEnvwe - 统一的空间电荷包络计算入口函数
%
% 修改重点：
% 1. 中心计算：保持使用高阶传输矩阵 (R矩阵 + Tijk二阶项)。
% 2. Q铁包络：k1 直接修正为 k1 / (1 + delta_slice)。
% 3. Cavity：保持原样。

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

me_c2 = 0.511e6;

if nargin < 3
    % ================= 单一包络模式 (保持不变) =================
    init_params_single = [betax, betay, alphax, alphay, emitx, emity, sigz1, Q_total, E0];
    [s_array, sigx_array, sigy_array, k1_values] = SCEqga(beamline, init_params_single, 1);
    
    gamma_array = zeros(size(s_array));
    current_E = E0;
    gamma_array(1) = current_E / me_c2;
    
    current_s = 0;
    point_idx = 1;
    
    for i = 1:length(beamline)
        element = beamline(i);
        if strcmpi(element.type, 'cavity')
            if isfield(element, 'num_cells') && isfield(element, 'frequency')
                element.length = element.num_cells * 0.5 * 3e8 / (element.frequency * 1e6);
            end
            A = 0; if isfield(element, 'gradient'), A = element.gradient * 1e6; end
            phi = 0; if isfield(element, 'phase'), phi = element.phase * pi/180; end
            energy_gain = A * element.length * cos(phi);
            current_E = current_E + energy_gain;
        end
        s_end = current_s + element.length;
        while point_idx <= length(s_array) && s_array(point_idx) <= s_end
            gamma_array(point_idx) = current_E / me_c2;
            point_idx = point_idx + 1;
        end
        current_s = s_end;
    end
    gamma_array(point_idx:end) = current_E / me_c2;
    
    results = struct();
    results.s_array = s_array;
    results.sigx_array = sigx_array;
    results.sigy_array = sigy_array;
    results.sigpx_array = []; results.sigpy_array = [];
    results.x_centroid = []; results.y_centroid = [];
    results.px_centroid = []; results.py_centroid = [];
    results.gamma_array = gamma_array;
    results.k1_values = k1_values;
    results.emitx_proj_array = []; results.emity_proj_array = [];
    return;
else
    % ================= 基于切片的多包络计算模式 =================
    beamfile = beamfile1;
    
    try
        slice_data_matrix = load(beamfile);
    catch
        error('无法读取文件: %s', beamfile);
    end
    
    % 解析切片数据
    total_particles = slice_data_matrix(1, 1);
    z_length = 1*slice_data_matrix(end, 1);
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
    % z_centroid_slices = slice_info(:, 14); 
    delta_centroid_slices = slice_info(:, 15); % 读取 delta
    
    % 读取二阶矩 <x*pz> 和 <px*pz> 用于中心计算
    if size(slice_info, 2) >= 17
        x_pz_slices = slice_info(:, 16);  
        px_pz_slices = slice_info(:, 17); 
    else
        x_pz_slices = x_centroid_slices .* delta_centroid_slices;
        px_pz_slices = px_centroid_slices .* delta_centroid_slices;
    end
    
    % 物理常数
    rc = 2.8179e-15;
    e = 1.602e-19;
    % c = 2.998e8; % 此处不需要c，除非用于梯度换算，但新逻辑直接用k1
    
    gamma0 = E0 / me_c2;
    z_val = z_length / num_slice;
    
    if calc_k1_flag, k1_values = []; end
    
    options = odeset('RelTol', 5e-6, 'AbsTol', 5e-8, 'MaxStep', 0.03, 'Stats', 'off');
    
    % 预筛选有效切片
    valid_slices = slice_particle_counts > 0 & ...
                   emitx_geo_slices > 1e-15 & ...
                   emity_geo_slices > 1e-15 & ...
                   betax_slices > 1e-10 & ...
                   betay_slices > 1e-10;
    
    valid_slice_indices = find(valid_slices);

    % === 1. 中心轨迹计算：保持高阶传输矩阵方法 ===
    [s_transport, T126_evol, T116_evol, T216_evol, T226_evol, ...
     T346_evol, T336_evol, T436_evol, T446_evol] = ...
        track_second_order_evolution(beamline, 0.001);
    
    num_elements = length(beamline);
    total_length = sum([beamline.length]);
    estimated_points = ceil(total_length / 0.02) + num_elements * 10;
    
    all_slice_envelope(num_slice) = struct('s', [], 'sigx', [], 'sigpx', [], ...
                                           'sigy', [], 'sigpy', [], 'E', [], ...
                                           'norm_emitx', 0, 'norm_emity', 0, ...
                                           'particle_count', 0, 'valid', false);
    
    all_slice_centroid(num_slice) = struct('s_exit', [], 'x_exit', [], 'px_exit', [], ...
                                           'y_exit', [], 'py_exit', [], ...
                                           'particle_count', 0, 'valid', false);
    
    s_envelope_common = []; 
    
    % 预计算元素类型
    element_types = zeros(num_elements, 1);
    for i = 1:num_elements
        switch lower(beamline(i).type)
            case 'drift', element_types(i) = 1;
            case 'quadrupole', element_types(i) = 2;
            case 'cavity', element_types(i) = 3;
            case {'sextupole', 'marker'}, element_types(i) = 0;
            otherwise, element_types(i) = 4;
        end
    end
    
    % ========== 逐切片计算 ==========
    for idx = 1:length(valid_slice_indices)
        slice_idx = valid_slice_indices(idx);
        
        % 提取切片参数
        betax_slice = betax_slices(slice_idx); alphax_slice = alphax_slices(slice_idx);
        betay_slice = betay_slices(slice_idx); alphay_slice = alphay_slices(slice_idx);
        emitx_geo_slice = emitx_geo_slices(slice_idx); emity_geo_slice = emity_geo_slices(slice_idx);
        slice_particles = slice_particle_counts(slice_idx);
        
        x0_slice = x_centroid_slices(slice_idx); px0_slice = px_centroid_slices(slice_idx);
        y0_slice = y_centroid_slices(slice_idx); py0_slice = py_centroid_slices(slice_idx);
        delta0_slice = delta_centroid_slices(slice_idx); % 获取 delta
        
        x_pz_slice = x_pz_slices(slice_idx); px_pz_slice = px_pz_slices(slice_idx);
        
        norm_emitx_slice = emitx_geo_slice * gamma0;
        norm_emity_slice = emity_geo_slice * gamma0;
        
        % 初始包络
        sig0x_slice = sqrt(betax_slice * emitx_geo_slice);
        sigp0x_slice = -alphax_slice * sqrt(emitx_geo_slice / betax_slice);
        sig0y_slice = sqrt(betay_slice * emity_geo_slice);
        sigp0y_slice = -alphay_slice * sqrt(emity_geo_slice / betay_slice);
        
        if ~all(isfinite([sig0x_slice, sigp0x_slice, sig0y_slice, sigp0y_slice])) || sig0x_slice <= 0 || sig0y_slice <= 0
            continue;
        end
        
        Q_slice = Q_total * slice_particles / total_particles;
        N_slice_real = abs(Q_slice) / e;
        
        % 状态向量：[sigx, sigpx, sigy, sigpy, Energy]
        % 注意：这里初始能量设为 E0，因为 delta 的影响已经体现在 k1 的修正中，
        % 或者体现在中心轨迹的 T 矩阵中。如果需要包络本身的 gamma 随 delta 变化，
        % 可以设为 E0 * (1 + delta0_slice)，但根据您的要求“k1直接根据输入值”，
        % 我们主要通过修改 k1 来体现色散对聚焦的影响。
        current_y_slice = [sig0x_slice, sigp0x_slice, sig0y_slice, sigp0y_slice, E0];
        current_s_slice = 0;
        
        % === 2. 中心轨迹计算 (保持不变) ===
        num_exit_points = length(s_transport);
        x_centroid_slice = zeros(num_exit_points, 1);
        px_centroid_slice = zeros(num_exit_points, 1);
        y_centroid_slice = zeros(num_exit_points, 1);
        py_centroid_slice = zeros(num_exit_points, 1);
        
        x_centroid_slice(1) = x0_slice; px_centroid_slice(1) = px0_slice;
        y_centroid_slice(1) = y0_slice; py_centroid_slice(1) = py0_slice;
        
        for i = 2:num_exit_points
            current_beamline_segment = beamline(1:i-1);
            R = eye(6);
            for j = 1:length(current_beamline_segment)
                R = getElementMatrix(current_beamline_segment(j)) * R;
            end
            
            x_centroid_slice(i) = R(1,1)*x0_slice + R(1,2)*px0_slice;
            px_centroid_slice(i) = R(2,1)*x0_slice + R(2,2)*px0_slice;
            y_centroid_slice(i) = R(3,3)*y0_slice + R(3,4)*py0_slice;
            py_centroid_slice(i) = R(4,3)*y0_slice + R(4,4)*py0_slice;
            
            % 二阶修正 (色散项)
            x_centroid_slice(i) = x_centroid_slice(i) + T116_evol(i) * x_pz_slice + T126_evol(i) * px_pz_slice;
            px_centroid_slice(i) = px_centroid_slice(i) + T216_evol(i) * x_pz_slice + T226_evol(i) * px_pz_slice;
        end
        
        % === 3. 包络计算 (ODE) ===
        s_envelope_slice = zeros(estimated_points, 1);
        sigx_envelope_slice = zeros(estimated_points, 1); sigpx_envelope_slice = zeros(estimated_points, 1);
        sigy_envelope_slice = zeros(estimated_points, 1); sigpy_envelope_slice = zeros(estimated_points, 1);
        E_envelope_slice = zeros(estimated_points, 1);
        
        current_envelope_idx = 0;
        ode_failed = false;
        
        for i = 1:num_elements
            element = beamline(i);
            if element.length == 0 || element_types(i) == 0, continue; end
            
            s_start = current_s_slice;
            s_end = current_s_slice + element.length;
            
            if any(~isfinite(current_y_slice)) || current_y_slice(1) <= 0 || current_y_slice(3) <= 0
                ode_failed = true; break;
            end
            
            switch element_types(i)
                case 1 % Drift
                    [z_seg, y_seg] = ode45(@(z,y) envelopeODE_sliced(z, y, 0, 0, 0, N_slice_real, z_val, rc, sigz, norm_emitx_slice, norm_emity_slice, me_c2), ...
                                           [s_start, s_end], current_y_slice, options);
                    
                case 2 % Quadrupole
                    % === 修改核心：k1 直接修正 ===
                    k1_nominal = element.k1;
                    % 公式：k_eff = k_nom / (1 + delta)
                    k1_effective = k1_nominal / (1 + delta0_slice);
                    
                    % 记录 k1 值 (仅用于输出显示)
                    if calc_k1_flag && slice_idx == valid_slice_indices(1)
                        k1_values = [k1_values, k1_nominal]; % 记录名义值或有效值均可，此处记录名义值
                    end
                    
                    % 使用通用的 envelopeODE_sliced，传入修正后的 k1_effective
                    % 注意：这里不再使用 envelopeODE_quad_sliced，因为那个函数内部会根据能量反算k1
                    [z_seg, y_seg] = ode45(@(z,y) envelopeODE_sliced(z, y, k1_effective, 0, 0, N_slice_real, z_val, rc, sigz, norm_emitx_slice, norm_emity_slice, me_c2), ...
                                           [s_start, s_end], current_y_slice, options);
                    
                case 3 % Cavity (保持原样，不做色散修正)
                    if isfield(element, 'num_cells') && isfield(element, 'frequency')
                        element.length = element.num_cells * 0.5 * 3e8 / (element.frequency * 1e6);
                    end
                    A = 0; if isfield(element, 'gradient'), A = element.gradient * 1e6; end
                    phi = 0; if isfield(element, 'phase'), phi = element.phase * pi/180; end
                    
                    energy_gain = A * element.length * cos(phi);
                    
                    % Entrance kick
                    if energy_gain > 0
                        dgamma_over_gamma_entrance = A*cos(phi) / current_y_slice(5);
                        M21_entrance = -dgamma_over_gamma_entrance / 2;
                        current_y_slice(2) = current_y_slice(2) + M21_entrance * current_y_slice(1);
                        current_y_slice(4) = current_y_slice(4) + M21_entrance * current_y_slice(3);
                    end
                    
                    % ODE 求解 (Cavity)
                    [z_seg, y_seg] = ode45(@(z,y) envelopeODE_sliced(z, y, 0, A, phi, N_slice_real, z_val, rc, sigz, norm_emitx_slice, norm_emity_slice, me_c2), ...
                                           [s_start, s_end], current_y_slice, options);
                    
                    if isempty(z_seg) || any(any(~isfinite(y_seg))) || y_seg(end,1) <= 0 || y_seg(end,3) <= 0
                        ode_failed = true; break;
                    end
                    
                    % Exit kick
                    current_y_slice = y_seg(end,:)';
                    if energy_gain > 0
                        dgamma_over_gamma_exit = A*cos(phi) / current_y_slice(5);
                        M21_exit = dgamma_over_gamma_exit / 2;
                        current_y_slice(2) = current_y_slice(2) + M21_exit * current_y_slice(1);
                        current_y_slice(4) = current_y_slice(4) + M21_exit * current_y_slice(3);
                    end
                    
                case 4 % 其他
                    [z_seg, y_seg] = ode45(@(z,y) envelopeODE_sliced(z, y, 0, 0, 0, N_slice_real, z_val, rc, sigz, norm_emitx_slice, norm_emity_slice, me_c2), ...
                                           [s_start, s_end], current_y_slice, options);
            end
            
            if isempty(z_seg) || any(any(~isfinite(y_seg))) || y_seg(end,1) <= 0 || y_seg(end,3) <= 0
                ode_failed = true; break;
            end
            
            if element_types(i) ~= 3
                current_y_slice = y_seg(end,:)';
            end
            
            % 存储数据
            n_new = length(z_seg);
            if current_envelope_idx == 0, idx_range = 1:n_new;
            else, idx_range = current_envelope_idx+1:current_envelope_idx+n_new-1; z_seg = z_seg(2:end); y_seg = y_seg(2:end, :); end
            
            s_envelope_slice(idx_range) = z_seg;
            sigx_envelope_slice(idx_range) = y_seg(:,1); sigpx_envelope_slice(idx_range) = y_seg(:,2);
            sigy_envelope_slice(idx_range) = y_seg(:,3); sigpy_envelope_slice(idx_range) = y_seg(:,4);
            E_envelope_slice(idx_range) = y_seg(:,5);
            
            current_envelope_idx = idx_range(end);
            current_s_slice = s_end;
        end
        
        if ode_failed, continue; end
        
        % 裁剪并保存
        all_slice_envelope(slice_idx).s = s_envelope_slice(1:current_envelope_idx);
        all_slice_envelope(slice_idx).sigx = sigx_envelope_slice(1:current_envelope_idx);
        all_slice_envelope(slice_idx).sigpx = sigpx_envelope_slice(1:current_envelope_idx);
        all_slice_envelope(slice_idx).sigy = sigy_envelope_slice(1:current_envelope_idx);
        all_slice_envelope(slice_idx).sigpy = sigpy_envelope_slice(1:current_envelope_idx);
        all_slice_envelope(slice_idx).E = E_envelope_slice(1:current_envelope_idx);
        all_slice_envelope(slice_idx).norm_emitx = norm_emitx_slice;
        all_slice_envelope(slice_idx).norm_emity = norm_emity_slice;
        all_slice_envelope(slice_idx).particle_count = slice_particles;
        all_slice_envelope(slice_idx).valid = true;
        
        all_slice_centroid(slice_idx).s_exit = s_transport;
        all_slice_centroid(slice_idx).x_exit = x_centroid_slice;
        all_slice_centroid(slice_idx).px_exit = px_centroid_slice;
        all_slice_centroid(slice_idx).y_exit = y_centroid_slice;
        all_slice_centroid(slice_idx).py_exit = py_centroid_slice;
        all_slice_centroid(slice_idx).particle_count = slice_particles;
        all_slice_centroid(slice_idx).valid = true;
        
        if isempty(s_envelope_common), s_envelope_common = s_envelope_slice(1:current_envelope_idx); end
    end
    
    % ================= 后处理：加权平均 (保持不变) =================
    if isempty(s_envelope_common), s_envelope_common = s_transport; end
    num_dense_points = length(s_envelope_common);
    
    sum_weights = zeros(num_dense_points, 1);
    sum_x_centroid = zeros(num_dense_points, 1); sum_y_centroid = zeros(num_dense_points, 1);
    sum_px_centroid = zeros(num_dense_points, 1); sum_py_centroid = zeros(num_dense_points, 1);
    sum_gamma = zeros(num_dense_points, 1);
    
    % 1. 计算加权中心
    for idx = 1:length(valid_slice_indices)
        slice_idx = valid_slice_indices(idx);
        weight = all_slice_centroid(slice_idx).particle_count;
        
        x_interp = interp1(all_slice_centroid(slice_idx).s_exit, all_slice_centroid(slice_idx).x_exit, s_envelope_common, 'linear', 'extrap');
        y_interp = interp1(all_slice_centroid(slice_idx).s_exit, all_slice_centroid(slice_idx).y_exit, s_envelope_common, 'linear', 'extrap');
        px_interp = interp1(all_slice_centroid(slice_idx).s_exit, all_slice_centroid(slice_idx).px_exit, s_envelope_common, 'linear', 'extrap');
        py_interp = interp1(all_slice_centroid(slice_idx).s_exit, all_slice_centroid(slice_idx).py_exit, s_envelope_common, 'linear', 'extrap');
        E_interp = interp1(all_slice_envelope(slice_idx).s, all_slice_envelope(slice_idx).E, s_envelope_common, 'linear', 'extrap');
        
        sum_weights = sum_weights + weight;
        sum_x_centroid = sum_x_centroid + x_interp * weight;
        sum_y_centroid = sum_y_centroid + y_interp * weight;
        sum_px_centroid = sum_px_centroid + px_interp * weight;
        sum_py_centroid = sum_py_centroid + py_interp * weight;
        sum_gamma = sum_gamma + (E_interp / me_c2) * weight;
    end
    
    mean_x_centroid = sum_x_centroid ./ sum_weights; mean_y_centroid = sum_y_centroid ./ sum_weights;
    mean_px_centroid = sum_px_centroid ./ sum_weights; mean_py_centroid = sum_py_centroid ./ sum_weights;
    mean_gamma = sum_gamma ./ sum_weights;
    
    % 2. 计算二阶矩
    sum_x2 = zeros(num_dense_points, 1); sum_xxp = zeros(num_dense_points, 1); sum_xp2 = zeros(num_dense_points, 1);
    sum_y2 = zeros(num_dense_points, 1); sum_yyp = zeros(num_dense_points, 1); sum_yp2 = zeros(num_dense_points, 1);
    
    for idx = 1:length(valid_slice_indices)
        slice_idx = valid_slice_indices(idx);
        weight = all_slice_centroid(slice_idx).particle_count;
        
        sigx_interp = interp1(all_slice_envelope(slice_idx).s, all_slice_envelope(slice_idx).sigx, s_envelope_common, 'linear', 'extrap');
        sigpx_interp = interp1(all_slice_envelope(slice_idx).s, all_slice_envelope(slice_idx).sigpx, s_envelope_common, 'linear', 'extrap');
        sigy_interp = interp1(all_slice_envelope(slice_idx).s, all_slice_envelope(slice_idx).sigy, s_envelope_common, 'linear', 'extrap');
        sigpy_interp = interp1(all_slice_envelope(slice_idx).s, all_slice_envelope(slice_idx).sigpy, s_envelope_common, 'linear', 'extrap');
        E_interp = interp1(all_slice_envelope(slice_idx).s, all_slice_envelope(slice_idx).E, s_envelope_common, 'linear', 'extrap');
        
        x_c = interp1(all_slice_centroid(slice_idx).s_exit, all_slice_centroid(slice_idx).x_exit, s_envelope_common, 'linear', 'extrap');
        px_c = interp1(all_slice_centroid(slice_idx).s_exit, all_slice_centroid(slice_idx).px_exit, s_envelope_common, 'linear', 'extrap');
        y_c = interp1(all_slice_centroid(slice_idx).s_exit, all_slice_centroid(slice_idx).y_exit, s_envelope_common, 'linear', 'extrap');
        py_c = interp1(all_slice_centroid(slice_idx).s_exit, all_slice_centroid(slice_idx).py_exit, s_envelope_common, 'linear', 'extrap');
        
        dx = x_c - mean_x_centroid; dpx = px_c - mean_px_centroid;
        dy = y_c - mean_y_centroid; dpy = py_c - mean_py_centroid;
        
        gamma_val = E_interp / me_c2;
        emitx_geo = all_slice_envelope(slice_idx).norm_emitx ./ gamma_val;
        emity_geo = all_slice_envelope(slice_idx).norm_emity ./ gamma_val;
        
        x2_slice = sigx_interp.^2 + dx.^2;
        y2_slice = sigy_interp.^2 + dy.^2;
        xxp_slice = (sigx_interp .* sigpx_interp) + dx .* dpx;
        yyp_slice = (sigy_interp .* sigpy_interp) + dy .* dpy;
        
        xp2_slice_intrinsic = (emitx_geo.^2 + (sigx_interp .* sigpx_interp).^2) ./ (sigx_interp.^2);
        yp2_slice_intrinsic = (emity_geo.^2 + (sigy_interp .* sigpy_interp).^2) ./ (sigy_interp.^2);
        
        xp2_slice = xp2_slice_intrinsic + dpx.^2;
        yp2_slice = yp2_slice_intrinsic + dpy.^2;
        
        sum_x2 = sum_x2 + x2_slice * weight; sum_xxp = sum_xxp + xxp_slice * weight; sum_xp2 = sum_xp2 + xp2_slice * weight;
        sum_y2 = sum_y2 + y2_slice * weight; sum_yyp = sum_yyp + yyp_slice * weight; sum_yp2 = sum_yp2 + yp2_slice * weight;
    end
    
    Sigma_x = sqrt(sum_x2 ./ sum_weights); Sigma_y = sqrt(sum_y2 ./ sum_weights);
    Sigma_px = (sum_xxp ./ sum_weights) ./ Sigma_x; Sigma_py = (sum_yyp ./ sum_weights) ./ Sigma_y;
    
    term_x = (sum_x2 .* sum_xp2) ./ (sum_weights.^2) - (sum_xxp ./ sum_weights).^2;
    Emit_x_proj = sqrt(max(term_x, 0));
    term_y = (sum_y2 .* sum_yp2) ./ (sum_weights.^2) - (sum_yyp ./ sum_weights).^2;
    Emit_y_proj = sqrt(max(term_y, 0));
    
    results = struct();
    results.s_array = s_envelope_common;
    results.sigx_array = Sigma_x; results.sigy_array = Sigma_y;
    results.sigpx_array = Sigma_px; results.sigpy_array = Sigma_py;
    results.x_centroid = mean_x_centroid; results.y_centroid = mean_y_centroid;
    results.px_centroid = mean_px_centroid; results.py_centroid = mean_py_centroid;
    results.gamma_array = mean_gamma;
    results.emitx_proj_array = Emit_x_proj; results.emity_proj_array = Emit_y_proj;
    if calc_k1_flag && ~exist('k1_values', 'var'), k1_values = []; end
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
    dydt = [0; 0; 0; 0; 0]; return;
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
    sc_x = 0; sc_y = 0;
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

% 注意：envelopeODE_quad_sliced 已被移除，因为现在直接使用 envelopeODE_sliced 并传入修正后的 k1
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