close all;
clear all;
cons=[];

emitx = 5e-7;      % 归一化发射度 x (m·rad)
emity = 5e-7;      % 归一化发射度 y (m·rad)
betax = 23;        % beta函数 x (m)
betay = 143;        % beta函数 y (m)
alphax = 72;        % alpha函数 x
alphay = -4.3;        % alpha函数 y
E_target = 10.8e6;   % 目标能量 10 MeV (eV)
Q=5e-11;

% 物理常数
me_c2 = 0.511e6;   % 电子静止能量 (eV)

% 计算gamma因子
gamma = E_target / me_c2;

% 计算几何发射度
emitx_geo = emitx / gamma;  % 几何发射度 x
emity_geo = emity / gamma;  % 几何发射度 y

% 计算目标出口处的束斑尺寸 (σ)
sigx_exit = sqrt(betax * emitx_geo);
sigy_exit = sqrt(betay * emity_geo);

% 计算目标出口处的包络导数 (σ')
dsigx_exit = -alphax * sqrt(emitx_geo / betax);
dsigy_exit = -alphay * sqrt(emity_geo / betay);

% Drift 元件
ELE.MATCH_D0 = createElement('drift', 0.5, 0, 'MATCH_D0');
ELE.MATCH_D1 = createElement('drift', 0.3, 0, 'MATCH_D1');
ELE.MATCH_D2 = createElement('drift', 0.3, 0, 'MATCH_D2');
ELE.MATCH_D3 = createElement('drift', 0.3, 0, 'MATCH_D3');
ELE.MATCH_D4 = createElement('drift', 2.5, 0, 'MATCH_D4');
ELE.MATCH_D5 = createElement('drift', 0.32, 0, 'MATCH_D5');
ELE.MATCH_D6 = createElement('drift', 0.32, 0, 'MATCH_D6');
ELE.MATCH_D7 = createElement('drift', 2.7, 0, 'MATCH_D7');
ELE.D= createElement('drift', 0.1, 0, 'MATCH_D7');

% Quadrupole 元件 (长度, k1, 名称)
% k1 值放大10倍
ELE.MATCH_Q1 = createElement('quadrupole', 0.12, -9.33688275, 'MATCH_Q1');
ELE.MATCH_Q2 = createElement('quadrupole', 0.12, 13.64206447, 'MATCH_Q2');
ELE.MATCH_Q3 = createElement('quadrupole', 0.12, 5.28057923, 'MATCH_Q3');
ELE.MATCH_Q4 = createElement('quadrupole', 0.11, -5.89991838, 'MATCH_Q4');
ELE.MATCH_Q5 = createElement('quadrupole', 0.12, 7.86813968, 'MATCH_Q5');
ELE.MATCH_Q6 = createElement('quadrupole', 0.12, -3.76615007, 'MATCH_Q6');
ELE.MATCH_Q7 = createElement('quadrupole', 0.12, -11.92418163, 'MATCH_Q7');

% ========== Linac Section (加速段) ==========
% Drift 元件
ELE.LINAC_D0 = createElement('drift', 0.173, 0, 'LINAC_D0');
ELE.LINAC_D1 = createElement('drift', 0.285, 0, 'LINAC_D1');
ELE.LINAC_D2 = createElement('drift', 0.663, 0, 'LINAC_D2');

% Quadrupole 元件 (长度, k1, 名称)
% k1 值放大10倍，第三个重复第一个
ELE.LINAC_Q1 = createElement('quadrupole', 0.2, .5906, 'LINAC_Q1');
ELE.LINAC_Q2 = createElement('quadrupole', 0.2, -.13735, 'LINAC_Q2');
ELE.LINAC_Q3 = createElement('quadrupole', 0.2, 0.208, 'LINAC_Q3');

% Cavity 元件 (长度, [梯度MV/m, CELL数, 频率MHz, 相位deg], 名称)
% 梯度: 15.15 MV/m, 1.038 m 长度, 1300 MHz, -10 deg
ELE.CAV = createElement('cavity', 1.038, [15.2, 9, 1300, 8], 'CAV');

% ========== 构建完整 beamline ==========
% Match Section (7个四极铁 + drifts)
MATCH_SEQ = [ELE.MATCH_D0, ELE.MATCH_Q1, ELE.MATCH_D1, ...
    ELE.MATCH_Q2, ELE.MATCH_D2, ...
    ELE.MATCH_Q3, ELE.MATCH_D3, ...
    ELE.MATCH_Q4, ELE.MATCH_D4, ...
     ELE.MATCH_Q5, ELE.MATCH_D5, ...
     ELE.MATCH_Q6, ELE.MATCH_D6, ...
     ELE.MATCH_Q7, ELE.MATCH_D7
 ];

% Cavity Cluster 1 (8个 cavity)
CAV_CLUSTER_1 = [
    ELE.LINAC_D0, ELE.CAV, ELE.LINAC_D0, ...
     ELE.LINAC_D0, ELE.CAV, ELE.LINAC_D0, ...
     ELE.LINAC_D0, ELE.CAV, ELE.LINAC_D0, ...
     ELE.LINAC_D0, ELE.CAV, ELE.LINAC_D0, ...
     ELE.LINAC_D0, ELE.CAV, ELE.LINAC_D0, ...
     ELE.LINAC_D0, ELE.CAV, ELE.LINAC_D0, ...
     ELE.LINAC_D0, ELE.CAV, ELE.LINAC_D0, ...
     ELE.LINAC_D0, ELE.CAV, ELE.LINAC_D0];

% Cavity Cluster 2 (8个 cavity)
CAV_CLUSTER_2 = CAV_CLUSTER_1;

% Cavity Cluster 3 (8个 cavity)
CAV_CLUSTER_3 = CAV_CLUSTER_1;

% Cavity Cluster 4 (8个 cavity)
CAV_CLUSTER_4 = CAV_CLUSTER_1;

% 完整 Linac Section
LINAC_SEQ = [
    CAV_CLUSTER_1, ELE.LINAC_D1, ELE.LINAC_Q1, ELE.LINAC_D2];%, ...
    % CAV_CLUSTER_2, ELE.LINAC_D1, ELE.LINAC_Q2, ELE.LINAC_D2, ...
    % CAV_CLUSTER_3, ELE.LINAC_D1, ELE.LINAC_Q1, ELE.LINAC_D2, ...
    % CAV_CLUSTER_4
%];
data1 = load('TT_RF_Track.txt');
% data2 = load('sy.line_dat');
% data3 = load('norx.line_dat');
% data4 = load('nory.line_dat');
% 完整 beamline
beamline = [MATCH_SEQ, LINAC_SEQ];
% 初始参数
init_params = [betax, betay, alphax, alphay, emitx, emity, 1e-3,Q, E_target, 0.001];


[results] = SCEnvwe( beamline, init_params, 'fvv9p_slice_params.txt');
%[s_array, sigx_array, sigy_array, k1_values, emitx_proj_array, emity_proj_array, x_centroid, y_centroid]=SCEnv(ELE.MATCH_D0, init_params, 'fvv9mer_slice_params.txt');
%[s_array, sigx_array, sigy_array, k1_values] = SCEnv(beamline, init_params);%, 'beam_11_slice_params.txt');
figure;
% 第一张图：x方向束流包络(sigma_x)演化
subplot(2, 1, 1);
plot(results.s_array,results.sigx_array * 1e3, 'r-', 'LineWidth', 1.5);
hold on;
plot(results.s_array,results.sigy_array * 1e3, 'b-', 'LineWidth', 1.5);
hold on;

scatter(data1(:,1),100*data1(:,3),'red','filled','o','LineWidth',.05);
scatter(data1(:,1),100*data1(:,3),'blue','filled','o','LineWidth',.05);
xlim([0,20]);
xlabel('s (m)', 'FontSize', 12);
ylabel('\sigma (mm)', 'FontSize', 12);
title('\sigma_{x&y}', 'FontSize', 14);


% 第三张图：x方向质心动量(centpx)演化
subplot(2, 1, 2);
plot(results.s_array,results.gamma_array.* results.emitx_proj_array * 1e6, 'r-', 'LineWidth', 1.5);
hold on;


plot(results.s_array,results.gamma_array.* results.emity_proj_array * 1e6, 'b-', 'LineWidth', 1.5);
scatter(data1(:,1),1e6*data1(:,3),'red','filled','o','LineWidth',.05);
scatter(data1(:,1),1e6*data1(:,3),'blue','filled','o','LineWidth',.05);
xlabel('s (m)', 'FontSize', 12);
ylabel('\epsilon (mm·mrad)', 'FontSize', 15);
xlim([0,20]);
title('Emit', 'FontSize', 14);



% figure;
% 
% % 第一张图：x方向束流包络(sigma_x)演化
% subplot(4, 1, 1);
% plot(results.s_array, results.sigy_array * 1e3, 'b-', 'LineWidth', 1.5);
% xlabel('s (m)', 'FontSize', 12);
% ylabel('\sigma_x (mm)', 'FontSize', 12);
% title('Y sig', 'FontSize', 14);
% grid on;
% 
% % 第二张图：x方向质心位置(centx)演化
% subplot(4, 1, 2);
% plot(results.s_array, results.y_centroid * 1e3, 'r-', 'LineWidth', 1.5);
% xlabel('s (m)', 'FontSize', 12);
% ylabel('x_{cent} (mm)', 'FontSize', 12);
% title('Y cent', 'FontSize', 14);
% grid on;
% 
% % 第三张图：x方向质心动量(centpx)演化
% subplot(4, 1, 3);
% plot(results.s_array, results.py_centroid * 1e3, 'g-', 'LineWidth', 1.5);
% xlabel('s (m)', 'FontSize', 12);
% ylabel('p_{x,cent} (mrad)', 'FontSize', 12);
% title('PY cent', 'FontSize', 14);
% grid on;
% 
% % 第四张图：x方向发射度(emittance)演化
% subplot(4, 1, 4);
% plot(results.s_array, results.gamma_array.*results.emity_proj_array * 1e6, 'm-', 'LineWidth', 1.5);
% xlabel('s (m)', 'FontSize', 12);
% ylabel('\epsilon_x (mm·mrad)', 'FontSize', 12);
% %ylim([0,0.6]);
% title('Y emit', 'FontSize', 14);
% grid on;