% 清除工作区和关闭所有图表
clear all;
close all;

% 基本参数
sigma = 30e-5;                         % 束长 m
N = 600;                              % 计算点数
rho = 10;                             % 半径 m
bunch_charge = 100;                   % 电子电荷 pC
h =.02;                             % 真空盒直径 m
n_imag = 5;                          % 镜像电荷数
total_angle = pi/9;                  % 总偏转角度
L = total_angle * rho;                % 总路径长度 m
max_phi =pi/9;%( 24*sigma/rho)^(1/3);                      % 最大φ角度
x=1.3;%下游距离 单位m, 稳态则设x=0 
dx=0.02;

x_init=.1;
k=8.987551787e9;
c = 299792458;                        % 光速 m/s

all_results = zeros(N, n_imag+1);     % 预分配空间，存储n=0到n=n_imag的结果
bunch_charge0 = bunch_charge*1e-12;   % 从pC转换为C
z_range = linspace(-5*sigma, 5*sigma, N);  % z的范围

e =1.602176634e-19;  % 电子电荷
% 定义高斯分布作为粒子线密度函数（全局）
lambda = @(z_prime) bunch_charge0/e*exp(-(z_prime)^2/(2*sigma^2)) / (sigma * sqrt(2*pi));
gam=1e9*e;
% 电子电量常数 (C)



if x==0

%% 第一部分：计算各阶镜像电荷的尾场

for n = 0:n_imag
    
    
    l = n*h;       % 屏蔽参数 l

    
    
    % 初始化CSR Wakefield数组
    wakefield_results = zeros(size(z_range));
    
    % 新的phi与z-z'的关系
    z_diff = @(phi, rho, l) phi*rho - sqrt(l^2 + 4*(phi/2 - phi^3/48)^2*rho^2);
    
    % 使用解析雅可比行列式表达式
    jacobian = @(phi, rho, l) (-rho + (4*(1/2 - phi^2/16)*(phi/2 - phi^3/48)*rho^2)/sqrt(l^2 + 4*(phi/2 - phi^3/48)^2*rho^2));
    
    % 直接计算尾场函数
    wake_func = @(phi, rho, l, e) direct_wake_calculation(phi, rho, l, e, c,k,0);
    
    % 循环计算每个z点的CSR Wakefield
    for i = 1:length(z_range)
        current_z = z_range(i);
        wake_sum = 0;
        lambda_z = lambda(current_z);

        if n > 0
            fac = 2;  % 镜像电荷因子
        else
            fac = 1;  % 实际电荷
        end
        
        % 积分上下限设置
        phi_upper = max_phi;
        phi_lower = 1E-17;
        
        % 创建积分点
        phi_points = 1000;
        phi_values = linspace(phi_upper, phi_lower, phi_points);
        dphi = phi_values(2) - phi_values(1);
        
        % 对phi积分
        for j = 1:length(phi_values)-1
            phi = phi_values(j);
            
            % 计算z-z'和z'
            z_diff_val = z_diff(phi, rho, l);
            z_prime = current_z - z_diff_val;
            dphi0=1e-5;
            % 计算被积函数
            wake_term = wake_func(phi, rho, l, e);
            lambda_value = lambda(z_prime);
            jac1=z_diff(phi, rho, l);
              jac2=z_diff(phi+dphi0, rho, l);
            %jac = jacobian(phi, rho, l);
            jac=-(jac2-jac1)/dphi0;
            % 积分计算，包含局部项修正
            integrand = wake_term * (lambda_value - lambda_z) * jac;
            
            % 考虑镜像电荷影响
            integrand = integrand * ((-1)^n);
            
            % 累加积分
            wake_sum = wake_sum + integrand * dphi;
        end
        
        % 保存结果
        wakefield_results(i) = fac*wake_sum;
        
        % 显示进度
        if mod(i, 20) == 0
            fprintf('已计算 n=%d, %d/%d 个点\n', n, i, length(z_range));
        end
    end
    
    % 将当前n值的结果保存到数组中
    all_results(:, n+1) = wakefield_results;
end

%% 第二部分：计算能量损失和绘制图形

% 计算叠加尾场
total_wakefields = zeros(N, n_imag+1);
for n = 0:n_imag
    if n == 0
        total_wakefields(:, n+1) = all_results(:, n+1);
    else
        total_wakefields(:, n+1) = total_wakefields(:, n) + all_results(:, n+1);
    end
end

% 找到n=0情况的最大值用于归一化
max_n0 =gam/L;%max(abs(all_results(:, 1)));
all_results_normalized = all_results / max_n0;
total_wakefields_normalized = total_wakefields / max_n0;

% 计算线密度

dz = z_range(2) - z_range(1);

% 计算各阶n对应的束团能量损失
bunch_energy_loss_eV = zeros(n_imag+1, 1);
for n = 0:n_imag

    ene_sum=0;
    for i=1:N
    lambdaz=lambda(z_range(i));
    % 正确公式：∫λ(z) * (dW(z)/ds) * dz * L
    integrand = lambdaz * total_wakefields(i, n+1);
    energy_loss_joules = integrand*dz * L;
    ene_sum=ene_sum+energy_loss_joules;
    % 转换为电子伏特
    
    end
bunch_energy_loss_eV(n+1) =ene_sum / e/(bunch_charge0/e)/1e6;

    fprintf('计算n=%d时的平均能量损失: %.6e MeV\n', n, bunch_energy_loss_eV(n+1));
end

% 创建从红[.6,0,0]到蓝[0,0,1]的渐变颜色
colors = zeros(n_imag+1, 3);
for i = 1:n_imag+1
    t = (i-1)/(n_imag);  % 从0到1的归一化参数
    colors(i,:) = [.8*(1-t), 0, t];  % 从[.6,0,0]渐变到[0,0,1]
end

% 绘制所有n值的归一化结果在一个图上
figure('Position', [100, 100, 900, 600]);
hold on;


for n = 0:n_imag
    plot(z_range/sigma, all_results_normalized(:, n+1), 'LineWidth', 2, 'Color', colors(n+1,:));
end

% 添加黑色的总尾场曲线
plot(z_range/sigma, total_wakefields_normalized(:, end), 'k-', 'LineWidth', 3);

hold off;
grid off;
box on;
ax = gca;

% 设置轴的线宽，加粗边框和轴线
ax.LineWidth = 2; 
xlabel('z/\sigma', 'FontName', 'Cambria Math', 'FontSize', 16);
ylabel('norm.CSR wake', 'FontName', 'Cambria Math', 'FontSize', 16);
title('CSR wake with different n and total', 'FontName', 'Cambria Math', 'FontSize', 18);

% 添加图例
legend_strs = cell(n_imag+2, 1);

for n = 0:n_imag
    legend_strs{n+1} = ['n = ', num2str(n), ', l = ', num2str(n*h)];
end
legend_strs{n_imag+2} = 'Total wakefield';
legend(legend_strs, 'Location', 'best', 'FontName', 'Cambria Math', 'FontSize', 12);
set(gca, 'FontName', 'Cambria Math', 'FontSize', 14);

% figure('Position', [100, 100, 300, 220]);
% 
% plot(z_range/sigma, total_wakefields_normalized(:, end), 'k-', 'LineWidth', 3);
% grid off;
% box on;
% ax = gca;
% ax.LineWidth = 2;

% 创建束团能量损失随n变化的散点连线图
figure('Position', [100, 100, 900, 600]);
plot(0:n_imag, bunch_energy_loss_eV, 'o-', 'LineWidth', 3, 'MarkerSize', 10, 'Color', [0.7 0.3 0.3], 'MarkerFaceColor', [0.7 0.3 0.3]);
grid off;
box on;
ax = gca;
ax.LineWidth = 2;
xlabel('n', 'FontName', 'Cambria Math', 'FontSize', 16);
ylabel('Electron Energy Loss (MeV)', 'FontName', 'Cambria Math', 'FontSize', 16);

set(gca, 'FontName', 'Cambria Math', 'FontSize', 14);

% 显示最终结果


else

all_imag_results=zeros(length(z_range),n_imag+2);
for n = 0:n_imag

    x_range = x_init:dx:x+x_init;  % 使用您已指定的dx和x值

% 为每个x值初始化结果数组
all_wakefield_results1 = zeros(length(z_range), length(x_range));
all_wakefield_results2 = zeros(length(z_range), length(x_range));
all_wakefield_results3 = zeros(length(z_range), length(x_range));
all_wakefield_results4 = zeros(length(z_range), length(x_range));

% 循环遍历不同的x值
for x_idx = 1:length(x_range)
    x_i = x_range(x_idx);
    fprintf('计算第%d个镜像电荷 x = %.3f (%d/%d)\n',n, x_i, x_idx, length(x_range));
    
    l = n*h;  % 屏蔽参数 l
    
    % 初始化四种不同的CSR Wakefield数组
    wakefield_results1 = zeros(size(z_range));  % 用于 lambda_value - lambda_z
    wakefield_results2 = zeros(size(z_range));  % 用于 -lambda_z
    wakefield_results3 = zeros(size(z_range));  % 用于 lambda_value
    wakefield_results4 = zeros(size(z_range));  % 用于 特殊项
    
    % 新的phi与z-z'的关系
    z_diff = @(phi, rho, l, x_i) x+rho*phi-sqrt(4 * rho^2 * sin(phi/2)^2 + x^2 + 4 * rho * x * sin(phi/2) * cos(phi/2)+l^2);
   % z_diff = @(phi, rho, l, x_i) rho*phi^3/24*((rho*phi+4*x)/(rho*phi+x));
    % % 使用解析雅可比行列式表达式
   jacc = @(phi, rho, l, x_i) rho*phi^2/8*(rho*phi+2*x)^2/(rho*phi+x)^2;
    % 直接计算尾场函数
    wake_func = @(phi, rho, l, e, c, k, x_i) direct_wake_calculation(phi, rho, l, e, c, k, x_i);
    
    % 预先计算积分点，避免在内循环中重复计算
    phi_upper = max_phi;
    phi_lower = 1E-17;
    phi_points = 1000;
    phi_values = linspace(phi_upper, phi_lower, phi_points);  % 注意从小到大排序
    dphi = phi_values(2) - phi_values(1);
    
    % 循环计算每个z点的CSR Wakefield
    for i = 1:length(z_range)
        current_z = z_range(i);
        wake_sum1 = 0;
        wake_sum2 = 0;
        wake_sum3 = 0;
        wake_sum4 = 0;
        lambda_z = lambda(current_z);
        
        if n > 0
            fac = 2;  % 镜像电荷因子
        else
            fac = 1;  % 实际电荷
        end
        
        % 对phi积分
        for j = 1:length(phi_values)
            phi = phi_values(j);
            
            % 计算z-z'和z'
            z_diff_val = z_diff(phi, rho, l, x_i);
            z_prime = current_z - z_diff_val;
            dphi0=1e-7;
            % 计算被积函数
            wake_term = wake_func(phi, rho, l, e, c, k, x_i);
            lambda_value = lambda(z_prime);
            jac1 = z_diff(phi, rho, l, x_i);
            jac2 = z_diff(phi+dphi0, rho, l, x_i);

            jac=-(jac2-jac1)/dphi0;
            %jac=-1*jacc(phi, rho, l,x_i);
            
            % 计算四种不同的积分
            integrand1 = 1*wake_term * (lambda_value)*jac;
            % integrand2 = 0*wake_term * (-lambda_z) * jac;
            % integrand3 = wake_term * jac;
            % integrand4 = 0;  % 这项暂时设为0，如果需要可以取消注释下一行
             %integrand4 = wake_term * -lambda_z*((current_z * rho * phi^3) / (6 * sigma^2) - (current_z * rho^2 * phi^4) / (8 * x_i * sigma^2)) * jac;
            
            % 考虑镜像电荷影响
            %integrand1 = integrand1 * ((-1)^n);
            % integrand2 = integrand2 * ((-1)^n);
            % integrand3 = integrand3 * ((-1)^n);
            % integrand4 = integrand4 * ((-1)^n);
            
            % 累加积分
            wake_sum1 = wake_sum1 + integrand1 * dphi;
            % wake_sum2 = wake_sum2 + integrand2 * dphi;
            % wake_sum3 = wake_sum3 + integrand3 * dphi;
            % wake_sum4 = wake_sum4 + integrand4 * dphi;
        end
        
        % 保存当前z点的结果
        wakefield_results1(i) = (-1)^n*fac * wake_sum1-4*k* e^2*lambda_z/(2*x_i);
        % wakefield_results2(i) = fac * wake_sum2;
        % wakefield_results3(i) = fac * wake_sum3;
        % wakefield_results4(i) = fac * wake_sum4;
        % 
        % 显示进度

    end
    
    % 保存当前x值的结果
    all_wakefield_results1(:, x_idx) = wakefield_results1';
    % all_wakefield_results2(:, x_idx) = wakefield_results2;
    % all_wakefield_results3(:, x_idx) = wakefield_results3;
    % all_wakefield_results4(:, x_idx) = wakefield_results4;
end
% 归一化
all_wakefield_results1 = all_wakefield_results1 * dx;
% all_wakefield_results2 = all_wakefield_results2 * dx;
% all_wakefield_results3 = all_wakefield_results3 * dx;
% all_wakefield_results4 = all_wakefield_results4 * dx;

wakefield_normalized1 = sum(all_wakefield_results1, 2);  % 2表示对每一行的元素求和
% wakefield_normalized2 = sum(all_wakefield_results2, 2);
% wakefield_normalized3 = sum(all_wakefield_results3, 2);
% wakefield_normalized4 = sum(all_wakefield_results4, 2);

all_imag_results(:,n+1)=wakefield_normalized1(:);




end


% 创建从红[.6,0,0]到蓝[0,0,1]的渐变颜色
colors = zeros(n_imag+1, 3);
for i = 1:n_imag+1
    t = (i-1)/(n_imag);  % 从0到1的归一化参数
    colors(i,:) = [.8*(1-t), 0, t];  % 从[.6,0,0]渐变到[0,0,1]
end

% 绘制所有n值的归一化结果在一个图上
figure('Position', [100, 100, 900, 600]);
hold on;


for n = 0:n_imag
    plot(z_range/sigma, all_imag_results(:, n+1)/gam, 'LineWidth', 2, 'Color', colors(n+1,:));
end
all_imag_results(:,end)=sum(all_imag_results(:,1:n_imag+1),2);
% 添加黑色的总尾场曲线
plot(z_range/sigma, all_imag_results(:, end)/gam, 'k-', 'LineWidth', 3);

hold off;
grid off;
box on;
ax = gca;

% 设置轴的线宽，加粗边框和轴线
ax.LineWidth = 2; 
xlabel('z/\sigma', 'FontName', 'Cambria Math', 'FontSize', 16);
ylabel('norm.CSR wake', 'FontName', 'Cambria Math', 'FontSize', 16);
title('CSR wake with different n and total', 'FontName', 'Cambria Math', 'FontSize', 18);

% 添加图例
legend_strs = cell(n_imag+2, 1);

for n = 0:n_imag
    legend_strs{n+1} = ['n = ', num2str(n), ', l = ', num2str(n*h)];
end
legend_strs{n_imag+2} = 'Total wakefield';
legend(legend_strs, 'Location', 'best', 'FontName', 'Cambria Math', 'FontSize', 12);
set(gca, 'FontName', 'Cambria Math', 'FontSize', 14);

% figure('Position', [100, 100, 300, 220]);
% 
% plot(z_range/sigma, total_wakefields_normalized(:, end), 'k-', 'LineWidth', 3);
% grid off;
% box on;
% ax = gca;
% ax.LineWidth = 2;



% 显示最终结果

end



% 直接计算dW/ds的函数
function wake = direct_wake_calculation(phi, rho, l, e, c,k,x)
    % 从基本定义计算尾场
    % 考虑相对论带电粒子在弯曲轨道中的辐射
    if x==0
    % β = 1（相对论极限）
    beta = 1;
    
    % 计算重要几何量
    L_gamma_r = 2 * rho * sin(phi/2);
    L_gamma = sqrt(L_gamma_r^2 + l^2);
    R = L_gamma;
    cos_theta2 = L_gamma_r / L_gamma;
    
    % 计算向量点积
    n_dot_beta = beta * cos(phi/2) * cos_theta2;
    n_dot_beta_dot_prime = (beta^2 * c/rho) * sin(phi/2) * cos_theta2;
    beta_dot_beta_dot_prime = (beta^3 * c/rho) * sin(phi);
    n_dot_n = 1;
    beta_dot_beta_prime = beta^2 * cos(phi);
    n_dot_beta_prime = beta * cos(phi/2) * cos_theta2;
    
    % 计算分子部分
    numerator = n_dot_beta * n_dot_beta_dot_prime -... 
                beta_dot_beta_dot_prime * n_dot_n -... 
                beta_dot_beta_prime * n_dot_beta_dot_prime + ...
                beta_dot_beta_dot_prime * n_dot_beta_prime;
    
    % 计算分母部分
    denominator = c * R * (1 - n_dot_beta_prime)^3;
    
    % 组合计算dW/ds
    if denominator <= 0
        wake = 0; % 避免除以零或负数
    else
        wake =k* e^2 * (numerator / denominator);
    end
    
    % 处理可能的无穷大或NaN
    if isinf(wake) || isnan(wake)
        wake = 0;
    end

    else



    
    
       % wake =4 *k* e^2 *rho / ( (rho*phi + 2*x)^2);
  beta = 1; % 相对论极限

        % 计算R = L_gamma（使用提供的公式）
        R = sqrt(4 * rho^2 * sin(phi/2)^2 + x^2 + 4 * rho * x * sin(phi/2) * cos(phi/2));

        % 计算theta角度
       % 使用余弦定理计算theta
        % 余弦定理: cos(theta) = (a^2 + b^2 - c^2)/(2*a*b)
        % 其中a和b是临边，c是对边
         L_arc = 2 * rho * sin(phi/2);
        cos_theta = (R^2 + L_arc^2 - x^2) / (2 * R * L_arc);

        % 确保cos_theta在[-1,1]范围内（处理数值误差）
        cos_theta = max(min(cos_theta, 1), -1);

        % 计算theta
        theta = acos(cos_theta);
        R_tot=sqrt(R^2+l^2);
        cos_a=R/R_tot;
        % 使用提供的公式计算向量点积
        n_dot_beta = beta * cos(theta - phi/2)*cos_a;
        n_dot_beta_dot_prime = (beta^2 * c / rho) * cos(pi/2-(theta + phi/2))*cos_a;
        beta_dot_beta_dot_prime = (beta^3 * c / rho) * sin(phi);
        n_dot_n = 1;
        beta_dot_beta_prime = beta^2 * cos(phi);
        n_dot_beta_prime = beta * cos(theta + phi/2)*cos_a;

        % 计算分子部分
        numerator = n_dot_beta * n_dot_beta_dot_prime -... 
                    beta_dot_beta_dot_prime * n_dot_n -... 
                    beta_dot_beta_prime * n_dot_beta_dot_prime + ...
                    beta_dot_beta_dot_prime * n_dot_beta_prime;

        % 计算分母部分
        denominator = c * R_tot * (1 - n_dot_beta_prime)^3;

        % 组合计算dW/ds
        if denominator <= 0
            wake = 0; % 避免除以零或负数
        else
            wake = k * e^2 * (numerator / denominator);
        end

        % 处理可能的无穷大或NaN
        if isinf(wake) || isnan(wake)
            wake = 0;
        end

   % 处理可能的无穷大或NaN

    end
end