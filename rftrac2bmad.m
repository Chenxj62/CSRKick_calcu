% 修改束团纵向切片能散的脚本
clear all; clc;

% 设置参数
infile = 'Output_V5_B6d_200K.txt';   % 输入文件名
outputfile = 'fv12v5200kin.dat';    % 输出文件名
nc=50e-12;
c=2.99792e8;

% 设置Cambridge字体
set(0, 'DefaultAxesFontName', 'Cambria');
set(0, 'DefaultTextFontName', 'Cambria');
set(0, 'DefaultAxesFontSize', 12);
set(0, 'DefaultTextFontSize', 10);

%% 读取文件
% 读取文件头信息
 [x, xp, y, yp, z, p]=textread(infile,'%n%n%n%n%n%n','headerlines',0);


np=size(x);

Ene=mean(p);

fprintf('束团能量为:%f MeV',Ene);

%% 创建纵向切片
x=x/1000;
xp=xp/1000;
y=y/1000;
yp=yp/1000;
z=z-mean(z);
z=-1*z/1000;

p=(p-mean(p))/mean(p);







%% 写入修改后的文件
fid = fopen(outputfile, 'w');
a= {'!' 'ASCII::3';'0' '   ! ix_ele'; '1' '   ! n_bunch';np '   ! n_particle';[] 'BEGIN_BUNCH';[] 'Electron';nc '   ! bunch_charge_tot';-0 '        ! z_center'; 0 '        ! t_center'};

for i=1:1:3
    fprintf(fid, '%s', a{i,1});
    fprintf(fid, '%s\n', a{i,2});
end

for i=4:1:9
    fprintf(fid, '%d         ', a{i,1});
    fprintf(fid, '%s\n', a{i,2});
end

% 写入粒子数据（使用修改后的动量）
for j = 1:np
    fprintf(fid, '%.15f                 ', [x(j), xp(j), y(j), yp(j), z(j), p(j)]);
    fprintf(fid, '\n');
end

fprintf(fid, 'END_BUNCH\n');
fclose(fid);

fprintf('\n修改后的束流数据已保存到: %s\n', outputfile);



