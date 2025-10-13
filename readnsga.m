% 给定文本内容
clear all;
close all;
inputFileName = 'populations.txt';
outputFileName = 'pop.txt';

plotis=0;

fidIn = fopen(inputFileName, 'r');
fidOut = fopen(outputFileName, 'w');

% 逐行检查并写入不包含字母的行
tline = fgetl(fidIn);
while ischar(tline)
    % 修改正则表达式，排除所有字母和#号
    if isempty(regexp(tline,  '[a-df-zA-DF-Z#]', 'once')) && ~isempty(strtrim(tline))
        fprintf(fidOut, '%s\n', tline);
    end
    tline = fgetl(fidIn);
end

% 关闭文件句柄
fclose(fidIn);
fclose(fidOut);

% 读取过滤后的文件
inputfile = 'pop.txt';  % 使用过滤后的文件
data = dlmread(inputfile);
  % filteredData2 =data2(data2(:,end-1)<.8&data2(:,end)<0.1, :);
  k=size(data,2);
  %data(:,k+1)=(data(:,k-1)+data(:,k))/0.6;
   % filteredData1 =data(data(:,end-1)<0.0032&data(:,end)<.01, :);
   
 filteredData1 =data(end-200:end, :);
% filteredData1 =filteredData1(filteredData1(:,end-1)<.01&filteredData1(:,end)<.12, :);
  kend=size(data,2);
    sorted_matrix = sortrows(filteredData1, kend);


%%%%%%%%%%%%%%%%%%%%%write file
 x=sorted_matrix(1,:);
 yout=filteredData1;

 outindex=1;

 outy=yout(outindex,:);

 











     if plotis==1
%    sorted_matrix=sorted_matrix(5:end,:);
    scatter(sorted_matrix(2:end,end-1),sorted_matrix(2:end,end), 30,'r','filled');
   hold on;
%   scatter(sorted_matrix(2,kend-1)/100,sorted_matrix(2, kend)/12, 50,'k','filled');
%     hold on;
%   scatter(sorted_matrix(end,kend-1)/100,sorted_matrix(end, kend)/12, 50,'k','filled');
%scatter(0.6,0.6, 50,'k','filled');
%   text(0.605,0.59, 'B','FontSize',15, 'VerticalAlignment','bottom', 'HorizontalAlignment','left');
 %  scatter(0.492,0.5194, 50,'k','filled');
 %  text(0.492,0.51947, 'A','FontSize',15, 'VerticalAlignment','top', 'HorizontalAlignment','left');
%    text(sorted_matrix(2,kend-1)/100,sorted_matrix(2, kend)/12, 'A','FontSize',12, 'VerticalAlignment','bottom', 'HorizontalAlignment','left');
% 
% 
     ylabel('\epsilon_{y}[\mumrad]');
     xlabel('\epsilon_{x}[\mumrad]');
 box on;
hold off;
ax = gca;

% 设置轴的线宽，加粗边框和轴线
ax.LineWidth = 2; % 你可以根据需要调整这个值
set(gca, 'FontSize', 12); % 设置字体大小
set(gca, 'FontWeight', 'bold'); % 设置字体粗细
  

    end



















%  pd=[filteredData1(:,10),filteredData1(:,11)/100,filteredData1(:,11)/10];
% % pd=[filteredData1(:,7),filteredData1(:,8)];
% p=parallelplot(pd);
% 
% %   p.Color={'k','r'};
% %  binEdges = [.06 .062];
% %  bins = {'1'};
% %  groupHeight = discretize(pd(:,2),binEdges,'categorical',bins);
% %   p.LineWidth = [0.5 3.5]
% %  p.GroupData = groupHeight;
%   p.CoordinateTickLabels = {'f1','f2','f3','f4'};
%    set(gcf,'position',[0,0,1500,500])