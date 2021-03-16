function [ r_xy, X, Y ] = r_xy_near( I, n, direction )
% 参数说明
% I -- 输入的图像灰度矩阵
% n -- 在I中按照一定要求随机抽取的像素点的个数
% direction -- 将direction指定方向的像素作为相邻像素
%              Left: 水平向左的像素
%              Right: 水平向右的像素
%              Down: 垂直向下的像素
%              Up; 垂直向上的像素
%              Left_Up: 对角线左上角的像素
%              Left_Down: 对角线左下角的像素
%              Right_Up: 对角线右上角的像素
%              Right_Down: 对角线右下角的像素
% r_xy -- 最终计算得到的两个长度为n的向量X和Y之间的相关系数
% X -- 随机抽取的像素点的灰度值（向量）
% Y -- 随机抽取的像素点的相邻点的灰度值（向量）

I = double(I);
[H, W] = size(I);

% 根据direction参数确定随机取样点的横纵坐标范围
switch direction
    case 'Left'
        R = floor( unifrnd(1, H, [1, n]) );
        C = floor( unifrnd(2, W, [1, n]) );
        delta_r = 0; delta_c = -1;
    case 'Right'
        R = floor( unifrnd(1, H, [1, n]) );
        C = floor( unifrnd(1, W-1, [1, n]) );
        delta_r = 0; delta_c = 1;
    case 'Up'
        R = floor( unifrnd(2, H, [1, n]) );
        C = floor( unifrnd(1, W, [1, n]) );
        delta_r = -1; delta_c = 0;
    case 'Down'
        R = floor( unifrnd(1, H-1, [1, n]) );
        C = floor( unifrnd(1, W, [1, n]) );
        delta_r = 1; delta_c = 0;
    case 'Left_Up'
        R = floor( unifrnd(2, H, [1, n]) );
        C = floor( unifrnd(2, W, [1, n]) );
        delta_r = -1; delta_c = -1;
    case 'Left_Down'
        R = floor( unifrnd(1, H-1, [1, n]) );
        C = floor( unifrnd(2, W, [1, n]) );
        delta_r = 1; delta_c = -1;
    case 'Right_Up'
        R = floor( unifrnd(2, H, [1, n]) );
        C = floor( unifrnd(1, W-1, [1, n]) );
        delta_r = -1; delta_c = 1;
    case 'Right_Down'
        R = floor( unifrnd(1, H-1, [1, n]) );
        C = floor( unifrnd(1, W-1, [1, n]) );
        delta_r = 1; delta_c = 1;
end

X = []; Y = [];  % 初始化随机取样点处的灰度值向量
for i = 1 : n
    X(i) = I( R(i), C(i) );
    Y(i) = I( R(i) + delta_r, C(i) + delta_c );
end

% 计算向量X与向量Y的相关系数
mx = mean(X); my = mean(Y);
cov_xy = mean( (X - mx) .* (Y - my) );
cov_x = mean( (X - mx) .* (X - mx) );
cov_y = mean( (Y - my) .* (Y - my) );
r_xy = cov_xy / sqrt( cov_x * cov_y );

end