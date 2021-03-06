function [x0, y0] = henon_mappinginv_once( x1, y1, M, N, a, c)
%% 说明
%   功能
%   对输入的二元组x1, y1作一次Henon逆映射得到x0, y0, 映射的表达式为
%   x0 = y1 - c mod N
%   y0 = x1 - 1 + a * (x0 ^ 2) mod M

%   注: 这里的"Henon逆映射"并不是严格数学意义上的逆映射, 应该理解为逆操作
%       举个例子, 显然地, 当M < N时, Henon映射不是一个满射, 自然就不是一个双射, 不可逆

%   参数
%   x1, y1: 输入的一个二元组, 代表逆映射前的(图像 / 比特)矩阵的某一像素的二维坐标
%   x0, y0: 输出的一个二元组, 代表原坐标(x1, y1)经逆向映射后在新(图像 / 比特)矩阵中的位置
%   M, N: 分别代表图像(比特)矩阵的行数和列数
%   a, c: Henon逆映射系统的参数. 要求a是大于0的实数

%% 执行逆映射
x0 = mod( y1 - c, N );
y0 = mod( x1 - 1 + a * (x0 ^ 2), M );

end
