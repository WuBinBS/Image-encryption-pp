%%  m文件说明
%   对《混沌映射与比特重组的图像加密》一文的图像加密算法和解密算法的性能分析
%   以下的每个分析都需要调用到独立的自定义函数文件cmbr_encryption.m, cmbr_decryption.m

%% 明文敏感性分析
%   分析方法
%   在明文图像中随机选取一个像素, 将该像素的值随机调为0到255的另外一个
%   不同的值, 再进行加密, 计算相应的指标R和U. 这样的操作过程重复若干次(例如100次)
%   需要额外调用到本m文件下方的函数pixel_change_rate( )和
%   normalized_mean_change_intensity( )

%   初始化, 设置参数
clear;  clc;  close all;
I = imread('fruit.png'); I = I(:, 1 : 323);  % 取方形区域
[H, W] = size(I);
exp_times = 100;  % 随机选取像素, 随机调整其灰度值的过程重复操作的次数
x0 = 0.234;  S = 1280;  % 加密密钥

%   按照上述分析方法重复实验
x_pos = 1 + floor( ( H - 1 ) * rand(1, exp_times) );
y_pos = 1 + floor( ( W - 1 ) * rand(1, exp_times) );
R = zeros(1, exp_times);  % 初始化R的结果序列
U = zeros(1, exp_times);  % 初始化U的结果序列

tic
for k = 1 : exp_times
    I1 = I;
    t = I( x_pos(k), y_pos(k) );
    while I1( x_pos(k), y_pos(k) ) == t
        I1( x_pos(k), y_pos(k) ) = floor( rand * 255 );
    end
    [ C, ~, ~ ] = cmbr_encryption( I, x0, S );
    [ C1, ~, ~ ] = cmbr_encryption( I1, x0, S );
    R(k) = pixel_change_rate( C, C1 );
    U(k) = normalized_mean_change_intensity( C, C1 );
    disp(['第', num2str(k), '次操作完成, 剩余', num2str(exp_times - k), '次']);
end
toc

%   可视化抽取的像素位置分布, 可视化序列R和序列U
figure,  scatter(x_pos, y_pos);  title([num2str(exp_times), '次抽取的像素位置分布']);
figure,  histogram(R);  title([num2str(exp_times), '次实验所得R直方图']);
figure,  histogram(U);  title([num2str(exp_times), '次实验所得U直方图']);

%%  密文敏感性分析 
%   分析方法
%   在密文图像中随机选取一个像素, 将该像素的值随机调为0到255的另外一个
%   不同的值, 再进行解密, 计算相应的指标R和U. 这样的操作过程重复若干次(例如100次)
%   需要额外调用到本m文件下方的函数pixel_change_rate( )和
%   normalized_mean_change_intensity( )

%   初始化, 设置参数
clear;  clc;  close all;
I = imread('fruit.png');  I = I(:, 1 : 323);  
[H, W] = size(I);
exp_times = 100;  % 随机选取像素, 随机调整其灰度值的过程重复操作的次数
x0 = 0.234;  S = 1280;  % 加密密钥
[ C, ~, s ] = cmbr_encryption( I, x0, S );

%   按照上述分析方法重复实验
x_pos = 1 + floor( ( H - 1 ) * rand(1, exp_times) );
y_pos = 1 + floor( ( W - 1 ) * rand(1, exp_times) );
R = zeros(1, exp_times);  % 初始化R的结果序列
U = zeros(1, exp_times);  % 初始化U的结果序列

tic
for k = 1 : exp_times
    C1 = C;
    t = C( x_pos(k), y_pos(k) );
    while C1( x_pos(k), y_pos(k) ) == t
        C1( x_pos(k), y_pos(k) ) = floor( rand * 255 );
    end
    [ I0, ~ ] = cmbr_decryption( C, x0, s, S );
    [ I1, ~ ] = cmbr_decryption( C1, x0, s, S );
    R(k) = pixel_change_rate( I0, I1 );
    U(k) = normalized_mean_change_intensity( I0, I1 );
    disp(['第', num2str(k), '次操作完成, 剩余', num2str(exp_times - k), '次']);
end
toc

%   可视化抽取的像素位置分布, 可视化序列R和序列U
figure,  scatter(x_pos, y_pos);  title([num2str(exp_times), '次抽取的像素位置分布']);
figure,  histogram(R);  title([num2str(exp_times), '次实验所得R直方图']);
figure,  histogram(U);  title([num2str(exp_times), '次实验所得U直方图']);

%%  解, 密钥敏感性分析
%   需要额外调用到本m文件下方的函数pixel_change_rate( )和
%   normalized_mean_change_intensity( )

clear;  clc;  close all;
I = imread('fruit.png');  I = I(:, 1 : 323);
x1 = 0.234;  x2 = 0.234 + 1e-10;  
S1 = 1280;  S2 = 1280 + 1;
[ C1, ~, s1 ] = cmbr_encryption( I, x1, S1 );  s2 = s1 + 1;

%%   单独针对密钥x0
clc;  close all;
[ C1, ~, ~ ] = cmbr_encryption( I, x1, S1 );  
figure,  imshow(C1);  title('密钥x0 = 0.234, S = 1280');
[ C2, ~, ~ ] = cmbr_encryption( I, x2, S1 );
figure,  imshow(C2);  title('密钥x0 = 0.234 + 1e-10, S = 1280');
disp(['同一幅明文图像fruit.png在密钥S = 1280, 而x0分别为0.234, 0.234 + 1e-10下', ...
    '加密得到的两幅密文图像的像素变化率R和归一化平均变化强度U分别为:']);
disp(['R = ', num2str( pixel_change_rate( C1, C2 ) )]);
disp(['U = ', num2str( normalized_mean_change_intensity( C1, C2 ) )]);

%%   单独针对密钥S
clc;  close all;
[ C1, ~, ~ ] = cmbr_encryption( I, x1, S1 );  
figure,  imshow(C1);  title('密钥x0 = 0.234, S = 1280');
[ C2, ~, ~ ] = cmbr_encryption( I, x1, S2 );
figure,  imshow(C2);  title('密钥x0 = 0.234, S = 1280 + 1');
disp(['同一幅明文图像fruit.png在密钥x0为0.234, 而S分别为1280, 1281下', ...
    '加密得到的两幅密文图像的像素变化率R和归一化平均变化强度U分别为:']);
disp(['R = ', num2str( pixel_change_rate( C1, C2 ) )]);
disp(['U = ', num2str( normalized_mean_change_intensity( C1, C2 ) )]);

%%   单独针对解钥x0
clc;  close all;
[ I1, ~ ] = cmbr_decryption( C1, x1, s1, S1 );
figure,  imshow(I1);  title('解钥x0 = 0.234, S = 1280, s = 13122840');
[ I2, ~ ] = cmbr_decryption( C1, x2, s1, S1 );
figure,  imshow(I2);  title('解钥x0 = 0.234 + 1e-10, S = 1280, s = 13122840');
disp(['同一幅密文图像在解钥S = 1280, s = 7290671, 而x0分别为0.234, 0.234 + 1e-10', ...
    '下解密得到的两幅明文图像的像素变化率R和归一化平均变化强度U分别为:']);
disp(['R = ', num2str( pixel_change_rate( I1, I2 ) )]);
disp(['U = ', num2str( normalized_mean_change_intensity( I1, I2 ) )]);

%%   单独针对解钥s
clc;  close all;
[ I1, ~ ] = cmbr_decryption( C1, x1, s1, S1 );
figure,  imshow(I1); title('解钥s = 13122840, S = 1280, x0 = 0.234');
[ I2, ~ ] = cmbr_decryption( C1, x1, s2, S1 );
figure,  imshow(I2); title('解钥s = 13122840 + 1, S = 1280, x0 = 0.234');
disp(['同一幅密文图像在解钥x0 = 0.234, S = 1280, 而s分别为7290671, 7290671 + 1', ...
    '下解密得到的两幅明文图像的像素变化率R和归一化平均变化强度U分别为:']);
disp(['R = ', num2str( pixel_change_rate( I1, I2 ) )]);
disp(['U = ', num2str( normalized_mean_change_intensity( I1, I2 ) )]);

%%   单独针对解钥S
clc;  close all;
[ I1, ~ ] = cmbr_decryption( C1, x1, s1, S1 );
figure,  imshow(I1);  title('解钥x0 = 0.234, S = 1280, s = 13122840');
[ I2, ~ ] = cmbr_decryption( C1, x1, s1, S2 );
figure,  imshow(I2);  title('解钥x0 = 0.234, S = 1280 + 1, s = 13122840');
disp(['同一幅密文图像在解钥x0 = 0.234, s = 7290671, 而S分别为1280, 1280 + 1', ...
    '下解密得到的两幅明文图像的像素变化率R和归一化平均变化强度U分别为:']);
disp(['R = ', num2str( pixel_change_rate( I1, I2 ) )]);
disp(['U = ', num2str( normalized_mean_change_intensity( I1, I2 ) )]);

%%  相邻像素的相关性分析
%   需要额外调用函数文件 r_xy_near.m

clear;  clc;  close all;
I = imread('fruit.png');  I = I(:, 1 : 323);
x0 = 0.234;  S = 1280;  n = 20000;
% direction = 'Left';
% direction = 'Right';
% direction = 'Up';
% direction = 'Down';
% direction = 'Left_Up';
% direction = 'Left_Down';
% direction = 'Right_Up';
direction = 'Right_Down';

[ r_xy_I, X_I, Y_I ] = r_xy_near( I, n, direction );
[ C, ~, ~ ] = cmbr_encryption( I, x0, S );
[ r_xy_C, X_C, Y_C ] = r_xy_near( C, n, direction );

figure,  scatter( X_I, Y_I );  title(['明文图像', num2str(n), '对相邻像素的灰度值分布']);
xlabel('(x,y)的像素值');  ylabel('(x,y)的相邻像素值');
disp(['明文图像随机取', num2str(n), '对相邻像素计算得相关系数: ', num2str(r_xy_I)]);
figure,  scatter( X_C, Y_C );  title(['密文图像', num2str(n), '对相邻像素的灰度值分布']);
xlabel('(x,y)的像素值');  ylabel('(x,y)的相邻像素值');
disp(['密文图像随机取', num2str(n), '对相邻像素计算得相关系数: ', num2str(r_xy_C)]);

%%  信息熵分析
%   对于一幅密文图像, 它的理想的信息熵的值是8.
%   调用Matlab自带函数entropy( )计算信息熵

clear;  clc;  close all;
I = imread('fruit.png');  I = I(:, 1 : 323);
x0 = 0.234;  S = 1280;
[C, ~, ~] = cmbr_encryption( I, x0, S );
disp(['密文图像的信息熵为: ', num2str( entropy(C) )]);

%%  加解密时间随着尺寸的变化

clear;  clc;  close all;
x0 = 0.234;  S = 1280;
size_sequence = 50 : 50 : 1000;
time_en = [];  time_de = [];
for k = size_sequence
    I = floor( rand(k, k) * 255 );
    
    tic
    [ C, ~, s ] = cmbr_encryption( I, x0, S );
    toc
    time_en = [time_en, toc];
    
    tic
    [ ~, ~ ] = cmbr_decryption( C, x0, s, S );
    toc
    time_de = [time_de, toc];
    
    disp(['图像大小为: ', num2str(k), ' * ', num2str(k), '的加解密已完成']);
end

figure,  plot(size_sequence, time_en, 'r');
title('加密算法用时与明文图像大小的关系曲线');
xlabel('明文图像大小');  ylabel('加密算法用时');

figure,  plot(size_sequence, time_de, 'b');
title('解密算法用时与密文图像大小的关系曲线');
xlabel('密文图像大小');  ylabel('解密算法用时');

%%  计算像素变化率R的函数, R越接近1 - 2 ^ (-8)(大约是0.9961), 明文敏感性越好
function R = pixel_change_rate( I1, I2 )
%   I1, I2是待比较的两个图像矩阵, 要求尺寸相同
[H1, W1] = size(I1);  [H2, W2] = size(I2);
if H1 ~= H2 || W1 ~= W2
    error('输入的两个图像矩阵规格不统一');
end

I1 = double(I1);  I2 = double(I2);
R = sum( sum( I1~=I2 ) ) / (H1 * W1);
end

%%  计算归一化平均变化强度U的函数, U的理想值为0.3446
function U = normalized_mean_change_intensity( I1, I2 )
%   I1, I2是待比较的两个图像矩阵, 要求尺寸相同
[H1, W1] = size(I1);  [H2, W2] = size(I2);
if H1 ~= H2 || W1 ~= W2
    error('输入的两个图像矩阵规格不统一');
end

I1 = double(I1);  I2 = double(I2);
U = sum( sum( abs(I1 - I2) ) ) / ( H1 * W1 * 255 );
end
