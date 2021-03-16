%%  m文件说明
%   对《混沌映射与比特重组的图像加密》一文的图像加密算法和解密算法的仿真实验
%   需要调用自定义函数文件cmbr_encryption.m, cmbr_decryption.m

%%  初始化与设置参数
clear;  clc;  close all;

% I0 = imread('fruit.png');  % 与原文相似的蔬果图, 大小为323 * 500
% I0 = I0(:, 1:323);  % 截取方形区域

% I0 = imread('cat.png');  % 招财猫图, 大小为651 * 500

% I0 = zeros(256, 256);  % 结果: 仿真能够进行

I0 = ones(256, 256) .* 255;  % 结果: 仿真能够进行

x0 = 0.234;  S = 1280;

%%  加密部分

%   展示原始明文图像I0和其直方图
figure,  imshow(I0);  title('明文图像');
figure,  imhist(I0);  title('明文图像直方图');

%   用密钥x0对明文图像I0进行加密, 得到中间密文图像C0, 密文图像C和解钥s, 并计时 
tic
[ C, C0, s ] = cmbr_encryption( I0, x0, S );
toc

%   展示中间密文图像C0, 最终密文图像C, 和各自的直方图
figure,  imshow(C0);  title('中间密文图像');
figure,  imhist(C0);  title('中间密文图像直方图');

figure,  imshow(C);  title('最终密文图像');
figure,  imhist(C);  title('最终密文图像直方图');

% close all
%%  解密部分

%   展示初始密文图像C和其直方图
figure,  imshow(C);  title('初始密文图像');
figure,  imhist(C);  title('初始密文图像直方图');

%   用解钥x0和s对密文图像C进行解密, 得到中间密文图像C1和明文图像I1, 并计时
tic
[I1, C1] = cmbr_decryption( C, x0, s, S );
toc

%   展示中间密文图像C1, 明文图像I1以及各自的直方图
figure,  imshow(C1);  title('中间密文图像');
figure,  imhist(C1);  title('中间密文图像直方图');

figure,  imshow(I1);  title('解密所得最终图像');
figure,  imhist(I1);  title('解密所得最终图像直方图');

if isequal(C0, C1)
    disp('解密过程得到的中间密文图像与加密过程的一致, 无损');
else
    disp('解密过程得到的中间密文图像与加密过程的不完全一致, 有损');
end

if isequal(I0, I1)
    disp('解密过程得到的解密图像与原明文图像一致, 无损');
else
    disp('解密过程得到的解密图像与原明文图像不完全一致, 有损');
end

% close all