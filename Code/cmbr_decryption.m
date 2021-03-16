function [ I, C0 ] = cmbr_decryption( C, x0, s, S )
%% 说明
%   功能
%	cmbr是chaotic mapping & bit recombination的缩写
%   先对初始密文图像C进行灰度加密逆操作得到中间密文图像C0
%   再对C0进行两个阶段的置乱逆操作得到解密图像I
%   注意程序编写过程中默认图像是256个灰度级的

%   参数
%   C: 初始密文图像
%   x0: 输入Tent系统进行迭代的初值, 是解密系统的解钥之一
%   s: 原明文图像的灰度总和, 是解密系统的解钥之一
%   S: 解密系统的解钥之一
%   C0: 解密过程得到的中间密文图像
%   I: 解密过程得到的最终的解密图像

%% 输入参数的检验

% 要求初值x0是一个大于0且小于1的实数
if ~isreal(x0) || x0 <= 0 || x0 >= 1
    error('Tent映射初值x0必须是一个大于0且小于1的实数');
end

%% 利用解钥 s 还原出Tent混沌系统的控制参数mju和初始迭代次数times
% 然后一次性生成解密用的混沌序列E, F, R

% 初始化, 生成Tent混沌系统参数
C = double(C);
[M ,N] = size(C);
mju = 2 ^ ( s / (M * N * 255) );
times = 10 ^ 3 + mod(s, 10 ^ 3);

% 迭代times次消除初态效应的影响
x_iter = [];
for i = 1 : times
    x_iter = tent_mapping_once(x0, mju);
    x0 = x_iter;
end

% 生成混沌序列E
E = zeros(1, M);
for i = 1 : M
    E(i) = tent_mapping_once(x0, mju);
    x0 = E(i);
end

% 生成混沌序列F
F = zeros(1, 8 * N);
for i = 1 : 8 * N
    F(i) = tent_mapping_once(x0, mju);
    x0 = F(i);
end

% 生成混沌序列R
R = zeros(1, M * N);
for i = 1 : M * N
    R(i) = tent_mapping_once(x0, mju);
    x0 = R(i);
end

%% 第一阶段: 由初始密文图像C还原出中间密文图像C0
C0 = zeros(M, N);
for i = 1 : M
    for j = 1 : N
        R_index = N * (i - 1) + j;
        D = mod( ceil( R(R_index) * (2 ^ 48) ), 2 ^ 8 );  % 用到混沌序列R中对应的值
        if i == 1 || j == 1
            T = S;
        else
            t2 = mod(N * (i - 1) + j - 1, N);  
            t1 = 1 + ( N * (i - 1) + j - 1 - t2 ) / N;
            T = C( t1, t2 );
        end
        C0(i, j) = mod( C(i, j) - D - T, 2 ^ 8 );
%         while C0(i, j) < 0
%             C0(i, j) = C0(i, j) + 2 ^ 8;
%         end
    end
end

%% 第二阶段: 对中间密文图像C0分解成8个Bit矩阵, 对每个矩阵实行henon逆映射

II = zeros( M, 8 * N );
for i = 1 : M
    for j = 1 : N
        t = dec2bin( C0(i, j) );
        for k = 1 : numel(t)
            II(i, j + N * (8 - k)) = str2double( t(end - k + 1) );
        end
    end
end

% 求出每个Bit面对应的Henon映射系统参数a, c, iter_times
a = zeros(1, 8);  c = zeros(1, 8);  iter_times = zeros(1, 8);
for bit = 1 : 8
    F_index = 1 + 8 * N * (bit - 1) / 8;
    iter_times(bit) = 1 + mod(  ceil( F(F_index) * (10 ^ 14) ), 5  );
    a(bit) = mod(  ceil( F(F_index) * (10 ^ 14) ), 2 ^ 8  );
    c(bit) = mod(  ceil( (F(F_index) ^ 2) * (10 ^ 14) ), 2 ^ 8  );
end

% 对每个Bit面开始划分方形区域并在每个方形区域进行Henon置乱逆操作
t = zeros(M, 8 * N);
mark_h = 0;  mark_w = 0;
min_hw = min(M, N);  max_hw = max(M, N);
if max_hw == M
    delta_h = min_hw;  delta_w = 0;
else
    delta_h = 0;  delta_w = min_hw;
end

while 1
    r = mod(max_hw, min_hw);  q = ( max_hw - r ) / min_hw;
    for bit = 1 : 8
        for i = 0 : min_hw - 1
            for j = 0 : min_hw - 1
                x_pos2 = i;  y_pos2 = j;
                for iter = 1 : iter_times(bit)
                    [x_pos1, y_pos1] = ...
                        henon_mappinginv_once( x_pos2, y_pos2, ...
                            min_hw, min_hw, a(bit), c(bit) );
                    x_pos2 = x_pos1;  y_pos2 = y_pos1;
                end
                for num = 1 : q
                    t( mark_h + (num - 1) * delta_h + x_pos1 + 1, ...
                       (bit - 1) * N + mark_w + (num - 1) * delta_w + y_pos1 + 1) = ...
                    II( mark_h + (num - 1) * delta_h + i + 1, ...
                       (bit - 1) * N + mark_w + (num - 1) * delta_w + j + 1);
                end
            end
        end
    end
    
    if r == 0
        break
    else
        max_hw = min_hw;  min_hw = r;
        mark_h = mark_h + q * delta_h;  mark_w = mark_w + q * delta_w;
        if delta_h == 0
            delta_h = min_hw;  delta_w = 0;
        else
            delta_h = 0;  delta_w = min_hw;
        end
    end
    
end
II = t;

%% 第三阶段: 对进行henon逆映射后的数字图像II再进行先整列后整行的逆置乱操作
% 然后将8个Bit矩阵合并成十进制像素灰度值的明文图像矩阵I

% 整列逆置乱
[~, P_F] = sort(F);
II( :, P_F ) = II;

% 整行逆置乱
[~, P_E] = sort(E);
II( P_E, : ) = II;

% Bit矩阵合并
I = zeros(M, N);
for i = 1 : M
    for j = 1 : N
        t = '';
        for pos = 1 : 8
            t = [ t, num2str( II( i, 8 * (j - 1) + pos ) ) ];
        end
        I(i ,j) = bin2dec(t);
    end
end

C0 = uint8(C0);
I = uint8(I);

end