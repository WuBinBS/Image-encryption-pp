function [ C, C0, s ] = cmbr_encryption( I, x0, S )
%% 说明
%	cmbr是chaotic mapping & bit recombination的缩写
%   功能
%   1. 计算原图像I的像素灰度值总和s, 作为解密程序所用的解钥之一
%   2. 先利用混沌映射和比特重组对明文图像I进行两阶段的置乱加密得到中间密文图像C0
%       再对中间密文图像C0进行灰度(扩散)加密得到最终密文图像C
%   注意程序编写过程中默认图像是256个灰度级的

%   参数
%   x0: 输入Tent系统进行迭代的初值, 相当于加密系统的密钥
%   I: 输入的明文图像
%   S: 输入的密钥之一, 供第三阶段的灰度加密用
%   s: 输出的明文图像像素灰度值总和, 作为解密系统的解钥之一
%   C0: 经过两阶段置乱加密得到的中间密文图像
%   C: C0再经过扩散加密得到的最终密文图像

%% 输入参数的检验

% 要求初值x0是一个大于0且小于1的实数
if ~isreal(x0) || x0 <= 0 || x0 >= 1
    error('Tent映射初值x0必须是一个大于0且小于1的实数');
end

%% 第一阶段: 全局置乱加密

% 初始化参数
I = double(I);  % 数据类型变为double类
[M, N] = size(I);  % 计算图像尺寸
s = sum( sum(I) );  % 计算图像I中像素灰度值的总和
mju = 2 ^ ( s / (M * N * 255) );  % 设置Tent混沌系统的控制参数mju
times = 10 ^ 3 + mod(s, 10 ^ 3);  % 设置Tent混沌系统初始迭代的次数times

% 对Tent混沌系统输入初始密钥x0, 并根据控制参数mju进行k次迭代消除初态效应的影响
x_iter = [];
for i = 1 : times
    x_iter = tent_mapping_once(x0, mju);
    x0 = x_iter;
end

% 将像素矩阵中每个像素灰度值转换成8位二进制数，
% 原矩阵中每个像素(二进制)灰度值的每一位都作为新矩阵II的一个元素
% 因此新矩阵II的规模是M * (8 * N)
II = [];
for i = 1 : M
    II_row = [];
    for j = 1 : N
        II_row = [II_row, zeros(1, 8)];
        t = dec2bin( I(i, j) );  % 将灰度值I(i,j)转化为二进制, 格式是string
        
        % 从string格式的t中以double格式提取每一位, 赋到相应的位置
        for k = 1 : numel(t)
            II_row(end - k + 1) = str2double( t(end - k + 1) );
        end
    end
    II = [II ; II_row];
end

% Tent混沌系统继续迭代M次, 产生长度为M的混沌序列E
% 将E按照升序排序，得到位置向量P_E, 利用生成的位置向量P_E
% 对已经转成比特的数字图像矩阵II进行整行置乱
E = zeros(1, M);
for i = 1 : M
    E(i) = tent_mapping_once( x0, mju );
    x0 = E(i);
end

[~, P_E] = sort(E);
II = II( P_E, : );

% Tent混沌系统继续迭代 8 * N 次, 产生长度为 8 * N 的混沌序列F
% 将F按照升序排序, 得到位置向量P_F, 利用生成的位置向量P_F
% 对数字图像矩阵II进行整列置乱
F = zeros(1, 8 * N);
for i = 1 : 8 * N
    F(i) = tent_mapping_once( x0, mju );
    x0 = F(i);
end

[~, P_F] = sort(F);
II = II( :, P_F );

%% 第二阶段: 各Bit面置乱加密

% 将第一阶段得到的置乱矩阵II从左到右均分成8个 M * N 的比特矩阵
% 对8个矩阵分别使用Henon映射进行置乱, 具体操作是, 每个矩阵分割出若干个方形区域,
% 对每个方形区域进行Henon映射置乱

% 求出每个Bit面对应的Henon映射系统参数a, c, iter_times
a = zeros(1, 8);  c = zeros(1, 8);  iter_times = zeros(1, 8);
for bit = 1 : 8
    F_index = 1 + 8 * N * (bit - 1) / 8;
    iter_times(bit) = 1 + mod(  ceil( F(F_index) * (10 ^ 14) ), 5  );
    a(bit) = mod(  ceil( F(F_index) * (10 ^ 14) ), 2 ^ 8  );
    c(bit) = mod(  ceil( (F(F_index) ^ 2) * (10 ^ 14) ), 2 ^ 8  );
end

% 对每个Bit面开始划分方形区域并在每个方形区域进行Henon映射置乱
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
                x_pos1 = i;  y_pos1 = j;
                for iter = 1 : iter_times(bit)
                    [x_pos2, y_pos2] = ...
                        henon_mapping_once( x_pos1, y_pos1, ...
                            min_hw, min_hw, a(bit), c(bit) );
                    x_pos1 = x_pos2;  y_pos1 = y_pos2;
                end
                for num = 1 : q
                    t( mark_h + (num - 1) * delta_h + x_pos2 + 1, ...
                       (bit - 1) * N + mark_w + (num - 1) * delta_w + y_pos2 + 1) = ...
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

% 8个比特矩阵各自迭代完成后, 将这8个比特矩阵按垂直于Bit平面的方向合并, 将比特转化为
% 十进制像素值, 得到中间密文图像C0
C0 = zeros(M, N);
for i = 1 : M
    for j = 1 : N
        t = '';
        for bit = 1 : 8
            t = [ t, num2str( II( i, j + N * (bit - 1) ) ) ];
        end
        C0(i ,j) = bin2dec(t);
    end
end

%% 第三阶段: 中间密文图像C0再加密 -- 灰度(扩散)加密

% Tent混沌系统继续迭代 M * N 次, 由此产生长度为 M * N 的混沌序列R
R = zeros(1, M * N);
for i = 1 : M * N
    R(i) = tent_mapping_once( x0, mju );
    x0 = R(i);
end

% 进行简单的扩散加密, 将混沌序列R中的M * N个元素与C0中的M * N个元素一一对应起来
% 加密后的矩阵C的每个像素(i, j)处的灰度值的计算公式为
% C(i, j) = ( D(i, j) + C0(i, j) ) mod 2 ^ 8 = 256
% 其中 D(i, j) = ceil( R(i, j) * 2 ^ 48 ) mod 2 ^ 8 = 256
C = zeros(M, N);
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
        C(i, j) = mod( D + C0(i, j) + T, 2 ^ 8 );
    end
end
C0 = uint8(C0);
C = uint8(C);

end