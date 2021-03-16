function y = tent_mapping_once( x, mju )
%% 说明
%   功能
%   该函数对输入的x作一次Tent映射, 得到y, 映射的表达式为
%   if 0 < x <= 0.5, y = mju * x
%   if 0.5 < x < 1, y = mju * (1 - x)
%   即 y = mju * min{x, (1-x)}

%   参数
%   x: 输入映射系统的值, 要求是大于0且小于1的实数
%   y: 映射系统输出的值, 一定是一个大于0小于1的实数
%   mju: Tent映射系统本身的参数值

%% 执行映射
if x <= 0.5
    y = mju * x;
else
    y = mju * (1 - x);
end

end