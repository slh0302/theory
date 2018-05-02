% line search Function
function [F, g] = fpS303(X_value, m)
% 返回线搜索的所需的X点处的函数值
% 调用
% [F,g] = Func Func(X_value, func, nOfvar)
% Input
%     X_value: 需要计算的点坐标，n维数组
%     func：已经定义的符号函数，例如 syms x1,x2; f = x1^2 + x2^2;，默认变量名字为x1,x2,x3...
%     nOfvar: 变量的个数
%
% output
%     F: X_value处的函数值
%     g: X_value处的导数值
% Create:   2018.04.17
% Coder:    Su LiHui
    numOfvar = length(X_value);
    if nargin == 1
        m = length(X_value);
    end
    f1_t = 0;
    f2_t = 0;
    g = zeros(1, numOfvar);
    var_x = sym('x',[1, numOfvar]);
    parfor i = 1: m
        f1 = var_x(i)^2;
        f2 = i * (var_x(i)) *0.5;
        f1_t =  f1_t  + double(subs(f1, {var_x(i)}, {X_value(i)}));
        f2_t =  f2_t  + double(subs(f2, {var_x(i)}, {X_value(i)}));
    end
    F = f1_t + f2_t^2 + f2_t^4;
    parfor i = 1:m
        g(i) = X_value(i)*2 + 2 * f2_t* i * 0.5 +  4 * f2_t^3 * i * 0.5;
    end
end
