% line search Function
function [F, g] = Func(X_value, func, nOfvar, g_func)
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
    if nargin == 3
        g_func = jacobian(func);
    end
    F = double(subs(func, num2cell(sym('x',[1, nOfvar])), num2cell(X_value)));
    g = double(subs(g_func, num2cell(sym('x',[1, nOfvar])), num2cell(X_value)));
end
