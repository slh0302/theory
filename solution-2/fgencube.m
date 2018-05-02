% line search Function
function [F, g] = fgencube(X_value, m)
% �����������������X�㴦�ĺ���ֵ
% ����
% [F,g] = Func Func(X_value, func, nOfvar)
% Input
%     X_value: ��Ҫ����ĵ����꣬nά����
%     func���Ѿ�����ķ��ź��������� syms x1,x2; f = x1^2 + x2^2;��Ĭ�ϱ�������Ϊx1,x2,x3...
%     nOfvar: �����ĸ���
%
% output
%     F: X_value���ĺ���ֵ
%     g: X_value���ĵ���ֵ
% Create:   2018.04.17
% Coder:    Su LiHui
    numOfvar = length(X_value);
    if nargin == 1
        m = length(X_value);
    end
    var_x = sym('x',[1, numOfvar]);
    f = (var_x(1) - 1) ^2 ;
    before_x = f;
    gf = sym('x',[1, numOfvar]);
    g = zeros(1, numOfvar);
    F = double(subs(f, {var_x(1)}, {X_value(1)}));
    for i=2:m
        f_tmp = 100 * (var_x(i) - var_x(i-1)^3)^2;
        g_f = before_x + f_tmp;
        gf(i-1) = diff(g_f, var_x(i-1));
        before_x =  f_tmp;
        
        % ����ֵ�͵���ֵ����
        F = F + double(subs(f_tmp, {var_x(i-1), var_x(i)}, {X_value(i-1), X_value(i)}));
        if i == 2
            g(i-1) = double(subs(gf(i-1), {var_x(i-1), var_x(i)}, {X_value(i-1), X_value(i)}));
        else
            g(i-1) = double(subs(gf(i-1), {var_x(i-2), var_x(i-1), var_x(i)}, {X_value(i-2), X_value(i-1), X_value(i)}));
        end
    end 
    g(numOfvar) = double(subs(diff(before_x, var_x(numOfvar)), {var_x(numOfvar-1), var_x(numOfvar)}, {X_value(numOfvar-1), X_value(numOfvar)}));
end
