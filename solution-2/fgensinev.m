% line search Function
function [F, g] = fgensinev(X_value, m)
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
    c1 = 1e-4;
    c2 = 4;
    F =0;
    g = zeros(1, numOfvar);
    var_x = sym('x',[1, numOfvar]);
    for i = 1: m-1   
        f1 = (var_x(i+1) - sin(var_x(i)))^2 / c1;
        f2 = var_x(i)^2/c2;
        if i == 1
            gf_tmp = 0;
        else
            gf_tmp =  (var_x(i) - sin(var_x(i-1)))^2 / c1;
        end
        f = f1 + f2;
        gf = gf_tmp + f;
        gf_f = diff(gf, var_x(i));
        F = F + double(subs(f, {var_x(i), var_x(i+1)}, {X_value(i), X_value(i+1)}));
        if i==1
            g(i) = double(subs(gf_f, {var_x(i), var_x(i+1)}, {X_value(i), X_value(i+1)}));
        else
             g(i) = double(subs(gf_f, {var_x(i-1), var_x(i), var_x(i+1)}, {X_value(i-1), X_value(i), X_value(i+1)}));
        end
        
        if i == m-1
            gf_f = diff(f, var_x(i+1));
            g(i+1) = double(subs(gf_f, {var_x(i), var_x(i+1)}, {X_value(i), X_value(i+1)}));
        end
    end

end
