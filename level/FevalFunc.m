% line search Function
function [F, g] = Func(X_value, func, nOfvar)
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
    [F, g] = func(X_value, nOfvar);
end
