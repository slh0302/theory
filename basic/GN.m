function [y, reInfo] = GN(f, Rx, theta, X)
% �ú���ʹ�� BFGS �㷨���x*.
%
% ����
%  [y, reInfo] = BFGS(f, line_method, theta, X, @Func, f, numOfvar)
%
% Input 
% f:        �Ѿ�����ķ��ź��������� syms x1,x2; f = x1^2 + x2^2;
% line_method:
%       line_method.ctr:        ���������õ�������������
%                                   ����ѡ��ֵΪ@boarmgld, @bowlf, @bostwlf
%       line_method.mthd:       ������ʹ�õĲ�ֵ����
%                                   ����ѡ��ֵΪ@bointrplt22, @bointrplt33
%       line_method.max_iter:  ����������������Ĭ��ֵ10
%       line_method.opt:  �Ƿ�ʹ�þ�ȷ�������� 
%                                           0ֵ��ʾ��ȷ������; 1ֵ��ʾ�Ǿ�ȷ
%       line_method.inextract:  �Ƿ�ʹ�ý��˷���ȡ��������������ޣ� 
%                                           0ֵ��ʾʹ�ý��˷�;
%                                          ����1��ֵ��ʾ������������ޣ����� 10
%       line_method.step:       ʹ�ý��˷�ʱ���ã�Ϊ���˷�������
%                                          ����ֵΪ0.01����0.001�����帴�ֲ������ĵ�
% theta:    �����ľ��ȣ�Ĭ��ֵ1e-8
% X:          �������ͣ� Ϊ��ʼ�㣬�ͷ��ź�����˳��һһ��Ӧ:
%                    [x1,x2,x3]->[0,10,20]
% lineFunc: �̶�ֵ:    @Func, Ҳ���Բ���Func�Լ�����һ��������������ĺ���ֵ�͵���ֵ�ú���
%
% Output
% y:   ���ŵ���ĺ���ֵ
% reinfo:     
%         reInfo.all�� �ú����ĵ��ô���
%         reInfo.iter �� �������ĵ�������;
%         reInfo.feva_num �� �����������ĵ��ô���;

% Create:   2018.04.17
% Coder:    Su LiHui
numOfvalue  = length(X);
var_x = num2cell(sym('x',[1, numOfvalue]));
cell_x =  num2cell(X);

% ����һ�׵���
Jrx = jacobian(Rx);
g = jacobian(f);

% ����
gkp1 = double(subs(g, var_x, cell_x));
xk = X;
k = 1;
while norm(gkp1, 2) > theta
    Jk = double(subs(Jrx, var_x, num2cell(xk)));
    Rk = double(subs(Rx, var_x, num2cell(xk)));
    dk = - (Jk' * Jk) \ (Jk'* Rk');
    fprintf('Iter %d:  \n', k );
    disp('Jk= ');
    disp(Jk);
    disp('dk= ');
    disp(dk);
    xkp1 = xk + dk';
    disp('Xk+1= ');
    disp(xkp1);
    
    gkp1 = double(subs(g, var_x, num2cell(xkp1)));
    xk = xkp1;
    k = k + 1;
end

    y =  double(subs(f, var_x, num2cell(xk)));
    reInfo = 0;
end