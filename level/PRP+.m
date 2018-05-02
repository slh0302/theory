function [y, reInfo] = PRP+(f, theta, X)
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
g = jacobian(f);

x0 = X;
g0 = double(subs(g, var_x, cell_x));
d0 = -g0;
alph = minValue(f, x0, d0);

x1 = x0 + alph * d0;
g1 = double(subs(g, var_x, num2cell(x1)));

% ����
gk = g1;
gk_1 = g0;
xk_1 = x1;
dk_1 = d0;
k = 1;
while norm(gk_1, 2) > theta
    fprintf('Iter %d:   Alph_k = %f,\n', k, alph );
    disp('Xk= ');
    disp(double(xk_1));
    disp('gk= ');
    disp(double(gk));
    betak_1 = (gk*gk') / (gk_1*gk_1');
    disp('beta_k-1= ');
    disp(double(betak_1));

    dk = -gk + betak_1 * dk_1;
    
    alph = minValue(f, xk_1, dk);

%     syms A B C D;
%     f_alph = -(2*A*B+8*C*D-4*B-8*D)/(2*B^2 + 8*D^2);
%     input = [xk_1(1), dk(1), xk_1(2), dk(2)];
%     alph = double(subs(f_alph, {'A','B','C','D'}, input);

    % ����xk ��xk_1
    xk = xk_1 + alph * dk;
    xk_1 = xk; % ����x2
    % ��һ�μ���
    gk_1 = gk;
    dk_1 = dk;
    gk = double(subs(g, var_x, num2cell(xk)));
    
    k = k + 1;
end

    y =  double(subs(f, var_x, num2cell(xk)));
    reInfo = 0;
end