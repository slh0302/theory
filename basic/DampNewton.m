%% Newton
function [y, reInfo] = DampNewton(f, line_method, theta, X, lineFunc, varargin)
% �ú���ʹ�� ����ţ�� �㷨���x*.
%
% ����
%  [y, reInfo] = DampNewton(f, line_method, theta, X, @Func, f, numOfvar)
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

g = jacobian(f);
G =  jacobian(g);
% ��ʼ��
var_x = num2cell(sym('x',[1, numOfvalue]));
cell_x =  num2cell(X);
y0 = double(subs(f, var_x, cell_x));
g0 = double(subs(g, var_x, cell_x));
G0 = double(subs(G, var_x, cell_x));
% ��ʼ�������
if det(G0) == 0
    disp("Gk ����, ʹ�ø��ݶȷ���");
    dk0 =  -g0' ;
else
    dk0 =  - G0 \ g0' ;
end
disp(dk0);

% ��ʼ��������
Point = X;
if ~line_method.inextract
    [~, right] = JinTui(f, Point, dk0, line_method.step);
else
    right = line_method.inextract;
end
BeginPoint = Point;
Rule.crtr = line_method.crtr;
Rule.mthd = line_method.mthd;
if ~line_method.opt
    Rule.opt = [line_method.opt, right , line_method.max_iter, 0.05];
else
    Rule.opt = [line_method.opt, right , line_method.max_iter, 0.95, 0.05];
end
Step = dk0';
[StepSize, info_search, ~] = bolinesearch(lineFunc, BeginPoint, Step, Rule, varargin{:});

% taylor calc
xk = Point + StepSize * dk0';
xk_1 = Point;
k = 0;
yk_1 = y0;
yk = double(subs(f, var_x, num2cell(xk)));
dk = dk0;
gk = g0;

% �����������
iter_num = info_search(2);
feva_num = info_search(3);
while abs(yk - yk_1)>theta && norm(gk, 2) > theta
    fprintf('Iter %d:   Alph_k = %f,  ||dyk||=%.12f  yk-1: %.12f  yk=%.12f ||gk||=%f  xk= \n', k, StepSize, abs(yk_1-yk), yk_1, yk, norm(gk, 2) );
    disp(xk);
    yk_1 = yk;
    xk_1 = xk;
    % ����dk 
    gk = double(subs(g, var_x,  num2cell(xk_1)));
    Gk = double(subs(G, var_x, num2cell(xk_1)));
    % �ж�Gk�Ƿ�����
    if abs(det(Gk)) == 0
        disp("Gk ����");
        dk =  - gk' ;
    else
        dk =  - Gk \ gk' ;
    end
    % ������
    Point = xk_1;
    if ~line_method.inextract
        [~, right] = JinTui(f, Point, dk, line_method.step);
    else
        right = line_method.inextract;
    end
    BeginPoint = Point;
    if ~line_method.opt
        Rule.opt = [line_method.opt, right , line_method.max_iter, 0.05];
    else
        Rule.opt = [line_method.opt, right , line_method.max_iter, 0.95, 0.05];
    end
    Step = dk';
    [StepSize, info_search, ~] = bolinesearch(lineFunc, BeginPoint, Step, Rule, varargin{:});
    disp(StepSize); 
    if isnan(StepSize)
        StepSize = 1;
    end
    % ��һ�ε�������
    xk = xk_1 + StepSize * dk';
    yk = double(subs(f, var_x, num2cell(xk)));
   
    % ��Ϣ����
    iter_num = iter_num + info_search(2);
    feva_num = feva_num + info_search(3);
    k = k +1;
end
% ���һ�����
fprintf('Iter %d:   Alph_k = %f,  ||dyk||=%.12f  yk-1: %.12f  yk=%.12f ||gk||=%f  xk= \n', k, StepSize, abs(yk_1-yk), yk_1, yk, norm(gk, 2) );
disp(xk);

y = yk;
reInfo.all = k;
reInfo.iter = iter_num;
reInfo.feva_num = feva_num;

end