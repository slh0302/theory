%% SR1
function [y, reInfo] = SR1(f, line_method, theta, X, lineFunc, varargin)
% �ú���ʹ�� SR1 �㷨���x*.
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
%
% Create:   2018.04.17
% Coder:    Su LiHui

% ��ʼ������
numOfvalue  = length(X);
var_x = num2cell(sym('x',[1, numOfvalue]));
cell_x =  num2cell(X);

% ����һ�׵���
g = jacobian(f);
G =  jacobian(g);

y0 = double(subs(f, var_x, cell_x));
g0 = double(subs(g, var_x, cell_x));
G0 = double(subs(G, var_x, cell_x));
[~,p] = chol(G0);
if p ~= 0
    disp('G0��������')
    % DFP ����:  ��ʼ��H0 = I
    [G0,  ~] = FixLM(G0); 
    H0 =  G0 \ eye(numOfvalue);
else
    H0 = G0 \ eye(numOfvalue);
end

% �����ݶȷ���
dk0 = - H0 * g0';
disp(dk0);

% calc alphk -> step size, extract line search
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

% xk���³�ʼ��
xk = Point + StepSize * dk0';
xk_1 = Point;

% ��ʼ����
k = 0;
yk_1 = y0;
yk = double(subs(f, var_x, num2cell(xk)));
% ��ʼ��
dk = dk0;
gk_1 = g0;
Hk_1= H0;

% �����������
iter_num = info_search(2);
feva_num = info_search(3);
while abs(yk - yk_1)>theta && norm(gk_1, 2) > theta
    % ��ӡ�����Ϣ
    fprintf('Iter %d:   Alph_k = %f,  ||dyk||=%.12f  yk-1: %.12f  yk=%.12f ||gk_1||=%f  xk= \n', k, StepSize, abs(yk_1-yk), yk_1, yk, norm(gk_1, 2) );
    disp(xk);
    
    % ���¿�ʼ����
    % gk���� 
    gk = double(subs(g, var_x,  num2cell(xk)));
    % SR1 ����
    sk = xk - xk_1;
    delt_gk  = gk - gk_1;
    Hk = double(Hk_1 + (sk' - Hk_1 * delt_gk') * (sk' - Hk_1 * delt_gk')' / ((sk' - Hk_1 * delt_gk')' * delt_gk'));
    [~,p] = chol(Hk);
    if p ~= 0
        disp('HK��������');
        Hk = eye(numOfvalue);
    end
    % �������ݶȼ���
    dk =  - Hk * gk' ;
    
    % �µıȽϣ����¸�ֵ
    yk_1 = yk;
    xk_1 = xk;
    gk_1 = gk;
    Hk_1 = Hk;
    
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
        Rule.opt = [line_method.opt, right , line_method.max_iter, 0.85, 0.15];
    end
    Step = dk';
    [StepSize, info_search, ~] = bolinesearch(lineFunc, BeginPoint, Step, Rule, varargin{:});
    
    % ��һ��XK�ļ���
    xk = xk_1 + StepSize * dk';
    yk = double(subs(f, var_x, num2cell(xk)));
    
    % ��Ϣ����
    iter_num = iter_num + info_search(2);
    feva_num = feva_num + info_search(3);
    k = k +1;
end   
y = yk;
reInfo.all = k;
reInfo.iter = iter_num;
reInfo.feva_num = feva_num;

end