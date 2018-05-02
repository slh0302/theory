function [y, reInfo] = PRP(f, line_method, theta, X, varargin)
% �ú���ʹ�� FR �㷨���x*.
%
% ����
%  [y, reInfo] = PRP(f, line_method, theta, X, @Func, f, numOfvar)
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

% ����һ�׵���
x0 = X;
[~, g0] = f(x0, varargin{:});
d0 = -g0;

% ����alph-k
Point = X;
if ~line_method.inextract
    [~, right] = JinTui(f, Point, d0, line_method.step);
else
    right = line_method.inextract;
end
BeginPoint = Point;
Rule.crtr = line_method.crtr;
Rule.mthd = line_method.mthd;
if ~line_method.opt
    Rule.opt = [line_method.opt, right , line_method.max_iter, 0.05];
else
    Rule.opt = [line_method.opt, right , line_method.max_iter, 0.85, 0.15];
end
Step = d0;
[StepSize, info_search, perf] = bolinesearch(f, BeginPoint, Step, Rule, varargin{:});
% alph = minValue(f, x0, d0);
alph = StepSize;
% ��������
iter_num = info_search(2);
feva_num = info_search(3);

x1 = perf.x;
g1 = perf.g;

% ����
gk = g1;
gk_1 = g0;
xk_1 = x1;
dk_1 = d0;
yk_1 = perf.F;
k = 1;
while norm(gk_1, 2) > theta
    fprintf('Iter %d:   Alph_k = %f,\n', k, alph );
    disp('Xk= ');
    disp(double(xk_1));
    disp('gk= ');
    disp(double(gk));
    % PRP �����ĸ��½��
    betak_1 = (gk*(gk' - gk_1')) / (gk_1*gk_1');
    disp('beta_k-1= ');
    disp(double(betak_1));

    dk = -gk + betak_1 * dk_1;
    
    % ����
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
    Step = dk;
    [StepSize, info_search, perf] = bolinesearch(f, BeginPoint, Step, Rule, varargin{:});
%     alph = minValue(f, xk_1, dk);
    alph = StepSize;

%     syms A B C D;
%     f_alph = -(2*A*B+8*C*D-4*B-8*D)/(2*B^2 + 8*D^2);
%     input = [xk_1(1), dk(1), xk_1(2), dk(2)];
%     alph = double(subs(f_alph, {'A','B','C','D'}, input);

    % ����xk ��xk_1
    xk = perf.x;
    xk_1 = xk; % ����x2
    yk = perf.F;
      
    % ��һ�μ���
    gk_1 = gk;
    dk_1 = dk;
    gk = perf.g;   
    
    % ��Ϣ����
    iter_num = iter_num + info_search(2);
    feva_num = feva_num + info_search(3);
    k = k + 1;
    if abs(yk - yk_1) < theta
        disp('Over');
        break;
    else
        yk_1 = yk;
    end

end

y =  yk_1;
reInfo.all = k;
reInfo.iter = iter_num;
reInfo.feva_num = feva_num;
end