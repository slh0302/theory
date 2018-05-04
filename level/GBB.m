function [y, reInfo] = GBB(f, line_method, theta, X, varargin)
% �ú���ʹ�� FR �㷨���x*.
%
% ����
%  [y, reInfo] = BB(f, line_method, theta, X, @Func, f, numOfvar)
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

% ����һ�׵���
x0 = X;
[y0, g0] = f(x0, varargin{:});

% ����alph-k
kexi = line_method.others(1);
theta1 = line_method.others(2);
theta2 = line_method.others(3);
gama = line_method.others(4);
M = line_method.others(5);

% ����
gk_1 = g0;
xk_1 = x1;
yk_1 = y0;
alphk_1 = alph0;
before_f = zeros(1, M);
before_f(1) = y0;
k = 0;
iter_num = 0;
feva_num = 0;
while norm(gk_1, 2) > theta
    fprintf('Iter %d:   Alph_k = %f,\n', k, alphk_1 );
    disp('Xk= ');
    disp(double(xk_1));
    disp('gk= ');
    disp(double(gk_1));
    disp('yk= ');
    disp(double(yk_1));
    % ��������
    % alph ����
    alphk_1 = 1;
    if alphk_1 <= kexi || alphk_1 >= (1/ kexi)
        alphk_1 = chooseValue(gk_1);
    end
    lamdak_1 = 1/ alphk_1;
    % func, xk, gk, lamda, M
    [StepSize, info_search, perf] = nlineSearch(f, xk_1, gk_1, lamdak_1, k, M, before_f, gama, theta1, theta2);
    % alph = minValue(f, x0, d0);
    xk = xk_1 - StepSize * gk_1;
    yk = perf.F;
    gk = perf.g;
    alphk_1 = -(gk * yk') / (StepSize * (gk * gk'));
    gk_1 = gk;
    xk_1 = xk;
    % ��Ϣ����
    iter_num = iter_num + info_search(2);
    feva_num = feva_num + info_search(3);
    k = k + 1;
    if abs(yk - yk_1) < theta
        disp('Over');
        break;
    else
        yk_1 = yk;
        before_f( mod(k, m) + 1 ) = yk;
    end

end

y =  yk_1;
reInfo.all = k;
reInfo.iter = iter_num;
reInfo.feva_num = feva_num;
end

function beta = chooseValue(gk)
    norm_gk = norm(gk,2 );
    if norm_gk > 1
        beta = 1;
    elseif norm_gk <= 1 && norm_gk >= 1e-5
        beta = 1/ norm_gk;
    else
        beta = 1e5;
    end
end

function [stepsize, info, perf]  = nlineSearch(func, xk, gk, lamda, k, M, value_f, gama, theta1, theta2, max_iter)
    max_j = min(k, M);
    iter = 0;
    while iter < max_iter
        [flamda, ~] = func(xk - lamda * gk);
        info(2) = info(2) + 1;
        info(3) = info(3) + 1;
        for j = 0:max_j
            index = mod(k - j, M) + 1;
            fk_j  = value_f(index);
            if flamda > (fk_j - gama*lamda*(gk*gk'))
                break
            end
        end
        [fph0, gph0] = func(xk - theta1 * lamda* gk);
        [fph1, ~] = func(xk - theta2 * lamda* gk);
         info(3) = info(3) + 2;
        new_theta = -0.5 * gph0 * (theta2 - theta1)^2 / (fph1 - fph0 - gph0 * (theta2 - theta1) );
        if new_theta < 0
            new_theta = theta1;
        end
        lamda = new_theta * lamda;
        iter = iter + 1;
    end
    stepsize = lamda;
    [fp, gp] = func(xk -  lamda* gk);
    perf.F = fp;
    perf.g = gp;
    perf.x = xk - lamda*gk;
end