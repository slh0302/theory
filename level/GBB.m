function [y, reInfo] = GBB(f, line_method, theta, X, varargin)
% 该函数使用 FR 算法求解x*.
%
% 调用
%  [y, reInfo] = BB(f, line_method, theta, X, @Func, f, numOfvar)
%
% Input 
% f:        已经定义的符号函数，例如 syms x1,x2; f = x1^2 + x2^2;
% line_method:
%       line_method.ctr:        计算所采用的线搜索方法，
%                                   可以选的值为@boarmgld, @bowlf, @bostwlf
%       line_method.mthd:       线搜索使用的插值外卖
%                                   可以选的值为@bointrplt22, @bointrplt33
%       line_method.max_iter:  线搜索迭代次数，默认值10
%       line_method.opt:  是否使用精确线搜索： 
%                                           0值表示精确线搜索; 1值表示非精确
%       line_method.inextract:  是否使用进退法获取线搜索区间的上限： 
%                                           0值表示使用进退法;
%                                          大于1的值表示搜索区间的上限，例如 10
%       line_method.step:       使用进退法时设置，为进退法参数：
%                                          建议值为0.01或者0.001，具体复现参数见文档
% theta:    收敛的精度，默认值1e-8
% X:          数组类型， 为初始点，和符号函数的顺序一一对应:
%                    [x1,x2,x3]->[0,10,20]
% lineFunc: 固定值:    @Func, 也可以参照Func自己定义一个计算线搜索点的函数值和导数值得函数
%
% Output
% y:   最优点出的函数值
% reinfo:     
%         reInfo.all： 该函数的调用次数
%         reInfo.iter ： 线搜索的迭代次数;
%         reInfo.feva_num ： 线搜索函数的调用次数;

% Create:   2018.04.17
% Coder:    Su LiHui

% 计算一阶导数
x0 = X;
[y0, g0] = f(x0, varargin{:});

% 计算alph-k
kexi = line_method.others(1);
theta1 = line_method.others(2);
theta2 = line_method.others(3);
gama = line_method.others(4);
M = line_method.others(5);

% 迭代
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
    % 迭代计算
    % alph 计算
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
    % 信息计算
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