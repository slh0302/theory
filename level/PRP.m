function [y, reInfo] = PRP(f, line_method, theta, X, varargin)
% 该函数使用 FR 算法求解x*.
%
% 调用
%  [y, reInfo] = PRP(f, line_method, theta, X, @Func, f, numOfvar)
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
numOfvalue  = length(X);
var_x = num2cell(sym('x',[1, numOfvalue]));

% 计算一阶导数
x0 = X;
[~, g0] = f(x0, varargin{:});
d0 = -g0;

% 计算alph-k
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
% 迭代计算
iter_num = info_search(2);
feva_num = info_search(3);

x1 = perf.x;
g1 = perf.g;

% 迭代
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
    % PRP 方法的更新结果
    betak_1 = (gk*(gk' - gk_1')) / (gk_1*gk_1');
    disp('beta_k-1= ');
    disp(double(betak_1));

    dk = -gk + betak_1 * dk_1;
    
    % 计算
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

    % 保存xk 到xk_1
    xk = perf.x;
    xk_1 = xk; % 保存x2
    yk = perf.F;
      
    % 下一次计算
    gk_1 = gk;
    dk_1 = dk;
    gk = perf.g;   
    
    % 信息计算
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