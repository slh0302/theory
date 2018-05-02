%% Newton
function [y, reInfo] = DampNewton(f, line_method, theta, X, lineFunc, varargin)
% 该函数使用 阻尼牛顿 算法求解x*.
%
% 调用
%  [y, reInfo] = DampNewton(f, line_method, theta, X, @Func, f, numOfvar)
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

g = jacobian(f);
G =  jacobian(g);
% 初始化
var_x = num2cell(sym('x',[1, numOfvalue]));
cell_x =  num2cell(X);
y0 = double(subs(f, var_x, cell_x));
g0 = double(subs(g, var_x, cell_x));
G0 = double(subs(G, var_x, cell_x));
% 初始方向计算
if det(G0) == 0
    disp("Gk 奇异, 使用负梯度方向");
    dk0 =  -g0' ;
else
    dk0 =  - G0 \ g0' ;
end
disp(dk0);

% 初始的线搜索
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

% 迭代计算次数
iter_num = info_search(2);
feva_num = info_search(3);
while abs(yk - yk_1)>theta && norm(gk, 2) > theta
    fprintf('Iter %d:   Alph_k = %f,  ||dyk||=%.12f  yk-1: %.12f  yk=%.12f ||gk||=%f  xk= \n', k, StepSize, abs(yk_1-yk), yk_1, yk, norm(gk, 2) );
    disp(xk);
    yk_1 = yk;
    xk_1 = xk;
    % 计算dk 
    gk = double(subs(g, var_x,  num2cell(xk_1)));
    Gk = double(subs(G, var_x, num2cell(xk_1)));
    % 判断Gk是否奇异
    if abs(det(Gk)) == 0
        disp("Gk 奇异");
        dk =  - gk' ;
    else
        dk =  - Gk \ gk' ;
    end
    % 线搜索
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
    % 下一次迭代计算
    xk = xk_1 + StepSize * dk';
    yk = double(subs(f, var_x, num2cell(xk)));
   
    % 信息计算
    iter_num = iter_num + info_search(2);
    feva_num = feva_num + info_search(3);
    k = k +1;
end
% 最后一次输出
fprintf('Iter %d:   Alph_k = %f,  ||dyk||=%.12f  yk-1: %.12f  yk=%.12f ||gk||=%f  xk= \n', k, StepSize, abs(yk_1-yk), yk_1, yk, norm(gk, 2) );
disp(xk);

y = yk;
reInfo.all = k;
reInfo.iter = iter_num;
reInfo.feva_num = feva_num;

end