function [y, reInfo] = GN(f, Rx, theta, X)
% 该函数使用 BFGS 算法求解x*.
%
% 调用
%  [y, reInfo] = BFGS(f, line_method, theta, X, @Func, f, numOfvar)
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
cell_x =  num2cell(X);

% 计算一阶导数
Jrx = jacobian(Rx);
g = jacobian(f);

% 迭代
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