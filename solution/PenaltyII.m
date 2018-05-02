%%  脚本说明
%Penalty II n = 2,4,6,8,10
% 可以设置的参数
%           numOfvar: 
 %                问题的规模，可以取得值为2,4,6,8,10
%           CalcutationFunc:
%                 求解方法：可以选的参数为：
%                   @DampNewton, @LM, @MixedNT, @DFP, @BFGS, @SR1    
%           line_method:
%                   line_method.ctr:        计算所采用的线搜索方法，
%                                               可以选的值为@boarmgld, @bowlf, @bostwlf
%                   line_method.mthd:       线搜索使用的插值外卖
%                                               可以选的值为@bointrplt22, @bointrplt33
%                   line_method.max_iter:  线搜索迭代次数，默认值10
%                   line_method.opt:  是否使用精确线搜索： 
%                                                       0值表示精确线搜索; 1值表示非精确
%                   line_method.inextract:  是否使用进退法获取线搜索区间的上限： 
%                                                       0值表示使用进退法;
%                                                      大于1的值表示搜索区间的上限，例如 10
%                   line_method.step:       使用进退法时设置，为进退法参数：
%                                                      建议值为0.01或者0.001，具体复现参数见文档
% 输出结果：打印到控制台
%
% Create:   2018.04.17
% Coder:    Su LiHui

clear;
% 定义Penalty II函数
numOfvar = 10;
var_x = sym('x',[1, numOfvar]);
f = 0;
a = 1e-5;
numOfvalue = numOfvar;
for i=1:numOfvalue
    if i == 1
        f_tmp = var_x(i) - 0.2;
    else
        value = exp(i/10) + exp((i-1)/10);
        f_tmp = sqrt(a) *(exp(var_x(i)/10) + exp(var_x(i - 1)/10) -  value);
    end
    f = f + f_tmp * f_tmp';
    disp(f_tmp);
end
% n -> 2n
for i=numOfvalue+1:2*numOfvalue
    if i == 2* numOfvalue
        ft = 0;
        for j = 1: numOfvalue
            ft = ft + (numOfvalue - j + 1) * (var_x(j)^2);
        end
        f_tmp = ft - 1;
    else
        f_tmp = sqrt(a) *(exp((var_x(i - numOfvalue + 1))/10) - exp(-1/10));
    end
    f = f + f_tmp * f_tmp';
     disp(f_tmp);
end

% X0， 初始点
X = ones(1, numOfvar) /2;

% Newton jintui 0.001
CalcutationFunc = @SR1;
line_method.crtr = @bostwlf;                  %  boarmgld, bowlf, bostwlf
line_method.mthd = @bointrplt33;       %  bointrplt22, bointrplt33
line_method.opt = 1;                           %  0 extract line search; 1 inextract
line_method.max_iter = 10;
line_method.inextract = 0;
line_method.step = 0.001;
theta = 1e-8;
[y, info_Num] = CalcutationFunc(f, line_method, theta, X, @Func, f, numOfvar);
fprintf('Final Result,  f=%.12f  \n', y);
fprintf('func: %d, iter: %d, feva: %d  \n', info_Num.all, info_Num.iter, info_Num.feva_num);
disp('done');


