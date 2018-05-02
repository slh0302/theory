%% 脚本说明
% B3d function. m = 3,5,10,15,20
% 可以设置的参数
%           m: 
 %                问题的规模，可以取得值为3,5,10,15,20
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

% 定义 B3d function
numOfvar = 3;
var_x = sym('x',[1, numOfvar]);
f = 0;
m =3;
for i=1:m
    f_tmp = exp(-0.1 * i * var_x(1)) - exp(-0.1 * i * var_x(2)) - var_x(3) * (exp(-0.1*i) - exp(-i));
    f = f + f_tmp^2;
    disp(f_tmp);
end

% 算法调用参数设置
CalcutationFunc = @DampNewton;
line_method.crtr = @boarmgld;                 %  boarmgld, bowlf, bostwlf
line_method.mthd = @bointrplt33;       %  bointrplt22, bointrplt33
line_method.opt = 0;                            %  0 extract line search; 1 inextract
line_method.max_iter = 10;
line_method.inextract = 0;
line_method.step = 0.01;                  
theta = 1e-8;
X = [0, 10 , 20];
[y, info_Num] = CalcutationFunc(f, line_method, theta, X, @Func, f, numOfvar);
fprintf('Final Result,  f=%f  \n', y);
fprintf('func: %d, iter: %d, feva: %d  \n', info_Num.all, info_Num.iter, info_Num.feva_num);
disp('done');



