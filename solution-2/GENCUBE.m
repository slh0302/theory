clear;
numOfvar = 10e3;
m = numOfvar;
var_x = sym('x',[1, numOfvar]);
f = (var_x(1) - 1) ^2 ;
before_x = f;
before_i = 1;
gf = sym('x',[1, numOfvar]);
for i=2:m
    f_tmp = 100 * (var_x(i) - var_x(i-1)^3)^2;
    g_f = before_x + f_tmp;
    gf(before_i) = diff(g_f, var_x(before_i));
    before_x =  f_tmp;
    before_i = i;
    f = f + f_tmp;
end
disp('func init done');
% 
% 
% ≥ı º÷µ
X = ones(1, numOfvar);

CalcutationFunc = @FR;
line_method.crtr = @boarmgld;                 %  boarmgld, bowlf, bostwlf
line_method.mthd = @bointrplt33;       %  bointrplt22, bointrplt33
line_method.opt = 0;                            %  0 extract line search; 1 inextract
line_method.max_iter = 10;
line_method.inextract = 0;
line_method.step = 0.01;                  
theta = 1e-8;
[y, info_Num] = CalcutationFunc(f, gf, line_method, theta, X, @Func, f, numOfvar, gf);
fprintf('Final Result,  f=%f  \n', y);
fprintf('func: %d, iter: %d, feva: %d  \n', info_Num.all, info_Num.iter, info_Num.feva_num);
disp('done');