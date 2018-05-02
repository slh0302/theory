clear;
% numOfvar  = 100;
% m = numOfvar;
% var_x = sym('x',[1, numOfvar]);
% % func
% c1 = 1e-4;
% c2 = 4;
% f = 0;
% for i = 1: m-1   
%     f1 = (var_x(i+1) - sin(var_x(i)))^2 / c1;
%     f2 = var_x(i)^2/c2;
%     f = f + f1 + f2;
% end
% ≥ı º÷µ
numOfvar = 20;
X = zeros(1, numOfvar) + -1;
X(1) = 4.712389;
% [f,G] = fgensinev(X);
CalcutationFunc = @FR;
line_method.crtr = @bostwlf;                 %  boarmgld, bowlf, bostwlf
line_method.mthd = @bointrplt33;       %  bointrplt22, bointrplt33
line_method.opt = 0;                            %  0 extract line search; 1 inextract
line_method.max_iter = 10;
line_method.inextract = 0;
line_method.step = 0.000002;   
func = @fgensinev;
theta = 1e-8;
% [y, info_Num] = CalcutationFunc(f, line_method, theta, X, @Func, f, numOfvar);
[y, info_Num] = CalcutationFunc(func, line_method, theta, X);
fprintf('Final Result,  f=%f  \n', y);
fprintf('func: %d, iter: %d, feva: %d  \n', info_Num.all, info_Num.iter, info_Num.feva_num);
disp('done');