clear;
numOfvar  =1000;
% m = numOfvar;
% var_x = sym('x',[1, numOfvar]);
% % func
% f = 0;
% f1_t = 0;
% f2_t = 0;
% for i = 1: m
%         f1 = var_x(i)^2;
%         f2 = i * (var_x(i)) *0.5;
%         f1_t = f1_t + f1;
%         f2_t = f2_t + f2;
%         disp(f2_t);
% end
% f = f1_t + f2_t^2 + f2_t^4;
% X = ones(1, numOfvar) * 0.1;
% y = double(subs(f, var_x, X));
% kkk = jacobian(f);
% for i = 1:m
%     disp(double(subs(diff(f, var_x(i)), var_x, X)));
% end
% g = double(subs(jacobian(f), var_x, X));
% ≥ı º÷µ
% numOfvar = 7;
X = ones(1, numOfvar) * 0.1;
 
% [f1,G] = fS303(X);
CalcutationFunc = @FR;
line_method.crtr = @bostwlf;                 %  boarmgld, bowlf, bostwlf
line_method.mthd = @bointrplt33;       %  bointrplt22, bointrplt33
line_method.opt = 0;                            %  0 extract line search; 1 inextract
line_method.max_iter = 10;
line_method.inextract = 0;
line_method.step = 0.00002;   
func = @fS303;
theta = 1e-8;
% [y, info_Num] = CalcutationFunc(f, line_method, theta, X, @Func, f, numOfvar);
[y, info_Num] = CalcutationFunc(func, line_method, theta, X);
fprintf('Final Result,  f=%f  \n', y);
fprintf('func: %d, iter: %d, feva: %d  \n', info_Num.all, info_Num.iter, info_Num.feva_num);
disp('done');