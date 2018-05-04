clear;
% ≥ı º÷µ
numOfvar = 10;
X = ones(1, numOfvar) * 2;
CalcutationFunc = @FR;
line_method.crtr = @bostwlf;                 %  boarmgld, bowlf, bostwlf
line_method.mthd = @bointrplt33;       %  bointrplt22, bointrplt33
line_method.opt = 0;                            %  0 extract line search; 1 inextract
line_method.max_iter = 10;
line_method.inextract = 0;
line_method.step = 0.00002;   
kexi = 1e-10;
theta1 = 0.1;
theta2 = 0.5;
beta= 1e-4;
M = 10;
line_method.others = [kexi, theta1, theta2, beta, M ];
func = @fgencube;
theta = 1e-8;
[y, info_Num] = CalcutationFunc(func, line_method, theta, X);
fprintf('Final Result,  f=%f  \n', y);
fprintf('func: %d, iter: %d, feva: %d  \n', info_Num.all, info_Num.iter, info_Num.feva_num);
disp('done');