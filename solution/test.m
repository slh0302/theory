clear;
G= [2,1;0,8];
norm(G,2)
syms x1 x2;
r1x = x1.^3 - x2 - 1;
r2x =x1.^2 - x2;
fr = [r1x, r2x];
fx = (x1.^3 - x2 - 1)*(x1.^3 - x2 - 1)+ (x1.^2 - x2) * (x1.^2 - x2);
X = [1.5, 2.25];
GN(fx, fr, 1e-8, X);
disp(Y2);
% 
% 
% 初始值
% X = [0,0];
% [y,~] = FR(f, 1e-8, X);

% syms A B C D;
% f = (2*A*B+8*C*D-4*B-8*D)/(2*B^2 + 8*D^2);
% Y0 = double(subs(f, {'A','B','C','D'}, [4, 153.6, 8, 259.2]));
% disp(Y0);

% % 算法调用参数设置
% CalcutationFunc = @BFGS;
% line_method.crtr = @boarmgld;                 %  boarmgld, bowlf, bostwlf
% line_method.mthd = @bointrplt33;       %  bointrplt22, bointrplt33
% line_method.opt = 0;                            %  0 extract line search; 1 inextract
% line_method.max_iter = 10;
% line_method.inextract = 2;
% line_method.step = 0.01;                  
% theta = 1e-8;
% X = [0, 0];
% [y, info_Num] = CalcutationFunc(f, line_method, theta, X, @Func, f, numOfvar);
% fprintf('Final Result,  f=%f  \n', y);
% fprintf('func: %d, iter: %d, feva: %d  \n', info_Num.all, info_Num.iter, info_Num.feva_num);
% disp('done');
