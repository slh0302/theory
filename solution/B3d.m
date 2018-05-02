%% �ű�˵��
% B3d function. m = 3,5,10,15,20
% �������õĲ���
%           m: 
 %                ����Ĺ�ģ������ȡ��ֵΪ3,5,10,15,20
%           CalcutationFunc:
%                 ��ⷽ��������ѡ�Ĳ���Ϊ��
%                   @DampNewton, @LM, @MixedNT, @DFP, @BFGS, @SR1    
%           line_method:
%                   line_method.ctr:        ���������õ�������������
%                                               ����ѡ��ֵΪ@boarmgld, @bowlf, @bostwlf
%                   line_method.mthd:       ������ʹ�õĲ�ֵ����
%                                               ����ѡ��ֵΪ@bointrplt22, @bointrplt33
%                   line_method.max_iter:  ����������������Ĭ��ֵ10
%                   line_method.opt:  �Ƿ�ʹ�þ�ȷ�������� 
%                                                       0ֵ��ʾ��ȷ������; 1ֵ��ʾ�Ǿ�ȷ
%                   line_method.inextract:  �Ƿ�ʹ�ý��˷���ȡ��������������ޣ� 
%                                                       0ֵ��ʾʹ�ý��˷�;
%                                                      ����1��ֵ��ʾ������������ޣ����� 10
%                   line_method.step:       ʹ�ý��˷�ʱ���ã�Ϊ���˷�������
%                                                      ����ֵΪ0.01����0.001�����帴�ֲ������ĵ�
% ����������ӡ������̨
%
% Create:   2018.04.17
% Coder:    Su LiHui

% ���� B3d function
numOfvar = 3;
var_x = sym('x',[1, numOfvar]);
f = 0;
m =3;
for i=1:m
    f_tmp = exp(-0.1 * i * var_x(1)) - exp(-0.1 * i * var_x(2)) - var_x(3) * (exp(-0.1*i) - exp(-i));
    f = f + f_tmp^2;
    disp(f_tmp);
end

% �㷨���ò�������
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



