%%  ���⣺ 
%       Least squares Problems:
%           3.7 Analysis of an Enzyme Reaction
%%  �ű�˵��
% ���� Analysis of an Enzyme Reaction ����
% m = 11, n= 4
% �������õĲ���
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

clear;
numOfvar = 4;
data_ui = [4.0e0, 2.0e0, 1.0e0, 5.0e-1, ...
                2.5e-1, 1.67e-1, 1.25e-1,1.0e-1, ...
                8.33e-2, 7.14e-2, 6.25e-2];
data_yi = [1.957e-1, 1.947e-1, 1.735e-1, 1.6e-1, ...
                8.44e-2, 6.27e-2, 4.56e-2, 3.42e-2, ...
                3.23e-2, 2.35e-2, 2.46e-2];
             
var_x = sym('x',[1, numOfvar]);
f = 0;
numOfvalue = 11;
for i=1:numOfvalue
    tmp1 =  data_ui(i) * (data_ui(i) +  var_x(2));
    tmp2 = data_ui(i)  * (data_ui(i) +  var_x(3)) +  var_x(4);
    f_tmp = data_yi(i) - var_x(1) * tmp1 / tmp2;
    f = f + f_tmp * f_tmp';
    disp(f_tmp);
end

% Xs�� ��ʼ��
X = [2.5e-1,  3.9e-1, 4.15e-1, 3.9e-1];

% Newton jintui 0.001
CalcutationFunc = @SR1;      %  DampNewton, LM, MixedNT, DFP, BFGS, SR1    
line_method.crtr = @boarmgld;                %  boarmgld, bowlf, bostwlf
line_method.mthd = @bointrplt33;       %  bointrplt22, bointrplt33
line_method.opt = 1;                            %  0 extract line search; 1 inextract
line_method.max_iter = 10;                   %  Default 10
line_method.inextract = 0;
line_method.step = 0.01;
theta = 1e-8;
[y, info_Num] = CalcutationFunc(f, line_method, theta, X, @Func, f, numOfvar);
fprintf('Final Result,  f=%.12f  \n', y);
fprintf('func: %d, iter: %d, feva: %d  \n', info_Num.all, info_Num.iter, info_Num.feva_num);
disp('done');


