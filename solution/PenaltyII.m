%%  �ű�˵��
%Penalty II n = 2,4,6,8,10
% �������õĲ���
%           numOfvar: 
 %                ����Ĺ�ģ������ȡ��ֵΪ2,4,6,8,10
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
% ����Penalty II����
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

% X0�� ��ʼ��
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


