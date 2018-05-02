% %% ���˷�
% function [left, right] = JinTui(f, point, dk, step)
% % ���˷�������ȷ���µ�������
% %
% % ����
% % [left, right] = JinTui(f, point, dk)
% % [left, right] = JinTui(f, point, dk, step)
% % input
% % f:       �Ѿ�����ķ��ź��������� syms x1,x2; f = x1^2 + x2^2;
% % point: ���˷�����ʼ�㣬 nά����
% % dk:     ���˷���ʼ����ķ��� nά������
% % step:  ���˷������󲽳���Ĭ��ֵ0.001
% %
% % output
% % left:     ���˷��ĳ�ʼ��߽�
% % right:    ���˷��ĳ�ʼ�ұ߽�
% %
% % Create:   2018.04.17
% % Coder:    Su LiHui
% 
% if nargin==3
%     % ��ֻ����������ʱ��Ĭ�����ò���Ϊ0.001
%     step=0.001;
% end
% t = 2;
% alph = 0;
% [y0, ~] = f(point + alph * dk);
% alph = alph + step;
% [y1, ~] = f(point + alph * dk);
% total_step = 0;
% 
% % �ж��Ƿ�С��0
% if y1 - y0 > 0
%     step = - step;
%     alph =  0 + step;
%     [y1, ~] = f(point + alph * dk);
% end
% 
% % ��ʼ��alph���
% left_alph = 0;
% right_alph = alph;
% 
% % �жϵ�������
% while y1 < y0
%     left_alph = alph;
%     step = step * t;
%     alph = alph + step;
%     right_alph = alph;
%     y0 = y1;
%     [y1, ~] = f(point + alph * dk);
%     total_step = total_step + 1;
% end
% left = min(left_alph, right_alph);
% right = max(left_alph, right_alph); 
% 
% end
% 
% 
%% ���˷�
function [left, right] = JinTui(f, point, dk, step)
% ���˷�������ȷ���µ�������
%
% ����
% [left, right] = JinTui(f, point, dk)
% [left, right] = JinTui(f, point, dk, step)
% input
% f:       �Ѿ�����ķ��ź��������� syms x1,x2; f = x1^2 + x2^2;
% point: ���˷�����ʼ�㣬 nά����
% dk:     ���˷���ʼ����ķ��� nά������
% step:  ���˷������󲽳���Ĭ��ֵ0.001
%
% output
% left:     ���˷��ĳ�ʼ��߽�
% right:    ���˷��ĳ�ʼ�ұ߽�
%
% Create:   2018.04.17
% Coder:    Su LiHui

if nargin==3
    % ��ֻ����������ʱ��Ĭ�����ò���Ϊ0.001
    step=0.001;
end
t = 2;
alph = 0;
numOfvar = length(point);
var_x = num2cell(sym('x',[1, numOfvar]));
y0 = double(subs(f, var_x, num2cell(point + alph * dk')));
alph = alph + step;
y1 = double(subs(f, var_x, num2cell(point + alph * dk')));
total_step = 0;

% �ж��Ƿ�С��0
if y1 - y0 > 0
    step = - step;
    alph =  0 + step;
    y1 = double(subs(f, var_x, num2cell(point + alph * dk')));
end

% ��ʼ��alph���
left_alph = 0;
right_alph = alph;

% �жϵ�������
while y1 < y0
    left_alph = alph;
    step = step * t;
    alph = alph + step;
    right_alph = alph;
    y0 = y1;
    y1 = double(subs(f, var_x, num2cell(point + alph * dk')));
    total_step = total_step + 1;
end
left = min(left_alph, right_alph);
right = max(left_alph, right_alph); 

end
