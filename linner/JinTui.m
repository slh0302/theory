% %% 进退法
% function [left, right] = JinTui(f, point, dk, step)
% % 进退法，用于确定下单峰区间
% %
% % 调用
% % [left, right] = JinTui(f, point, dk)
% % [left, right] = JinTui(f, point, dk, step)
% % input
% % f:       已经定义的符号函数，例如 syms x1,x2; f = x1^2 + x2^2;
% % point: 进退法的起始点， n维数组
% % dk:     进退法初始计算的方向， n维度数组
% % step:  进退法的增大步长，默认值0.001
% %
% % output
% % left:     进退法的初始左边界
% % right:    进退发的初始右边界
% %
% % Create:   2018.04.17
% % Coder:    Su LiHui
% 
% if nargin==3
%     % 当只有三个参数时，默认设置步长为0.001
%     step=0.001;
% end
% t = 2;
% alph = 0;
% [y0, ~] = f(point + alph * dk);
% alph = alph + step;
% [y1, ~] = f(point + alph * dk);
% total_step = 0;
% 
% % 判断是否小于0
% if y1 - y0 > 0
%     step = - step;
%     alph =  0 + step;
%     [y1, ~] = f(point + alph * dk);
% end
% 
% % 初始化alph结果
% left_alph = 0;
% right_alph = alph;
% 
% % 判断单调区间
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
%% 进退法
function [left, right] = JinTui(f, point, dk, step)
% 进退法，用于确定下单峰区间
%
% 调用
% [left, right] = JinTui(f, point, dk)
% [left, right] = JinTui(f, point, dk, step)
% input
% f:       已经定义的符号函数，例如 syms x1,x2; f = x1^2 + x2^2;
% point: 进退法的起始点， n维数组
% dk:     进退法初始计算的方向， n维度数组
% step:  进退法的增大步长，默认值0.001
%
% output
% left:     进退法的初始左边界
% right:    进退发的初始右边界
%
% Create:   2018.04.17
% Coder:    Su LiHui

if nargin==3
    % 当只有三个参数时，默认设置步长为0.001
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

% 判断是否小于0
if y1 - y0 > 0
    step = - step;
    alph =  0 + step;
    y1 = double(subs(f, var_x, num2cell(point + alph * dk')));
end

% 初始化alph结果
left_alph = 0;
right_alph = alph;

% 判断单调区间
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
