%% LM 使得G保持正定
function  [newG, judge] = FixLM(G, vk)
% 该函数使用 DFP 算法求解x*.
%
% 调用
%  [newG, judge] = FixLM(G)
%  [newG, judge] = FixLM(G, vk)
%
% Input 
% G:        初始的海森矩阵，若正定则不修正直接返回该值
% vk:       增量系数，默认值0.01
%
% Output
% newG:   修正过后的海森矩阵
% judge： 初始的海森矩阵是否正定

% Create:   2018.04.18
% Coder:    Su LiHui
    if nargin == 1
        vk = 0.01;
    end
    I = eye(length(G));
    [~, p] = chol(G);
    judge = p;
    while p ~= 0
        G = G + vk * I;
        vk = vk * 2;
        [~, p] = chol(G);
    end
    newG = G;
end