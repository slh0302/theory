%% LM ʹ��G��������
function  [newG, judge] = FixLM(G, vk)
% �ú���ʹ�� DFP �㷨���x*.
%
% ����
%  [newG, judge] = FixLM(G)
%  [newG, judge] = FixLM(G, vk)
%
% Input 
% G:        ��ʼ�ĺ�ɭ����������������ֱ�ӷ��ظ�ֵ
% vk:       ����ϵ����Ĭ��ֵ0.01
%
% Output
% newG:   ��������ĺ�ɭ����
% judge�� ��ʼ�ĺ�ɭ�����Ƿ�����

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