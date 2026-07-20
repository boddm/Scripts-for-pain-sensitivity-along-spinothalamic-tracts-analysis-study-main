function [X, meanX, stdX] = myPLS_norm(X, grouping, mode)
% =============================
% 功能：
% 对输入数据X进行归一化处理，支持全体zscore、组内zscore、标准差归一化等多种方式。
% 用于PLS分析前的数据标准化，保证不同变量/组间可比性。
%
% =============================
% 输入：
% X        ：N×V矩阵，原始数据
%            N：受试者数量，V：变量数量
% grouping ：N×1向量，受试者分组信息
% mode     ：归一化方式
%            0 = 不归一化
%            1 = 跨所有受试者z评分（均值0，方差1）
%            2 = 组内z评分（每组单独z评分，默认）
%            3 = 跨所有受试者标准归一化（无中心化，方差1）
%            4 = 组内标准归一化（无中心化，方差1）
%
% =============================
% 输出：
% X        ：归一化后的数据矩阵
% meanX    ：均值（每组或全体）
% stdX     ：标准差（每组或全体）
%
% =============================
% 示例：
% [X, meanX, stdX] = myPLS_norm(X, grouping, mode);
%
% =============================
% 历史：
% 2024-06-09  初版，标准化注释与格式，详细参数说明。
% =============================

%% ========== 获取分组信息 ========== %%
groupIDs = unique(grouping); % 所有组编号
nGroups = length(groupIDs);  % 组数

%% ========== 归一化主流程 ========== %%
switch mode
    case 1
        %% 跨所有受试者z评分
        meanX = mean(X); % 均值
        stdX = std(X);   % 标准差
        X = zscore(X);   % 归一化
    case 2
        %% 组内z评分
        for iG = 1:nGroups
            idx = find(grouping == groupIDs(iG)); % 本组索引
            meanX(iG, :) = mean(X(idx, :));       % 本组均值
            stdX(iG, :) = std(X(idx, :));         % 本组标准差
            X(idx, :) = zscore(X(idx, :));        % 本组归一化
        end
    case 3
        %% 跨所有受试者标准归一化（无中心化）
        meanX = mean(X); % 均值
        stdX = sqrt(mean(X.^2, 1)); % 标准差（无中心化）
        X2 = stdX;
        X = X ./ repmat(X2, [size(X, 1) 1]); % 归一化
    case 4
        %% 组内标准归一化（无中心化）
        for iG = 1:nGroups
            idx = find(grouping == groupIDs(iG)); % 本组索引
            meanX(iG, :) = mean(X(idx, :));       % 本组均值
            stdX(iG, :) = sqrt(mean(X(idx,:).^2, 1)); % 本组标准差（无中心化）
            X2 = stdX(iG, :);
            X(idx, :) = X(idx, :) ./ repmat(X2, [size(X(idx, :), 1) 1]); % 本组归一化
        end
end

end