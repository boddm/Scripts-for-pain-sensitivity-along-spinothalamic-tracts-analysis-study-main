function [Y0, design_names, nDesignScores] = myPLS_getY(behavMode, behaviorData, grouping, group_names, behav_names)
% =============================
% 功能：
%   根据分析模式（行为/对比/交互）生成PLS分析所需的设计矩阵Y0及其变量名。
%   支持2组/3组对比、行为与交互项自动生成。
%
% =============================
% 输入：
%   behavMode    ：字符串，行为分析类型
%                 'behavior'：仅行为变量
%                 'contrast'：组间对比
%                 'contrastBehav'：行为+对比
%                 'contrastBehavInteract'：行为+对比+交互
%   behaviorData ：N×B矩阵，原始行为/设计数据
%                 N：受试者数量，B：行为/设计变量数量
%   grouping     ：N×1向量，受试者分组信息
%   group_names  ：1×G元胞，组名称（G为组数）
%   behav_names  ：1×B元胞，行为变量名称
%
% =============================
% 输出：
%   Y0           ：N×D矩阵，PLS分析用设计矩阵
%                 D：设计变量数量，依赖于模式
%   design_names ：1×D元胞，设计变量名称
%   nDesignScores：整数，设计变量数量
%
% =============================
% 示例：
%   [Y0, design_names, nDesignScores] = myPLS_getY(behavMode, behaviorData, grouping, group_names, behav_names);
%
% =============================
% 历史：
%   2024-06-09  初版，标准化注释与格式，详细参数说明。
%   2024-06-10  增加详细中文注释，细化每个变量、流程、分支、异常处理说明。
% =============================

%% ========== 获取分组信息 ========== %%
groupIDs = unique(grouping);         % groupIDs：所有分组的唯一ID（如[1 2]）
nGroups  = length(groupIDs);         % nGroups：分组数
nBehav   = size(behaviorData, 2);    % nBehav：行为变量数
nSubj    = size(behaviorData, 1);    % nSubj：受试者数

% 保证名称为行向量，便于后续拼接
if size(behav_names, 1) > size(behav_names, 2)
    behav_names = behav_names';
end
if size(group_names, 1) > size(group_names, 2)
    group_names = group_names';
end

%% ========== 生成对比向量（如需） ========== %%
% 若分析模式包含对比（contrast），则生成组间对比向量c及其名称
if contains(behavMode, 'contrast')
    if nGroups == 1
        % 仅有一个组无法做对比分析，报错
        error('仅有一个组，无法进行对比分析');
    elseif nGroups == 2
        % 两组对比，生成一个对比向量c
        group1ID = grouping == groupIDs(1); nGroup1 = nnz(group1ID); % 第一组索引及人数
        group2ID = grouping == groupIDs(2); nGroup2 = nnz(group2ID); % 第二组索引及人数
        c = zeros(nSubj, 1);              % 初始化对比向量
        c(group1ID) = -1 / nGroup1;       % 第一组赋负权重
        c(group2ID) = 1 / nGroup2;        % 第二组赋正权重
        c_orig = c;                       % 保留原始对比向量
        c = orth(c_orig);                 % 正交化，保证正交性
        if sign(c(1)) ~= sign(c_orig(1)); c = -c; end % 保证方向一致
        contrastName = {['对比（' group_names{2} ' > ' group_names{1} '）']};
    elseif nGroups == 3
        % 三组对比，生成两个对比向量c(:,1), c(:,2)
        group1ID = grouping == groupIDs(1); nGroup1 = nnz(group1ID);
        group2ID = grouping == groupIDs(2); nGroup2 = nnz(group2ID);
        group3ID = grouping == groupIDs(3); nGroup3 = nnz(group3ID);
        c = zeros(nSubj, 2);              % 初始化对比矩阵
        % 第一对比：第3组 vs 第1/2组
        c(group1ID, 1) = -1 / (nGroup1 + nGroup2);
        c(group2ID, 1) = -1 / (nGroup1 + nGroup2);
        c(group3ID, 1) = 1 / nGroup3;
        % 第二对比：第2组 vs 第1组
        c(group1ID, 2) = -1 / nGroup1;
        c(group2ID, 2) = 1 / nGroup2;
        c_orig = c;                       % 保留原始对比矩阵
        c(:, 1) = orth(c_orig(:, 1));
        c(:, 2) = orth(c_orig(:, 2));
        if sign(c(1, 1)) ~= sign(c_orig(1, 1))
            c(:, 1) = -c(:, 1);
        end
        if sign(c(1, 2)) ~= sign(c_orig(1, 2))
            c(:, 2) = -c(:, 2);
        end
        contrastName = {['对比（' group_names{3} '>' group_names{1} '/' group_names{2} '）'], ['对比（' group_names{2} '>' group_names{1} '）']};
    else
        % 仅支持2组或3组对比，超出报错
        error('仅支持2组或3组对比');
    end
end

%% ========== 根据模式生成Y0 ========== %%
switch behavMode
    case 'behavior'
        % 仅行为变量，直接输出原始行为数据
        Y0 = behaviorData;
        design_names = behav_names;
    case 'contrast'
        % 仅对比，输出对比向量
        Y0 = c;
        design_names = contrastName;
    case 'contrastBehav'
        % 行为+对比，拼接对比向量和行为数据
        Y0 = [c, behaviorData];
        design_names = [contrastName, behav_names];
    case 'contrastBehavInteract'
        % 行为+对比+交互项
        if nGroups == 2
            % 两组时，每个行为变量生成主效应和交互项
            Y0 = [];
            for iBeh = 1:nBehav
                Y0(:, 2 * iBeh - 1) = behaviorData(:, iBeh);      % 主效应
                Y0(:, 2 * iBeh)     = zscore(Y0(:, 2 * iBeh - 1));% 标准化主效应
                Y0(:, 2 * iBeh)     = Y0(:, 2 * iBeh) .* c(:, 1); % 交互项：主效应*对比
                design_names{2 * iBeh - 1} = behav_names{iBeh};
                design_names{2 * iBeh}     = [behav_names{iBeh} '*' contrastName{1}];
            end
        elseif nGroups == 3
            % 三组时，每个行为变量生成主效应和两个交互项
            Y0 = [];
            for iBeh = 1:length(behav_names)
                Y0(:, 3 * iBeh - 2) = behaviorData(:, iBeh);         % 主效应
                Y0(:, 3 * iBeh - 1) = -zscore(Y0(:, 3 * iBeh - 2));  % 标准化主效应（负）
                Y0(:, 3 * iBeh)     = zscore(Y0(:, 3 * iBeh - 2));   % 标准化主效应（正）
                Y0(:, 3 * iBeh - 1) = Y0(:, 3 * iBeh - 1) .* c(:, 1);% 交互项1
                Y0(:, 3 * iBeh)     = Y0(:, 3 * iBeh) .* c(:, 2);    % 交互项2
                design_names{3 * iBeh - 2} = behav_names{iBeh};
                design_names{3 * iBeh - 1} = [behav_names{iBeh} '*' contrastName{1}];
                design_names{3 * iBeh}     = [behav_names{iBeh} '*' contrastName{2}];
            end
        end
        % 拼接对比向量和交互项
        Y0 = [c, Y0];
        design_names = [contrastName, design_names];
    otherwise
        % 未知行为类型，报错
        error('未定义的行为类型')
end

%% ========== 获取设计变量数量 ========== %%
nDesignScores = size(Y0, 2); % 设计变量数量

end