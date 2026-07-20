function all_boot_orders = myPLS_get_boot_orders(nBootstraps, grouping, grouped_boot)

% =============================
% 功能：
% 生成PLS引导抽样（bootstrap）时的重采样顺序矩阵。支持组内/组间抽样。
%
% =============================
% 输入：
% nBootstraps ：整数，引导样本数
% grouping    ：N×1向量，受试者分组信息
%               例如[1,1,2]表示前两名为组1，第三名为组2
% grouped_boot：逻辑值/0或1，引导抽样是否组内抽样
%               0 = 忽略分组，整体抽样
%               1 = 组内独立抽样
%
% =============================
% 输出：
% all_boot_orders ：N×nBootstraps矩阵，每列为一次bootstrap的受试者索引
%                   N：受试者数量
%
% =============================
% 示例：
% all_boot_orders = myPLS_get_boot_orders(nBootstraps, grouping, grouped_boot);
%
% =============================
% 历史：
% 2024-06-09  初版，标准化注释与格式，详细参数说明。
% =============================

%% ========== 检查是否组内抽样 ========== %%
if ~grouped_boot
    grouping = ones(size(grouping)); % 若非组内抽样，所有受试者视为一组
end

%% ========== 获取分组信息 ========== %%
groupIDs = unique(grouping); % 所有组的编号
nGroups = length(groupIDs);  % 组数
nSubj = length(grouping);    % 受试者总数

all_boot_orders = nan(nSubj, nBootstraps); % 初始化

%% ========== 组内/组间重采样主循环 ========== %%
for iG = 1:nGroups
    groupID = find(grouping == groupIDs(iG)); % 本组受试者索引
    num_subj_group = length(groupID);         % 本组受试者数
    [boot_order_tmp, ~] = rri_boot_order(num_subj_group, 1, nBootstraps); % 本组重采样顺序
    all_boot_orders(groupID, :) = groupID(boot_order_tmp); % 填充到总矩阵
end

end