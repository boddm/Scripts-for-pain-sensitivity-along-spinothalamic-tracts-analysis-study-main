function res = myPLS_analysis(input, pls_opts)
% =============================
% 功能：
% 主函数，执行PLS分析。包括数据归一化、协方差矩阵计算、奇异值分解、PLS得分与载荷计算、排列检验与引导检验，输出所有分析结果。
% 支持分组、行为/对比/交互分析、各类归一化与统计检验。
%
% =============================
% 输入：
% input    ：结构体，包含分析输入数据：
%   .brain_data   ：N×M矩阵，原始影像数据
%                   N：受试者数量
%                   M：成像变量数量
%   .behav_data   ：N×B矩阵，原始行为/设计数据
%                   B：行为/设计变量数量
%   .grouping     ：N×1向量，受试者分组信息
%                   例如[1,1,2]表示前两名为组1，第三名为组2
%   .group_names  ：1×G元胞，组名称（可选，G为组数）
%   .behav_names  ：1×B元胞，行为变量名称（可选）
%   .img_names    ：1×M元胞，影像变量名称（可选）
%
% pls_opts ：结构体，PLS分析参数，需包含如下字段：
%   .nPerms              ：整数，排列检验次数
%   .nBootstraps         ：整数，引导样本数
%   .normalization_img   ：整数，影像数据归一化方式
%                         0 = 不归一化
%                         1 = 跨所有受试者z评分（均值0，方差1）
%                         2 = 组内z评分（每组单独z评分，默认）
%                         3 = 跨所有受试者标准归一化（无中心化，方差1）
%                         4 = 组内标准归一化（无中心化，方差1）
%   .normalization_behav ：整数，行为/设计数据归一化方式
%                         0 = 不归一化
%                         1 = 跨所有受试者z评分（均值0，方差1）
%                         2 = 组内z评分（每组单独z评分，默认）
%                         3 = 跨所有受试者标准归一化（无中心化，方差1）
%                         4 = 组内标准归一化（无中心化，方差1）
%   .grouped_PLS         ：逻辑值/0或1，是否分组PLS
%                         0 = 所有受试者整体计算PLS
%                         1 = 组间协方差矩阵拼接（行为PLS）
%   .grouped_perm        ：逻辑值/0或1，排列检验是否组内排列
%                         0 = 忽略分组排列
%                         1 = 组内排列
%   .grouped_boot        ：逻辑值/0或1，引导抽样是否组内抽样
%                         0 = 忽略分组引导
%                         1 = 组内引导
%   .boot_procrustes_mod ：整数，引导过程Procrustes旋转模式
%                         1 = 仅U旋转（默认）
%                         2 = U和V平均旋转
%   .save_boot_resampling：逻辑值/0或1，是否保存引导重采样数据
%                         0 = 不保存
%                         1 = 保存
%   .behav_type          ：字符串，行为分析类型
%                         'behavior'：标准行为PLS（默认）
%                         'contrast'：仅计算两个组对比
%                         'contrastBehav'：结合对比和行为测量
%                         'contrastBehavInteract'：考虑组-行为交互效应
%
% =============================
% 输出：
% res ：结构体，包含所有PLS分析结果：
%   .X0, .Y0          ：未归一化输入矩阵
%   .meanX0, .stdX0   ：X0均值、标准差
%   .meanY0, .stdY0   ：Y0均值、标准差
%   .X, .Y            ：归一化输入矩阵
%   .design_names     ：设计变量名称（如行为/对比/交互）
%   .grouping         ：分组信息
%   .group_names      ：组名称
%   .img_names        ：影像变量名称（如有）
%   .R                ：L×M矩阵，brain-behavior协方差矩阵
%   .U                ：B×L矩阵，行为显著性
%   .V                ：M×L矩阵，影像显著性
%   .S                ：L×L对角矩阵，奇异值
%   .explCovLC        ：每个LC解释的协方差比例
%   .LC_pvals         ：每个LC的p值（排列检验）
%   .Lx               ：N×L矩阵，脑得分
%   .Ly               ：N×(L×组数)矩阵，行为/设计得分
%   .LC_img_loadings  ：M×L矩阵，corr(Lx,X)
%   .LC_behav_loadings：B×(L×组数)矩阵，corr(Ly,Y)
%   .Sp_vect          ：排列检验零分布
%   .boot_results     ：结构体，引导检验结果，包含：
%       .Ub_vect    ：B×L×P矩阵，引导样本U
%       .Vb_vect    ：M×L×P矩阵，引导样本V
%       .Lxb        ：N×L×P矩阵，引导样本脑得分
%       .Lyb        ：N×(L×组数)×P矩阵，引导样本行为得分
%       .LC_img_loadings_boot ：M×L×P矩阵，引导样本影像载荷
%       .LC_behav_loadings_boot：B×(L×组数)×P矩阵，引导样本行为载荷
%       .*_mean     ：各统计量均值
%       .*_std      ：各统计量标准差
%       .*_lB       ：95%置信区间下限
%       .*_uB       ：95%置信区间上限
%       .U_stad     ：U的标准误差比
%       .V_stad     ：V的标准误差比
%
% =============================
% 示例：
% res = myPLS_analysis(input, pls_opts);
%
% =============================
% 历史：
% 2024-06-09  初版，标准化注释与格式，详细参数说明。
% =============================

%% ========== 常量与输入信息 ========== %%
% 受试者数量
nSubj = size(input.brain_data, 1); % N
% 行为/设计变量数量
nBehav = size(input.behav_data, 2); % B
% 成像变量数量
nImg = size(input.brain_data, 2);   % M

%% ========== 获取成像矩阵X0 ========== %%
X0 = input.brain_data; % 原始影像数据

%% ========== 生成行为/对比/交互矩阵Y0 ========== %%
[Y0, ~, nDesignScores] = myPLS_getY(pls_opts.behav_type, input.behav_data, input.grouping, input.group_names, input.behav_names); % 生成设计矩阵

disp('... 输入数据信息 ...')
disp(['受试者数量: ' num2str(nSubj)]);
disp(['成像变量数量 (体素/连接): ' num2str(nImg)]);
disp(['设计变量数量 (行为/对比): ' num2str(nDesignScores)]);
disp(' ')

%% ========== 归一化输入数据矩阵X和Y ========== %%
% (如果Y中有组对比，则组不会被考虑)
[X, ~, ~] = myPLS_norm(X0, input.grouping, pls_opts.normalization_img); % 影像数据归一化
[Y, ~, ~] = myPLS_norm(Y0, input.grouping, pls_opts.normalization_behav); % 行为数据归一化

%% ========== 交叉协方差矩阵 ========== %%
R = myPLS_cov(X, Y, input.grouping, pls_opts.grouped_PLS); % 计算交叉协方差

%% ========== 奇异值分解 ========== %%
[U, S, V] = svd(R, 'econ'); % SVD分解

% 潜在成分数（LCs）
nLC = min(size(S)); 

%% ICA惯例：将LCs转换为最大值为正
for iLC = 1:nLC
    [~, maxID] = max(abs(V(:, iLC)));
    if sign(V(maxID, iLC)) < 0
        V(:, iLC) = -V(:, iLC);
        U(:, iLC) = -U(:, iLC);
    end
end

% 每个LC解释的协方差仅用于完整PLSC报告，CPM后续不使用。
% explCovLC = (diag(S).^2) / sum(diag(S.^2));

%% ========== CPM流程不需要PLS得分、载荷和排列检验，故关闭这些计算 ========== %%
% [Lx, Ly, corr_Lx_X, corr_Ly_Y, corr_Lx_Y, corr_Ly_X] = myPLS_get_PLS_scores_loadings(X, Y, V, U, input.grouping, pls_opts);
% corr_Lx_Ly = diag(corr(Lx, Ly));
% [Sp_vect, corr_Lx_Ly_vect] = myPLS_permutations_xiu(X, Y, U, V, input.grouping, pls_opts);
% LC_pvals = myPLS_get_LC_pvals(Sp_vect, S, pls_opts);
% corr_Lx_Ly_pvals = myPLS_get_corr_pvals(corr_Lx_Ly_vect, corr_Lx_Ly, pls_opts);

%% ========== 引导测试以确定PLS载荷稳定性 ========== %%
boot_results = myPLS_bootstrapping(X0, Y0, U, V, S, input.grouping, pls_opts); % 引导检验

%% ========== 汇总所有结果变量到结构体 ========== %%
res.X = X;           % 归一化影像数据
res.Y = Y;           % 归一化行为数据
res.V = V;           % 影像显著性；PLSC_CPMmask中Inf/NaN修正需要使用
% CPM后续不使用以下完整PLSC报告字段，故不再保存。
% res.X0 = X0;
% res.meanX0 = meanX0;
% res.stdX0 = stdX0;
% res.Y0 = Y0;
% res.meanY0 = meanY0;
% res.stdY0 = stdY0;
% res.design_names = design_names;
% res.grouping = input.grouping;
% res.group_names = input.group_names;
% if isfield(input, 'img_names')
%     res.img_names = input.img_names;
% end
% res.R = R;
% res.U = U;
% res.S = S;
% res.explCovLC = explCovLC;
% res.LC_pvals = LC_pvals;
% res.Lx = Lx;
% res.Ly = Ly;
% res.corr_Lx_Ly = corr_Lx_Ly;
% res.corr_Lx_Ly_pvals = corr_Lx_Ly_pvals;
% res.LC_img_loadings = corr_Lx_X;
% res.LC_behav_loadings = corr_Ly_Y;
% res.Sp_vect = Sp_vect;
res.boot_results = boot_results; % 引导检验结果

end
