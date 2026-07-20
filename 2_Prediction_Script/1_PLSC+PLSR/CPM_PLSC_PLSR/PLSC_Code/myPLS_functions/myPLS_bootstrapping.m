function boot_results = myPLS_bootstrapping(X0, Y0, U, V, S, grouping, pls_opts)
% =============================
% 功能：
% 对X和Y进行有放回引导抽样（bootstrap），可选择是否考虑分组。每个引导样本均重新计算PLS分析，输出各类统计量及得分。
% 支持组内/组间抽样、Procrustes旋转、保存引导样本等多种选项。
% 
% =============================
% 输入：
% X0             ：N×M矩阵，原始影像数据（未归一化）
%                  N：受试者数量
%                  M：影像变量数量
% Y0             ：N×B矩阵，原始行为/设计数据（未归一化）
%                  B：行为/设计变量数量
% U              ：B×L矩阵，行为/设计显著性
%                  L：潜在成分（LCs）数量
% V              ：M×L矩阵，影像显著性
% S              ：L×L对角矩阵，奇异值
% grouping       ：N×1向量，受试者分组信息
%                  例如：[1,1,2] 表示前两名为组1，第三名为组2
% pls_opts       ：结构体，PLS分析参数，需包含如下字段：
%   .nBootstraps         ：整数，引导样本数（如1000）
%   .grouped_PLS         ：逻辑值/0或1，是否分组PLS
%                         0 = 在所有受试者上计算PLS
%                         1 = 组间协方差矩阵拼接（传统行为PLS）
%   .grouped_boot        ：逻辑值/0或1，引导抽样是否组内抽样
%                         0 = 忽略组别，整体抽样
%                         1 = 组内独立抽样
%   .boot_procrustes_mod ：整数，Procrustes旋转模式
%                         1 = 仅在U上计算旋转（只旋转U）
%                         2 = U和V分别旋转，取平均旋转
%   .save_boot_resampling：逻辑值/0或1，是否保存引导样本原始向量
%                         0 = 不保存
%                         1 = 保存
%   .normalization_img   ：整数，影像数据归一化方式
%                         0 = 不归一化
%                         1 = 跨所有受试者z评分（均值为0，方差为1）
%                         2 = 组内z评分（默认，每组单独z评分）
%                         3 = 跨所有受试者标准归一化（无中心化，方差为1）
%                         4 = 组内标准归一化（无中心化，方差为1）
%   .normalization_behav ：整数，行为/设计数据归一化方式
%                         0 = 不归一化
%                         1 = 跨所有受试者z评分（均值为0，方差为1）
%                         2 = 组内z评分（默认，每组单独z评分）
%                         3 = 跨所有受试者标准归一化（无中心化，方差为1）
%                         4 = 组内标准归一化（无中心化，方差为1）
% 
% =============================
% 输出：
% boot_results ：结构体，包含所有引导抽样结果及统计量：
%   .Ub_vect    ：B×L×P矩阵，bootstrap行为显著性向量
%                 B：行为/设计变量数，L：潜在成分数，P：引导样本数
%   .Vb_vect    ：M×L×P矩阵，bootstrap影像显著性向量
%                 M：影像变量数，L：潜在成分数，P：引导样本数
%   .Lxb        ：N×L×P矩阵，bootstrap影像得分
%                 N：受试者数，L：潜在成分数，P：引导样本数
%   .Lyb        ：N×(L×组数)×P矩阵，bootstrap行为得分
%                 N：受试者数，L：潜在成分数，组数：分组数量，P：引导样本数
%   .LC_img_loadings_boot ：M×L×P矩阵，bootstrap影像载荷
%   .LC_behav_loadings_boot：B×(L×组数)×P矩阵，bootstrap行为载荷
%   .*_mean     ：各统计量均值（如Ub_mean, Vb_mean等）
%   .*_std      ：各统计量标准差
%   .*_lB       ：95%置信区间下限
%   .*_uB       ：95%置信区间上限
%   .U_stad     ：U的标准误差比（U ./ 标准误）
%   .V_stad     ：V的标准误差比（V ./ 标准误）
% 
% =============================
% 示例：
% boot_results = myPLS_bootstrapping(X0, Y0, U, V, S, grouping, pls_opts);
% 
% =============================
% 历史：
% 2024-06-09  初版，标准化注释与格式，详细参数说明。
% =============================

%% ========== 设置随机数生成器，保证可复现 ========== %%
rng(1);

%% ========== 检查输入维度 ========== %%
if(size(X0, 1) ~= size(Y0, 1))
    error('输入参数X和Y的行数应相同');
end

%% ========== 引导抽样主循环 ========== %%
disp('... 引导抽样 ...') % 显示进程

% 计算引导抽样顺序（是否组内抽样由grouped_boot决定）
all_boot_orders = myPLS_get_boot_orders(pls_opts.nBootstraps, grouping, pls_opts.grouped_boot);

% 初始化存储变量
Ub_vect = [];
Vb_vect = [];
boot_results = struct();

for iB = 1:pls_opts.nBootstraps
    %% 每50次显示一次进度
    if mod(iB, 100) == 0
        disp(num2str(iB));
    end

    %% 重采样X和Y并归一化
    Xb = X0(all_boot_orders(:, iB), :); % 按照引导顺序采样
    Xb = myPLS_norm(Xb, grouping, pls_opts.normalization_img); % 归一化
    Yb = Y0(all_boot_orders(:, iB), :);
    Yb = myPLS_norm(Yb, grouping, pls_opts.normalization_behav);

    %% 生成重采样X和Y之间的交叉协方差矩阵
    Rb = myPLS_cov(Xb, Yb, grouping, pls_opts.grouped_PLS);

    %% Rb的奇异值分解
    [Ub, Sb, Vb] = svd(Rb, 'econ');

    %% Procrustes变换（校正轴旋转/反射）
    switch pls_opts.boot_procrustes_mod
        case 1
            % 仅在U上计算旋转
            rotatemat_full = rri_bootprocrust(U, Ub);
            
            % 旋转和重新缩放Ub和Vb
            Vb = Vb * Sb * rotatemat_full;
            Ub = Ub * Sb * rotatemat_full;
            Vb = Vb./repmat(diag(S)', size(Vb, 1), 1);
            Ub = Ub./repmat(diag(S)', size(Ub, 1), 1);
        case 2
            % U和V分别旋转，取平均
            rotatemat1 = rri_bootprocrust(U, Ub);
            rotatemat2 = rri_bootprocrust(V, Vb);
            rotatemat_full = (rotatemat1 + rotatemat2)/2;
            Vb = Vb * rotatemat_full;
            Ub = Ub * rotatemat_full;
        otherwise
            error('pls_opts.boot_procrustes_mod中存在无效值！');
    end

    %% 存储所有引导样本的向量
    Ub_vect(:, :, iB) = Ub;
    Vb_vect(:, :, iB) = Vb;

    %% CPM流程只使用V_stad；跳过bootstrap得分、载荷和相关矩阵
    % [Lxb, Lyb, corr_Lxb_Xb, corr_Lyb_Yb, corr_Lxb_Yb, corr_Lyb_Xb] = myPLS_get_PLS_scores_loadings(Xb, Yb, Vb, Ub, grouping, pls_opts);
    % boot_results.Lxb(:, :, iB) = Lxb;
    % boot_results.Lyb(:, :, iB) = Lyb;
    % boot_results.LC_img_loadings_boot(:, :, iB) = corr_Lxb_Xb;
    % boot_results.LC_behav_loadings_boot(:, :, iB) = corr_Lyb_Yb;
    % boot_results.corr_Lxb_Yb(:, :, iB) = corr_Lxb_Yb;
    % boot_results.corr_Lyb_Xb(:, :, iB) = corr_Lyb_Xb;
    % boot_results.corr_Lxb_Lyb(:, :, iB) = corr(Lxb, Lyb)';
end

%% ========== 计算引导统计量 ========== %%
boot_stats = myPLS_bootstrap_stats(Ub_vect, Vb_vect);

% 将所有统计字段保存到boot_results中，便于后续扩展
fN = fieldnames(boot_stats);
for iF = 1:length(fN)
    boot_results.(fN{iF}) = boot_stats.(fN{iF});
end

% 如果要求保存引导抽样数据，则保存
if pls_opts.save_boot_resampling
    boot_results.Ub_vect = Ub_vect;
    boot_results.Vb_vect = Vb_vect;
end


%% ========== U和V稳定性 ========== %%
boot_results.U_stad = U ./ sqrt(sum((Ub_vect - boot_stats.Ub_mean).^2, 3) ./ (pls_opts.nBootstraps-1));
boot_results.V_stad = V ./ sqrt(sum((Vb_vect - boot_stats.Vb_mean).^2, 3) ./ (pls_opts.nBootstraps-1));

disp(' ')

end
