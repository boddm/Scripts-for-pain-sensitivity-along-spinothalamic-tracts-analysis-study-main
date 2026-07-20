function boot_stats = myPLS_bootstrap_stats(Ub_vect, Vb_vect)

% =============================
% 功能：
% 计算PLS bootstrap重采样下的显著性和载荷的均值、标准差、置信区间等统计量。
% 输出所有bootstrap分布的统计描述，用于后续显著性判断和可视化。
%
% =============================
% 输入：
% Ub_vect      ：B×L×P矩阵，bootstrap行为/设计显著性
%                B：行为/设计变量数，L：潜在成分数，P：bootstrap样本数
% Vb_vect      ：M×L×P矩阵，bootstrap影像显著性
%                M：影像变量数，L：潜在成分数，P：bootstrap样本数
% boot_results ：结构体，包含bootstrap结果，需包含：
%   .Lxb                    ：N×L×P矩阵，bootstrap脑得分
%   .Lyb                    ：N×(L×组数)×P矩阵，bootstrap行为得分
%   .LC_img_loadings_boot   ：M×L×P矩阵，bootstrap影像载荷
%   .LC_behav_loadings_boot ：B×(L×组数)×P矩阵，bootstrap行为载荷
%
% =============================
% 输出：
% boot_stats ：结构体，包含所有bootstrap统计量：
%   .Ub_mean, .Ub_std, .Ub_lB, .Ub_uB                       ：B×L矩阵，行为/设计显著性均值、标准差、置信区间
%   .Vb_mean, .Vb_std, .Vb_lB, .Vb_uB                       ：M×L矩阵，影像显著性均值、标准差、置信区间
%   .LC_behav_loadings_mean, .LC_behav_loadings_std         ：B×(L×组数)矩阵，行为/设计载荷均值、标准差
%   .LC_behav_loadings_lB, .LC_behav_loadings_uB            ：B×(L×组数)矩阵，行为/设计载荷置信区间
%   .LC_img_loadings_mean, .LC_img_loadings_std             ：M×L矩阵，影像载荷均值、标准差
%   .LC_img_loadings_lB, .LC_img_loadings_uB                ：M×L矩阵，影像载荷置信区间
%
% =============================
% 示例：
% boot_stats = myPLS_bootstrap_stats(Ub_vect, Vb_vect, boot_results);
%
% =============================
% 历史：
% 2024-06-09  初版，标准化注释与格式，详细参数说明。
% =============================

%% ========== 行为/设计和影像显著性统计 ========== %%
% CPM流程只使用V_stad，因此这里只保留计算V_stad所需的Ub/Vb统计量。
boot_stats.Ub_mean = mean(Ub_vect, 3); % 行为/设计显著性均值
boot_stats.Ub_std = std(Ub_vect, [], 3); % 行为/设计显著性标准差
boot_stats.Ub_lB = prctile(Ub_vect, 2.5, 3); % 95%置信区间下界
boot_stats.Ub_uB = prctile(Ub_vect, 97.5, 3); % 95%置信区间上界

boot_stats.Vb_mean = mean(Vb_vect, 3); % 影像显著性均值
boot_stats.Vb_std = std(Vb_vect,[], 3); % 影像显著性标准差
boot_stats.Vb_lB = prctile(Vb_vect, 2.5, 3); % 95%置信区间下界
boot_stats.Vb_uB = prctile(Vb_vect, 97.5, 3); % 95%置信区间上界

%% ========== 以下统计量仅用于完整PLSC绘图/报告，CPM后续不使用，已关闭 ========== %%
% boot_stats.LC_behav_loadings_mean = mean(boot_results.LC_behav_loadings_boot, 3);
% boot_stats.LC_behav_loadings_std = std(boot_results.LC_behav_loadings_boot, [], 3);
% boot_stats.LC_behav_loadings_lB = prctile(boot_results.LC_behav_loadings_boot, 2.5, 3);
% boot_stats.LC_behav_loadings_uB = prctile(boot_results.LC_behav_loadings_boot, 97.5, 3);
% boot_stats.LC_img_loadings_mean = mean(boot_results.LC_img_loadings_boot, 3);
% boot_stats.LC_img_loadings_std = std(boot_results.LC_img_loadings_boot, [], 3);
% boot_stats.LC_img_loadings_lB = prctile(boot_results.LC_img_loadings_boot, 2.5, 3);
% boot_stats.LC_img_loadings_uB = prctile(boot_results.LC_img_loadings_boot, 97.5, 3);
% boot_stats.corr_Lxb_Yb_mean = mean(boot_results.corr_Lxb_Yb, 3);
% boot_stats.corr_Lxb_Yb_std = std(boot_results.corr_Lxb_Yb, [], 3);
% boot_stats.corr_Lxb_Yb_lB = prctile(boot_results.corr_Lxb_Yb, 2.5, 3);
% boot_stats.corr_Lxb_Yb_uB = prctile(boot_results.corr_Lxb_Yb, 97.5, 3);
% boot_stats.corr_Lyb_Xb_mean = mean(boot_results.corr_Lyb_Xb, 3);
% boot_stats.corr_Lyb_Xb_std = std(boot_results.corr_Lyb_Xb, [], 3);
% boot_stats.corr_Lyb_Xb_lB = prctile(boot_results.corr_Lyb_Xb, 2.5, 3);
% boot_stats.corr_Lyb_Xb_uB = prctile(boot_results.corr_Lyb_Xb, 97.5, 3);
% boot_stats.corr_Lxb_Lyb_mean = mean(boot_results.corr_Lxb_Lyb, 3);
% boot_stats.corr_Lxb_Lyb_std = std(boot_results.corr_Lxb_Lyb, [], 3);
% boot_stats.corr_Lxb_Lyb_lB = prctile(boot_results.corr_Lxb_Lyb, 2.5, 3);
% boot_stats.corr_Lxb_Lyb_uB = prctile(boot_results.corr_Lxb_Lyb, 97.5, 3);
% boot_stats.corr_Lxb_Lyb_pvals = myPLS_get_corr_pvals(boot_results.corr_Lxb_Lyb, corr_Lx_Ly, pls_opts);

end
