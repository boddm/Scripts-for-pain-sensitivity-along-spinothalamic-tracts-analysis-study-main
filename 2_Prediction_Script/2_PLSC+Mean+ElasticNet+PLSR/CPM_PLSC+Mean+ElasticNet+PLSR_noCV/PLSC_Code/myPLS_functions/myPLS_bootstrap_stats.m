function boot_stats = myPLS_bootstrap_stats(Ub_vect, Vb_vect)

% =============================
% Purpose:
% Compute statistics (mean, std, confidence intervals) for PLS bootstrap
% significance and loadings. Outputs descriptive statistics of all bootstrap
% distributions for subsequent significance testing and visualization.
%
% =============================
% Inputs:
% Ub_vect      : B×L×P matrix, bootstrap behavioral/design salience
%                B: number of behavioral/design variables, L: number of latent components, P: number of bootstrap samples
% Vb_vect      : M×L×P matrix, bootstrap imaging salience
%                M: number of imaging variables, L: number of latent components, P: number of bootstrap samples
% boot_results : structure containing bootstrap results, must include:
%   .Lxb                    : N×L×P matrix, bootstrap brain scores
%   .Lyb                    : N×(L×num_groups)×P matrix, bootstrap behavioral scores
%   .LC_img_loadings_boot   : M×L×P matrix, bootstrap imaging loadings
%   .LC_behav_loadings_boot : B×(L×num_groups)×P matrix, bootstrap behavioral loadings
%
% =============================
% Outputs:
% boot_stats : structure containing all bootstrap statistics:
%   .Ub_mean, .Ub_std, .Ub_lB, .Ub_uB                       : B×L matrices, behavioral/design salience mean, std, confidence intervals
%   .Vb_mean, .Vb_std, .Vb_lB, .Vb_uB                       : M×L matrices, imaging salience mean, std, confidence intervals
%   .LC_behav_loadings_mean, .LC_behav_loadings_std         : B×(L×num_groups) matrices, behavioral/design loadings mean, std
%   .LC_behav_loadings_lB, .LC_behav_loadings_uB            : B×(L×num_groups) matrices, behavioral/design loadings confidence intervals
%   .LC_img_loadings_mean, .LC_img_loadings_std             : M×L matrices, imaging loadings mean, std
%   .LC_img_loadings_lB, .LC_img_loadings_uB                : M×L matrices, imaging loadings confidence intervals
%
% =============================
% Example:
% boot_stats = myPLS_bootstrap_stats(Ub_vect, Vb_vect, boot_results);
%
% =============================

%% ========== Behavioral/Design and Imaging Saliency Statistics ========== %%
% CPM process only uses V_stad, so here we only retain Ub/Vb statistics needed for V_stad calculation.
boot_stats.Ub_mean = mean(Ub_vect, 3); % behavioral/design saliency mean
boot_stats.Ub_std = std(Ub_vect, [], 3); % behavioral/design saliency standard deviation
boot_stats.Ub_lB = prctile(Ub_vect, 2.5, 3); % 95% confidence interval lower bound
boot_stats.Ub_uB = prctile(Ub_vect, 97.5, 3); % 95% confidence interval upper bound

boot_stats.Vb_mean = mean(Vb_vect, 3); % imaging saliency mean
boot_stats.Vb_std = std(Vb_vect,[], 3); % imaging saliency standard deviation
boot_stats.Vb_lB = prctile(Vb_vect, 2.5, 3); % 95% confidence interval lower bound
boot_stats.Vb_uB = prctile(Vb_vect, 97.5, 3); % 95% confidence interval upper bound

%% ========== The following statistics are only used for full PLSC plotting/reporting ========== %%
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
