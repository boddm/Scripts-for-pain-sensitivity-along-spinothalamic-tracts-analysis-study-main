function res = myPLS_analysis(input, pls_opts)
% =============================
% Function:
% Main function to perform PLS analysis. Includes data normalization, covariance matrix computation,
% singular value decomposition, PLS scores and loadings calculation, permutation test and bootstrap
% test, and outputs all analysis results. Supports grouping, behavior/contrast/interaction analysis,
% various normalization methods and statistical tests.
%
% =============================
% Input:
% input    : Structure containing analysis input data:
%   .brain_data   : N×M matrix, raw imaging data
%                   N: number of subjects
%                   M: number of imaging variables
%   .behav_data   : N×B matrix, raw behavior/design data
%                   B: number of behavior/design variables
%   .grouping     : N×1 vector, subject grouping information
%                   e.g., [1,1,2] means first two subjects are group 1,
%                   third subject is group 2
%   .group_names  : 1×G cell, group names (optional, G is number of groups)
%   .behav_names  : 1×B cell, behavior variable names (optional)
%   .img_names    : 1×M cell, imaging variable names (optional)
%
% pls_opts : Structure, PLS analysis parameters, must contain the following fields:
%   .nPerms              : Integer, number of permutation tests
%   .nBootstraps         : Integer, number of bootstrap samples
%   .normalization_img   : Integer, imaging data normalization method
%                         0 = no normalization
%                         1 = z-score across all subjects (mean 0, variance 1)
%                         2 = z-score within groups (z-score per group, default)
%                         3 = standard normalization across all subjects (no centering, variance 1)
%                         4 = standard normalization within groups (no centering, variance 1)
%   .normalization_behav : Integer, behavior/design data normalization method
%                         0 = no normalization
%                         1 = z-score across all subjects (mean 0, variance 1)
%                         2 = z-score within groups (z-score per group, default)
%                         3 = standard normalization across all subjects (no centering, variance 1)
%                         4 = standard normalization within groups (no centering, variance 1)
%   .grouped_PLS         : Logical/0 or 1, whether to use grouped PLS
%                         0 = compute PLS across all subjects
%                         1 = concatenate covariance matrices between groups (behavior PLS)
%   .grouped_perm        : Logical/0 or 1, whether permutations are within groups
%                         0 = ignore grouping for permutations
%                         1 = permute within groups
%   .grouped_boot        : Logical/0 or 1, whether bootstrap sampling is within groups
%                         0 = ignore grouping for bootstrap
%                         1 = bootstrap within groups
%   .boot_procrustes_mod : Integer, Procrustes rotation mode during bootstrap
%                         1 = rotate U only (default)
%                         2 = average rotation of U and V
%   .save_boot_resampling: Logical/0 or 1, whether to save bootstrap resampling data
%                         0 = do not save
%                         1 = save
%   .behav_type          : String, type of behavior analysis
%                         'behavior': standard behavior PLS (default)
%                         'contrast': compute contrast between two groups only
%                         'contrastBehav': combine contrast and behavior measures
%                         'contrastBehavInteract': consider group-behavior interaction effects
%
% =============================
% Output:
% res : Structure containing all PLS analysis results:
%   .X0, .Y0          : Unnormalized input matrices
%   .meanX0, .stdX0   : Mean and standard deviation of X0
%   .meanY0, .stdY0   : Mean and standard deviation of Y0
%   .X, .Y            : Normalized input matrices
%   .design_names     : Design variable names (e.g., behavior/contrast/interaction)
%   .grouping         : Grouping information
%   .group_names      : Group names
%   .img_names        : Imaging variable names (if available)
%   .R                : L×M matrix, brain-behavior covariance matrix
%   .U                : B×L matrix, behavior salience
%   .V                : M×L matrix, imaging salience
%   .S                : L×L diagonal matrix, singular values
%   .explCovLC        : Proportion of covariance explained by each LC
%   .LC_pvals         : p-values for each LC (permutation test)
%   .Lx               : N×L matrix, brain scores
%   .Ly               : N×(L×number of groups) matrix, behavior/design scores
%   .LC_img_loadings  : M×L matrix, corr(Lx,X)
%   .LC_behav_loadings: B×(L×number of groups) matrix, corr(Ly,Y)
%   .Sp_vect          : Null distribution from permutation test
%   .boot_results     : Structure, bootstrap test results, containing:
%       .Ub_vect    : B×L×P matrix, bootstrap samples of U
%       .Vb_vect    : M×L×P matrix, bootstrap samples of V
%       .Lxb        : N×L×P matrix, bootstrap samples of brain scores
%       .Lyb        : N×(L×number of groups)×P matrix, bootstrap samples of behavior scores
%       .LC_img_loadings_boot : M×L×P matrix, bootstrap samples of imaging loadings
%       .LC_behav_loadings_boot: B×(L×number of groups)×P matrix, bootstrap samples of behavior loadings
%       .*_mean     : Mean of each statistic
%       .*_std      : Standard deviation of each statistic
%       .*_lB       : Lower bound of 95% confidence interval
%       .*_uB       : Upper bound of 95% confidence interval
%       .U_stad     : Standard error ratio of U
%       .V_stad     : Standard error ratio of V
%
% =============================
% Example:
% res = myPLS_analysis(input, pls_opts);
%
% =============================

%% ========== Constants and Input Information ========== %%
% Number of subjects
nSubj = size(input.brain_data, 1); % N
% Number of behavior/design variables
nBehav = size(input.behav_data, 2); % B
% Number of imaging variables
nImg = size(input.brain_data, 2);   % M

%% ========== Imaging Matrix X0 ========== %%
X0 = input.brain_data; % Original imaging data

%% ========== Design Matrix Y0 ========== %%
[Y0, ~, nDesignScores] = myPLS_getY(pls_opts.behav_type, input.behav_data, input.grouping, input.group_names, input.behav_names); % Generate design matrix

disp('... Input Data Information ...')
disp(['Number of subjects: ' num2str(nSubj)]);
disp(['Number of imaging variables (voxels/connections): ' num2str(nImg)]);
disp(['Number of design variables (behavior/contrast): ' num2str(nDesignScores)]);
disp(' ')

%% ========== Normalize Input Data Matrices X and Y ========== %%
% (如果Y中有组对比，则组不会被考虑)
[X, ~, ~] = myPLS_norm(X0, input.grouping, pls_opts.normalization_img); % Normalize imaging data
[Y, ~, ~] = myPLS_norm(Y0, input.grouping, pls_opts.normalization_behav); % Normalize behavior data

%% ========== Cross-Covariance Matrix ========== %%
R = myPLS_cov(X, Y, input.grouping, pls_opts.grouped_PLS); % Calculate cross covariance matrix

%% ========== Singular Value Decomposition ========== %%
[U, S, V] = svd(R, 'econ'); % SVD decomposition

% 潜在成分数（LCs）
nLC = min(size(S));

%% ICA Convention: Convert LCs to Positive Values
for iLC = 1:nLC
    [~, maxID] = max(abs(V(:, iLC)));
    if sign(V(maxID, iLC)) < 0
        V(:, iLC) = -V(:, iLC);
        U(:, iLC) = -U(:, iLC);
    end
end

% explCovLC = (diag(S).^2) / sum(diag(S.^2));

%% ========== Skip PLS scores, loadings and permutation tests ========== %%
% [Lx, Ly, corr_Lx_X, corr_Ly_Y, corr_Lx_Y, corr_Ly_X] = myPLS_get_PLS_scores_loadings(X, Y, V, U, input.grouping, pls_opts);
% corr_Lx_Ly = diag(corr(Lx, Ly));
% [Sp_vect, corr_Lx_Ly_vect] = myPLS_permutations(X, Y, U, V, input.grouping, pls_opts);
% LC_pvals = myPLS_get_LC_pvals(Sp_vect, S, pls_opts);
% corr_Lx_Ly_pvals = myPLS_get_corr_pvals(corr_Lx_Ly_vect, corr_Lx_Ly, pls_opts);

%% ========== Bootstrap Test to Determine PLS Loadings Stability and Significance ========== %%
boot_results = myPLS_bootstrapping(X0, Y0, U, V, S, input.grouping, pls_opts); % 引导检验

%% ========== Aggregate All Result Variables to Structure ========== %%
res.X = X;           % Normalized imaging data
res.Y = Y;           % Normalized behavior data
res.V = V;           % Imaging salience; needed for Inf/NaN correction in PLSC_CPMmask
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
res.boot_results = boot_results; % Bootstrap test results

end
