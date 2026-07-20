function res = myPLS_analysis(input, pls_opts)
% =============================
% Purpose:
% Main routine for PLS analysis. Performs data normalization, covariance-matrix computation, singular-value decomposition, PLS score & loading calculation, permutation and bootstrap testing, and returns all results.
% Supports grouping, behavior/contrast/interaction analysis, various normalizations and statistical tests.
%
% =============================
% Inputs:
% input    : structure containing analysis data:
%   .brain_data   : N×M matrix, raw imaging data
%                   N: number of subjects
%                   M: number of imaging variables
%   .behav_data   : N×B matrix, raw behavior/design data
%                   B: number of behavior/design variables
%   .grouping     : N×1 vector, subject group labels
%                   e.g. [1,1,2] means first two subjects in group 1, third in group 2
%   .group_names  : 1×G cell, group names (optional, G = number of groups)
%   .behav_names  : 1×B cell, behavior variable names (optional)
%   .img_names    : 1×M cell, imaging variable names (optional)
%
% pls_opts : structure, PLS options, must contain:
%   .nPerms              : integer, number of permutations
%   .nBootstraps         : integer, number of bootstrap samples
%   .normalization_img   : integer, imaging-data normalization
%                         0 = none
%                         1 = z-score across all subjects (mean 0, var 1)
%                         2 = z-score within each group (default)
%                         3 = standard normalization across all (no centering, var 1)
%                         4 = standard normalization within groups (no centering, var 1)
%   .normalization_behav : integer, behavior/design-data normalization
%                         0 = none
%                         1 = z-score across all subjects (mean 0, var 1)
%                         2 = z-score within each group (default)
%                         3 = standard normalization across all (no centering, var 1)
%                         4 = standard normalization within groups (no centering, var 1)
%   .grouped_PLS         : logical/0 or 1, perform grouped PLS
%                         0 = compute PLS on all subjects together
%                         1 = concatenate group-wise covariance matrices (behavior PLS)
%   .grouped_perm        : logical/0 or 1, permute within groups
%                         0 = ignore groups during permutation
%                         1 = permute within each group
%   .grouped_boot        : logical/0 or 1, bootstrap within groups
%                         0 = ignore groups during bootstrap
%                         1 = bootstrap within each group
%   .boot_procrustes_mod : integer, Procrustes rotation mode during bootstrap
%                         1 = rotate U only (default)
%                         2 = rotate both U and V
%   .save_boot_resampling: logical/0 or 1, save bootstrap resamples
%                         0 = do not save
%                         1 = save
%   .behav_type          : string, type of behavior analysis
%                         'behavior'        : standard behavior PLS (default)
%                         'contrast'        : two-group contrast only
%                         'contrastBehav'   : contrast plus behavior measure
%                         'contrastBehavInteract': group-behavior interaction
%
% =============================
% Outputs:
% res : structure containing all PLS results:
%   .X0, .Y0          : un-normalized input matrices
%   .meanX0, .stdX0   : mean and std of X0
%   .meanY0, .stdY0   : mean and std of Y0
%   .X, .Y            : normalized input matrices
%   .design_names     : design variable names (behavior/contrast/interaction)
%   .grouping         : group labels
%   .group_names      : group names
%   .img_names        : imaging variable names (if provided)
%   .R                : L×M matrix, brain-behavior covariance matrix
%   .U                : B×L matrix, behavior saliences
%   .V                : M×L matrix, imaging saliences
%   .S                : L×L diagonal matrix, singular values
%   .explCovLC        : proportion of covariance explained by each LC
%   .LC_pvals         : p-value for each LC (permutation test)
%   .Lx               : N×L matrix, brain scores
%   .Ly               : N×(L×#groups) matrix, behavior/design scores
%   .LC_img_loadings  : M×L matrix, corr(Lx,X)
%   .LC_behav_loadings: B×(L×#groups) matrix, corr(Ly,Y)
%   .Sp_vect          : permutation null distribution
%   .boot_results     : structure, bootstrap results, containing:
%       .Ub_vect    : B×L×P matrix, bootstrap U
%       .Vb_vect    : M×L×P matrix, bootstrap V
%       .Lxb        : N×L×P matrix, bootstrap brain scores
%       .Lyb        : N×(L×#groups)×P matrix, bootstrap behavior scores
%       .LC_img_loadings_boot : M×L×P matrix, bootstrap imaging loadings
%       .LC_behav_loadings_boot: B×(L×#groups)×P matrix, bootstrap behavior loadings
%       .*_mean     : mean of each statistic
%       .*_std      : std of each statistic
%       .*_lB       : 95% CI lower bound
%       .*_uB       : 95% CI upper bound
%       .U_stad     : standard error ratio for U
%       .V_stad     : standard error ratio for V
%
% =============================
% Example:
% res = myPLS_analysis(input, pls_opts);
% =============================

%% ========== Constants and input information ========== %%
% Number of subjects
nSubj = size(input.brain_data, 1); % N
% Number of behavior/design variables
nBehav = size(input.behav_data, 2); % B
% Number of imaging variables
nImg = size(input.brain_data, 2);   % M

%% ========== Acquire imaging matrix X0 ========== %%
X0 = input.brain_data; % Raw imaging data

%% ========== Generate behavior/contrast/interaction matrix Y0 ========== %%
[Y0, design_names, nDesignScores] = myPLS_getY(pls_opts.behav_type, input.behav_data, input.grouping, input.group_names, input.behav_names); % Generate design matrix

disp('... Input data information ...')
disp(['Number of subjects: ' num2str(nSubj)]);
disp(['Number of imaging variables (voxels/connections): ' num2str(nImg)]);
disp(['Number of design variables (behavior/contrast): ' num2str(nDesignScores)]);
disp(' ')

%% ========== Normalize input data matrices X and Y ========== %%
% (Groups are ignored if Y contains group contrasts)
[X, meanX0, stdX0] = myPLS_norm(X0, input.grouping, pls_opts.normalization_img); % Normalize imaging data
[Y, meanY0, stdY0] = myPLS_norm(Y0, input.grouping, pls_opts.normalization_behav); % Normalize behavior data

%% ========== Cross-covariance matrix ========== %%
R = myPLS_cov(X, Y, input.grouping, pls_opts.grouped_PLS); % Compute cross-covariance

%% ========== Singular value decomposition ========== %%
[U, S, V] = svd(R, 'econ'); % SVD decomposition

% Number of latent components (LCs)
nLC = min(size(S)); 

%% ICA convention: make LCs positive at maximum value
for iLC = 1:nLC
    [~, maxID] = max(abs(V(:, iLC)));
    if sign(V(maxID, iLC)) < 0
        V(:, iLC) = -V(:, iLC);
        U(:, iLC) = -U(:, iLC);
    end
end

% Covariance explained by each LC
explCovLC = (diag(S).^2) / sum(diag(S.^2));

%% ========== Compute PLS scores and loadings ========== %%
[Lx, Ly, corr_Lx_X, corr_Ly_Y, corr_Lx_Y, corr_Ly_X] = myPLS_get_PLS_scores_loadings(X, Y, V, U, input.grouping, pls_opts); % scores & loadings

%% ========== Permutation test to assess LC significance ========== %%
Sp_vect = myPLS_permutations(X, Y, U, input.grouping, pls_opts); % permutation test

% Compute p-values from permutation null distribution
LC_pvals = myPLS_get_LC_pvals(Sp_vect, S, pls_opts); % calculate p-values

%% ========== Bootstrapping to test stability of PLS loadings ========== %%
boot_results = myPLS_bootstrapping(X0, Y0, U, V, S, input.grouping, pls_opts); % bootstrap test

%% ========== Save all result variables in struct ========== %%
res.X0 = X0; % raw imaging data
res.meanX0 = meanX0; % imaging mean
res.stdX0 = stdX0;   % imaging std
res.Y0 = Y0;         % raw behavior data
res.meanY0 = meanY0; % behavior mean
res.stdY0 = stdY0;   % behavior std
res.X = X;           % normalized imaging data
res.Y = Y;           % normalized behavior data
res.design_names = design_names; % design variable names
res.grouping = input.grouping;   % grouping info
res.group_names = input.group_names; % group names
if isfield(input, 'img_names')
    res.img_names = input.img_names;
end
res.R = R;           % cross-covariance matrix
res.U = U;           % behavior saliences
res.S = S;           % singular values
res.V = V;           % imaging saliences
res.explCovLC = explCovLC; % covariance explained by each LC
res.LC_pvals = LC_pvals;   % p-values for LCs
res.Lx = Lx;         % brain scores
res.Ly = Ly;         % behavior/design scores
res.LC_img_loadings = corr_Lx_X; % imaging loadings
res.LC_behav_loadings = corr_Ly_Y; % behavior loadings
res.Sp_vect = Sp_vect;     % permutation null distribution
res.LC_pvals = LC_pvals;   % p-values for LCs
res.boot_results = boot_results; % bootstrap test results

end
