function boot_results = myPLS_bootstrapping(X0, Y0, U, V, S, grouping, pls_opts)
% =============================
% Function:
% Performs bootstrap resampling with replacement on X and Y, with optional grouping.
% Each bootstrap sample recomputes PLS analysis, outputting various statistics and scores.
% Supports within-group/between-group resampling, Procrustes rotation, saving bootstrap samples, and more.
%
% =============================
% Inputs:
% X0             : N×M matrix, raw imaging data (un-normalized)
%                  N: number of subjects
%                  M: number of imaging variables
% Y0             : N×B matrix, raw behavior/design data (un-normalized)
%                  B: number of behavior/design variables
% U              : B×L matrix, behavior/design salience
%                  L: number of latent components (LCs)
% V              : M×L matrix, imaging salience
% S              : L×L diagonal matrix, singular values
% grouping       : N×1 vector, subject group information
%                  e.g., [1,1,2] means first two subjects are group 1, third is group 2
% pls_opts       : struct, PLS analysis parameters, must contain the following fields:
%   .nBootstraps         : integer, number of bootstrap samples (e.g., 1000)
%   .grouped_PLS         : logical/0 or 1, whether to perform grouped PLS
%                         0 = compute PLS across all subjects
%                         1 = concatenate covariance matrices across groups (traditional behavioral PLS)
%   .grouped_boot        : logical/0 or 1, whether to resample within groups
%                         0 = ignore groups, resample from all subjects
%                         1 = independent resampling within each group
%   .boot_procrustes_mod : integer, Procrustes rotation mode
%                         1 = compute rotation on U only (rotate U only)
%                         2 = rotate U and V separately, then average
%   .save_boot_resampling: logical/0 or 1, whether to save raw bootstrap vectors
%                         0 = do not save
%                         1 = save
%   .normalization_img   : integer, imaging data normalization method
%                         0 = no normalization
%                         1 = z-score across all subjects (mean=0, variance=1)
%                         2 = z-score within groups (default, z-score per group)
%                         3 = standard normalization across all subjects (no centering, variance=1)
%                         4 = standard normalization within groups (no centering, variance=1)
%   .normalization_behav : integer, behavior/design data normalization method
%                         0 = no normalization
%                         1 = z-score across all subjects (mean=0, variance=1)
%                         2 = z-score within groups (default, z-score per group)
%                         3 = standard normalization across all subjects (no centering, variance=1)
%                         4 = standard normalization within groups (no centering, variance=1)
%
% =============================
% Outputs:
% boot_results : struct, contains all bootstrap results and statistics:
%   .Ub_vect    : B×L×P matrix, bootstrap behavior salience vectors
%                 B: number of behavior/design variables, L: number of LCs, P: number of bootstrap samples
%   .Vb_vect    : M×L×P matrix, bootstrap imaging salience vectors
%                 M: number of imaging variables, L: number of LCs, P: number of bootstrap samples
%   .Lxb        : N×L×P matrix, bootstrap imaging scores
%                 N: number of subjects, L: number of LCs, P: number of bootstrap samples
%   .Lyb        : N×(L×num_groups)×P matrix, bootstrap behavior scores
%                 N: number of subjects, L: number of LCs, num_groups: number of groups, P: number of bootstrap samples
%   .LC_img_loadings_boot : M×L×P matrix, bootstrap imaging loadings
%   .LC_behav_loadings_boot: B×(L×num_groups)×P matrix, bootstrap behavior loadings
%   .*_mean     : mean of each statistic (e.g., Ub_mean, Vb_mean, etc.)
%   .*_std      : standard deviation of each statistic
%   .*_lB       : lower bound of 95% confidence interval
%   .*_uB       : upper bound of 95% confidence interval
%   .U_stad     : ratio of U to standard error (U ./ standard error)
%   .V_stad     : ratio of V to standard error (V ./ standard error)
%
% =============================
% Example:
% boot_results = myPLS_bootstrapping(X0, Y0, U, V, S, grouping, pls_opts);
%
% =============================

%% ========== Set random number generator for reproducibility ========== %%
rng(1);

%% ========== Check input dimensions ========== %%
if(size(X0, 1) ~= size(Y0, 1))
    error('Input parameters X and Y must have the same number of rows');
end

%% ========== Bootstrap main loop ========== %%
disp('... Bootstrapping ...') % Display progress

% Compute bootstrap resampling order (grouped_boot determines whether to resample within groups)
all_boot_orders = myPLS_get_boot_orders(pls_opts.nBootstraps, grouping, pls_opts.grouped_boot);

% Initialize storage variables
Ub_vect = [];
Vb_vect = [];
boot_results = struct();

for iB = 1:pls_opts.nBootstraps
    if mod(iB, 100) == 0
        disp(num2str(iB));
    end

    %% Resample X and Y and normalize
    Xb = X0(all_boot_orders(:, iB), :); % Bootstrap resampling
    Xb = myPLS_norm(Xb, grouping, pls_opts.normalization_img); % Normalization
    Yb = Y0(all_boot_orders(:, iB), :);
    Yb = myPLS_norm(Yb, grouping, pls_opts.normalization_behav);

    %% Generate cross-covariance matrix between resampled X and Y
    Rb = myPLS_cov(Xb, Yb, grouping, pls_opts.grouped_PLS);

    %% Singular value decomposition of Rb
    [Ub, Sb, Vb] = svd(Rb, 'econ');

    %% Procrustes transformation (correcting axis rotation/reflection)
    switch pls_opts.boot_procrustes_mod
        case 1
            % Compute rotation on U only
            rotatemat_full = rri_bootprocrust(U, Ub);

            % Rotate and rescale Ub and Vb
            Vb = Vb * Sb * rotatemat_full;
            Ub = Ub * Sb * rotatemat_full;
            Vb = Vb./repmat(diag(S)', size(Vb, 1), 1);
            Ub = Ub./repmat(diag(S)', size(Ub, 1), 1);
        case 2
            % Rotate U and V separately, then average
            rotatemat1 = rri_bootprocrust(U, Ub);
            rotatemat2 = rri_bootprocrust(V, Vb);
            rotatemat_full = (rotatemat1 + rotatemat2)/2;
            Vb = Vb * rotatemat_full;
            Ub = Ub * rotatemat_full;
        otherwise
            error('Invalid value in pls_opts.boot_procrustes_mod!');
    end

    %% Store vectors for all bootstrap samples
    Ub_vect(:, :, iB) = Ub;
    Vb_vect(:, :, iB) = Vb;

    %% CPM process only uses V_stad; skip bootstrap scores, loadings, and correlation matrices
    % [Lxb, Lyb, corr_Lxb_Xb, corr_Lyb_Yb, corr_Lxb_Yb, corr_Lyb_Xb] = myPLS_get_PLS_scores_loadings(Xb, Yb, Vb, Ub, grouping, pls_opts);
    % boot_results.Lxb(:, :, iB) = Lxb;
    % boot_results.Lyb(:, :, iB) = Lyb;
    % boot_results.LC_img_loadings_boot(:, :, iB) = corr_Lxb_Xb;
    % boot_results.LC_behav_loadings_boot(:, :, iB) = corr_Lyb_Yb;
    % boot_results.corr_Lxb_Yb(:, :, iB) = corr_Lxb_Yb;
    % boot_results.corr_Lyb_Xb(:, :, iB) = corr_Lyb_Xb;
    % boot_results.corr_Lxb_Lyb(:, :, iB) = corr(Lxb, Lyb)';
end

%% ========== Compute Bootstrap Statistics ========== %%
boot_stats = myPLS_bootstrap_stats(Ub_vect, Vb_vect);

% Save all statistical fields to boot_results for easy extension
fN = fieldnames(boot_stats);
for iF = 1:length(fN)
    boot_results.(fN{iF}) = boot_stats.(fN{iF});
end

% Save bootstrap resampling data if requested
if pls_opts.save_boot_resampling
    boot_results.Ub_vect = Ub_vect;
    boot_results.Vb_vect = Vb_vect;
end


%% ========== U and V Stability ========== %%
boot_results.U_stad = U ./ sqrt(sum((Ub_vect - boot_stats.Ub_mean).^2, 3) ./ (pls_opts.nBootstraps-1));
boot_results.V_stad = V ./ sqrt(sum((Vb_vect - boot_stats.Vb_mean).^2, 3) ./ (pls_opts.nBootstraps-1));

disp(' ')

end
