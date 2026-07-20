function [feature_mask, y_mu, y_sigma, plsr_beta, selector_info] = PLSC_CPMmask(x_train, y_train, bsr_k, n_comp, bsr_thresh, selector_options)
% Function:
%   Initial feature screening using PLSC, combined with LASSO/ElasticNet and PLSR for regression modeling
% Inputs:
%   x_train - Training set connectivity matrix, dimensions: samples x features
%   y_train - Training set behavioral data
%   bsr_k - Latent variable index used in PLSC
%   n_comp - Number of PLSR components
%   bsr_thresh - Bootstrap Ratio threshold
%   selector_options - LASSO/ElasticNet secondary screening parameters
% Outputs:
%   feature_mask - Final feature selection mask
%   y_mu - Training set behavioral mean
%   y_sigma - Training set behavioral standard deviation
%   plsr_beta - PLSR regression coefficients
%   selector_info - Secondary screening information
% Example:
%   [mask, y_mu, y_sigma, beta, info] = PLSC_CPMmask(connMat, behavData, 1, 4, 2, selector_options)

%% Parameter settings
n_bootstraps = 5000;
norm_type = 1;         % Normalization type
subj_groups = ones(size(x_train, 1), 1);  % Subject grouping


%% Default secondary screening parameters
if nargin < 6 || isempty(selector_options)
    % Default to original PLSC process
    selector_options = struct('method', 'none');
end

%% 1. Set PLSC toolbox path
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'PLSC_Code')));


%% 2. Define all input parameters
fprintf('Define PLSC input parameters...\n');
input = struct();
input.brain_data = x_train;
input.behav_data = y_train;
input.grouping = subj_groups;

% --- PLS options ---
pls_opts = struct();
pls_opts.nBootstraps = n_bootstraps;
pls_opts.normalization_img = norm_type;
pls_opts.normalization_behav = norm_type;
pls_opts.grouped_PLS = 0;
pls_opts.grouped_perm = 0;
pls_opts.grouped_boot = 0;
pls_opts.boot_procrustes_mod = 2;
pls_opts.save_boot_resampling = 0;
pls_opts.behav_type = 'behavior';

%% 3. Check all input parameters
% !!! Check input parameters before running PLS analysis to ensure validity !!!
[input, pls_opts] = myPLS_initialize(input, pls_opts);


%% 4. Run PLSC analysis (including updated Bootstrap process)
fprintf('Run PLSC analysis...\n');
res = myPLS_analysis(input, pls_opts);


%% 5. Extract normalized data and Bootstrap ratio
fprintf('Extract PLSC results...\n');
x_train_norm = res.X;
y_train_norm = res.Y;

right_bsr = res.boot_results.V_stad;

% Bootstrap ratio may be Inf/NaN when standard error is close to 0
% Use old program's handling method
inf_right_idx = find(~isfinite(right_bsr));
for iter_inf = 1:numel(inf_right_idx)
    right_bsr(inf_right_idx(iter_inf)) = res.V(inf_right_idx(iter_inf));
end
clear inf_right_idx iter_inf

if bsr_k > size(right_bsr, 2)
    error('PLSC_CPMmask:InvalidBSRK', 'bsr_k=%d超过当前PLSC成分数%d。', bsr_k, size(right_bsr, 2));
end

%% 6. Feature selection
fprintf('Feature selection based on Bootstrap Ratio...\n');
feature_mask = zeros(size(right_bsr, 1), 1);
pos_feat_idx = right_bsr(:, bsr_k) > bsr_thresh;
neg_feat_idx = right_bsr(:, bsr_k) < -bsr_thresh;
feature_mask(pos_feat_idx) = 1;
feature_mask(neg_feat_idx) = -1;

selected_feat_idx = (feature_mask ~= 0);
if ~any(selected_feat_idx)
    error('PLSC_CPMmask:NoPLSCFeatures', '当前BSR阈值下PLSC没有筛出任何特征。');
end


%% 7. Generate target for secondary screening
y_selector_target = mean(y_train_norm, 2);


%% 8. Fit secondary screening model
[feature_mask, selected_features, selector_info] = fit_plsc_second_stage_selector( ...
    feature_mask, x_train_norm, y_selector_target, selector_options);


%% 9. Normalize behavioral data
[y_train_mean, y_mu, y_sigma] = normalize_behav(y_train);


%% 10. PLS regression prediction
fprintf('PLS regression...\n');

effective_n_comp = min([n_comp, size(selected_features, 2), size(selected_features, 1) - 1]);
if effective_n_comp < 1
    error('PLSC_CPMmask:NotEnoughFeatures', 'Insufficient number of features available for PLSR after secondary screening.');
end

selector_info.requested_n_comp = n_comp;
selector_info.effective_n_comp = effective_n_comp;

if effective_n_comp < n_comp
    % Display PLSR component number auto-adjustment process
    fprintf('PLS components adjusted from %d to %d to match the number of features after screening.\n', n_comp, effective_n_comp);
end

[~, ~, ~, ~, plsr_beta, ~, ~, ~] = plsregress(selected_features, y_train_mean, effective_n_comp); % PLSR
end
