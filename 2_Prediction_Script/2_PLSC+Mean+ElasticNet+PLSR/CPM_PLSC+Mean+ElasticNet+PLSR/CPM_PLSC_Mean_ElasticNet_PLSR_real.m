%% Main Script (Optimized Version)
% Function: Perform CPM-PLSR analysis, including cross-validation and performance evaluation
% Input: None
% Output: Save results to mat file
% Example: Run the script directly

clear; close all; clc;

%% Set reproducibility
rng(42);  % Fixed random seed, ensure reproducibility

%% Environment Setup
scripts_path = fileparts(mfilename('fullpath'));
if isempty(scripts_path)
    scripts_path = pwd;
end
addpath(genpath(scripts_path));

path = 'path';
base_file = fullfile(path, 'Processing');
output_file = fullfile(path, 'Results_PLSR');

if ~exist(output_file, 'dir')
    mkdir(output_file)
end


%% Load data
Data = load(fullfile(base_file, 'data.mat'));

x0 = Data.brain_data';
y0 = Data.beh_data;

%% Generate cross-validation parameters
fprintf('===== Generate cross-validation parameters =====\n');

n_subjs = size(x0, 2);
n_permutations = 1000;
n_folds = 51;

%% Generate perm_indices
fprintf('Generating perm_indices...\n');
perm_idxs = zeros(n_permutations, n_subjs);
for np = 1:n_permutations
    perm_idxs(np, :) = randperm(n_subjs);
end

%% Generate fold_sizes and fold_bounds (shared across all permutations)
fprintf('Generating fold_sizes...\n');
fold_sizes = repmat(floor(n_subjs / n_folds), n_folds, 1);
fold_sizes(1:mod(n_subjs, n_folds)) = fold_sizes(1:mod(n_subjs, n_folds)) + 1;

fold_bounds = zeros(n_folds, 2);
start_idx = 1;
for nf = 1:n_folds
    end_idx = start_idx + fold_sizes(nf) - 1;
    fold_bounds(nf, :) = [start_idx, end_idx];
    start_idx = end_idx + 1;
end

%% Real label cross-validation for the first permutation of subjects

% Use real label cross-validation for the first permutation of subjects
real_perm_idx = 1;

%% Parameters
bsr_k = 1;
bsr_thresh = [1.96 2.58 3.29];
n_comp = 5;

%% Generate second-stage feature selection parameters
second_stage_selector.method = 'elasticnet';   % Optional: 'none', 'lasso', 'elasticnet'
second_stage_selector.alpha_grid = [0.01 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 ...
    0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 0.99]; % ElasticNet will search for optimal Alpha; LASSO fixed at 1
second_stage_selector.use_cv = true;
second_stage_selector.alpha_index = 6;
second_stage_selector.cv_folds = 5;
second_stage_selector.num_lambda = 100;
second_stage_selector.lambda_values = logspace(-1, 1, 100);
second_stage_selector.lambda_index = 50;
second_stage_selector.coef_tol = 1e-4;
second_stage_selector.use_1se = false;

%% Initialize storage for results
behav_cols = 1:size(y0, 2);
mode_label = 'all';
y_mode = y0;

% Show current behavior mode progress
fprintf('===== Use all behavior data, behavior columns: %s =====\n', mat2str(behav_cols));

for nbsr = 1:numel(bsr_thresh)
    for nc = 1:5
        % Clear or reset variables before each inner loop
        clear perf_results beta feature_mask
        perf_results.behav_mode = mode_label;
        perf_results.behav_cols = behav_cols;
        perf_results.second_stage_selector = second_stage_selector;
        y_test_norm = zeros(size(y_mode, 1), 1);  % Test set Z-score normalized behavior scores
        y_test_mean = zeros(size(y_mode, 1), 1);  % Test set standardized behavior scores
        y_test_pred = zeros(size(y_mode, 1), 1);  % Test set predicted behavior scores

        for nf = 1:n_folds
            % Show current cross-validation progress
            fprintf('\t%s real label cross-validation - Fold %d/%d\n', mode_label, nf, n_folds);

            %% Split training set / test set
            test_subj_idx = perm_idxs(real_perm_idx, fold_bounds(nf, 1):fold_bounds(nf, 2));
            train_subj_idx = setdiff(1:size(y_mode, 1), test_subj_idx);

            %% Get connection matrix
            x_train = x0(:, train_subj_idx);
            x_test = x0(:, test_subj_idx);

            %% Behavioral data (without Z-score normalization)
            y_train = y_mode(train_subj_idx, :);
            y_test = y_mode(test_subj_idx, :);

            %% CPM prediction
            [feature_mask, y_test_pred(test_subj_idx, :), perf_results.y_mu{nf}, perf_results.y_sigma{nf}, ...
                perf_results.beta{nf}, perf_results.selector_info{nf}] = ...
                cpm_lr_D1D2(x_train, x_test, y_train, [], bsr_k, nc, bsr_thresh(nbsr), second_stage_selector);

            %% Test set Z-score normalized behavior scores (using training set parameters, no rescaling)
            [y_test_norm(test_subj_idx, :), ~, ~] = normalize_behav( ...
                y_test, perf_results.y_mu{nf}, perf_results.y_sigma{nf});

            %% Test set Z-score normalized behavior scores (using training set parameters, no rescaling)
            y_test_mean(test_subj_idx, :) = y_test_norm(test_subj_idx, :);

            %% Store feature mask
            perf_results.mask{nf} = feature_mask;
        end

        %% Calculate real performance
        fprintf('===== Calculate real performance =====\n');

        for t1 = 1:1
            for t2 = 1:n_selected_behav
                [perf_results.R(t1, t2), perf_results.P(t1, t2)] = corr(y_test_mean(:, t1), y_test_pred(:, t2));
                deno = y_test_mean(:, t1) - 0;
                nume = y_test_mean(:, t1) - y_test_pred(:, t2);
                perf_results.q2(t1, t2) = 1 - mean(nume.^2) / mean(deno.^2);
            end
        end

        %% Save results
        fprintf('===== Save results =====\n');

        save_name = sprintf('result_real_%s_%s_lv%d_bsr%.2f_ncomp%d_fold%d.mat', mode_label, second_stage_selector.method, bsr_k, bsr_thresh(nbsr), nc, n_folds);

        save(fullfile(output_file, save_name), 'perf_results', 'y_test_pred', 'y_test_norm', 'y_test_mean');
        fprintf('Results saved to: %s\n', fullfile(output_file, save_name));

    end
end

rmpath(genpath(scripts_path));