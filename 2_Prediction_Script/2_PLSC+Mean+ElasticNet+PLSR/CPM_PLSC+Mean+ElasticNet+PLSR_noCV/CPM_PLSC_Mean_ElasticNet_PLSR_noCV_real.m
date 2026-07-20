%% Main Script - Final Model Version
% Purpose: Perform CPM-PLSR analysis, create final prediction model (no cross-validation)
% Input: None
% Output: Save final model to mat file
% Example: Run the script directly
clear; close all; clc;

%% Set random seed (ensure reproducibility)
rng(42);  % Fixed random seed, ensure reproducibility of results across runs

%% Environment setup and path configuration
scripts_path = fileparts(mfilename('fullpath'));
if isempty(scripts_path)
    scripts_path = pwd;
end
addpath(genpath(scripts_path));

path = 'Path';
base_file = fullfile(path, 'Processing');
output_file = fullfile(path, 'Results_PLSR');

if ~exist(output_file, 'dir')
    mkdir(output_file)
end

%% Load data
Data = load(fullfile(base_file, 'data.mat'));

x0 = Data.brain_data;
y0 = Data.beh_data;

%% Final Model Parameters Setting
fprintf('===== Create Final Prediction Model =====\n');

%% Parameter Setting
bsr_k = 1;
bsr_thresh = [1.96 2.58 3.29];
n_comp = 5;

%% Secondary Feature Selection Parameters Setting
second_stage_selector.method = 'elasticnet';   % Optional: 'none', 'lasso', 'elasticnet'
second_stage_selector.alpha_grid = [0.01 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 ...
    0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 0.99]; % ElasticNet will search for optimal alpha in these values; LASSO is fixed to 1
second_stage_selector.use_cv = true;
second_stage_selector.alpha_index = 11;
second_stage_selector.cv_folds = 5;
second_stage_selector.num_lambda = 100;
second_stage_selector.lambda_values = logspace(-2, 2, 100);
second_stage_selector.lambda_index = 50;
second_stage_selector.coef_tol = 1e-4;
second_stage_selector.use_1se = false;

%% Create Final Prediction Model with Different Parameter Combinations
behav_cols = 1:size(y0, 2);
mode_label = 'all';
y_mode = y0;
n_selected_behav = size(y_mode, 2);

fprintf('===== Use All Behavioral Data, Behavioral Columns: %s =====\n', mat2str(behav_cols));

for nbsr = 1:numel(bsr_thresh)
    for nc = 1:n_comp
        fprintf('\tCreate Final Prediction Model - %s, BSR Threshold %.2f, PLS Components %d\n', mode_label, bsr_thresh(nbsr), nc);
        clear model_params

        %% Train Model with All Data
        x_train = x0;
        y_train = y_mode;

        %% Calculate Image Normalization Parameters
        x_mean = mean(x_train);
        x_std = std(x_train) + 1e-8;

        %% Normalize Image Data
        x_train_norm = (x_train - x_mean) ./ x_std;

        %% CPM Model Training with All Data
        [feature_mask, y_mu, y_sigma, reg_beta, selector_info] = ...
            PLSC_CPMmask(x_train, y_train, bsr_k, nc, bsr_thresh(nbsr), second_stage_selector);

        %% Get selected feature indices
        selected_feature_indices = feature_mask ~= 0;

        %% Save model parameters
        model_params.feature_mask = feature_mask;
        model_params.selected_feature_indices = selected_feature_indices;
        model_params.reg_beta = reg_beta;
        model_params.bsr_k = bsr_k;
        model_params.bsr_thresh = bsr_thresh(nbsr);
        model_params.n_comp = nc;
        model_params.n_selected_behav = n_selected_behav;
        model_params.behav_mode = mode_label;
        model_params.behav_cols = behav_cols;
        model_params.second_stage_selector = second_stage_selector;
        model_params.selector_info = selector_info;

        %% Save normalization parameters
        model_params.x_mean = x_mean;
        model_params.x_std = x_std;
        model_params.y_mu = y_mu;
        model_params.y_sigma = y_sigma;

        %% Calculate model performance on training set
        selected_train_features = x_train_norm(:, selected_feature_indices);
        y_train_pred = [ones(size(selected_train_features, 1), 1), selected_train_features] * reg_beta;

        [y_train_mean, ~, ~] = normalize_behav(y_train, y_mu, y_sigma);

        % Calculate performance metrics
        [model_params.R, model_params.P] = corr(y_train_mean, y_train_pred);
        deno = y_train_mean - 0;
        nume = y_train_mean - y_train_pred;
        model_params.q2 = 1 - mean(nume.^2) / mean(deno.^2);

        model_params.y_train_pred = y_train_pred;
        model_params.y_train_mean = y_train_mean;
        model_params.y_train_score = y_train_mean;

        % Display current final model performance
        fprintf('Model Performance: q2 = %.4f, r = %.4f, p = %.4f\n', model_params.q2, model_params.R, model_params.P);
        %% Save Final Model
        fprintf('===== Save Final Model =====\n');

        save_name = sprintf('final_model_%s_%s_lv%d_bsr%.2f_ncomp%d.mat', mode_label, second_stage_selector.method, bsr_k, bsr_thresh(nbsr), nc);

        save(fullfile(output_file, save_name), 'model_params');
        fprintf('Final model saved to: %s\n', fullfile(output_file, save_name));
    end
end

rmpath(genpath(scripts_path));