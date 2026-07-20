%% Test Set Prediction Script
% Function: Load the trained final model and make predictions on test set data
% Input: None
% Output: Save prediction results to mat file
% Example: Run the script directly

clear; close all; clc;

%% Environment setup
scripts_path = fileparts(mfilename('fullpath'));
if isempty(scripts_path)
    scripts_path = pwd;
end
addpath(genpath(scripts_path));

% Test Data Path
test_data_path = 'Path_test';
test_base_file = fullfile(test_data_path, 'Processing_Test');
model_path = fullfile(test_data_path, 'Results_PLSR');
output_file = fullfile(test_data_path, 'Test_Predictions');

if ~exist(output_file, 'dir')
    mkdir(output_file)
end

%% Load test data
fprintf('===== Loading test data =====\n');
test_data = load(fullfile(test_base_file, 'data.mat'));
x_test = test_data.brain_data;
y_test = test_data.beh_data;

% Check test set column count
fprintf('Test set dimensions: %d rows x %d columns\n', size(y_test, 1), size(y_test, 2));

%% Model parameter settings
bsr_k = 1;
bsr_thresh = [1.96 2.58 3.29];
n_comp = 5;

%% Second-stage feature selection method settings
second_stage_selector.method = 'elasticnet';   % Must be consistent with the final model script: 'none', 'lasso', 'elasticnet'

%% Predict on each model
mode_label = 'all';

% Display current test process
fprintf('===== Using all behavioral data for prediction =====\n');

for nbsr = 1:numel(bsr_thresh)
    for nc = 1:n_comp
        % Display current model prediction process
        fprintf('\tUsing model for prediction - %s, BSR Threshold %.2f, PLS Components %d\n', mode_label, bsr_thresh(nbsr), nc);

        %% Load model
        model_name = sprintf('final_model_%s_%s_lv%d_bsr%.2f_ncomp%d.mat', mode_label, second_stage_selector.method, bsr_k, bsr_thresh(nbsr), nc);
        model_file = fullfile(model_path, model_name);

        if ~exist(model_file, 'file')
            fprintf('Warning: Model file not found: %s\n', model_file);
            continue;
        end

        load(model_file, 'model_params');

        %% Predict on test set
        y_pred = predict_with_final_model(x_test, model_params);

        %% Calculate prediction performance
        % Use test set parameters to normalize test set labels
        [y_test_eval, y_test_mu, y_test_sigma] = normalize_behav(y_test);

        % Calculate performance metrics
        [R, P] = corr(y_test_eval, y_pred);
        deno = y_test_eval - 0;
        nume = y_test_eval - y_pred;
        q2 = 1 - mean(nume.^2) / mean(deno.^2);

        % Display performance metrics for current test set model
        fprintf('Test Set Performance: q2 = %.4f, r = %.4f, p = %.4f\n', q2, R, P);

        %% Save prediction results
        prediction_results.y_pred = y_pred;
        prediction_results.y_true = y_test;
        prediction_results.y_true_eval = y_test_eval;
        prediction_results.y_true_eval_mu = y_test_mu;
        prediction_results.y_true_eval_sigma = y_test_sigma;
        prediction_results.R = R;
        prediction_results.P = P;
        prediction_results.q2 = q2;
        prediction_results.model_params = model_params;

        save_name = sprintf('test_prediction_newmethod_muly_%s_%s_lv%d_bsr%.2f_ncomp%d.mat', mode_label, second_stage_selector.method, bsr_k, bsr_thresh(nbsr), nc);
        save(fullfile(output_file, save_name), 'prediction_results');
        fprintf('Prediction results saved to: %s\n', fullfile(output_file, save_name));

        % Clear variables for next model
        clear model_params prediction_results y_pred y_test_eval y_test_mu y_test_sigma R P q2
    end
end

fprintf('===== All models predictions completed =====\n');

rmpath(genpath(scripts_path));
