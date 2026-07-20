function [feature_mask, y_test_pred, y_mu, y_sigma, reg_beta, selector_info] = ...
    cpm_lr_D1D2(x_train, x_test, y_train, ~, bsr_k, n_comp, bsr_thresh, selector_options)
% Function:
%   Predict test set behavioral scores using PLSC, secondary feature selection, and PLSR model
% Inputs:
%   x_train - Training set connectivity matrix, dimensions: features x samples
%   x_test - Test set connectivity matrix, dimensions: features x samples
%   y_train - Training set behavioral data
%   bsr_k - Latent variable number used in PLSC
%   n_comp - Number of PLSR components
%   bsr_thresh - Bootstrap Ratio threshold
%   selector_options - LASSO/ElasticNet secondary selection parameters
% Outputs:
%   feature_mask - Final feature selection mask
%   y_test_pred - Test set predicted scores
%   y_mu - Training set behavioral mean
%   y_sigma - Training set behavioral standard deviation
%   reg_beta - PLSR regression coefficients
%   selector_info - Secondary selection information
% Example:
%   [mask, pred] = cpm_lr_D1D2(trainMat, testMat, trainBehav, [], 1, 4, 2, selector_options)

%% Transpose connectivity matrix (features x samples -> samples x features)
x_train = x_train';
x_test = x_test';

%% Set default secondary selection parameters
if nargin < 8 || isempty(selector_options)
    % Keep original PLSC+PLSR process when parameters are not provided
    selector_options = struct('method', 'none');
end

%% Get feature mask, behavior normalization parameters, and regression coefficients
[feature_mask, y_mu, y_sigma, reg_beta, selector_info] = ...
    PLSC_CPMmask(x_train, y_train, bsr_k, n_comp, bsr_thresh, selector_options);

%% Get non-zero feature indices
selected_feature_indices = feature_mask ~= 0;

%% Normalize test set data
x_test_norm = (x_test - mean(x_train)) ./ (std(x_train) + 1e-8);

%% Predict using selected features
selected_test_features = x_test_norm(:, selected_feature_indices);
y_test_pred = [ones(size(selected_test_features, 1), 1), selected_test_features] * reg_beta;
end
