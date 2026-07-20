function [feature_mask, selected_features, selector_info] = fit_plsc_second_stage_selector(feature_mask, x_train_norm, y_selector_target, selector_options)
% Purpose:
%   Perform second-stage sparse feature selection using LASSO or ElasticNet on PLSC candidate features
% Inputs:
%   feature_mask - Initial feature mask obtained from PLSC
%   x_train_norm - Normalized training set feature matrix
%   y_selector_target - PLSC-standardized training set behavior selection target
%   selector_options - Second-stage selection parameters
% Outputs:
%   feature_mask - Final feature mask after second-stage selection
%   selected_features - Training set feature matrix after second-stage selection
%   selector_info - Information about the second-stage selection process
% Example:
%   [mask, selected_features, info] = fit_plsc_second_stage_selector(mask, x_train_norm, y_selector_target, selector_options)

%% Normalize second-stage selection parameters
selector_options = normalize_selector_options(selector_options, size(x_train_norm, 1));


%% Initialize second-stage selection information
selector_info = struct();
selector_info.method = selector_options.method;
selector_info.use_cv = selector_options.use_cv;
selector_info.used = false;
selector_info.n_plsc_features = nnz(feature_mask ~= 0);
selector_info.n_selected_features = selector_info.n_plsc_features;
selector_info.alpha = NaN;
selector_info.lambda = NaN;
selector_info.alpha_index = NaN;
selector_info.lambda_index = NaN;
selector_info.mse = NaN;
selector_info.fallback = '';
selector_info.plsc_feature_indices = find(feature_mask ~= 0);
selector_info.selected_feature_indices = selector_info.plsc_feature_indices;
selector_info.coefficients = [];


%% Extract PLSC candidate features
plsc_selected_idx = (feature_mask ~= 0);
selected_features = x_train_norm(:, plsc_selected_idx);

if strcmpi(selector_options.method, 'none')
    return
end


%% Check PLSC candidate features
if isempty(selected_features)
    error('fit_plsc_second_stage_selector:NoPLSCFeatures', ...
        'PLSC/BSR did not select any features, cannot perform %s second-stage selection.', selector_options.method);
end


%% Display second-stage selection process
fprintf('Executing %s second-stage selection...\n', upper(selector_options.method));


%% Initialize best model record
best_mse = inf;
best_alpha = NaN;
best_lambda = NaN;
best_coef = [];
best_keep = [];


%% Execute LASSO/ElasticNet selection based on CV usage
if selector_options.use_cv

    %% Iterate over Alpha candidates and find optimal models
    for alpha_idx = 1:numel(selector_options.alpha_grid)
        % Current Alpha controls the regularization ratio between L1 and L2 penalties
        current_alpha = selector_options.alpha_grid(alpha_idx);

        [coef_mat, fit_info] = lasso(selected_features, y_selector_target, ...
            'Alpha', current_alpha, ...
            'CV', selector_options.cv_folds, ...
            'Standardize', false, ...
            'Lambda', selector_options.lambda_values);

        % Default to use the Lambda with the minimum cross-validation MSE
        lambda_idx = fit_info.IndexMinMSE;
        if selector_options.use_1se && isfield(fit_info, 'Index1SE')
            lambda_idx = fit_info.Index1SE;
        end

        % Extract non-zero coefficient features
        coef_vec = coef_mat(:, lambda_idx);
        keep_idx = abs(coef_vec) > selector_options.coef_tol;

        % If current Lambda has no non-zero coefficients, fallback to the Lambda with the minimum MSE that has non-zero coefficients
        if ~any(keep_idx)
            nonzero_counts = sum(abs(coef_mat) > selector_options.coef_tol, 1);
            candidate_idx = find(nonzero_counts > 0);
            if ~isempty(candidate_idx)
                [~, best_candidate_pos] = min(fit_info.MSE(candidate_idx));
                lambda_idx = candidate_idx(best_candidate_pos);
                coef_vec = coef_mat(:, lambda_idx);
                keep_idx = abs(coef_vec) > selector_options.coef_tol;
            end
        end

        % Record the model with the minimum cross-validation MSE for each Alpha candidate
        current_mse = fit_info.MSE(lambda_idx);
        if current_mse < best_mse
            best_mse = current_mse;
            best_alpha = current_alpha;
            best_lambda = fit_info.Lambda(lambda_idx);
            best_coef = coef_vec;
            best_keep = keep_idx;
            best_alpha_index = alpha_idx;
            best_lambda_index = lambda_idx;
        end
    end

else

    %% Execute LASSO/ElasticNet selection with fixed Alpha and Lambda indices
    best_alpha_index = selector_options.alpha_index;
    best_alpha = selector_options.alpha_grid(best_alpha_index);

    [coef_mat, fit_info] = lasso(selected_features, y_selector_target, ...
        'Alpha', best_alpha, ...
        'Standardize', false, ...
        'Lambda', selector_options.lambda_values);

    best_lambda_index = min(max(1, selector_options.lambda_index), size(coef_mat, 2));
    best_lambda = fit_info.Lambda(best_lambda_index);
    best_coef = coef_mat(:, best_lambda_index);
    best_keep = abs(best_coef) > selector_options.coef_tol;
    best_mse = NaN;
end


%% Save optimal second-stage selection information
selector_info.used = true;
selector_info.alpha = best_alpha;
selector_info.lambda = best_lambda;
selector_info.alpha_index = best_alpha_index;
selector_info.lambda_index = best_lambda_index;
selector_info.mse = best_mse;
selector_info.coefficients = best_coef;


%% Handle no non-zero coefficients case
if isempty(best_keep) || ~any(best_keep)
    selector_info.fallback = 'no_nonzero_coefficients_keep_plsc_features';
    warning('fit_plsc_second_stage_selector:NoNonzeroCoefficients', ...
        '%s did not select any non-zero coefficient features, keeping original PLSC features.', upper(selector_options.method));
    return
end


%% Update final feature mask based on non-zero coefficients
plsc_feature_indices = find(plsc_selected_idx);
feature_mask(plsc_feature_indices(~best_keep)) = 0;
selected_features = x_train_norm(:, feature_mask ~= 0);


%% Save final second-stage selection results
selector_info.n_selected_features = nnz(feature_mask ~= 0);
selector_info.selected_feature_indices = find(feature_mask ~= 0);


%% Display second-stage selection results
fprintf('%s second-stage selection: PLSC features %d -> Non-zero coefficient features %d (Alpha=%.2f, Lambda=%.4g)\n', ...
    upper(selector_options.method), selector_info.n_plsc_features, ...
    selector_info.n_selected_features, selector_info.alpha, selector_info.lambda);
end


function selector_options = normalize_selector_options(selector_options, n_samples)
% Function:
%   Normalize LASSO/ElasticNet secondary screening parameters
% Inputs:
%   selector_options - Original secondary screening parameters
%   n_samples - Number of training samples
% Outputs:
%   selector_options - Normalized secondary screening parameters
% Example:
%   selector_options = normalize_selector_options(selector_options, size(x_train_norm, 1))

%% Set default parameters
if nargin < 1 || isempty(selector_options)
    selector_options = struct();
end

%% Check screening method type
if ~isfield(selector_options, 'method') || isempty(selector_options.method)
    selector_options.method = 'none';
end

%% Check screening method type
selector_options.method = lower(selector_options.method);
if ~ismember(selector_options.method, {'none', 'lasso', 'elasticnet'})
    error('fit_plsc_second_stage_selector:InvalidMethod', ...
        'second_stage_selector.method must be "none", "lasso", or "elasticnet".');
end

%% Check training sample count
if ~strcmpi(selector_options.method, 'none') && n_samples < 2
    error('fit_plsc_second_stage_selector:NotEnoughSamples', ...
        'LASSO/ElasticNet secondary screening requires at least 2 training samples.');
end

%% Set Alpha search range
if strcmpi(selector_options.method, 'lasso')
    selector_options.alpha_grid = 1;
elseif ~isfield(selector_options, 'alpha_grid') || isempty(selector_options.alpha_grid)
    selector_options.alpha_grid = [0.1 0.3 0.5 0.7 0.9 1.0];
end

%% Check if cross-validation is used
if ~isfield(selector_options, 'use_cv') || isempty(selector_options.use_cv)
    selector_options.use_cv = false;
end

%% Set fixed Alpha index
if ~isfield(selector_options, 'alpha_index') || isempty(selector_options.alpha_index)
    selector_options.alpha_index = ceil(numel(selector_options.alpha_grid) / 2);
end
selector_options.alpha_index = min(max(1, selector_options.alpha_index), numel(selector_options.alpha_grid));

%% Set cross-validation folds count
if ~isfield(selector_options, 'cv_folds') || isempty(selector_options.cv_folds)
    selector_options.cv_folds = min(10, n_samples);
end
selector_options.cv_folds = min(max(2, selector_options.cv_folds), n_samples);

%% Set Lambda candidates count
if ~isfield(selector_options, 'num_lambda') || isempty(selector_options.num_lambda)
    selector_options.num_lambda = 100;
end

if ~isfield(selector_options, 'lambda_values') || isempty(selector_options.lambda_values)
    selector_options.lambda_values = logspace(-4, 0, selector_options.num_lambda);
end

%% Set fixed Lambda index
if ~isfield(selector_options, 'lambda_index') || isempty(selector_options.lambda_index)
    selector_options.lambda_index = ceil(numel(selector_options.lambda_values) / 2);
end
selector_options.lambda_index = min(max(1, selector_options.lambda_index), numel(selector_options.lambda_values));

%% Set non-zero coefficient tolerance threshold
if ~isfield(selector_options, 'coef_tol') || isempty(selector_options.coef_tol)
    selector_options.coef_tol = 1e-8;
end

%% Check if 1-SE rule is used
if ~isfield(selector_options, 'use_1se') || isempty(selector_options.use_1se)
    selector_options.use_1se = false;
end
end
