function y_pred_norm = predict_with_final_model(x_test, model_params)
% Function:
%   Predict test set behavioral scores using the trained final model
% Inputs:
%   x_test - Test set connectivity matrix, dimensions: samples x features
%   model_params - Final model parameter structure
% Outputs:
%   y_pred_norm - Predicted scores for test set
% Example:
%   y_pred = predict_with_final_model(test_data, model_params)

%% Normalize test set data (using training set parameters)
x_test_norm = (x_test - model_params.x_mean) ./ model_params.x_std;

%% Predict using selected features
selected_test_features = x_test_norm(:, model_params.selected_feature_indices);

% Calculate predictions using PLSR coefficients
y_pred_norm = [ones(size(selected_test_features, 1), 1), selected_test_features] * model_params.reg_beta;

end
