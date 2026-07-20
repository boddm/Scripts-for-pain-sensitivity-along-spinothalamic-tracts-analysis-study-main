function [y_score, y_mu, y_sigma] = normalize_behav(y, y_mu, y_sigma)
% Function:
%   Normalize behavior scores to Z-score and take the mean
% Input:
%   y - Behavior data matrix
%   y_mu - Behavior mean; fitted using current data if empty
%   y_sigma - Behavior standard deviation; fitted using current data if empty
% Output:
%   y_score - Total behavior score after Z-score and mean
%   y_mu - Behavior mean
%   y_sigma - Behavior standard deviation
% Example:
%   [y_score, y_mu, y_sigma] = normalize_behav(y)

%% Check input behavior data
if isempty(y) || size(y, 2) < 1
    error('normalize_behav:NoBehaviors', 'At least 1 behavior behavior is required.');
end

%% Fit or load behavior mean
if nargin < 2 || isempty(y_mu)
    y_mu = mean(y, 1);
end

%% Fit or load behavior standard deviation
if nargin < 3 || isempty(y_sigma)
    y_sigma = std(y, 0, 1);
end

%% Handle zero standard deviation columns
zero_sigma_idx = (y_sigma == 0);
y_sigma(zero_sigma_idx) = 1;

%% Calculate Z-score behavior scores
y_z = (y - y_mu) ./ y_sigma;
y_z(:, zero_sigma_idx) = 0;

y_score = mean(y_z, 2);

end
