function [y_score, y_mu, y_sigma] = normalize_behav(y, y_mu, y_sigma)
% Purpose:
%   Z-score normalize behavioral indicators column-wise, then take mean
% Inputs:
%   y - behavioral data matrix
%   y_mu - behavioral mean; fitted from current data if empty
%   y_sigma - behavioral standard deviation; fitted from current data if empty
% Outputs:
%   y_score - behavioral total score (mean of Z-scores)
%   y_mu - behavioral mean
%   y_sigma - behavioral standard deviation
% Example:
%   [y_score, y_mu, y_sigma] = normalize_behav(y)

%% Check input behavioral data
if isempty(y) || size(y, 2) < 1
    error('normalize_behav:NoBehaviors', 'At least 1 behavioral indicator is required.');
end


%% Fit or reuse behavioral normalization parameters if not provided
if nargin < 2 || isempty(y_mu) || nargin < 3 || isempty(y_sigma)
    [y_z, y_mu, y_sigma] = zscore(y, 0, 1);
else
    y_z = (y - y_mu) ./ y_sigma;
end


%% Handle zero-standard deviation columns in behavioral data matrix
zero_sigma_idx = (y_sigma == 0);
y_sigma(zero_sigma_idx) = 1;
y_z(:, zero_sigma_idx) = 0;

y_score = mean(y_z, 2);

end
