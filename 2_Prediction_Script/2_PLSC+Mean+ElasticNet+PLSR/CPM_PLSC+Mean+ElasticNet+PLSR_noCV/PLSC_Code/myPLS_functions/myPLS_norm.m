function [X, meanX, stdX] = myPLS_norm(X, grouping, mode)
% =============================
% Function:
% Performs normalization on input data X, supporting various methods including
% global z-score, within-group z-score, and standard deviation normalization.
% Used for data standardization before PLS analysis to ensure comparability
% across different variables/groups.
%
% =============================
% Inputs:
% X        : N×V matrix, raw data
%            N: number of subjects, V: number of variables
% grouping : N×1 vector, subject group information
% mode     : normalization method
%            0 = no normalization
%            1 = z-score across all subjects (mean 0, variance 1)
%            2 = within-group z-score (z-score per group, default)
%            3 = standard normalization across all subjects (no centering, variance 1)
%            4 = within-group standard normalization (no centering, variance 1)
%
% =============================
% Outputs:
% X        : normalized data matrix
% meanX    : mean (per group or global)
% stdX     : standard deviation (per group or global)
%
% =============================
% Example:
% [X, meanX, stdX] = myPLS_norm(X, grouping, mode);
%
% =============================

%% ========== Get group information ========== %%
groupIDs = unique(grouping); % All group IDs
nGroups = length(groupIDs);  % Number of groups

%% ========== Normalization main process ========== %%
switch mode
    case 1
        %% Z-score across all subjects
        meanX = mean(X); % Mean
        stdX = std(X);   % Standard deviation
        X = zscore(X);   % Normalization
    case 2
        %% Within-group z-score
        for iG = 1:nGroups
            idx = find(grouping == groupIDs(iG)); % Index of current group
            meanX(iG, :) = mean(X(idx, :));       % Group mean
            stdX(iG, :) = std(X(idx, :));         % Group standard deviation
            X(idx, :) = zscore(X(idx, :));        % Group normalization
        end
    case 3
        %% Standard normalization across all subjects (no centering)
        meanX = mean(X); % Mean
        stdX = sqrt(mean(X.^2, 1)); % Standard deviation (no centering)
        X2 = stdX;
        X = X ./ repmat(X2, [size(X, 1) 1]); % Normalization
    case 4
        %% Within-group standard normalization (no centering)
        for iG = 1:nGroups
            idx = find(grouping == groupIDs(iG)); % Index of current group
            meanX(iG, :) = mean(X(idx, :));       % Group mean
            stdX(iG, :) = sqrt(mean(X(idx,:).^2, 1)); % Group standard deviation (no centering)
            X2 = stdX(iG, :);
            X(idx, :) = X(idx, :) ./ repmat(X2, [size(X(idx, :), 1) 1]); % Group normalization
        end
end
end