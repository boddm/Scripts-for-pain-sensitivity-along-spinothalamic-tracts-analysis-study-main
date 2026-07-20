function [X, meanX, stdX] = myPLS_norm(X, grouping, mode)
% This function applies normalization on data before running PLS
%
% Inputs:
% - X            : N x V matrix, N is #subjects, V is #variables
% - grouping     : N x 1 vector, subject grouping (e.g. diagnosis)
%                  e.g. [1,1,2] = subjects 1 and 2 belong to group 1,
%                  subject 3 belongs to group 2
% - mode         : normalization option 
%                  0 = no normalization
%                  1 = zscore across all subjects
%                  2 = zscore within groups (default)
%                  3 = std normalization (no centering) across subjects 
%                  4 = std normalization (no centering) within groups 
%
% Outputs:
% - X            : normalized X matrix
%
%
% Function based on myPLS_norm by D. Van De Ville
% Modifications by D. Zoeller (Aug 2019): 
%   - adapted to work also for different group labels than 1,2,3,... (e.g. 2,3,4)
%   - removed group number input (can be derived from grouping vector)

%% ========== Number and IDs of groups ========== %%
groupIDs = unique(grouping); % all group labels
nGroups = length(groupIDs);  % number of groups

%% ========== main normalization pipeline ========== %%
switch mode
    case 1
        %% z-score across all subjects
        meanX = mean(X); % mean
        stdX = std(X);   % standard deviation
        X = zscore(X);   % normalization
    case 2
        %% z-score within groups
        for iG = 1:nGroups
            idx = find(grouping == groupIDs(iG)); % indices for this group
            meanX(iG, :) = mean(X(idx, :));       % group mean
            stdX(iG, :) = std(X(idx, :));         % group std
            X(idx, :) = zscore(X(idx, :));        % normalize within group
        end
    case 3
        %% standard normalization across all subjects (no centering)
        meanX = mean(X); % mean
        stdX = sqrt(mean(X.^2, 1)); % std (no centering)
        X2 = stdX;
        X = X ./ repmat(X2, [size(X, 1) 1]); % normalize
    case 4
        %% standard normalization within groups (no centering)
        for iG = 1:nGroups
            idx = find(grouping == groupIDs(iG)); % indices for this group
            meanX(iG, :) = mean(X(idx, :));       % group mean
            stdX(iG, :) = sqrt(mean(X(idx,:).^2, 1)); % group std (no centering)
            X2 = stdX(iG, :);
            X(idx, :) = X(idx, :) ./ repmat(X2, [size(X(idx, :), 1) 1]); % normalize within group
        end
end
end