function R = myPLS_cov(X, Y, grouping, grouped_PLS)
% =============================
% Function:
% Generate cross-covariance matrix R for PLS analysis. Supports grouped and
% non-grouped modes. In grouped mode, concatenates covariance matrices by rows.
%
% =============================
% Inputs:
% X           : N×M matrix, imaging data
%               N: number of subjects
%               M: number of imaging variables
% Y           : N×B matrix, behavioral/design data
%               B: number of behavioral/design variables
% grouping    : N×1 vector, subject group information
%               e.g., [1,1,2] means first two are group 1, third is group 2
% grouped_PLS : logical/0 or 1, whether to use grouped PLS
%               0 = ignore grouping, compute covariance on whole data
%               1 = compute covariance for each group and concatenate
%
% =============================
% Output:
% R           : cross-covariance matrix
%               B×M for non-grouped, (G×B)×M for grouped (G = number of groups)
%
% =============================
% Example:
% R = myPLS_cov(X, Y, grouping, grouped_PLS);
%
% =============================

%% ========== if grouped PLS ========== %%
if ~grouped_PLS
    grouping = ones(size(grouping)); % If not grouped, treat all subjects as one group
end

%% ========== get group information ========== %%
groupIDs = unique(grouping); % Group IDs
groups = length(groupIDs);   % Number of groups

%% ========== calculate cross-covariance matrices for each group and concatenate ========== %%
for iG = 1:groups
    Ysel = Y(grouping == groupIDs(iG), :); % Behavioral data for this group
    Xsel = X(grouping == groupIDs(iG), :); % Imaging data for this group
    R0 = Ysel.'*Xsel; % Cross-covariance for this group
    if ~exist('R', 'var')
        R = R0; % Assign directly for first group
    else
        R = [R; R0]; % Concatenate by rows for remaining groups
    end
end

end