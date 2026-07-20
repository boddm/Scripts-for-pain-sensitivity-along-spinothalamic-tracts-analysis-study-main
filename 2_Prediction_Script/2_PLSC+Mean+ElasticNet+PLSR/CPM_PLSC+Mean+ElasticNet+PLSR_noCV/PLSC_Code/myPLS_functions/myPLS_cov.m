function R = myPLS_cov(X, Y, grouping, grouped_PLS)
% =============================
% Function:
% Generate cross-covariance matrix R required for PLS analysis. Supports both
% grouped and non-grouped modes; in grouped mode, covariance matrices from
% each group are concatenated by rows.
%
% =============================
% Inputs:
% X           : N×M matrix, imaging data
%               N: number of subjects
%               M: number of imaging variables
% Y           : N×B matrix, behavioral/design data
%               B: number of behavioral/design variables
% grouping    : N×1 vector, subject grouping information
%               e.g., [1,1,2] means first two are group 1, third is group 2
% grouped_PLS : logical/0 or 1, whether to use grouped PLS
%               0 = ignore grouping, compute covariance for all data
%               1 = compute covariance separately for each group and concatenate
%
% =============================
% Output:
% R           : cross-covariance matrix
%               B×M matrix for non-grouped, (G×B)×M matrix for grouped (G = number of groups)
%
% =============================
% Example:
% R = myPLS_cov(X, Y, grouping, grouped_PLS);
%
% =============================

%% ========== Check if grouped PLS ========== %%
if ~grouped_PLS
    grouping = ones(size(grouping)); % If not grouped, treat all subjects as one group
end

%% ========== Get grouping information ========== %%
groupIDs = unique(grouping); % All group IDs
groups = length(groupIDs);   % Number of groups

%% ========== Compute cross-covariance for each group and concatenate ========== %%
for iG = 1:groups
    Ysel = Y(grouping == groupIDs(iG), :); % Behavioral data for this group
    Xsel = X(grouping == groupIDs(iG), :); % Imaging data for this group
    R0 = Ysel.'*Xsel; % Cross-covariance for this group
    if ~exist('R', 'var')
        R = R0; % First group: direct assignment
    else
        R = [R; R0]; % Other groups: concatenate by rows
    end
end

end