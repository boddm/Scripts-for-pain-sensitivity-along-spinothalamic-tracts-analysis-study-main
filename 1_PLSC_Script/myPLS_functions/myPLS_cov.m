function R = myPLS_cov(X, Y, grouping, grouped_PLS)
% This function generates a cross-covariance matrix (stacked) for the PLS
% analysis
% 
% Inputs:
% - X            : N x M matrix, N is #subjects, M is #imaging variables
% - Y            : N x B matrix, B is #behaviors
% - grouping     : N x 1 vector, subject grouping (e.g. diagnosis)
%                  e.g. [1,1,2] = subjects 1 and 2 belong to group 1,
%                  subject 3 belongs to group 2
% - grouped_PLS  : binary variable indicating if groups should be 
%                  considered when computing R
%              0 = PLS will computed over all subjects (ignoring grouping)
%              1 = R will be constructed by concatenating group-wise
%                  covariance matrices (as in conventional behavior PLS)
%
% Outputs:
% - R            : cross-covariance matrix
%
%
% Function based on myPLS_cov by Dimitri Van De Ville
% Modifications by D. Zoeller (Aug 2019): 
%   - adapted to work also for different group labels than 1,2,3,... (e.g. 2,3,4)
%   - removed group number input (can be derived from grouping vector)

%% ========== Check whether grouped PLS ========== %%
if ~grouped_PLS
    grouping = ones(size(grouping)); % If not grouped, treat all subjects as one group
end

%% ========== Number and IDs of groups ========== %%
groupIDs = unique(grouping); % IDs of all groups
groups = length(groupIDs);   % Number of groups

%% ========== Compute and concatenate group-wise cross-covariance ========== %%
for iG = 1:groups
    Ysel = Y(grouping == groupIDs(iG), :); % Behavioral data for this group
    Xsel = X(grouping == groupIDs(iG), :); % Imaging data for this group
    R0 = Ysel.'*Xsel; % Cross-covariance for this group
    if ~exist('R', 'var')
        R = R0; % Assign directly for the first group
    else
        R = [R; R0]; % Concatenate row-wise for remaining groups
    end
end

end