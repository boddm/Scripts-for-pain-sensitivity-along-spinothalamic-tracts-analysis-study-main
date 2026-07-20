function [Lx, Ly, corr_Lx_X, corr_Ly_Y, corr_Lx_Y, corr_Ly_X] = myPLS_get_PLS_scores_loadings(X, Y, V, U, grouping, pls_opts)
% This function computes PLS scores & loadings.
% Individual-specific PLS scores are subjects' expression of imaging and 
% behavior/design saliences (covariance patterns).
% Component-specific PLS loadings express the contribution of original 
% imaging or behavior/design variables to the latent components (LCs).
%
% Inputs:
% - X           : N x M, N is #subjects, M is #imaging variables, imaging data
% - Y           : N x B, B is #behavior/design variables, behavior/design data
% - V           : M x L, L is #LCs, imaging saliences
% - U           : B x (L x #groups), behavior/design saliences
% - grouping    : N x 1 vector, subject grouping (e.g. diagnosis)
%                 e.g. [1,1,2] = subjects 1 and 2 belong to group 1,
%                 subject 3 belongs to group 2
% - pls_opts    : options for the PLS analysis
% 
% Outputs:
% - Lx          : N x L, PLS imaging scores
% - Ly          : N x (L x #groups), PLS behavior/design scores
% - corr_Lx_X   : M x L, Pearson's correlation between imaging data and 
%                 PLS imaging scores (i.e., structure coefficients)
% - corr_Ly_Y   : B x (L x #groups), Pearson's correlation between behavior/design data
%                 and PLS behavior/design scores (i.e. structure coefficients)
% - corr_Lx_Y   : M x L, Pearson's correlation between 
%                 behavior/design data and PLS imaging scores
% - corr_Ly_X   : B x (L x #groups), Pearson's correlation between 
%                 behavior/design data and PLS imaging scores

%% ========== Check if grouped PLS is enabled ========== %%
if ~pls_opts.grouped_PLS
    grouping = ones(size(grouping)); % if not grouped PLS, treat all subjects as one group
end

%% ========== Extract group information ========== %%
groupIDs = unique(grouping); % group labels
nGroups = length(groupIDs);  % number of groups
nBehav = size(Y, 2);         % number of behavioral variables
nLC = size(V, 2);            % number of latent components

%% ========== Compute PLS scores ========== %%
% Imaging scores Lx = X * V
Lx = X * V; % N×L

% Compute behavioral/design scores Ly for each group
iter = 1;
for iG = 1:nGroups
    % Select U for this group
    idx = iter:iter + nBehav - 1;
    Usel = U(idx, :); % behavioral/design weights for this group
    
    % Select Y for this group and compute scores
    for iG2 = 1:nGroups
        this_groupID = find(grouping == groupIDs(iG2)); % subject indices for this group
        Ysel = Y(this_groupID, :); % behavioral data for this group
        Lyy(iG, this_groupID, :) = Ysel * Usel; % compute scores for this group
    end        
    iter = iter + nBehav;
end

% Merge behavioral/design scores across groups
for iG = 1:nGroups
    this_groupID = find(grouping == groupIDs(iG));
    if nLC == 1
        Ly(this_groupID, :) = Lyy(iG, this_groupID, :)';
    else
        Ly(this_groupID, :) = Lyy(iG, this_groupID, :);
    end
end

%% ========== Compute PLS loadings ========== %%
% Correlation between imaging scores and imaging data (structure coefficients)
corr_Lx_X = corr(Lx, X)'; % M×L

% Compute correlation between behavioral/design scores and behavioral/design data for each group
iter = 1;
for iG = 1:nGroups
    idx = iter:iter + nBehav - 1;
    this_groupID = find(grouping == groupIDs(iG));
    corr_Ly_Y(idx, :) = corr(Ly(this_groupID, :), Y(this_groupID, :))';
    iter = iter + nBehav;
end

% Compute correlation between imaging scores and behavioral/design data for each group
iter = 1;
for iG = 1:nGroups
    idx = iter:iter + nBehav - 1;
    this_groupID = find(grouping == groupIDs(iG));
    corr_Lx_Y(idx, :) = corr(Lx(this_groupID, :), Y(this_groupID, :))';
    iter = iter + nBehav;
end

% Correlation between behavioral/design scores and imaging data
corr_Ly_X = corr(Ly, X)'; % B×(L×number of groups)

% Progress display: PLS scores and loadings computed
disp('PLS scores and loadings computation completed.');

end