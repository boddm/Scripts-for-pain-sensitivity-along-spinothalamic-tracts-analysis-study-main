function all_boot_orders = myPLS_get_boot_orders(nBootstraps, grouping, grouped_boot)
% =============================
% Function:
% Generate resampling order matrix for PLS bootstrap. Supports within-group/between-group sampling.
%
% =============================
% Inputs:
% nBootstraps : Integer, number of bootstrap samples
% grouping    : N×1 vector, subject grouping information
%               e.g., [1,1,2] means first two subjects are group 1, third is group 2
% grouped_boot: Logical/0 or 1, whether bootstrap sampling is within-group
%               0 = ignore grouping, sample from entire dataset
%               1 = independent sampling within each group
%
% =============================
% Output:
% all_boot_orders : N×nBootstraps matrix, each column is subject indices for one bootstrap
%                   N: number of subjects
%
% =============================
% Example:
% all_boot_orders = myPLS_get_boot_orders(nBootstraps, grouping, grouped_boot);
%
% =============================

%% ========== Check if within-group sampling ========== %%
if ~grouped_boot
    grouping = ones(size(grouping)); % If not within-group sampling, treat all subjects as one group
end

%% ========== Get grouping information ========== %%
groupIDs = unique(grouping); % All group IDs
nGroups = length(groupIDs);  % Number of groups
nSubj = length(grouping);    % Total number of subjects

all_boot_orders = nan(nSubj, nBootstraps); % Initialize

%% ========== Within-group/between-group resampling main loop ========== %%
for iG = 1:nGroups
    groupID = find(grouping == groupIDs(iG)); % Subject indices in this group
    num_subj_group = length(groupID);         % Number of subjects in this group
    [boot_order_tmp, ~] = rri_boot_order(num_subj_group, 1, nBootstraps); % Resampling order for this group
    all_boot_orders(groupID, :) = groupID(boot_order_tmp); % Fill into the total matrix
end

end