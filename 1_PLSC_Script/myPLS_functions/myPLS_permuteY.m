function Yp = myPLS_permuteY(Y, grouping, grouped_perm)
% This function permutes the Y matrix during permutation testing
%
% Inputs:
% - Y               : N x B vector, N is #subjects, B is #behavior/design
%                     scores, behavior/design data (normalized)
% - grouping        : N x 1 vector, subject grouping (e.g. diagnosis)
%                     e.g. [1,1,2] = subjects 1 and 2 belong to group 1,
%                     subject 3 belongs to group 2
% - grouped_perm    : binary variable indicating if groups should be
%                     considered during permutation of Y
%                     0 = permutations ignoring grouping
%                     1 = permutations within group
%
% Outputs:
% - Yp              : N x B vector, permuted behavior/design data

%% ========== Check if permuting within group ========== %%
if ~grouped_perm
    grouping = ones(size(grouping)); % if not permuting within group, treat all subjects as one group
end

%% ========== Get grouping information ========== %%
groupIDs = unique(grouping); % all group labels
nGroups = length(groupIDs);  % number of groups

%% ========== Main loop: permute within each group ========== %%
Yp = zeros(size(Y)); % initialize output
for iG = 1:nGroups
    groupID = find(grouping == groupIDs(iG));   % indices of subjects in this group

    % permute subjects within this group
    thisY = Y(groupID, :); % data for this group

    for iY = 1:size(Yp, 2)
        perm_order = randperm(size(thisY, 1)); % generate permutation order
        thisYp = thisY(perm_order, iY); % permuted data

        % store permuted data into output matrix
        Yp(groupID, iY) = thisYp;
    end
end
end