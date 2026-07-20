function [Y0, design_names, nDesignScores] = myPLS_getY(behavMode, behaviorData, grouping, group_names, behav_names)
% Generate design matrix Y0 and variable names required for PLS analysis according to the analysis mode (behavior / contrast / interaction).
% Supports automatic generation of 2-group / 3-group contrasts, behavioral and interaction terms.
%
% =============================
% Input:
%   behavMode    : string, type of behavioral analysis
%                 'behavior'              : behavioral variables only
%                 'contrast'              : group contrast
%                 'contrastBehav'         : behavior + contrast
%                 'contrastBehavInteract' : behavior + contrast + interaction
%   behaviorData : N×B matrix, raw behavioral / design data
%                 N: number of subjects, B: number of behavioral / design variables
%   grouping     : N×1 vector, subject grouping information
%   group_names  : 1×G cell, group names (G: number of groups)
%   behav_names  : 1×B cell, behavioral variable names
%
% =============================
% Output:
%   Y0           : N×D matrix, design matrix for PLS analysis
%                 D: number of design variables, depends on mode
%   design_names : 1×D cell, design variable names
%   nDesignScores: integer, number of design variables

%% ========== Retrieve grouping information ========== %%
% Number and IDs of groups
groupIDs = unique(grouping);
nGroups = length(groupIDs);

% Number of behavior scores
nBehav = size(behaviorData,2);

% Number of subjects
nSubj = size(behaviorData,1); 

% Ensure behav_names and group_names are row vectors (1xN)
if size(behav_names, 1) > size(behav_names, 2)
    behav_names = behav_names';
end
if size(group_names, 1) > size(group_names, 2)
    group_names = group_names';
end

%% ========== Compute contrast (if necessary) ========== %%
% Contrasts are implemented for up to 3 groups, for 3-group analysis see
% Zoller et al., Schizophrenia Research 2018

if contains(behavMode, 'contrast')
    if nGroups == 1
        error('There is only one group - contrasts are only implemented for 2 or 3 groups');
    elseif nGroups == 2
        group1ID = grouping == groupIDs(1); nGroup1 = nnz(group1ID); % indices and size of group 1
        group2ID = grouping == groupIDs(2); nGroup2 = nnz(group2ID); % indices and size of group 2
        c = zeros(nSubj, 1);              % initialize contrast vector
        c(group1ID) = -1 / nGroup1;       % assign negative weight to group 1
        c(group2ID) = 1 / nGroup2;        % assign positive weight to group 2
        c_orig = c;                       % keep original contrast vector
        c = orth(c_orig);                 % orthogonalize to ensure orthogonality
        if sign(c(1)) ~= sign(c_orig(1)); c = -c; end % ensure consistent direction
        contrastName = {['contrast (' group_names{2} ' > ' group_names{1} ')']};
    elseif nGroups == 3
        % three-group contrast, generate two contrast vectors c(:,1), c(:,2)
        group1ID = grouping == groupIDs(1); nGroup1 = nnz(group1ID);
        group2ID = grouping == groupIDs(2); nGroup2 = nnz(group2ID);
        group3ID = grouping == groupIDs(3); nGroup3 = nnz(group3ID);
        c = zeros(nSubj, 2);              % initialize contrast matrix
        % first contrast: group 3 vs groups 1&2
        c(group1ID, 1) = -1 / (nGroup1 + nGroup2);
        c(group2ID, 1) = -1 / (nGroup1 + nGroup2);
        c(group3ID, 1) = 1 / nGroup3;
        % second contrast: group 2 vs group 1
        c(group1ID, 2) = -1 / nGroup1;
        c(group2ID, 2) = 1 / nGroup2;
        c_orig = c;                       % keep original contrast matrix
        c(:, 1) = orth(c_orig(:, 1));
        c(:, 2) = orth(c_orig(:, 2));
        if sign(c(1, 1)) ~= sign(c_orig(1, 1))
            c(:, 1) = -c(:, 1);
        end
        if sign(c(1, 2)) ~= sign(c_orig(1, 2))
            c(:, 2) = -c(:, 2);
        end
        contrastName = {['contrast (' group_names{3} '>' group_names{1} '/' group_names{2} ')'], ...
                        ['contrast (' group_names{2} '>' group_names{1} ')']};
    else
        % only 2 or 3 groups supported, throw error otherwise
        error('Only 2 or 3 group contrasts are supported');
    end
end

%% ========== Specify Y0 according to Mode ========== %%
switch behavMode
    case 'behavior'
        % behavior variables only, output original behavioral data directly
        Y0 = behaviorData;
        design_names = behav_names;
    case 'contrast'
        % contrast only, output contrast vector
        Y0 = c;
        design_names = contrastName;
    case 'contrastBehav'
        % behavior + contrast, concatenate contrast vector and behavioral data
        Y0 = [c, behaviorData];
        design_names = [contrastName, behav_names];
    case 'contrastBehavInteract'
        % behavior + contrast + interaction terms
        if nGroups == 2
            % for two groups, generate main effect and interaction for each behavioral variable
            Y0 = [];
            for iBeh = 1:nBehav
                Y0(:, 2 * iBeh - 1) = behaviorData(:, iBeh);      % main effect
                Y0(:, 2 * iBeh)     = zscore(Y0(:, 2 * iBeh - 1));% standardize main effect
                Y0(:, 2 * iBeh)     = Y0(:, 2 * iBeh) .* c(:, 1); % interaction: main effect * contrast
                design_names{2 * iBeh - 1} = behav_names{iBeh};
                design_names{2 * iBeh}     = [behav_names{iBeh} '*' contrastName{1}];
            end
        elseif nGroups == 3
            % for three groups, generate main effect and two interactions for each behavioral variable
            Y0 = [];
            for iBeh = 1:length(behav_names)
                Y0(:, 3 * iBeh - 2) = behaviorData(:, iBeh);         % main effect
                Y0(:, 3 * iBeh - 1) = -zscore(Y0(:, 3 * iBeh - 2));  % standardized main effect (negative)
                Y0(:, 3 * iBeh)     = zscore(Y0(:, 3 * iBeh - 2));   % standardized main effect (positive)
                Y0(:, 3 * iBeh - 1) = Y0(:, 3 * iBeh - 1) .* c(:, 1);% interaction 1
                Y0(:, 3 * iBeh)     = Y0(:, 3 * iBeh) .* c(:, 2);    % interaction 2
                design_names{3 * iBeh - 2} = behav_names{iBeh};
                design_names{3 * iBeh - 1} = [behav_names{iBeh} '*' contrastName{1}];
                design_names{3 * iBeh}     = [behav_names{iBeh} '*' contrastName{2}];
            end
        end
        % concatenate contrast vector and interaction terms
        Y0 = [c, Y0];
        design_names = [contrastName, design_names];
    otherwise
        % undefined behavior type, throw error
        error('Undefined behavior type')
end

%% ========== Get number of design scores ========== %%
nDesignScores = size(Y0, 2);
end