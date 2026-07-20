function [Y0, design_names, nDesignScores] = myPLS_getY(behavMode, behaviorData, grouping, group_names, behav_names)
% =============================
% Function:
%   Generate design matrix Y0 and variable names for PLS analysis based on analysis mode (behavior/contrast/interaction).
%   Supports 2-group/3-group contrasts, automatic generation of behavior and interaction terms.
%
% =============================
% Inputs:
%   behavMode    : String, behavior analysis type
%                 'behavior': behavior variables only
%                 'contrast': group contrast
%                 'contrastBehav': behavior + contrast
%                 'contrastBehavInteract': behavior + contrast + interaction
%   behaviorData : N×B matrix, raw behavior/design data
%                 N: number of subjects, B: number of behavior/design variables
%   grouping     : N×1 vector, subject grouping information
%   group_names  : 1×G cell, group names (G is number of groups)
%   behav_names  : 1×B cell, behavior variable names
%
% =============================
% Outputs:
%   Y0           : N×D matrix, design matrix for PLS analysis
%                 D: number of design variables, depends on mode
%   design_names : 1×D cell, design variable names
%   nDesignScores: Integer, number of design variables
%
% =============================
% Example:
%   [Y0, design_names, nDesignScores] = myPLS_getY(behavMode, behaviorData, grouping, group_names, behav_names);
%
% =============================

%% ========== Get group information ========== %%
groupIDs = unique(grouping);         % groupIDs: unique IDs for all groups (e.g., [1 2])
nGroups  = length(groupIDs);         % nGroups: number of groups
nBehav   = size(behaviorData, 2);    % nBehav: number of behavioral variables
nSubj    = size(behaviorData, 1);    % nSubj: number of subjects

% Ensure names are row vectors for concatenation
if size(behav_names, 1) > size(behav_names, 2)
    behav_names = behav_names';
end
if size(group_names, 1) > size(group_names, 2)
    group_names = group_names';
end

%% ========== Generate contrast vector ========== %%
% If the analysis mode includes contrast, generate group contrast vector c and its name
if contains(behavMode, 'contrast')
    if nGroups == 1
        % Only one group, cannot perform contrast analysis
        error('Only one group, cannot perform contrast analysis');
    elseif nGroups == 2
        % Two-group contrast, generate one contrast vector c
        group1ID = grouping == groupIDs(1); nGroup1 = nnz(group1ID); % First group index and count
        group2ID = grouping == groupIDs(2); nGroup2 = nnz(group2ID); % Second group index and count
        c = zeros(nSubj, 1);              % Initialize contrast vector
        c(group1ID) = -1 / nGroup1;       % Assign negative weight to first group
        c(group2ID) = 1 / nGroup2;        % Assign positive weight to second group
        c_orig = c;                       % Keep original contrast vector
        c = orth(c_orig);                 % Orthogonalize to ensure orthogonality
        if sign(c(1)) ~= sign(c_orig(1)); c = -c; end % Ensure consistent direction
        contrastName = {['Contrast (' group_names{2} ' > ' group_names{1} ')']};
    elseif nGroups == 3
        % Three-group contrast, generate two contrast vectors c(:,1), c(:,2)
        group1ID = grouping == groupIDs(1); nGroup1 = nnz(group1ID);
        group2ID = grouping == groupIDs(2); nGroup2 = nnz(group2ID);
        group3ID = grouping == groupIDs(3); nGroup3 = nnz(group3ID);
        c = zeros(nSubj, 2);              % Initialize contrast matrix
        % First contrast: Group 3 vs Group 1/2
        c(group1ID, 1) = -1 / (nGroup1 + nGroup2);
        c(group2ID, 1) = -1 / (nGroup1 + nGroup2);
        c(group3ID, 1) = 1 / nGroup3;
        % Second contrast: Group 2 vs Group 1
        c(group1ID, 2) = -1 / nGroup1;
        c(group2ID, 2) = 1 / nGroup2;
        c_orig = c;                       % Keep original contrast matrix
        c(:, 1) = orth(c_orig(:, 1));
        c(:, 2) = orth(c_orig(:, 2));
        if sign(c(1, 1)) ~= sign(c_orig(1, 1))
            c(:, 1) = -c(:, 1);
        end
        if sign(c(1, 2)) ~= sign(c_orig(1, 2))
            c(:, 2) = -c(:, 2);
        end
        contrastName = {['Contrast (' group_names{3} '>' group_names{1} '/' group_names{2} ')'], ['Contrast (' group_names{2} '>' group_names{1} ')']};
    else
        % Only support 2-group or 3-group contrast, error if exceeded
        error('Only support 2-group or 3-group contrast');
    end
end

%% ========== Generate Y0 based on mode ========== %%
switch behavMode
    case 'behavior'
        % Behavior variables only, output raw behavior data directly
        Y0 = behaviorData;
        design_names = behav_names;
    case 'contrast'
        % Contrast only, output contrast vector
        Y0 = c;
        design_names = contrastName;
    case 'contrastBehav'
        % Behavior + contrast, concatenate contrast vector and behavior data
        Y0 = [c, behaviorData];
        design_names = [contrastName, behav_names];
    case 'contrastBehavInteract'
        % Behavior + contrast + interaction terms
        if nGroups == 2
            % For two groups, generate main effect and interaction term for each behavior variable
            Y0 = [];
            for iBeh = 1:nBehav
                Y0(:, 2 * iBeh - 1) = behaviorData(:, iBeh);      % Main effect
                Y0(:, 2 * iBeh)     = zscore(Y0(:, 2 * iBeh - 1));% Standardized main effect
                Y0(:, 2 * iBeh)     = Y0(:, 2 * iBeh) .* c(:, 1); % Interaction term: main effect * contrast
                design_names{2 * iBeh - 1} = behav_names{iBeh};
                design_names{2 * iBeh}     = [behav_names{iBeh} '*' contrastName{1}];
            end
        elseif nGroups == 3
            % For three groups, generate main effect and two interaction terms for each behavior variable
            Y0 = [];
            for iBeh = 1:length(behav_names)
                Y0(:, 3 * iBeh - 2) = behaviorData(:, iBeh);         % Main effect
                Y0(:, 3 * iBeh - 1) = -zscore(Y0(:, 3 * iBeh - 2));  % Standardized main effect (negative)
                Y0(:, 3 * iBeh)     = zscore(Y0(:, 3 * iBeh - 2));   % Standardized main effect (positive)
                Y0(:, 3 * iBeh - 1) = Y0(:, 3 * iBeh - 1) .* c(:, 1);% Interaction term 1
                Y0(:, 3 * iBeh)     = Y0(:, 3 * iBeh) .* c(:, 2);    % Interaction term 2
                design_names{3 * iBeh - 2} = behav_names{iBeh};
                design_names{3 * iBeh - 1} = [behav_names{iBeh} '*' contrastName{1}];
                design_names{3 * iBeh}     = [behav_names{iBeh} '*' contrastName{2}];
            end
        end
        % Concatenate contrast vector and interaction terms
        Y0 = [c, Y0];
        design_names = [contrastName, design_names];
    otherwise
        % Undefined behavior type, throw error
        error('Undefined behavior type')
end

%% ========== Get number of design scores ========== %%
nDesignScores = size(Y0, 2); % Number of design scores

end