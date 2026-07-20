function [input, pls_opts, save_opts] = myPLS_initialize(input, pls_opts, save_opts)
% =============================
% Function:
%   Initialize input, parameters, and save options for PLS analysis. Set default values,
%   check input validity, and auto-complete missing fields to ensure smooth subsequent analysis.
%   Supports multiple input formats, grouping, behavior/contrast/interaction analysis,
%   normalization, plotting, and save parameters.
%
% =============================
% Inputs:
%   input    : Structure containing analysis input data, must include:
%       .brain_data/X0   : N×M matrix, raw imaging data
%       .behav_data/Y0   : N×B matrix, raw behavioral/design data
%       .grouping        : N×1 vector, subject grouping information
%       .group_names     : 1×G cell, group names (optional)
%       .behav_names     : 1×B cell, behavioral variable names (optional)
%       .img_names       : 1×M cell, imaging variable names (optional)
%   pls_opts : Structure, PLS analysis parameters, must include:
%       .behav_type           : Behavior analysis type (behavior/contrast/contrastBehav/contrastBehavInteract)
%       .grouped_PLS          : Whether to use grouped PLS (0/1)
%       .grouped_perm         : Whether permutation test is within-group (0/1)
%       .grouped_boot         : Whether bootstrap sampling is within-group (0/1)
%       .boot_procrustes_mod  : Bootstrap Procrustes rotation mode (1/2)
%       .save_boot_resampling : Whether to save bootstrap samples (0/1)
%       .normalization_img    : Imaging data normalization method (0~4, see myPLS_bootstrapping.m)
%       .normalization_behav  : Behavioral data normalization method (0~4, see myPLS_bootstrapping.m)
%   save_opts : Structure, save and plot parameters, must include:
%       .prefix        : String, result file name prefix
%       .output_path   : String, output directory
%       .alpha         : Significance level
%       .img_type      : Imaging type (volume/corrMat/barPlot)
%       .mask_file     : Mask file (required for volume)
%       .struct_file   : Structural image file (required for volume)
%       .BSR_thres     : BSR threshold
%       .load_thres    : Loading threshold
%       .grouped_plots : Whether to plot by group (0/1)
%       .plot_boot_samples : Whether to plot bootstrap samples (0/1)
%
% =============================
% Outputs:
%   input     : Structure, completed and standardized input data
%   pls_opts  : Structure, completed and standardized PLS parameters
%   save_opts : Structure, completed and standardized save and plot parameters
%
% =============================
% Example:
%   [input, pls_opts, save_opts] = myPLS_initialize(input, pls_opts, save_opts);
%
% =============================

% This function sets default values and checks the validity of inputs for myPLS_analysis

if nargin < 3 || isempty(save_opts)
    save_opts = struct();
end

% =============================
% Main program flow begins
% =============================

disp('... Initializing input parameters ...')

%% ========== 1. Input data validation and completion ========== %%
% Compatible with X0 and Y0 inputs, allowing field names X0/Y0 or brain_data/behav_data in input struct
% Reason: Some users may name input fields as X0/Y0, auto-completion needed for compatibility
if ~isfield(input, 'brain_data') && isfield(input, 'X0')
    % If brain_data not provided but X0 exists, auto-complete
    input.brain_data = input.X0;
elseif ~isfield(input, 'brain_data') && ~isfield(input, 'X0')
    % Reason: If neither brain_data nor X0 exists, imaging data is missing, cannot proceed
    error('Missing imaging data input (brain_data or X0 field)');
end

% Reason: Some users may name input fields as Y0, auto-completion needed for compatibility
if ~isfield(input, 'behav_data') && isfield(input, 'Y0')
    % If behav_data not provided but Y0 exists, auto-complete
    input.behav_data = input.Y0;
elseif ~isfield(input, 'behav_data') && ~isfield(input, 'Y0')
    % Reason: If neither behav_data nor Y0 exists, behavioral data is missing, cannot proceed
    error('Missing behavioral data input (behav_data or Y0 field)');
end

% Reason: PLS analysis requires consistent sample sizes for X and Y, otherwise data mismatch
if (size(input.brain_data, 1) ~= size(input.behav_data, 1))
    error('Sample size mismatch between imaging data (X) and behavioral data (Y), please check input!');
end

% ========== Get number of groups, behavioral variables, and imaging variables ========== %
% groupIDs: unique IDs for all groups (e.g., [1 2])
groupIDs = unique(input.grouping);    % Group ID vector
nGroups  = length(groupIDs);           % Number of groups
nBehav   = size(input.behav_data, 2);  % Number of behavioral variables
nImg     = size(input.brain_data, 2);  % Number of imaging variables

%% ========== Group names and behavioral variable names completion ========== %%
% Reason: If the number of group_names is less than the actual number of groups, information is incomplete and needs to be reset to standard naming
if isfield(input, 'group_names') && numel(input.group_names) < nGroups
    disp('!!! Number of group names is less than actual number of groups, will use standard naming automatically!');
    input.group_names = []; % Clear, will auto-complete later
end

% Reason: If group_names is not provided or empty, auto-generate standard group names to ensure subsequent workflow availability
if ~isfield(input, 'group_names') || isempty(input.group_names)
    input.group_names = cell(nGroups, 1);
    for iG = 1:nGroups
        input.group_names{iG} = ['Group' num2str(groupIDs(iG))];
    end
end

% Reason: If the number of group_names exceeds the number of groups, it is redundant and needs to be truncated
if numel(input.group_names) > nGroups
    disp('!!! Number of group names exceeds actual number of groups, only keeping the first few:');
    input.group_names = input.group_names(1:nGroups);
    for iG = 1:nGroups
        disp(['   Group' num2str(iG) ': "' input.group_names{iG} '"']);
    end
end

% Reason: If the number of behav_names is less than the actual number of behavioral variables, information is incomplete and needs to be reset to standard naming
if isfield(input, 'behav_names') && numel(input.behav_names) < nBehav
    disp('!!! Number of behavioral variable names is less than actual number of variables, will use standard naming automatically!');
    input.behav_names = [];
end

% Reason: If behav_names is not provided or empty, auto-generate standard behavioral variable names to ensure subsequent workflow availability
if ~isfield(input, 'behav_names') || isempty(input.behav_names)
    input.behav_names = cell(nBehav, 1);
    for iB = 1:nBehav
        input.behav_names{iB} = ['Behavior' num2str(iB)];
    end
end

% Reason: If the number of behav_names exceeds the number of behavioral variables, it is redundant and needs to be truncated
if numel(input.behav_names) > nBehav
    disp('!!! Number of behavioral variable names exceeds actual number of variables, only keeping the first few:');
    input.behav_names = input.behav_names(1:nBehav);
    for iB = 1:nBehav
        disp(['   Behavior' num2str(iB) ': "' input.behav_names{iB} '"']);
    end
end


%% ========== 2. PLS parameter completion and validation ========== %%
% Reason: If analysis type not specified, default to behavior PLS to ensure workflow availability
if ~isfield(pls_opts, 'behav_type') || isempty(pls_opts.behav_type)
    pls_opts.behav_type = 'behavior';
elseif ~(strcmp(pls_opts.behav_type, 'behavior') || strcmp(pls_opts.behav_type, 'contrast') || strcmp(pls_opts.behav_type, 'contrastBehav') || strcmp(pls_opts.behav_type, 'contrastBehavInteract'))
    % Reason: If analysis type not in allowed range, error to prevent subsequent errors
    error('Invalid behavior analysis type (behav_type), please check parameters!')
end

% Reason: If PLS grouping not specified, auto-set based on analysis type to ensure reasonable parameters
if ~isfield(pls_opts, 'grouped_PLS')
    disp('PLS grouping option not specified, will auto-set based on analysis type:')
    if contains(pls_opts.behav_type, 'contrast')
        pls_opts.grouped_PLS = 0;
        disp('   Analysis type is contrast, PLS will not consider grouping')
    else
        pls_opts.grouped_PLS = 1;
        disp('   Analysis type is behavior, PLS will proceed by group')
    end
end

% Reason: For contrast PLS, within-group normalization cannot be selected, and grouped PLS cannot be used, otherwise mathematical meaning is inconsistent
if contains(pls_opts.behav_type, 'contrast')
    if pls_opts.normalization_behav == 2 || pls_opts.normalization_behav == 4 || pls_opts.normalization_img == 2 || pls_opts.normalization_img == 4
        error('For contrast PLS analysis, within-group normalization cannot be selected, please set normalization method to all subjects!')
    end
    if pls_opts.grouped_PLS == 1
        error('For contrast PLS analysis, grouped PLS cannot be selected, please set grouped_PLS parameter to 0!')
    end
end

% Reason: If permutation test grouping not specified, auto-set based on analysis type to ensure reasonable parameters
if ~isfield(pls_opts, 'grouped_perm') || isempty(pls_opts.grouped_perm)
    disp('Permutation test grouping option not specified, will auto-set based on analysis type:')
    if contains(pls_opts.behav_type, 'contrast')
        pls_opts.grouped_perm = 0;
        disp('   Analysis type is contrast, permutation test will not consider grouping')
    else
        pls_opts.grouped_perm = 1;
        disp('   Analysis type is behavior, permutation test will proceed by group')
    end
end

% Reason: If bootstrap sampling grouping not specified, auto-set based on analysis type to ensure reasonable parameters
if ~isfield(pls_opts, 'grouped_boot') || isempty(pls_opts.grouped_boot)
    disp('Bootstrap sampling grouping option not specified, will auto-set based on analysis type:')
    if contains(pls_opts.behav_type, 'contrast')
        pls_opts.grouped_boot = 0;
        disp('   Analysis type is contrast, bootstrap sampling will not consider grouping')
    else
        pls_opts.grouped_boot = 1;
        disp('   Analysis type is behavior, bootstrap sampling will proceed by group')
    end
end

% Reason: When grouping options are inconsistent, warn user to avoid analysis logic confusion
if pls_opts.grouped_boot ~= pls_opts.grouped_perm
    disp('!!! Grouping options for permutation test and bootstrap sampling are inconsistent, please confirm analysis requirements!')
end
if pls_opts.grouped_PLS ~= pls_opts.grouped_perm
    disp('!!! Grouping options for PLS and permutation test are inconsistent, please confirm analysis requirements!')
end

% Reason: If Procrustes transformation mode not specified, default to 1 to ensure workflow availability
if ~isfield(pls_opts, 'boot_procrustes_mod') || isempty(pls_opts.boot_procrustes_mod)
    disp('Bootstrap Procrustes transformation mode not specified, using default value 1: only apply Procrustes transformation to U (behavior weights)')
    pls_opts.boot_procrustes_mod = 1;
end

% Reason: Procrustes mode only allows 1 or 2, other values will cause error
if pls_opts.boot_procrustes_mod ~= 1 && pls_opts.boot_procrustes_mod ~= 2
    error('Invalid bootstrap Procrustes transformation mode (boot_procrustes_mod), must be 1 or 2!');
end

% Reason: If whether to save bootstrap samples not specified, auto-set based on number of imaging variables to prevent storage pressure from large data
if ~isfield(pls_opts, 'save_boot_resampling') || isempty(pls_opts.save_boot_resampling)
    disp('Whether to save bootstrap samples not specified, will auto-set based on number of imaging variables:')
    if nImg > 1000
        pls_opts.save_boot_resampling = 0;
        disp('   Number of imaging variables greater than 1000, will not save bootstrap samples')
    else
        pls_opts.save_boot_resampling = 1;
        disp('   Number of imaging variables not greater than 1000, will save bootstrap samples')
    end
end

% Reason: If user forces to save bootstrap samples, number of imaging variables must be less than 1000 to prevent file size from growing too large
if pls_opts.save_boot_resampling && nImg > 1000
    disp('!!! Number of imaging variables greater than 1000, may cause file size to grow too large!')
end

disp(' ')

% =============================
% Program main process end
% =============================

end