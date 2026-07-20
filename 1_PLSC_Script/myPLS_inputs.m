% Script for defining all input parameters for PLS medical image processing toolbox
% 
% Parameters to be defined in this script:
%
%   - input : structure containing input data for analysis
%       - .X             : N x M matrix, N = number of subjects, M = number of imaging variables
%       - .Y             : N x B matrix, B = number of behavioral variables
%       - .grouping      : N x 1 vector, subject grouping for PLS analysis
%                               e.g. [1,1,2] = subjects 1 and 2 belong to group 1, subject 3 to group 2
%       - [.group_names] : group names (optional)
%       - [.behav_names] : behavioral variable names (optional)
%   - pls_opts : PLS analysis options
%       - .nPerms              : number of permutations
%       - .nBootstraps         : number of bootstrap samples
%       - .normalization_img   : imaging data normalization option
%       - .normalization_behav : behavioral data normalization option
%              0 = no normalization
%              1 = z-score standardization across all subjects
%              2 = z-score standardization within groups (default)
%              3 = standard-deviation normalization across all subjects (no centering)
%              4 = standard-deviation normalization within groups (no centering)
%       - .grouped_PLS : whether to account for groups when computing R (binary)
%                                0 = compute PLS across all subjects
%                                1 = concatenate covariance matrices by group (standard behavioral PLS)
%       - .grouped_perm : whether to account for groups during permutations (binary)
%                  0 = ignore groups during permutations
%                  1 = permute within groups
%       - .grouped_boot : whether to account for groups during bootstrapping (binary)
%              0 = ignore groups during bootstrapping
%              1 = bootstrap within groups
%       - .boot_procrustes_mod : bootstrapping Procrustes transformation mode
%              1 = standard (rotate U only)
%              2 = average rotation of U and V
%       - .save_boot_resampling : whether to save bootstrap resampling data (binary)
%              0 = do not save
%              1 = save
%       - .behav_type        : behavioral analysis type
%              'behavior'  standard behavioral PLS
%              'contrast'  contrast between two groups only
%              'contrastBehav'  contrast combined with behavioral variables
%              'contrastBehavInteract'  also includes group-by-behavior interaction
%   - save_opts: result saving and plotting options
%       - .output_path   : path for saving results
%       - .prefix        : prefix for all result files (optional)
%       - .img_type      : result plotting style
%              'volume'  voxel-wise NIfTI data, display bootstrap ratio on brain maps
%              'corrMat' ROI correlation matrix, display bootstrap ratio on correlation matrix
%              'barPlot' vectorized brain data, display bootstrap ratio on bar plots
%       - .mask_file     : gray-matter mask, required only when img_type='volume'
%       - .BSR_thres     : 2×1 vector, bootstrap ratio visualization thresholds, required only when img_type='volume'
%       - .struct_file   : structural image filename, required only when img_type='volume'
%       - .load_thres    : 2×1 vector, loadings visualization thresholds, required only when img_type='volume'
%       - .grouped_plots : whether to account for groups when plotting (binary)
%              0 = ignore groups when plotting
%              1 = account for groups when plotting
%       - .alpha         : significance level for latent components
%       - .plot_boot_samples : whether to plot bootstrap samples in bar plots (binary)
%       - .errorbar_mode : 'std' plots standard deviation, 'CI' plots 95% confidence interval
%       - .hl_stable     : whether to highlight stable bootstrap scores (binary)


%% ---------- Input Data ----------
%%%%%%%% Load your data here %%%%%%%%
load('data.mat'); % Load your data file

X0 = brain_data;
Y0behav = beh_data;

% If you want to regress out confounding variables, process them here
% --- Brain imaging data ---
% X0 is typically the brain-imaging data matrix with dimensions: subjects (rows) × imaging features (columns)
input.brain_data = X0;

% --- Behavioral data ---
% Y0 is the behavioral data matrix with dimensions: subjects (rows) × behavioral features (columns)
% Construction of Y0 depends on pls_opts.behav_type:
%   * If behav_type = 'contrast', behavData can be empty; Y0 relies only on grouping information
%   * If behav_type = 'behavior', Y0 is identical to behavData
%   * If behav_type = 'contrastBehav' or 'contrastBehavInteract', Y0 is built from both grouping and behavioral data
input.behav_data = Y0behav;

% --- Grouping data ---
% subj_grouping: grouping vector, binary variable, multiple groups allowed
input.grouping = diagnosis;

% --- Group names ---
% Specify group names here for plotting purposes
input.group_names = {'HC'};

% --- Behavioral variable names ---
input.behav_names = {'val1', 'val2', 'val3', 'val4', 'val5', 'val6', 'val7', 'val8'};

% --- Imaging variable names ---
% Used only when save_opts.img_type = 'barPlot'
for i = 1:size(X0, 2)
    input.img_names{i, 1} = ['I ' num2str(i)];
end

clear diagnosis brain_data Fa_data i

%% ---------- PLS Options ----------

% --- Permutations & Bootstrapping ---
pls_opts.nPerms = 5000;
pls_opts.nBootstraps = 5000;

% --- Data normalization options ---
% 0: no normalization
% 1: z-score standardization across all subjects
% 2: z-score standardization within groups (default for grouped PLSC, see Krishnan et al., 2011)
% 3: standard-deviation normalization across all subjects (no centering)
% 4: standard-deviation normalization within groups (no centering)
pls_opts.normalization_img = 1;
pls_opts.normalization_behav = 1;

% --- PLS grouping option ---
% 0: compute PLS across all subjects
% 1: concatenate covariance matrices by group (standard behavioral PLS, see Krishnan et al., 2011)
pls_opts.grouped_PLS = 0; 

% --- Permutation grouping option ---
% 0: ignore groups during permutations
% 1: permute within groups
pls_opts.grouped_perm = 0;

% --- Bootstrapping grouping option ---
% 0: ignore groups during bootstrapping
% 1: bootstrap within groups
pls_opts.grouped_boot = 0;

% --- Bootstrapping Procrustes transformation mode ---
% In some cases, relying only on U rotation leads to extremely low standard errors and infinite bootstrap ratios
% Mode 2 computes transformation matrices for both U and V
% 1: standard
% 2: average rotation of U and V
pls_opts.boot_procrustes_mod = 2;

% --- Save bootstrap resampling data? ---
% Choose whether to save bootstrap resampling data (recommended only when imaging dimensions are small)
pls_opts.save_boot_resampling = 0;

% --- Behavioral analysis type ---
% 'behavior'  standard behavioral PLS
% 'contrast'  contrast between two groups only
% 'contrastBehav'  contrast combined with behavioral variables
% 'contrastBehavInteract'  also includes group-by-behavior interaction
pls_opts.behav_type = 'behavior';

%% ---------- Result Saving & Plotting Options ----------
% --- Output path for results ---
save_opts.output_path = '../Results_PLSC';

% --- Prefix for all result files ---
% If undefined, toolbox default prefix is used
save_opts.prefix = sprintf('myPLS_%s_norm%d-%d', pls_opts.behav_type, pls_opts.normalization_img, pls_opts.normalization_behav);

% --- Grouping option for plotting ---
% 0: ignore groups when plotting
% 1: account for groups when plotting
save_opts.grouped_plots = 0;

% --- Significance level for latent components ---
save_opts.alpha = 0.30;

% --- Brain data type ---
% Specify result plotting style
% 'volume' voxel-wise NIfTI data, display bootstrap ratio on brain maps
% 'corrMat' ROI correlation matrix, display bootstrap ratio on correlation matrix
% 'barPlot' vectorized brain data, display bootstrap ratio on bar plots

% % Uncomment below to view volume plotting example:
save_opts.img_type = 'barPlot'; 

% % Uncomment below to view correlation-matrix plotting example:
% input.brain_data = input.brain_data(:,1:30135); 
% save_opts.img_type = 'corrMat';

% Uncomment below to view bar-plot example:
%input.brain_data = input.brain_data(:, 1:200); 
save_opts.img_type = 'barPlot';
% save_opts.fig_pos_img = [440   606   560   192];
save_opts.fig_pos_img = [606, 400, 600, 400];

save_opts.fig_pos_img_1D = [606, 400, 2000, 500];


% --- Brain visualization thresholds ---
% (Bootstrap ratio and loadings thresholds, required only when img_type = 'volume' or 'corrMat')
save_opts.BSR_thres = [-2.3, 2.3]; % negative/positive thresholds for bootstrap ratio visualization
save_opts.load_thres = [-0.4, 0.4]; % negative/positive thresholds for loadings visualization

% --- Brain mask ---
% (gray-matter mask, required only when img_type = 'volume')
save_opts.mask_file = 'example_mask.nii'; % binary mask filename to constrain analysis range

% --- Structural template file for visualization ---
% (background structural image, required only when img_type='volume')
save_opts.struct_file = 'example_struct.nii';

% --- Slice orientation for volume plotting ---
% (required only when img_type='volume')
save_opts.volume_orientation = 'axial'; %'axial','coronal','sagittal'

% --- Bar-plot options ---
save_opts.plot_boot_samples = 0; % whether to plot bootstrap samples in bar plots (binary)
save_opts.errorbar_mode = 'CI'; % 'std' plots standard deviation; 'CI' plots 95% confidence interval
save_opts.hl_stable = 1; % whether to highlight stable bootstrap scores (binary)

% --- Custom dimensions for behavioral bar plots ---
save_opts.fig_pos_behav = [606, 400, 600, 400];
