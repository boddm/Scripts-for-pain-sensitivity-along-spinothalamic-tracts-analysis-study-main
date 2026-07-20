function myPLS_plot_loadings_2D(var, var_name, signif_LC, pos_threshold, neg_threshold, save_opts)
% This function plots 2-dimensional loadings (e.g. correlation matrix)
% of all significant latent components (LCs)
%
% Inputs:
% - var            : M x L matrix, M is #imaging values (when vectorized),
%                     L is #components
% - var_name       : name of the variable for filenames
% - signif_LC      : significant latent components to plot
% - pos_threshold  : threshold for loadings' visualization (positive values)
% - neg_threshold  : threshold for loadings' visualization (negative values)
% - save_opts      : Options for result saving and plotting
%       - .output_path : output directory where figures are saved
%       - .prefix      : prefix of all results files (optional)

%% ========== Automatically determine the number of ROIs (assuming input is vectorized upper-triangle) ========== %%
nImg = size(var, 1); % Number of variables (length of vector)

nRois = 1;
while cumsum(1:nRois) <= nImg
    nRois = nRois + 1;
end
% nRois is the dimension of the correlation matrix (number of regions), inferred automatically

%% ========== Loop to plot 2D loadings for each significant LC ========== %%
for iter_lc = 1:numel(signif_LC)
    this_lc = signif_LC(iter_lc); % current LC index
    
    % Build output file name
    file_name = fullfile(save_opts.output_path, [save_opts.prefix '_LC' num2str(this_lc) '_' var_name]);
    
    %% ========== Restore as symmetric correlation matrix ========== %%
    % var(:,this_lc) is the vectorized loadings for the current LC; jVecToSymmetricMat converts it to nRois×nRois symmetric matrix
    CM = jVecToSymmetricMat(var(:, this_lc), nRois, 1);
    
    %% ========== Thresholding ========== %%
    % Zero-out positive loadings below pos_threshold and negative loadings above neg_threshold to highlight significant connections
    CM(CM > 0 & CM < pos_threshold) = 0;
    CM(CM < 0 & CM > neg_threshold) = 0;
    
    %% ========== Plotting ========== %%
    figure;
    imagesc(CM);
    colormap('jet');
    colorbar;
    xlabel('ROIs');
    ylabel('ROIs');
    set(gca, 'TickLabelInterpreter', 'none', 'FontSize', 6, 'Box', 'off');
    set(gcf, 'Color', 'w');
    title(['LC', num2str(this_lc), ' - ', var_name]);
    
    %% ========== Save figure ========== %%
    % Display progress: indicate save path
    fprintf('Saving 2D loadings plot for LC%d to: %s.jpg\n', this_lc, file_name);
    saveas(gcf, [file_name, '.jpg']);
    close(gcf);
    
end