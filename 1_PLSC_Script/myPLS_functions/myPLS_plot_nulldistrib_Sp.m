function myPLS_plot_nulldistrib_Sp(S, Sp_vect, signif_LC, save_opts)
% This function plots the null distribution of singular values obtained
% with permutation testing
%
% Inputs:
% - S          : L x L matrix, singular values (diagonal matrix)
% - Sp_vect    : matrix with permuted singular values
% - signif_LC      : significant latent components to plot
% - save_opts      : Options for result saving and plotting
%       - .output_path : output directory where figures are saved
%       - .prefix      : prefix of all results files (optional)
%% ========== Extract singular values for all LCs ========== %%
S_vect = diag(S); % S_vect: L×1 vector of actual singular values for all LCs

%% ========== Loop to plot null distribution for each significant LC ========== %%
for iLC = 1:length(signif_LC)
    this_lc = signif_LC(iLC); % index of the current LC to plot

    this_S = S_vect(this_lc);           % actual singular value of the current LC
    this_Sp_vect = Sp_vect(this_lc, :); % permutation distribution for the current LC (1×Nperm)

    % Build output file name
    file_name = fullfile(save_opts.output_path, [save_opts.prefix, '_LC', num2str(this_lc), '_nullDistribPermutedSingVals']);

    %% ========== Plot histogram ========== %%
    figure('Position', save_opts.fig_pos_img, 'Visible', 'off'); % create new figure
    histogram(this_Sp_vect, 50, 'FaceColor', [0, 0.44 0.74]); % plot null distribution histogram
    title(['LC' num2str(this_lc) ' - Null distribution of singular values'], 'FontSize', 14); % figure title
    set(gca, 'FontSize', 12, 'Box', 'off'); % axis styling
    set(gcf, 'Color', 'w'); % white background

    % Set x-axis range: slightly wider on the left, extend to the right of the actual singular value
    a = std(this_Sp_vect); % standard deviation of null distribution
    xlim([min(this_Sp_vect) - a, this_S + a]);

    hold on
    % Draw vertical line for the actual singular value
    line([this_S, this_S], ylim, 'LineStyle', ':', 'LineWidth', 1, 'Color', 'r');
    hold off

    % Save figure as tif with 300 dpi resolution
    print(gcf, [file_name, '.tif'], '-dpng', '-r300');

    close gcf
end
end