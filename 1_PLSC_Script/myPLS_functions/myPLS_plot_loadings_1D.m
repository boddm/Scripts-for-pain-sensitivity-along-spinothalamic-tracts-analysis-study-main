function myPLS_plot_loadings_1D(var_type, var_PLS_type, vars_PLS, vars_b_vect, vars_mean, vars_std, vars_lB, vars_uB, var_names, signif_LC, nGroups, fig_pos, save_opts)
% This function plots one-dimensional loadings (e.g. behavior, 1D imaging
% metrics) of all significant latent components (LCs) as barcharts.
%
% Inputs:
% - var_type       : string, type of variable plotted ('Imaging', 'Design' or 'Behavior')
% - var_PLS_type   : name of the PLS variable type for y-axis and for filenames
% - vars_PLS       : B x L matrix, B is #variables, L is #components
% - vars_b_vect    : B x L x P matrix, P is #bootstrap samples
% - vars_mean      : B x L matrix, bootstrapping mean
% - vars_std       : B x L matrix, bootstrapping standard deviation
% - vars_lB        : B x L matrix, lower bound of 95% CIs
% - vars_uB        : B x L matrix, upper bound of 95% CIs
% - var_names      : string, names of variables
% - signif_LC      : significant LCs to plot (e.g. [1,2])
% - nGroups        : number of groups used for PLS analysis, will determine
%                    if each group has its own set of saliences
% - fig_pos        : position of the figure to plot
% - save_opts: Options for results saving and plotting
%       - .output_path   : path where to save the results
%       - .prefix        : prefix of all results files
%       - .plot_boot_samples : binary variable indicating if bootstrap
%                          samples should be plotted in bar plots
%       - .errorbar_mode : 'std' = plotting standard deviations
%                          'CI' = plotting 95% confidence intervals
%       - .hl_stable	 : binary variable indicating if stable bootstrap
%                          scores should be highlighted
%       - .grouped_plots : binary variable indicating if groups should be 
%                          considered during plotting
%              0 = plotting ignoring grouping
%              1 = plotting cosidering grouping

%% ========== Variable Initialization ========== %%
% nVars: total number of variables (e.g., behavioral variables, imaging metrics)
nVars = size(vars_PLS, 1);          
% nBootstraps: number of bootstrap samples
nBootstraps = size(vars_b_vect, 3); 
% nLCs: number of latent components (LCs)
nLCs = size(vars_PLS, 2);           
% nBehav: number of variables per group (e.g., behavioral variables grouped; for Imaging, all variables)
nBehav = nVars / nGroups;

%% ========== Set plot colours ========== %%
% Configure bar-chart colour scheme according to variable type and number of groups
if strcmp(var_type, 'Imaging')
    plot_col = [0, 0.5, 0.5]; % Single colour for Imaging
    if nGroups > 1
        error('Imaging data do not support multi-group display');
    end
end
if strcmp(var_type, 'Behavior') || strcmp(var_type, 'Design')
    plot_col =  [0, 0, 0.5;    % Group 1 – blue
                 1, 0.4, 0.4;  % Group 2 – red
                 0, 0.5, 0.5]; % Group 3 – cyan
    plot_col = plot_col(1:nGroups, :); % Keep only the first nGroups colours
    if nGroups > 3
        error('Multi-group plotting supports a maximum of 3 groups');
    end
end

%% ========== Reorder to group variables of the same type across groups ========== %%
% Purpose: place bars of the same variable from different groups adjacently for easier comparison
iter = 1;
for iB = 1:nBehav
    re_order(iter:iter+nGroups-1) = iB:nBehav:nVars;
    iter = iter+nGroups;
end

%% ========== Main loop: plot bar charts for each significant LC ========== %%
for iLC = 1:length(signif_LC)
    this_lc = signif_LC(iLC); % Current LC index
    % Build output file name, including prefix, LC number, variable type, and PLS type
    file_name = fullfile(save_opts.output_path, [save_opts.prefix, '_LC', num2str(this_lc), '_', var_type, '_', var_PLS_type]);
    figure('position', fig_pos, 'Visible', 'off'); % Create new figure window
    hold on;
    
    %%% Highlight stable bootstrap score regions (if CI does not cross zero)
    if save_opts.hl_stable
        % Calculate max y-axis range for highlighting background patches
        max1 = max(abs(vars_lB(:, this_lc)));
        max2 = max(abs(vars_uB(:, this_lc)));
        maxPatch = max(max1, max2) * 1.1;
        for iC = 1:nVars
            this_iC = find(re_order == iC); % Current variable index in re-ordered position
            % Check if CI is entirely positive or negative, decide whether to highlight
            if vars_mean(iC, this_lc) < 0
                if vars_uB(iC, this_lc) < 0
                    % Entire CI is negative, highlight
                    patch(this_iC + .5 * [1, -1, -1, 1], 2 * maxPatch * [-1, -1, 1, 1], ...
                        [1, 1, 0.8], 'edgecolor', 'none', 'LineStyle', 'none');
                end
            else
                if vars_lB(iC, this_lc) > 0
                    % Entire CI is positive, highlight
                    patch(this_iC + .5 * [1, -1, -1, 1], 2 * maxPatch * [-1, -1, 1, 1], ...
                        [1, 1, 0.8], 'edgecolor', 'none', 'LineStyle', 'none');
                end
            end
        end
        clearvars max*
        % Draw y = 0 reference line
        hl = refline(0, 0); % Horizontal line y=0
        hl.Color = 'k';
        hl.LineWidth = 1;
    end
    
    %% Create bar chart for main effect means
    % vars_mean(re_order, this_lc)：re-ordered mean for current LC
    b = bar(vars_mean(re_order, this_lc), 'FaceAlpha', .5); 
    b.FaceColor = 'flat';
    % Set bar colors by group
    for iG = 1:nGroups
        b.CData(iG:nGroups:end, :) = repmat(plot_col(iG, :), size(b.CData, 1) / nGroups, 1);
    end
    
    %% Create scatter plot for bootstrap samples (if required)
    if save_opts.plot_boot_samples
        % rand_x: generate jittered x-coordinates for each variable to avoid scatter overlap
        rand_x = (1:nVars) + 0.08 * randn(nBootstraps, 1); % add x-axis jitter to prevent overlap
        % vars_tmp: bootstrap distribution for all variables of current LC, reshaped to P×B
        vars_tmp = squeeze(vars_b_vect(re_order, this_lc, :))'; % reshape from B×P to P×B
        % plot scatter points, color assigned by group
        s = scatter(rand_x(:), vars_tmp(:), 12, repmat(repelem(plot_col, nBootstraps, 1), nBehav, 1), 'filled', 'MarkerFaceAlpha', .5);
    end
    
    %% Add error bars (standard deviation or confidence intervals)
    switch save_opts.errorbar_mode
        case 'CI'
            % CI error bars; ploterr allows separate upper/lower bounds
            h = ploterr(1:nVars, vars_mean(re_order, this_lc), [], {vars_lB(re_order, this_lc), vars_uB(re_order, this_lc)}, 'k.', 'abshhxy', 0.2);
        case 'std'
            % Standard-deviation error bars
            h = ploterr(1:nVars, vars_mean(re_order, this_lc), [], vars_std(re_order, this_lc), 'k.', 'abshhxy', 0.2);
    end
    set(h(1), 'marker', 'none'); % Remove marker symbols
    set(h(2), 'LineWidth', 0.5); % Error-bar line width
    
    %% Set axes, labels and title
    if exist('s', 'var')
        % If scatter points exist, auto-adjust y-axis range based on scatter
        y_lim = [min(s.YData) - 0.05, max(s.YData) + 0.05];
    else
        % Otherwise auto-adjust based on error bars
        y_dat = [h(1).YData h(2).YData];
        y_lim = [min(y_dat) - 0.05, max(y_dat) + 0.05];
    end
    
    xlims = get(gca, 'xlim');
    set(gca, 'Layer', 'Top', 'ylim', y_lim, 'xlim', [xlims(1) + 0.5, xlims(2) - 0.5]); 
    xticks(1:nVars);
    xticklabels(var_names);

    if strcmp(var_type, 'Imaging')
        set(gca, 'TickLabelInterpreter', 'none', 'FontSize', 8, 'Box', 'off');
    elseif strcmp(var_type, 'Behavior') || strcmp(var_type, 'Design')
        set(gca, 'TickLabelInterpreter', 'none', 'FontSize', 10, 'Box', 'off');
    end

    set(gcf, 'Color', 'w');
    xtickangle(60);
    xlabel([var_type, ' variables']);
    ylabel(var_PLS_type);
    title(['LC', num2str(this_lc), ' - ', var_type, ' ', var_PLS_type], 'FontSize', 16);
    
    %% Save figure
    % Save as tif format, 300 dpi resolution
    print(gcf, [file_name, '.tif'], '-dpng', '-r300');
    
    close gcf

    %% Write table
    % Call myPLS_table_loadings to write loadings, errors, etc. for current LC
    myPLS_table_loadings(var_type, var_PLS_type, vars_PLS, vars_std, vars_lB, vars_uB, save_opts, var_names, signif_LC)
end

end