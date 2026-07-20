function corr_LxLy = myPLS_plot_subjScores(Lx, Ly, names_groups, grouping, signif_LC, save_opts)
% This function plots the imaging & behavior scores of all significant 
% latent components (LCs). Subjects' scores are shown according to their 
% diagnostic group. Correlations between imaging & behavior scores is displayed.
%
% Inputs:
% - Lx            : N x L matrix, N is #subjects, L#components, imaging scores
% - Ly            : N x L matrix, behavioral  scores
% - names_groups  : string, names of diagnostic groups
% - grouping      : N x 1 vector, subject (diagnostic) grouping 
%                   e.g. [1,1,2] = subjects 1&2 belong to group 1,
%                   subject 3 belongs to group 2.
% - signif_LC     : significant latent components to plot, e.g. [1,2]
%
% Outputs:
% - corr_LxLy     : correlation between imaging and behavior scores
%% ========== Get grouping information ========== %%
groupIDs = unique(grouping); 
nGroups  = length(groupIDs); 

colors = {'b', 'r', 'c', 'g', 'm', 'y', 'w', 'k'}; % Matlab built-in colors
plot_colors = colors(1:nGroups);     % Select colors matching the number of groups

%% ========== Loop to plot scatter for each significant LC ========== %%
disp('Correlation between imaging and behavior scores');
for iter_lc = 1:length(signif_LC)
    this_lc = signif_LC(iter_lc); % current LC index to plot
    
    figure('position', save_opts.fig_pos_img, 'Visible', 'off'); % create new figure
    % plot scatter by group
    for iG = 1:nGroups
        % extract subject indices for current group
        idx = grouping == groupIDs(iG);
        % plot scatter for current group, x=imaging score, y=behavior score
        plot(Lx(idx, this_lc), Ly(idx, this_lc), [plot_colors{iG} '.'], 'MarkerSize', 12);
        hold on
    end
    hold off
    
    % set title, legend, axis labels
    title(['LC' num2str(this_lc) ' - Correlation between imaging and behavior scores'], 'FontSize', 16);
    legend(names_groups, 'Location', 'southeast');
    xlabel('Imaging score');
    ylabel('Behavior/design score');
    
    set(gcf, 'Color', 'w');
    set(gca, 'Box', 'off', 'FontSize', 14);
    
    % compute correlation between imaging and behavior scores for current LC
    corr_LxLy(iter_lc) = corr(Lx(:, this_lc), Ly(:, this_lc));
    disp(['LC' num2str(this_lc) ': r = ' num2str(corr_LxLy(iter_lc), '%0.2f')]);
    
    % save vector figure (eps), 300 dpi
    print(gcf, fullfile(save_opts.output_path, [save_opts.prefix '_LC' num2str(iter_lc) '_corrLxLy']), '-depsc2', '-vector', '-r300');
    
    close gcf
end

end