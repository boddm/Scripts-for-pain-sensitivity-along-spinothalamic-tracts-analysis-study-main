function myPLS_plot_screeplot_copy(S, Sp_vect, pvals_LC, save_opts)
% This function plots singular values (observed & surrogates obtained 
% with permutation testing), as well as the covariance explained by each
% latent component (LC).
%
% Inputs:
% - S           : L x L matrix, L is #LCs, singular values (diagonal matrix)
% - Sp_vect     : L x #permutations matrix, permuted singular values         
% - pvals_LC    : L x 1 vector with p-values for each LC
% - save_opts   : Options for result saving and plotting
%       - .output_path : output directory where figures are saved
%       - .prefix      : prefix of all results files (optional)

%% ========== Variable preparation and filename setup ========== %%
nLC = size(diag(S), 1); % nLC: number of latent components (LCs)
file_name = fullfile(save_opts.output_path, [save_opts.prefix '_LC_pvals_explainedCovariance']); % output filename

%% ========== Compute cumulative percentage of covariance explained by each LC ========== %%
cumul_explCov = cumsum(diag(S).^2 / sum(diag(S).^2)) * 100;

%% ========== Plot dual-y-axis figure using yyaxis ========== %%
figure('Position', save_opts.fig_pos_img); % create new figure and set position

% Left y-axis: singular values
yyaxis left
h1 = plot(1:nLC, diag(S), '-ok', 'MarkerFaceColor', 'w', 'LineWidth', 1); % black hollow circles line
ylabel('Singular values', 'FontSize', 14)
set(gca, 'YColor', 'k');
set(gca, 'Ytick', 0:100:ceil(max([diag(S); mean(Sp_vect(:))+2*std(Sp_vect(:))])));
hold on

% Right y-axis: cumulative percentage of explained covariance
yyaxis right
h2 = plot(1:nLC, cumul_explCov, '-o', 'Color', [0 0.4470 0.7410], 'MarkerFaceColor', [0 0.4470 0.7410], 'LineWidth', 1);
ylabel('Explained covariance (%)', 'FontSize', 14)
set(gca, 'YColor', [0 0.4470 0.7410]);
set(gca, 'YLim', [0 100]);
set(gca, 'Ytick', 0:20:100);

% Set x-axis
xlabel('Latent component', 'FontSize', 14);
set(gca, 'XTick', 1:nLC);
set(gcf, 'Color', 'w'); % set background to white

%% ========== Plot permutation distribution mean and standard deviation ========== %%
yyaxis left
h3 = errorbar(1:nLC, mean(Sp_vect, 2), std(Sp_vect,[], 2), 'r-', 'LineWidth', 1); % red error bars

%% ========== Legend ========== %%
legend([h1 h3], {'Observed', 'Permuted'}, 'FontSize', 14, 'Location', 'best');

%% ========== Display p-value for each LC and set x-axis labels ========== %%
str = cell(1, nLC);
for iLC = 1:nLC
    str{iLC} = sprintf('LC%d (p=%5.3f)', iLC, pvals_LC(iLC)); % generate label for each LC
    fprintf('%s\n', str{iLC}); % print p-values to console
end
set(gca, 'XTick', 1:nLC, 'XTickLabel', str, 'Box', 'off', 'TickDir', 'out'); % set x-axis labels

title('Covariance explained by each LC', 'FontSize', 16); % figure title

%% ========== Save figure ========== %%
print(gcf, [file_name '.tif'], '-dtiff', '-r300'); % save as tif with 300 dpi resolution

close(gcf)

end