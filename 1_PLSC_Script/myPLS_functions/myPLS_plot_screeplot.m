function myPLS_plot_screeplot(S, Sp_vect, pvals_LC, save_opts)
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

%% ========== Variable preparation and file name setting ========== %%
nLC = size(diag(S), 1); % nLC: number of latent components (LCs)
file_name = fullfile(save_opts.output_path, [save_opts.prefix '_LC_pvals_explainedCovariance']); % Output file name

%% ========== Compute cumulative percentage of covariance explained by each LC ========== %%
% cumul_explCov: cumulative percentage of covariance explained by each LC and all previous LCs
cumul_explCov = cumsum(diag(S).^2 / sum(diag(S).^2)) * 100;

%% ========== Plot singular values and explained covariance with dual y-axes ========== %%
figure('Position', save_opts.fig_pos_img, 'Visible', 'off'); % Create new figure window and set position
[AX, H1, H2] = plotyy(1:nLC, diag(S), 1:nLC, cumul_explCov, 'plot', 'plot'); % Left y-axis: singular values, right y-axis: cumulative explained covariance

%% ========== Set plot properties ========== %%
AX(2).YLim = [0 100]; % Right y-axis range set to 0~100%
set(H1, 'Marker', 'o'); % Add dots to singular value curve
set(H1, 'Color', 'k');  % Singular value curve color set to black
set(H2, 'Color', [0 0.4470 0.7410]); % Explained covariance curve color set to blue
set(AX, {'Ycolor'}, {'k'; [0 0.4470 0.7410]}) % Set left y-axis color to black, right y-axis color to blue
set(get(AX(1), 'XLabel'), 'String', 'Latent components', 'FontSize', 16); % x-axis label
set(get(AX(1), 'Ylabel'), 'String', 'Singular values', 'FontSize', 16)    % Left y-axis label
set(get(AX(2), 'Ylabel'), 'String', 'Explained covariance', 'FontSize', 16) % Right y-axis label
set(AX(1), 'Ytick', 0:5000:AX(1).YLim(2)); % Left y-axis tick positions
set(AX(2), 'Ytick', 0:20:100);             % Right y-axis tick positions
set(gcf, 'Color', 'w'); % Set background as white
hold on

%% ========== Add permuted singular value mean ± std error ========== %%
errorbar(1:nLC, mean(Sp_vect, 2), std(Sp_vect,[], 2), 'r-'); 
legend({'observed', 'surrogates'}, 'FontSize', 14); % legend

%% ========== Display p-value for each LC and set x-axis labels ========== %%
for iLC = 1:nLC
    str{iLC} = sprintf('LC%d (p=%5.3f)', iLC, pvals_LC(iLC)); % generate label for each LC
    fprintf('%s\n', str{iLC}); % display p-value in console
end
set(gca, 'XTick', 1:nLC, 'XTickLabel', str, 'Box', 'off', 'TickDir', 'out'); % set x-axis labels

title('Explained covariance by each LC', 'FontSize', 18); % figure title

%% ========== Save figure ========== %%
print(gcf, [file_name '.tif'], '-dtiff', '-r300'); % save figure as tif, 300 dpi

end