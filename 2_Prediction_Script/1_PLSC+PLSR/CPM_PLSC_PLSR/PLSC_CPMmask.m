function [mask, beta] = PLSC_CPMmask(inputConnMat, inputBehav, bsr_k, n_comp, bsr_thresh)
% 功能：使用PLSC方法进行特征选择和回归分析
% 输入：inputConnMat - 连接矩阵，inputBehav - 行为数据，bsrK - 特征选择参数，nComp - PLS成分数，bsrThresh - Bootstrap Ratio阈值
% 输出：mask - 特征选择掩码，beta - 回归系数
% 示例：[mask, beta] = PLSC_CPMmask(connMat, behavData, 1, 4, 2)
% 历史：2023-XX-XX 初始版本

%% 参数设置
n_bootstraps = 5000;
norm_type = 1;         % 归一化类型
subj_groups = ones(size(inputConnMat, 1), 1);  % 分组标签

%% 1. 设置PLSC工具箱路径
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'PLSC_Code')));

%% 2. 定义所有输入参数
% 输入来自CPM训练折中的函数参数。
fprintf('定义PLSC输入参数...\n');
input = struct();
input.brain_data = inputConnMat;
input.behav_data = inputBehav;
input.grouping = subj_groups;

% --- PLS选项 ---
pls_opts = struct();
pls_opts.nBootstraps = n_bootstraps;
pls_opts.normalization_img = norm_type;
pls_opts.normalization_behav = norm_type;
pls_opts.grouped_PLS = 0;
pls_opts.grouped_perm = 0;
pls_opts.grouped_boot = 0;
pls_opts.boot_procrustes_mod = 2;
pls_opts.save_boot_resampling = 0;
pls_opts.behav_type = 'behavior';

%% 3. 检查所有输入的有效性
% !!! 在运行PLS前，务必运行此函数以检查设置 !!!
[input, pls_opts] = myPLS_initialize(input, pls_opts);


%% 4. 运行PLSC分析（包括新版Bootstrap流程）
fprintf('运行PLSC分析...\n');
res = myPLS_analysis(input, pls_opts);


%% 5. 提取归一化数据和Bootstrap ratio
fprintf('提取PLSC结果...\n');
x_norm = res.X;
y_norm = res.Y;

rightBSR = res.boot_results.V_stad;

% 当标准误接近0时，Bootstrap ratio可能出现Inf/NaN；沿用旧程序的处理方式。
inf_right_idx = find(~isfinite(rightBSR));
for iter_inf = 1:numel(inf_right_idx)
    rightBSR(inf_right_idx(iter_inf)) = res.V(inf_right_idx(iter_inf));
end
clear inf_right_idx iter_inf

if bsr_k > size(rightBSR, 2)
    error('PLSC_CPMmask:InvalidBSRK', 'bsr_k=%d超过当前PLSC成分数%d。', bsr_k, size(rightBSR, 2));
end

%% 6. 特征选择
fprintf('执行特征选择...\n');
mask = zeros(size(rightBSR, 1), 1);
posFeatIdx = rightBSR(:, bsr_k) > bsr_thresh;
negFeatIdx = rightBSR(:, bsr_k) < -bsr_thresh;
mask(posFeatIdx) = 1;
mask(negFeatIdx) = -1;


%% 7. PLS回归预测
fprintf('执行PLS回归...\n');
selFeatIdx = mask ~= 0;

sel_Features = x_norm(:, selFeatIdx);
[~, ~, ~, ~, beta, ~, ~, ~] = plsregress(sel_Features, y_norm, n_comp); % PLSR

end
