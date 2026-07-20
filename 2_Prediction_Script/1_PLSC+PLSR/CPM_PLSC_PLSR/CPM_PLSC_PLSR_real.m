%% 主脚本优化版
% 功能：执行CPM-PLSR分析，包括交叉验证和性能评估
% 输入：无
% 输出：保存结果到mat文件
% 示例：直接运行脚本
% 历史：2024-05-20 优化版本

clear; close all; clc;

%% 设置随机种子（确保结果可重复）
rng(42);  % 固定随机种子，确保每次运行结果一致

%% 环境设置
scriptsPath = '/Users/boddm/Desktop/工作/2.Brain_Spinal_PAIN_吴_副本/Work2/Script/PLSR_Script_副本2/1_PLSC+PLSR/CPM_PLSC_PLSR'; % PLSC程序的路径，需要改
addpath(genpath(scriptsPath));

path = '/Users/boddm/Desktop/杂事/未命名文件夹/4.PLSC+PLSR/PLSC+PLSR_未回归/2_delete_33+26/hot40_hot100_vonfrey+FA_AD_MD_RD';
basefile = fullfile(path, 'Processing'); % 输入路径-需要改
outputFile = fullfile(path, 'Results_PLSR'); % 输出路径，需要改

mkdir(outputFile)

%%
Data = load(fullfile(basefile, 'data.mat')); % ibs的影像矩阵

connMatAll = Data.brain_data';  % n_edges x n_subj
behavDataAll = Data.beh_data;  % n_subj x 8_behaviors

%% 生成交叉验证参数
fprintf('===== 生成交叉验证参数 =====\n');

nSubjects = size(connMatAll, 2);
nPermutations = 1000;
nFolds = size(behavDataAll, 1);

%% 生成permIndices（所有排列共享）
fprintf('生成排列列表...\n');
permIndices = zeros(nPermutations, nSubjects);
for np = 1:nPermutations
    permIndices(np, :) = randperm(nSubjects);
end

%% 生成foldSizes和foldBounds（所有排列共享）
fprintf('生成交叉验证分组...\n');
foldSizes = repmat(floor(nSubjects / nFolds), nFolds, 1);
foldSizes(1:mod(nSubjects, nFolds)) = foldSizes(1:mod(nSubjects, nFolds)) + 1;

foldBounds = zeros(nFolds, 2);
startIdx = 1;
for nf = 1:nFolds
    endIdx = startIdx + foldSizes(nf) - 1;
    foldBounds(nf, :) = [startIdx, endIdx];
    startIdx = endIdx + 1;
end

%% 真实标签的10折交叉验证

% 使用第一次排列的真实标签
realPermIdx = 1;
nSelectedBehav = size(behavDataAll, 2); % 使用所有行为变量进行预测

%% 参数设置
bsrK = 1;           % 特征选择参数 (lv1)
bsrThresh = [1.96 2.58 3.29];      % Bootstrap Ratio阈值 (bsr2)
nComp = 4;          % PLS成分数 (ncomp4)

%% 初始化存储
for nbsr = 1
    for nc = 3%:nComp

        % 清除或重置内部循环前需要重置的变量
        clear perfResults beta mask
        behavTestNorm = zeros(size(behavDataAll, 1), nSelectedBehav);  % 重置归一化行为数据
        predTestScores = zeros(size(behavDataAll, 1), nSelectedBehav);  % 重置预测值

        for nf = 1:nFolds
            fprintf('\t真实标签交叉验证 - Fold %d/%d\n', nf, nFolds);

            %% 划分训练集/测试集
            testSubjIdx = permIndices(realPermIdx, foldBounds(nf, 1):foldBounds(nf, 2));
            trainSubjIdx = setdiff(1:size(behavDataAll, 1), testSubjIdx);

            %% 获取连接矩阵
            connMatTrain = connMatAll(:, trainSubjIdx);
            connMatTest = connMatAll(:, testSubjIdx);

            %% 行为数据 (不进行Z-score)
            behavTrain = behavDataAll(trainSubjIdx, 1:nSelectedBehav);
            behavTest = behavDataAll(testSubjIdx, 1:nSelectedBehav);

            %% 归一化测试集数据 (使用训练集参数)
            behavTestNorm(testSubjIdx, :) = (behavTest - mean(behavTrain)) ./ (std(behavTrain) + 1e-8);

            %% CPM预测
            [mask, predTestScores(testSubjIdx, :), perfResults.beta{nf}] = cpm_lr_D1D2(connMatTrain, connMatTest, behavTrain, [], bsrK, nc, bsrThresh(nbsr));

            %% 存储特征掩码
            perfResults.mask{nf} = mask;
        end

        %% 计算真实性能
        fprintf('===== 计算真实性能 =====\n');

        for t1 = 1:nSelectedBehav
            for t2 = t1 %1:nSelectedBehav
                [perfResults.R(t1, t2), perfResults.P(t1, t2)] = corr(behavTestNorm(:, t1), predTestScores(:, t2));
                deno = behavTestNorm(:, t1) - 0;
                nume = behavTestNorm(:, t1) - predTestScores(:, t2);
                perfResults.q2(t1, t2) = 1 - mean(nume.^2) / mean(deno.^2);
                % fprintf('行为变量 %d->%d:q2 = %.4f, r = %.4f, p = %.4f\n', t1, t2, perfResults.q2(t1, t2), perfResults.R(t1, t2), perfResults.P(t1, t2));
            end
        end

        %% 保存结果
        fprintf('===== 保存结果 =====\n');

        saveName = sprintf('result_real_allBehav_lv%d_bsr%d_ncomp%d_fold%d.mat', bsrK, bsrThresh(nbsr), nc, nFolds);

        save(fullfile(outputFile, saveName), 'perfResults', 'predTestScores', 'behavTestNorm');
        fprintf('结果已保存到: %s\n', fullfile(outputFile, saveName));

    end
end

rmpath(genpath(scriptsPath));
