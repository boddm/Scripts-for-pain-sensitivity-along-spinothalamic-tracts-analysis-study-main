function [mask, predTestScores, beta] = cpm_lr_D1D2(connMatTrain, connMatTest, behavTrain, ~, bsrK, nComp, bsrThresh)
    % 功能：使用CPM方法进行预测
    % 输入：connMatTrain - 训练集连接矩阵，connMatTest - 测试集连接矩阵，behavTrain - 训练集行为数据，bsrK - 特征选择参数，nComp - PLS成分数，bsrThresh - Bootstrap Ratio阈值
    % 输出：mask - 特征选择掩码，predTestScores - 测试集预测得分
    % 示例：[mask, pred] = cpm_lr_D1D2(trainMat, testMat, trainBehav, [], 1, 4, 2)
    % 历史：2023-XX-XX 初始版本

    %% 转置连接矩阵 (特征x样本 -> 样本x特征)
    connMatTrain = connMatTrain';
    connMatTest = connMatTest';

    %% 获取特征掩码和回归系数
    [mask, beta] = PLSC_CPMmask(connMatTrain, behavTrain, bsrK, nComp, bsrThresh);

    %% 获取非零特征索引
    selFeatIdx = mask ~= 0;

    %% 对测试集数据进行归一化
    connTestNorm = (connMatTest - mean(connMatTrain)) ./ (std(connMatTrain) + 1e-8);

    %% 使用选中的特征进行预测
    selTestFeatures = connTestNorm(:, selFeatIdx);
    predTestScores = [ones(size(selTestFeatures, 1), 1), selTestFeatures] * beta;

end