function [input, pls_opts, save_opts] = myPLS_initialize(input, pls_opts, save_opts)

% =============================
% 功能：
%   初始化PLS分析的输入、参数和保存选项。设置默认值、检查输入有效性、自动补全缺失字段，确保后续分析流程顺利。
%   支持多种输入格式、分组、行为/对比/交互分析、归一化、绘图与保存参数等。
%
% =============================
% 输入：
%   input    ：结构体，包含分析输入数据，需包含：
%       .brain_data/X0   ：N×M矩阵，原始影像数据
%       .behav_data/Y0   ：N×B矩阵，原始行为/设计数据
%       .grouping        ：N×1向量，受试者分组信息
%       .group_names     ：1×G元胞，组名称（可选）
%       .behav_names     ：1×B元胞，行为变量名称（可选）
%       .img_names       ：1×M元胞，影像变量名称（可选）
%   pls_opts ：结构体，PLS分析参数，需包含：
%       .behav_type           ：行为分析类型（behavior/contrast/contrastBehav/contrastBehavInteract）
%       .grouped_PLS          ：是否分组PLS（0/1）
%       .grouped_perm         ：排列检验是否组内排列（0/1）
%       .grouped_boot         ：引导抽样是否组内抽样（0/1）
%       .boot_procrustes_mod  ：引导Procrustes旋转模式（1/2）
%       .save_boot_resampling ：是否保存引导样本（0/1）
%       .normalization_img    ：影像数据归一化方式（0~4，详见myPLS_bootstrapping.m）
%       .normalization_behav  ：行为数据归一化方式（0~4，详见myPLS_bootstrapping.m）
%   save_opts ：结构体，保存与绘图参数，需包含：
%       .prefix        ：字符串，结果文件名前缀
%       .output_path   ：字符串，输出目录
%       .alpha         ：显著性水平
%       .img_type      ：影像类型（volume/corrMat/barPlot）
%       .mask_file     ：掩膜文件（volume时必需）
%       .struct_file   ：结构像文件（volume时必需）
%       .BSR_thres     ：BSR阈值
%       .load_thres    ：载荷阈值
%       .grouped_plots ：是否分组绘图（0/1）
%       .plot_boot_samples ：是否绘制引导样本（0/1）
%
% =============================
% 输出：
%   input     ：结构体，补全并标准化后的输入数据
%   pls_opts  ：结构体，补全并标准化后的PLS参数
%   save_opts ：结构体，补全并标准化后的保存与绘图参数
%
% =============================
% 示例：
%   [input, pls_opts, save_opts] = myPLS_initialize(input, pls_opts, save_opts);
%
% =============================
% 历史：
%   2024-06-09  初版，标准化注释与格式，详细参数说明。
%   2024-06-10  增加详细中文注释，细化每个变量、流程、分支、异常处理说明。
%   2024-06-10  为所有if判断添加缘由注释。
%   2024-06-10  全面中文化注释和输出。
% =============================

% 此函数设置默认值并检查myPLS_analysis输入的有效性

if nargin < 3 || isempty(save_opts)
    save_opts = struct();
end

% =============================
% 程序主流程开始
% =============================

disp('... 初始化输入参数 ...')

%% ========== 1. 输入数据检查与补全 ========== %%
% 兼容X0和Y0输入，允许input结构体中字段名为X0/Y0或brain_data/behav_data
% 缘由：有些用户可能用X0/Y0命名输入字段，为保证兼容性需自动补全
if ~isfield(input, 'brain_data') && isfield(input, 'X0')
    % 若未提供brain_data但有X0，则自动补全
    input.brain_data = input.X0;
elseif ~isfield(input, 'brain_data') && ~isfield(input, 'X0')
    % 缘由：若既无brain_data也无X0，说明影像数据缺失，无法继续分析
    error('缺少影像数据输入（brain_data 或 X0 字段）');
end

% 缘由：有些用户可能用Y0命名输入字段，为保证兼容性需自动补全
if ~isfield(input, 'behav_data') && isfield(input, 'Y0')
    % 若未提供behav_data但有Y0，则自动补全
    input.behav_data = input.Y0;
elseif ~isfield(input, 'behav_data') && ~isfield(input, 'Y0')
    % 缘由：若既无behav_data也无Y0，说明行为数据缺失，无法继续分析
    error('缺少行为数据输入（behav_data 或 Y0 字段）');
end

% 缘由：PLS分析要求X和Y的样本数一致，否则数据不匹配
if (size(input.brain_data, 1) ~= size(input.behav_data, 1))
    error('影像数据（X）与行为数据（Y）样本数不一致，请检查输入！');
end

% ========== 获取分组、行为变量、成像变量数量 ========== %
% groupIDs：所有分组的唯一ID（如[1 2]）
groupIDs = unique(input.grouping);    % 分组ID向量
nGroups  = length(groupIDs);           % 分组数
nBehav   = size(input.behav_data, 2);  % 行为变量数
nImg     = size(input.brain_data, 2);  % 成像变量数

%% ========== 组名与行为变量名补全 ========== %%
% 缘由：若group_names数量少于实际分组数，说明信息不全，需重置为标准命名
if isfield(input, 'group_names') && numel(input.group_names) < nGroups
    disp('！！！分组名称数量少于实际分组数，将自动使用标准命名！');
    input.group_names = []; % 清空，后续自动补全
end

% 缘由：若未提供group_names或为空，自动生成标准组名，保证后续流程可用
if ~isfield(input, 'group_names') || isempty(input.group_names)
    input.group_names = cell(nGroups, 1);
    for iG = 1:nGroups
        input.group_names{iG} = ['分组' num2str(groupIDs(iG))];
    end
end

% 缘由：若group_names数量多于分组数，说明多余，需截断
if numel(input.group_names) > nGroups
    disp('！！！分组名称数量多于实际分组数，仅保留前若干个：');
    input.group_names = input.group_names(1:nGroups);
    for iG = 1:nGroups
        disp(['   分组' num2str(iG) '："' input.group_names{iG} '"']);
    end
end

% 缘由：若behav_names数量少于实际行为变量数，说明信息不全，需重置为标准命名
if isfield(input, 'behav_names') && numel(input.behav_names) < nBehav
    disp('！！！行为变量名称数量少于实际变量数，将自动使用标准命名！');
    input.behav_names = [];
end

% 缘由：若未提供behav_names或为空，自动生成标准行为变量名，保证后续流程可用
if ~isfield(input, 'behav_names') || isempty(input.behav_names)
    input.behav_names = cell(nBehav, 1);
    for iB = 1:nBehav
        input.behav_names{iB} = ['行为' num2str(iB)];
    end
end

% 缘由：若behav_names数量多于行为变量数，说明多余，需截断
if numel(input.behav_names) > nBehav
    disp('！！！行为变量名称数量多于实际变量数，仅保留前若干个：');
    input.behav_names = input.behav_names(1:nBehav);
    for iB = 1:nBehav
        disp(['   行为' num2str(iB) '："' input.behav_names{iB} '"']);
    end
end


%% ========== 2. PLS参数补全与检查 ========== %%
% 缘由：若未指定分析类型，默认行为PLS，保证流程可用
if ~isfield(pls_opts, 'behav_type') || isempty(pls_opts.behav_type)
    pls_opts.behav_type = 'behavior';
elseif ~(strcmp(pls_opts.behav_type, 'behavior') || strcmp(pls_opts.behav_type, 'contrast') || strcmp(pls_opts.behav_type, 'contrastBehav') || strcmp(pls_opts.behav_type, 'contrastBehavInteract'))
    % 缘由：若分析类型不在允许范围，报错，防止后续出错
    error('行为分析类型（behav_type）无效，请检查参数！')
end

% 缘由：若未指定PLS分组，依据分析类型自动设置，保证参数合理
if ~isfield(pls_opts, 'grouped_PLS')
    disp('未指定PLS分组选项，将根据分析类型自动设置：')
    if contains(pls_opts.behav_type, 'contrast')
        pls_opts.grouped_PLS = 0;
        disp('   分析类型为对比，PLS不考虑分组')
    else
        pls_opts.grouped_PLS = 1;
        disp('   分析类型为行为，PLS按分组进行')
    end
end

% 缘由：对比PLS时，归一化不能选组内归一化，且不能分组PLS，否则数学含义不符
if contains(pls_opts.behav_type, 'contrast')
    if pls_opts.normalization_behav == 2 || pls_opts.normalization_behav == 4 || pls_opts.normalization_img == 2 || pls_opts.normalization_img == 4
        error('对比PLS分析时，不能选择组内归一化，请将归一化方式设置为全体受试者！')
    end
    if pls_opts.grouped_PLS == 1
        error('对比PLS分析时，不能选择分组PLS，请将grouped_PLS参数设为0！')
    end
end

% 缘由：若未指定排列检验分组，依据分析类型自动设置，保证参数合理
if ~isfield(pls_opts, 'grouped_perm') || isempty(pls_opts.grouped_perm)
    disp('未指定排列检验分组选项，将根据分析类型自动设置：')
    if contains(pls_opts.behav_type, 'contrast')
        pls_opts.grouped_perm = 0;
        disp('   分析类型为对比，排列检验不考虑分组')
    else
        pls_opts.grouped_perm = 1;
        disp('   分析类型为行为，排列检验按分组进行')
    end
end

% 缘由：若未指定引导抽样分组，依据分析类型自动设置，保证参数合理
if ~isfield(pls_opts, 'grouped_boot') || isempty(pls_opts.grouped_boot)
    disp('未指定引导抽样分组选项，将根据分析类型自动设置：')
    if contains(pls_opts.behav_type, 'contrast')
        pls_opts.grouped_boot = 0;
        disp('   分析类型为对比，引导抽样不考虑分组')
    else
        pls_opts.grouped_boot = 1;
        disp('   分析类型为行为，引导抽样按分组进行')
    end
end

% 缘由：分组选项不一致时提醒用户，避免分析逻辑混乱
if pls_opts.grouped_boot ~= pls_opts.grouped_perm
    disp('！！！排列检验与引导抽样的分组选项不一致，请确认分析需求！')
end
if pls_opts.grouped_PLS ~= pls_opts.grouped_perm
    disp('！！！PLS与排列检验的分组选项不一致，请确认分析需求！')
end

% 缘由：若未指定Procrustes变换模式，默认1，保证流程可用
if ~isfield(pls_opts, 'boot_procrustes_mod') || isempty(pls_opts.boot_procrustes_mod)
    disp('未指定引导Procrustes变换模式，采用默认值1：仅对U（行为权重）做Procrustes变换')
    pls_opts.boot_procrustes_mod = 1;
end

% 缘由：Procrustes模式只允许1或2，其他值报错
if pls_opts.boot_procrustes_mod ~= 1 && pls_opts.boot_procrustes_mod ~= 2
    error('引导Procrustes变换模式（boot_procrustes_mod）无效，只能为1或2！');
end

% 缘由：若未指定是否保存引导样本，依据成像变量数自动设置，防止大数据量导致存储压力
if ~isfield(pls_opts, 'save_boot_resampling') || isempty(pls_opts.save_boot_resampling)
    disp('未指定是否保存引导样本，将根据成像变量数自动设置：')
    if nImg > 1000
        pls_opts.save_boot_resampling = 0;
        disp('   成像变量数大于1000，不保存引导样本')
    else
        pls_opts.save_boot_resampling = 1;
        disp('   成像变量数不大于1000，将保存引导样本')
    end
end

% 缘由：若用户强制保存但变量数过多，提醒可能导致大文件
if pls_opts.save_boot_resampling && nImg > 1000
    disp('！！！已选择保存引导样本，但成像变量数大于1000，可能导致文件过大！')
end

disp(' ')

% =============================
% 程序主流程结束
% =============================

end