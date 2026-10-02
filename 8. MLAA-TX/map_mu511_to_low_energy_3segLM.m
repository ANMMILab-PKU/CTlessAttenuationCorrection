function [mu_low, deriv_eta] = map_mu511_to_low_energy_3segLM(mu_511, energy_type)
% MAP_MU511_TO_LOW_ENERGY_V2
% -------------------------------------------------------------------------
% 三段分段线性插值模型 (Three-Segment Piecewise Linear Mapping)
% 基于 NIST XCOM 数据，采用三个锚点:
%   (1) 水 (Water)          — 空气与软组织混合区基准
%   (2) 海绵骨 (Cancellous)  — 软组织与疏质骨过渡区
%   (3) 密质骨 (Cortical)    — 高密度骨骼区
%
% 相较于 V1 (两段模型) 的改进:
%   - 在骨骼区域引入中间锚点，减小光电效应非线性引起的映射偏差
%   - 软组织区域保持不变 (不影响软组织重建精度)
%
% 输入:
%   mu_511:      511 keV 下的衰减图 (列向量, mm^-1)
%   energy_type: 202 或 307
%
% 输出:
%   mu_low:    映射后的等效低能衰减系数 (列向量, mm^-1)
%   deriv_eta: 映射函数在各体素上的局部导数/斜率 (用于链式反投影)
%
% 使用建议:
%   中间锚点 (cancellous bone) 的默认值基于 NIST XCOM 数据估算，
%     1. 在重建图像中选取已知材料 ROI
%     2. 计算 ROI 内 mu_511 均值
%     3. 根据该材料在低能下的真实 mu 值反推各段斜率
% -------------------------------------------------------------------------

mu_511 = single(mu_511);

% =====================================================================
% 锚点定义 (Anchor Points)
% =====================================================================
% --- 锚点 1: 水 (Water, NIST XCOM) ---
mu_water_511 = single(0.00960);   % mm^-1

% --- 锚点 2: 海绵骨 (Cancellous Bone, 约50%骨基质+50%骨髓) ---
mu_cancel_511 = single(0.01370);  % mm^-1

% --- 锚点 3: 密质骨 (Cortical Bone, NIST XCOM ICRU) ---
mu_bone_511 = single(0.01780);    % mm^-1

% =====================================================================
% 各能量下的低能锚点值
% =====================================================================
if energy_type == 202
    mu_water_low  = single(0.01370);  % Water    @ 202 keV
    mu_cancel_low = single(0.02060);  % Cancel.  @ 202 keV (NIST 混合估计)
    mu_bone_low   = single(0.02745);  % Cortical @ 202 keV
elseif energy_type == 307
    mu_water_low  = single(0.01180);  % Water    @ 307 keV
    mu_cancel_low = single(0.01670);  % Cancel.  @ 307 keV (NIST 混合估计)
    mu_bone_low   = single(0.02150);  % Cortical @ 307 keV
else
    error('Unsupported energy type: %d. Must be 202 or 307.', energy_type);
end

% =====================================================================
% 计算各段斜率 (Piecewise Slopes)
% =====================================================================
% 段 1: [0, mu_water_511] — 空气至软组织 (Compton 散射主导)
k_soft = mu_water_low / mu_water_511;

% 段 2: [mu_water_511, mu_cancel_511] — 软组织至海绵骨 (光电效应逐渐增强)
k_cancel = (mu_cancel_low - mu_water_low) / (mu_cancel_511 - mu_water_511);

% 段 3: [mu_cancel_511, mu_bone_511+] — 海绵骨至密质骨及以上 (光电效应主导)
k_dense = (mu_bone_low - mu_cancel_low) / (mu_bone_511 - mu_cancel_511);

% =====================================================================
% 初始化输出
% =====================================================================
mu_low    = zeros(size(mu_511), 'single');
deriv_eta = zeros(size(mu_511), 'single');

% =====================================================================
% 分段掩模 (Segment Masks)
% =====================================================================
mask_soft   = (mu_511 <= mu_water_511);
mask_cancel = (mu_511 > mu_water_511) & (mu_511 <= mu_cancel_511);
mask_dense  = (mu_511 > mu_cancel_511);

% =====================================================================
% 段 1: 空气与软组织混合区 (Air-Soft Tissue Mixture)
% 此区域 Compton 散射主导，能量缩放近似线性
% η(μ_511) = k_soft * μ_511
% =====================================================================
mu_low(mask_soft)    = mu_511(mask_soft) .* k_soft;
deriv_eta(mask_soft) = k_soft;

% =====================================================================
% 段 2: 软组织至海绵骨过渡区 (Soft Tissue - Cancellous Bone Transition)
% 光电效应贡献逐渐增强，斜率介于 k_soft 和 k_dense 之间
% η(μ_511) = μ_water_low + k_cancel * (μ_511 - μ_water_511)
% =====================================================================
mu_low(mask_cancel)    = mu_water_low + k_cancel .* (mu_511(mask_cancel) - mu_water_511);
deriv_eta(mask_cancel) = k_cancel;

% =====================================================================
% 段 3: 海绵骨至密质骨区 (Cancellous - Cortical Bone)
% 光电效应完全主导，能量缩放比例最大
% η(μ_511) = μ_cancel_low + k_dense * (μ_511 - μ_cancel_511)
% =====================================================================
mu_low(mask_dense)    = mu_cancel_low + k_dense .* (mu_511(mask_dense) - mu_cancel_511);
deriv_eta(mask_dense) = k_dense;

end

