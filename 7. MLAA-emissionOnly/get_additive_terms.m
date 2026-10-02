function [r_em, s_202, s_307] = get_additive_terms(id1, id2, num_events)
% GET_ADDITIVE_TERMS 
% 获取当前 LOR 的加性背景项 (Scatter + Randoms)
%
% Inputs:
%   id1, id2:   Crystal IDs for the current chunk
%   num_events: Number of events in chunk
%
% Outputs:
%   r_em:  Emission 的 Randoms + Scatter
%   s_202: Transmission 202keV 的 Scatter
%   s_307: Transmission 307keV 的 Scatter
%
% 注意: 如果您的数据已经是 "Only True"，则此函数应返回 0，或者您可以利用此接口
% 引入预计算的 Scatter 分量以修正模型。

    % 示例: 假设没有加性项 (Only True 模式)
    % 您可以在这里读取外部文件，或者根据 id 计算
    
    r_em  = zeros(num_events, 1, 'single');
    s_202 = zeros(num_events, 1, 'single');
    s_307 = zeros(num_events, 1, 'single');
    
    % TODO: 
    % scatter_data = load_scatter_from_disk(id1, id2);
    % r_em = scatter_data;

end