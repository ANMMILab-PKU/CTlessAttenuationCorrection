function [h1_corr, h2_corr, K_corr] = calc_exact_derivatives_emission_260502(proj_mu, proj_act, y, r, s, NF)
% [Emission Sparse Correction] 计算精确的增量导数修正 (True - Base)

proj_mu = min(max(single(proj_mu), 0), 50);
proj_act = max(single(proj_act), 0);
y = single(y); r = single(r); s = single(s); NF = single(NF);

u = exp(-proj_mu) .* proj_act ./ NF;
u(isnan(u) | isinf(u)) = 0;

w = max(u + r + s, 1e-10);

% 1. 一阶修正: True (u - yu/w) - Base (u) = -yu/w
h1_corr = -u .* (y ./ w);

% 2. 二阶修正: True (-u + yu(r+s)/w^2) - Base (-u) = yu(r+s)/w^2
h2_corr = u .* y .* (r + s) ./ (w.^2);

% % 3. 三阶分离上限修正 (接受建议：Base 已经提供完整上界 u，修正为 0)
% K_corr = zeros(size(u), 'single');
% 3. 三阶分离上限修正 (绝对理论界限)
% 泊松发射似然的三阶导数上限精确对应于无本底的期望信号 u
K_corr = u;

end