function [h1, h2, K] = calc_exact_derivatives_transmission_260502(z, y, b, s)
% [Transmission] 计算精确的一阶、二阶导数及三阶导数上界

z = single(z); y = single(y); b = single(b); s = single(s);
z = min(max(z, 0), 50);
s_safe = max(s, 1e-6);

u = b .* exp(-z);
w = max(u + s_safe, 1e-10);

% 一阶导数
h1 = u .* (1 - y ./ w);

% 二阶导数 (修复代数漏洞：必须保留本底强凹性 -u)
h2 = -u + y .* s_safe .* u ./ (w.^2);

% 三阶导数解析上界
K = b + 0.1 .* y;   % 真实三阶导数上界（0.0962）

end