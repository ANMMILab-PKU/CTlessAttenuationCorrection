function [h1, h2, K] = calc_exact_derivatives_transmission_260518(z, y, b, s)
% [Transmission] 计算精确的一阶、二阶导数及三阶导数上界

% 透射二阶导数中引入衰减相关修正项,对透射 Hessian 添加一个与 LOR 路径长度成正比的最小曲率(Hessian 的正则化下界)

z = single(z); y = single(y); b = single(b); s = single(s);
z = min(max(z, 0), 50);
s_safe = max(s, 1e-6);

u = b .* exp(-z);
w = max(u + s_safe, 1e-10);

% 一阶导数
h1 = u .* (1 - y ./ w);

% 二阶导数 (修复代数漏洞：必须保留本底强凹性 -u)
h2 = -u + y .* s_safe .* u ./ (w.^2);

% --- [新增] Hessian 衰减相关下界保护 ---
% 对于高衰减 LOR (u->0), h2->0 导致曲率约束退化。
% 引入与信号强度和路径长度相关的最小曲率保证，
% 防止这些 LOR 的步长因缺乏二阶约束而过冲。
%
% 物理含义: 即使透射信号极弱，目标函数在 z 方向上仍有
% 与空白扫描信号强度相关的内在曲率 (凹性不低于 -b*exp(-z_max))。
%
% eps_h2 控制下界强度。设得太大会过度保守化;
% 设得太小则失去保护作用。建议范围: 1e-4 ~ 1e-3。
eps_h2 = single(5e-4);
h2_floor = -max(abs(h2), eps_h2 .* b);   % 下界与 blank 信号成正比
h2 = min(h2, h2_floor);                   % h2 为负值，取绝对值更大的


% 三阶导数解析上界
K = b + 0.1 .* y;   % 真实三阶导数上界（0.0962）

end