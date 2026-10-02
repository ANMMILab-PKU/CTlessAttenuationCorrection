function curv_optimal = calc_optimal_curvature_double_opt_0p001(z, y, y_bar, r)

% CALC_OPTIMAL_CURVATURE
% 引入局部双精度 (Local Double Precision) 解决浮点数抵消与溢出问题
% 兼顾解决:
% 1. z极小时 h(0)-h(z) 的灾难性截断误差 (背景极值)
% 2. z极大时 exp(z) 造成的 NaN 污染

% 全部LOR先赋予有保证的下界值，只对z足够大（数值精度稳定）的LOR才用一般公式覆盖

% =========================================================
% --- 局部双精度计算开始 ---
% =========================================================
z_d     = double(z);
y_d     = double(y);
y_bar_d = double(y_bar);
r_d     = double(r);

% 1. 恢复 phi (Unattenuated Signal)
% 由于使用了 double，极大拓展了数值下溢的边界。
% 为绝对防止人工极值或迭代中发散导致的 z > 700 引起 Inf 崩溃，对 z 加一个物理安全锁
% 现实中人体最大的线积分通常不超过 10~20，设为 50 是极其安全的物理边界。
z_safe_d = min(z_d, 50.0);
signal_d = y_bar_d - r_d;               % = b*exp(-z)
phi_d = signal_d .* exp(z_safe_d);      % = b（无衰减信号）

% ===== 第1步：所有LOR先赋予 z=0 极限值（曲率下界）=====
% 对应 Li Cheng 的 deri2_H0/deri2_h0，Fessler 的 ni_max
% r=0时退化为 phi_d = b_tr
curv_optimal_d = phi_d .* (1 - y_d .* r_d ./ (phi_d + r_d + 1e-300).^2);
curv_optimal_d = max(curv_optimal_d, 0);

% ===== 第2步：z>=0.1 的LOR用一般公式覆盖（Fessler阈值）=====
% 对应 Li Cheng 的 mask=(T>0) 配合 0.001 的数值安全余量
mask_general = (z_d >= 0.001);

if any(mask_general)
    zg   = z_d(mask_general);
    yg   = y_d(mask_general);
    rbg  = y_bar_d(mask_general);
    rg   = r_d(mask_general);
    sg   = signal_d(mask_general);
    phig = phi_d(mask_general);

    ybar0_g = phig + rg;
    h0_g = ybar0_g - yg .* log(ybar0_g + 1e-300);
    hz_g = rbg      - yg .* log(rbg    + 1e-300);
    hdot_g = (1 - yg ./ (rbg + 1e-300)) .* (-sg);

    tmp_g = 2 .* (h0_g - hz_g + zg .* hdot_g);
    c_gen = max(tmp_g, 0) ./ zg.^2;   % max(tmp,0) 对应 Fessler 的 max(tmp(i),0)

    curv_optimal_d(mask_general) = c_gen;
end

% =========================================================
% --- 局部双精度计算结束 ---
% =========================================================

% 7. 转回单精度输出 (极度重要，防止污染外部三维大矩阵的内存)
curv_optimal = single(curv_optimal_d);
% 保证曲率非负 (安全底线)
curv_optimal(curv_optimal < 0) = 0;
end

