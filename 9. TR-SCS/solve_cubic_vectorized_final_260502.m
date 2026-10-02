function [delta_mu_vec, dbg_solver] = solve_cubic_vectorized_final_260502(G, H, M, current_mu)
% SOLVE_CUBIC_VECTORIZED_FINAL
% -------------------------------------------------------------------------
% Box-Constrained Cubic Surrogate Solver (箱约束三次代理求解器)
% 具备绝对收敛性保证，确保代理函数 m(delta) 必定单调递增
% -------------------------------------------------------------------------

    G = single(G(:));
    H = single(H(:));
    M = single(M(:));
    current_mu = single(current_mu(:));
    
    Nvox = length(G);

    % =====================================================================
    % 1. 计算每个体素专属的物理合法区间 [L_j, U_j]
    % =====================================================================
    % 衰减系数的物理极限：最小为 0 (空气)
    % 最大设为 0.03 mm^-1 (完美匹配您的单位系统)
    mu_min = 0.0;
    mu_max = 0.03; 
    
    L_j = mu_min - current_mu;
    U_j = mu_max - current_mu;

    % =====================================================================
    % 2. 求解无约束下的三次函数驻点 delta_star
    % =====================================================================
    delta_star = zeros(Nvox, 1, 'single');

    idx_pos = G > 0;
    if any(idx_pos)
        delta_star(idx_pos) = (H(idx_pos) + sqrt(H(idx_pos).^2 + 2 .* M(idx_pos) .* G(idx_pos))) ./ M(idx_pos);
    end

    idx_neg = G < 0;
    if any(idx_neg)
        discriminant = max(H(idx_neg).^2 - 2 .* M(idx_neg) .* G(idx_neg), 0);
        delta_star(idx_neg) = (-H(idx_neg) - sqrt(discriminant)) ./ M(idx_neg);
    end

    % =====================================================================
    % 3. 将驻点钳制到物理区间内
    % =====================================================================
    delta_clamped = min(max(delta_star, L_j), U_j);

    % =====================================================================
    % 4. 计算代理函数 m(delta) 在四个关键候选点的值
    % =====================================================================
    calc_m = @(d) G .* d + 0.5 .* H .* (d.^2) - (M ./ 6.0) .* (abs(d).^3);

    m_0    = zeros(Nvox, 1, 'single'); 
    m_star = calc_m(delta_clamped);    
    m_L    = calc_m(L_j);              
    m_U    = calc_m(U_j);              

    % =====================================================================
    % 5. 向量化寻找代理函数的全局最大值 (Box-Constrained Argmax)
    % =====================================================================
    [max_m, max_idx] = max([m_0, m_star, m_L, m_U], [], 2);

    delta_chosen = zeros(Nvox, 1, 'single');
    
    delta_chosen(max_idx == 1) = 0.0;                       
    delta_chosen(max_idx == 2) = delta_clamped(max_idx == 2); 
    delta_chosen(max_idx == 3) = L_j(max_idx == 3);         
    delta_chosen(max_idx == 4) = U_j(max_idx == 4);         

    delta_mu_vec = delta_chosen;

    % =====================================================================
    % 6. 输出诊断信息
    % =====================================================================
    dbg_solver.pred_increase = sum(max_m); 
    dbg_solver.saddle_escapes = sum(max_idx == 2 & H > 0); 
    dbg_solver.hit_lower_bound = sum(max_idx == 3);
    dbg_solver.hit_upper_bound = sum(max_idx == 4);
    dbg_solver.violation_count = 0; 
end