function [reg_grad, reg_curv] = calc_regularization_quadratic(img, beta, img_dim, voxel_size)
% CALC_REGULARIZATION_QUADRATIC_LICHENG 
% 根据 Li Cheng 论文实现二次先验正则化 (Quadratic Prior)
%
% Inputs:
%   img:        当前的 mu-map (1D double/single)
%   beta:       正则化权重参数 (scalar)
%   img_dim:    图像尺寸 [nx, ny, nz]
%   voxel_size: 体素尺寸 [vx, vy, vz] (用于计算距离权重)
%
% Outputs:
%   reg_grad:   正则化项的一阶导数 (Gradient of Penalty), dR/dmu
%   reg_curv:   正则化项的 SQS 曲率 (Curvature of Penalty), D_reg
%
% Math:
%   R(mu) = 1/2 * sum_j sum_k w_jk (mu_j - mu_k)^2
%   Grad_j = beta * sum_k w_jk (mu_j - mu_k)
%   Curv_j = beta * 2 * sum_k w_jk  (Factor 2 comes from SQS derivation)

    nx = img_dim(1); ny = img_dim(2); nz = img_dim(3);
    img_3d = reshape(img, [nx, ny, nz]);
    
    reg_grad_3d = zeros(size(img_3d), 'single');
    reg_curv_3d = zeros(size(img_3d), 'single');
    
    % 定义 3D 26-邻域 (相对于中心 [0,0,0] 的偏移)
    % 如果想要严格的 6-邻域，只保留 distance=1 的偏移即可
    [x, y, z] = meshgrid(-1:1, -1:1, -1:1);
    shifts = [x(:) y(:) z(:)];
    shifts(14,:) = []; % 移除中心点 (0,0,0)
    
    % 预计算物理距离权重 w_jk = 1 / distance
    dists = sqrt( (shifts(:,1)*voxel_size(1)).^2 + ...
                  (shifts(:,2)*voxel_size(2)).^2 + ...
                  (shifts(:,3)*voxel_size(3)).^2 );
    weights = 1 ./ dists;
    
    % 遍历所有邻居
    for k = 1:size(shifts, 1)
        s = shifts(k, :);
        w = weights(k);
        
        % 获取邻居像素 (利用 circshift)
        img_neighbor = circshift(img_3d, s);
        
        % 边界处理: 
        % circshift 会循环移位，我们需要将移位后从另一侧“卷回来”的边界像素剔除
        % 创建一个 mask，标记有效的邻域连接
        mask = true(size(img_3d));
        if s(1) == 1, mask(1,:,:) = false; elseif s(1) == -1, mask(end,:,:) = false; end
        if s(2) == 1, mask(:,1,:) = false; elseif s(2) == -1, mask(:,end,:) = false; end
        if s(3) == 1, mask(:,:,1) = false; elseif s(3) == -1, mask(:,:,end) = false; end
        
        % 1. 梯度累加: sum( w_{jk} * (mu_j - mu_k) )
        diff = (img_3d - img_neighbor);
        term_grad = w * diff;
        
        % 仅在有效边界内累加
        reg_grad_3d(mask) = reg_grad_3d(mask) + term_grad(mask);
        
        % 2. 曲率累加: sum( w_{jk} ) * 2
        % 注意: Li Cheng 论文引用的 Fessler SQS 推导中，二次项曲率系数是 2
        term_curv = 2 * w; 
        reg_curv_3d(mask) = reg_curv_3d(mask) + term_curv;
    end
    
    % 应用 Beta 权重
    reg_grad = beta * reg_grad_3d(:);
    reg_curv = beta * reg_curv_3d(:);

end