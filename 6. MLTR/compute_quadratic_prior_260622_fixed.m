function [reg_grad, reg_curv] = compute_quadratic_prior_260622_fixed(image, beta)
% COMPUTE_QUADRATIC_PRIOR_260622
%
% Exact 6-neighborhood pairwise quadratic prior:
%
%   U(mu) = 1/2 * sum_{j~k} (mu_j - mu_k)^2
%
% Each undirected neighbor pair is counted once.
%
% Outputs:
%   reg_grad = beta * dU/dmu
%   reg_curv = beta * diagonal(d2U/dmu2)
%
% Boundary voxels use their true neighbor count:
%   interior = 6, face = 5, edge = 4, corner = 3.
%
% This version avoids imfilter/Image Processing Toolbox dependency and
% returns both outputs as single precision.

image = single(image);
beta = single(beta);

reg_grad = zeros(size(image), 'single');
reg_curv = zeros(size(image), 'single');

% Dimension 1 neighbor pairs
if size(image, 1) > 1
    d = image(1:end-1, :, :) - image(2:end, :, :);
    reg_grad(1:end-1, :, :) = reg_grad(1:end-1, :, :) + d;
    reg_grad(2:end,   :, :) = reg_grad(2:end,   :, :) - d;

    reg_curv(1:end-1, :, :) = reg_curv(1:end-1, :, :) + 1;
    reg_curv(2:end,   :, :) = reg_curv(2:end,   :, :) + 1;
end

% Dimension 2 neighbor pairs
if size(image, 2) > 1
    d = image(:, 1:end-1, :) - image(:, 2:end, :);
    reg_grad(:, 1:end-1, :) = reg_grad(:, 1:end-1, :) + d;
    reg_grad(:, 2:end,   :) = reg_grad(:, 2:end,   :) - d;

    reg_curv(:, 1:end-1, :) = reg_curv(:, 1:end-1, :) + 1;
    reg_curv(:, 2:end,   :) = reg_curv(:, 2:end,   :) + 1;
end

% Dimension 3 neighbor pairs
if size(image, 3) > 1
    d = image(:, :, 1:end-1) - image(:, :, 2:end);
    reg_grad(:, :, 1:end-1) = reg_grad(:, :, 1:end-1) + d;
    reg_grad(:, :, 2:end)   = reg_grad(:, :, 2:end)   - d;

    reg_curv(:, :, 1:end-1) = reg_curv(:, :, 1:end-1) + 1;
    reg_curv(:, :, 2:end)   = reg_curv(:, :, 2:end)   + 1;
end

reg_grad = beta .* reg_grad;
reg_curv = beta .* reg_curv;

end
