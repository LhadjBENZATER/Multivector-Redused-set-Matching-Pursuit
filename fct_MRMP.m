function [X_hat, Supp] = fct_MRMP(A, At, Y, K, bk, a)
% FCT_MRMP  Standard Multivector Reduced-set Matching Pursuit.
%
%   [X_hat, Supp] = FCT_MRMP(A, At, Y, K, bk, a)
%
% Fixed-threshold joint sparse recovery: at each iteration, candidate
% atoms are the top round(bk*K) entries of the normalized cross-channel
% correlation, and only those exceeding admission threshold `a` are
% added to the common support. The support is refit by least squares
% and hard-thresholded to K active bins per iteration.
%
% Inputs:
%   A, At - M x N sensing matrix and its N x M adjoint
%   Y     - M x L measurement matrix (L channels/SUs)
%   K     - target joint sparsity
%   bk    - scalar candidate-pool factor in (0,1]
%   a     - scalar admission threshold in (0,1)
%
% Outputs:
%   X_hat - N x L reconstructed matrix
%   Supp  - recovered common support (sorted bin indices)

    [M, L] = size(Y); %#ok<ASGLU>
    N = size(A, 2);
    r = Y;
    X_hat = zeros(N, L);
    Supp = [];

    for iter = 1:K
        corr = At * r;
        cr = abs(corr);
        cr = cr ./ (max(cr, [], 1) + eps);

        poolSize = max(round(bk * K), 1);
        [~, topLin] = maxk(cr(:), poolSize);
        [rowIdx, ~] = ind2sub(size(cr), topLin);

        admitted = rowIdx(cr(topLin) > a);
        Supp = unique([Supp; admitted]);

        if isempty(Supp)
            continue;
        end

        X_hat(Supp, :) = lsqminnorm(A(:, Supp), Y);

        [~, idx] = maxk(abs(X_hat), K, 1);
        mask = false(size(X_hat));
        mask(sub2ind(size(X_hat), idx(:), repelem(1:L, K)')) = true;
        X_hat(~mask) = 0;

        r = Y - A * X_hat;
        Supp = find(sum(X_hat ~= 0, 2));

        if norm(r, 'fro') < 1e-10 || numel(Supp) >= K
            break;
        end
    end
end
