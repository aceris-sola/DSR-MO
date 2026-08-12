function [Xf, Ff] = dsr_mo_density_filter(X, F, eps_frac, min_pts)
% dsr_mo_density_filter - DBSCAN-like density filtering of a candidate report set.
%   [Xf, Ff] = dsr_mo_density_filter(X, F, eps_frac, min_pts)
%
%   Mimics the "density-filtered reporting" idea of S-CARD-CMSA at the
%   reporting stage: group reported candidates by spatial density and keep a
%   single best (lowest-fitness) representative per dense cluster, thereby
%   removing redundant reports that share the same basin. This is a simplified
%   non-parametric analogue of S-CARD-CMSA's density-filtered secondary
%   insertion, applied as a plug-in post-processing module.
%
%   X        : (n x dim) reported candidate positions
%   F        : (n x 1)   objective values (minimization)
%   eps_frac : cluster radius as a fraction of the decision-space diagonal
%              (default 0.02). Used to build connected components.
%   min_pts  : min points to seed a cluster (default 1). With 1, every point
%              is a seed and clusters are transitive closures under eps.
%   Returns the filtered set (at most one representative per cluster).
%
%   No oracle information is used: the radius is derived only from the domain
%   diagonal, consistent with the black-box spirit of the pipeline.

    if nargin < 3 || isempty(eps_frac), eps_frac = 0.02; end
    if nargin < 4 || isempty(min_pts),  min_pts  = 1;    end

    [n, dim] = size(X);
    if n == 0
        Xf = X; Ff = F; return;
    end

    % Adaptive scale: median nearest-neighbour distance of the report set.
    % This reflects the intrinsic basin separation of the current landscape
    % (same principle as the no-oracle ecological-niche radius used in the
    % search stage), so merging is correct across problems with different
    % peak spacing. eps overrides it if provided non-empty.
    D0 = pdist2(X, X);
    D0(1:n+1:end) = inf;
    nn = min(D0, [], 2);
    adaptive = median(nn);
    if ~isfinite(adaptive) || adaptive < 1e-12
        adaptive = 1.0;
    end
    if nargin < 3 || isempty(eps_frac)
        eps = adaptive;   % adaptive default
    else
        diag = sqrt(sum((max(X, [], 1) - min(X, [], 1)).^2));
        if diag < 1e-12, diag = 1.0; end
        eps = eps_frac * diag;
    end

    % Sparse-ish adjacency via pairwise distance is fine for report sets (small n)
    D = pdist2(X, X);
    adj = D <= eps;

    % Union-Find connected components
    parent = 1:n;
    function r = find_p(i)
        while parent(i) ~= i
            i = parent(i);
        end
        r = i;
    end
    for i = 1:n
        for j = i+1:n
            if adj(i, j)
                ri = find_p(i); rj = find_p(j);
                if ri ~= rj
                    parent(ri) = rj;
                end
            end
        end
    end
    root = zeros(1, n);
    for i = 1:n
        root(i) = find_p(i);
    end

    % Keep the best point of each cluster; drop clusters smaller than min_pts
    [u, ~, ic] = unique(root);
    Xf = []; Ff = [];
    for k = 1:numel(u)
        members = find(ic == k);
        if numel(members) < min_pts
            continue;
        end
        [~, bidx] = min(F(members));
        Xf = [Xf; X(members(bidx), :)]; %#ok<AGROW>
        Ff = [Ff; F(members(bidx))];      %#ok<AGROW>
    end
end