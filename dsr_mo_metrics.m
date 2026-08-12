function m = dsr_mo_metrics(found_x, found_f, problem, accuracy)
    % M = dsr_mo_metrics(found_x, found_f, problem, accuracy)
    % 计算 CEC2026 单目标多模态优化的评测指标。
    %
    % found_x : (n x D) 算法报出的候选峰位置
    % found_f : (n x 1) 算法报出的目标值（最小化）
    % problem : ProblemMM 对象（含 globMinX / globMinF / nichRad / numGlobMin）
    % accuracy: 相对适应度容忍度（默认 1e-3）
    %
    % 返回 struct m：
    %   n_found   : 命中真实全局峰的数量
    %   pr        : Peak Ratio = n_found / numGlobMin
    %   precision : Precision  = n_found / n_reported
    %   f1        : F1 = 2*precision*pr/(precision+pr)
    %   coverage  : Peak Coverage = n_found / numGlobMin（与 pr 同义，用于汇总）
    %   sr        : success rate = 1 若全部命中
    %   dists     : 每个真实峰到最近命中点的距离 (numGlobMin x 1)

    if nargin < 4 || isempty(accuracy)
        accuracy = 1e-3;
    end

    true_peaks = problem.globMinX;
    n_min = problem.numGlobMin;
    glob_min_f = problem.globMinF;

    scale = abs(glob_min_f);
    if scale < 1e-12
        scale = 1.0;
    end
    tol = accuracy * scale;

    nich = problem.nichRad;
    peak_radius = min(nich(:));

    m = struct();
    m.dists = inf(n_min, 1);

    if isempty(found_x)
        m.n_found = 0; m.pr = 0; m.precision = 0; m.f1 = 0;
        m.coverage = 0; m.sr = 0;
        return;
    end

    found_f = found_f(:);
    dim = size(true_peaks, 2);
    found_x = reshape(found_x, [], dim);

    % 1) 精度过滤：只保留达到全局最优精度的点
    valid_mask = found_f <= glob_min_f + tol + 1e-9;
    valid_x = found_x(valid_mask, :);
    if isempty(valid_x)
        m.n_found = 0; m.pr = 0; m.precision = 0; m.f1 = 0;
        m.coverage = 0; m.sr = 0;
        return;
    end

    % 2) 去重：点之间距离 > peak_radius 才算不同峰
    valid_list = [];
    for i = 1:size(valid_x, 1)
        v = valid_x(i, :);
        dup = false;
        for j = 1:size(valid_list, 1)
            if norm(v - valid_list(j, :)) <= peak_radius
                dup = true;
                break;
            end
        end
        if ~dup
            valid_list = [valid_list; v]; %#ok<AGROW>
        end
    end

    % 3) 匹配每个真实峰
    n_found = 0;
    dists = inf(n_min, 1);
    for i = 1:n_min
        tp = true_peaks(i, :);
        dmin = inf;
        for j = 1:size(valid_list, 1)
            d = norm(tp - valid_list(j, :));
            if d < dmin, dmin = d; end
        end
        dists(i) = dmin;
        if dmin < peak_radius
            n_found = n_found + 1;
        end
    end

    n_rep = size(found_x, 1);
    pr = n_found / n_min;
    precision = n_found / n_rep;
    if precision + pr > 0
        f1 = 2 * precision * pr / (precision + pr);
    else
        f1 = 0;
    end

    m.n_found   = n_found;
    m.pr        = pr;
    m.precision = precision;
    m.f1        = f1;
    m.coverage  = pr;
    m.sr        = double(n_found >= n_min);
    m.dists     = dists;
end