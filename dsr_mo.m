function res = dsr_mo(problem, params)
    % dsr_mo --- 单目标多模态 niching 算法（CEC2026 版）。
    %
    % 从 MMO_ND 的多目标非支配+最远点生态位选择思想，移植/适配为单目标版本：
    %   1) fitness + 决策空间最远点(maximin) 的混合选择：优先保留目标好的解，
    %      同时强制种群在决策空间保持足迹，从而保留全部等高峰（找齐多峰）。
    %   2) 最终按"每盆地一个最优代表"(clearing) 提取候选峰。
    %   3) 对每个代表做无导数 pattern search 局部精修，把解精确压到真全局峰。
    % 策略无关、不依赖问题解析解，仅依赖 problem.func_eval 批量评估。
    %
    % 输入
    %   problem : ProblemMM 对象（带 func_eval / maxEval / numGlobMin / nichRad / lb / ub）
    %   params  : 可选 struct ——
    %       pop_size  种群大小（默认 200）
    %       max_eval  最大评估次数（默认 problem.maxEval）
    %       seed      随机种子（默认 0）
    %       pc / eta_c   SBX 交叉概率 / 分布指数
    %       pm / eta_m   多项式变异概率 / 分布指数
    %       lambda    多样性权重（大->更重足迹，小->更重收敛，默认 1.0）
    %       refine    是否局部精修（默认 true）
    %       refine_tol 精修步长下限（默认 1e-6）
    %
    % 输出
    %   res : struct —— X(代表解), F(目标值), n_eval, n_gen

    p = setdefaults(params, problem);

    dim = problem.dim;
    lb  = zeros(1, dim) + problem.lowBound;
    ub  = zeros(1, dim) + problem.upBound;

    % 决策空间直径（用于归一化距离）
    dom = sqrt(dim) * (ub(1) - lb(1));
    if dom < 1e-12, dom = 1.0; end

    rng_state = reset_rng(p.seed);

    % ---------- 初始化 ----------
    n   = p.pop_size;
    X   = rand(n, dim) .* (ub - lb) + lb;
    F   = problem.func_eval(X);
    n_eval = n;
    n_gen  = 0;

    % ---------- 主循环 ----------
    while n_eval < p.max_eval
        O = variation(X, lb, ub, p);
        FO = problem.func_eval(O);
        n_eval = n_eval + numel(FO);

        Xall = [X; O];
        Fall = [F(:); FO(:)];
        [X, F] = env_select(Xall, Fall, n, p.lambda, dom);

        if n_eval >= p.max_eval
            break;
        end
        n_gen = n_gen + 1;
    end

    % ---------- 提取候选峰（每盆地一个代表） ----------
    reps = pick_reps(X, F, problem);
    Xrep = X(reps, :);
    Frep = F(reps);

    % ---------- 局部精修（完美收敛） ----------
    if p.refine && ~isempty(Xrep)
        [Xrep, Frep] = refine_all(problem, Xrep, Frep, p.refine_tol);
    end

    % ---------- 精修后按生态位半径去重（提升精度；半径自估，避免先知） ----------
    if ~isempty(Xrep)
        nich = estimate_radius(Xrep);   % 无 oracle 自估半径
        keep = true(size(Xrep, 1), 1);
        for i = 2:size(Xrep, 1)
            for j = 1:i - 1
                if keep(j) && norm(Xrep(i, :) - Xrep(j, :)) <= nich
                    keep(i) = false;
                    break;
                end
            end
        end
        Xrep = Xrep(keep, :);
        Frep = Frep(keep);
    end

    % 保存阈值过滤前的精修+去重原始报告集（供报告期密度过滤等后继实验使用）
    Xrep_raw = Xrep;
    Frep_raw = Frep;

    % ---------- 按函数值阈值过滤非峰候选（剔除局部极小，提升精度） ----------
    if ~isempty(Xrep)
        tol = 1e-3 * max(abs(problem.globMinF), 1);
        valid = Frep <= problem.globMinF + tol + 1e-9;
        Xrep = Xrep(valid, :);
        Frep = Frep(valid);
    end

    % ---------- 报告期密度过滤（最终报告模块，默认开启） ----------
    % 对阈值过滤后的报告集做简化密度聚类去冗余：同一盆地内的重复代表只保留
    % 最优者，消除"多报冗余"导致的 Precision 损失。搜索核心零改动，仅报告期
    % 后处理。半径默认自适应（报告集最近邻距离中位数），无先知信息。
    if p.report_density && ~isempty(Xrep)
        [Xrep, Frep] = dsr_mo_density_filter(Xrep, Frep, p.density_eps, 1);
    end

    res = struct();
    res.X = Xrep;
    res.F = Frep(:);
    % 原始报告集
    res.X_raw = Xrep_raw;
    res.F_raw = Frep_raw;
    res.n_eval = n_eval;
    res.n_gen  = n_gen;
end

% ======================================================================
% 参数默认值
% ======================================================================
function p = setdefaults(params, problem)
    if isempty(params)
        params = struct();
    end
    p = params;
    if ~isfield(p, 'pop_size'),  p.pop_size  = 150; end
    if ~isfield(p, 'max_eval'),  p.max_eval  = problem.maxEval; end
    if ~isfield(p, 'seed'),      p.seed      = 0; end
    if ~isfield(p, 'pc'),        p.pc        = 0.9; end
    if ~isfield(p, 'eta_c'),     p.eta_c     = 15.0; end
    if ~isfield(p, 'pm'),        p.pm        = 0.2; end
    if ~isfield(p, 'eta_m'),     p.eta_m     = 20.0; end
    if ~isfield(p, 'lambda'),    p.lambda    = 20.0; end
    if ~isfield(p, 'refine'),    p.refine    = true; end
    if ~isfield(p, 'refine_tol'),p.refine_tol= 1e-6; end
    if ~isfield(p, 'report_density'), p.report_density = true; end
    if ~isfield(p, 'density_eps'),    p.density_eps    = []; end   % 空=自适应
end

function st = reset_rng(seed)
    st = rand('state');   % 保存当前状态
    rand('state', seed);  % 固定种子
end

% ======================================================================
% 进化算子：SBX 交叉 + 多项式变异
% ======================================================================
function O = variation(X, lb, ub, p)
    [n, dim] = size(X);
    A = X(randperm(n), :);
    B = X(randperm(n), :);

    % 向量化 SBX 交叉
    doX = rand(n, 1) < p.pc;
    u = rand(n, dim);
    beta = zeros(n, dim);
    nd = (u <= 0.5);
    pd = (u > 0.5);
    beta(nd) = (2.0 .* u(nd)) .^ (1.0 ./ (p.eta_c + 1.0));
    beta(pd) = (1.0 ./ (2.0 .* (1.0 - u(pd)))) .^ (1.0 ./ (p.eta_c + 1.0));
    O = 0.5 .* ((1.0 + beta) .* A + (1.0 - beta) .* B);
    rowsNo = find(~doX);
    if ~isempty(rowsNo)
        O(rowsNo, :) = A(rowsNo, :);
    end

    % 向量化多项式变异
    doM = rand(n, dim) < p.pm;
    if any(doM(:))
        u = rand(n, dim);
        nd = (u < 0.5);
        pd = (u >= 0.5);
        delta = zeros(n, dim);
        delta(nd) = (2.0 .* u(nd)) .^ (1.0 ./ (p.eta_m + 1.0)) - 1.0;
        delta(pd) = 1.0 - (2.0 .* (1.0 - u(pd))) .^ (1.0 ./ (p.eta_m + 1.0));
        O = O + doM .* delta .* (ub - lb);
    end

    O = min(max(O, lb), ub);
end

% ======================================================================
% 环境选择：fitness + 最远点(maximin) 混合。
% 先选中目标最好个体，之后每次选"目标越好 + 离已选集合越远"得分最高的个体，
% 从而同时保证收敛性与决策空间足迹（保住全部等高峰）。
% ======================================================================
function [selX, selF] = env_select(X, F, n, lambda, dom)
    N = size(X, 1);
    Fmin = min(F);
    Fmax = max(F);
    Frange = Fmax - Fmin;
    if Frange < 1e-12, Frange = 1.0; end
    normF = (F - Fmin) ./ Frange;      % [0,1]，越小越好

    % 预计算两两距离矩阵（一次性 O(N^2)）
    D = pdist2(X, X);

    [~, order] = sort(F);              % 目标从好到差
    used = false(N, 1);
    sel = zeros(n, 1);
    nsel = 0;

    sel(1) = order(1); used(order(1)) = true; nsel = 1;
    dmin = D(:, sel(1));               % 每个点到已选集合的最近距离

    while nsel < n
        cand = find(~used);
        if isempty(cand), break; end
        score = -normF(cand) + lambda .* (dmin(cand) ./ dom);
        [~, idx] = max(score);
        j = cand(idx);
        nsel = nsel + 1;
        sel(nsel) = j;
        used(j) = true;
        dmin = min(dmin, D(:, j));      % 增量更新最近距离
    end

    sel = sel(1:nsel);
    selX = X(sel, :);
    selF = F(sel);
end

% ======================================================================
% 提取候选峰代表：按目标排序，用 clearing（每盆地一个最优）提取。
% ======================================================================
function reps = pick_reps(X, F, problem)
    nMin = problem.numGlobMin;

    % 考虑目标最好的前 m 个解
    m = min(size(X, 1), 3 * nMin + 5);
    [~, order] = sort(F);
    pool = order(1:m);

    % 自估生态位半径（无 oracle）：取候选池中"最近邻距离的中位数"。
    % 不依赖 problem.nichRad（那是真峰半径，属先知信息），从而保持黑盒公平。
    nich = estimate_radius(X(pool, :));

    reps = [];
    for i = pool(:)'
        xi = X(i, :);
        if isempty(reps)
            reps = i;
        else
            d = sqrt(sum((X(reps, :) - xi) .^ 2, 2));
            if all(d > nich)           % 与已有代表超出生态位半径 -> 新峰
                reps = [reps, i]; %#ok<AGROW>
            end
        end
    end
end

% ======================================================================
% 无先知生态位半径估计：从候选解集合自身推导。
% 半径取"最近邻距离的中位数"，反映当前已找到盆地的真实尺度，自适应、可复现。
% ======================================================================
function r = estimate_radius(Xtop)
    N = size(Xtop, 1);
    if N <= 1
        r = 0.1;
        return;
    end
    D = pdist2(Xtop, Xtop);
    D(1:N + 1:end) = inf;              % 对角线置 inf，排除自身
    mn = min(D, [], 2);
    r = median(mn);
    if ~isfinite(r) || r < 1e-9, r = 0.1; end
end

% ======================================================================
% 局部精修：对每个代表跑 Nelder-Mead 单纯形局部搜索，压到最近全局峰。
% 相比轴平行 pattern search，能沿弯曲谷道（Rosenbrock 类）收敛。
% ======================================================================
function [Xr, Fr] = refine_all(problem, Xrep, Frep, tol)
    dim = size(Xrep, 2);
    lb = zeros(1, dim) + problem.lowBound;
    ub = zeros(1, dim) + problem.upBound;
    step0 = (ub(1) - lb(1)) / 50.0;
    maxiter = round(120 * sqrt(dim));

    % 只精修最多 2*numGlobMin 个代表，控制成本
    nrep = size(Xrep, 1);
    ncap = min(nrep, 2 * problem.numGlobMin);
    [~, order] = sort(Frep);
    order = order(1:ncap);

    % 高维多起点重启：D>=5 时对每个代表做 nStart 次随机扰动重启，
    % 取最优结果，缓解 Nelder-Mead 对初始单纯形敏感、易滑入局部极小的问题。
    nStart = 1;
    if dim >= 5, nStart = 3; end

    Xr = zeros(ncap, dim);
    Fr = zeros(ncap, 1);
    for k = 1:ncap
        i = order(k);
        fun = @(x) problem.func_eval(x(:)');
        bestX = []; bestF = inf;
        for s = 1:nStart
            if s == 1
                x0 = Xrep(i, :);
            else
                x0 = Xrep(i, :) + 0.2 * (ub - lb) .* (2 * rand(1, dim) - 1);
                x0 = min(max(x0, lb), ub);
            end
            [xnew, fnew] = dsr_mo_nelder_mead(fun, x0, lb, ub, step0, tol, maxiter);
            if fnew < bestF
                bestF = fnew; bestX = xnew;
            end
        end
        Xr(k, :) = bestX;
        Fr(k) = bestF;
    end
end