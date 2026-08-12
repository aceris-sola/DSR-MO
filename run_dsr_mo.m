function results = run_dsr_mo(dims, instances, pids, runs, params)
    % results = run_dsr_mo(dims, instances, pids, runs, params)
    % 在 CEC2026 多模态优化基准上运行单目标 niching 算法并评测。
    %
    % 参数（均可省略，默认跑快速子集）：
    %   dims      : 维度列表，默认 [2, 5]
    %   instances : instance 列表，默认 [1]
    %   pids      : problem ID 列表，默认 1:16
    %   runs      : 每个 problem 独立运行次数，默认 1
    %   params    : 传给 dsr_mo 的参数字典
    %
    % 返回 results : struct 数组，每元素一个 (pid, ins, dim) 的运行汇总。

    if nargin < 1 || isempty(dims),     dims = [2, 5]; end
    if nargin < 2 || isempty(instances),instances = [1]; end
    if nargin < 3 || isempty(pids),     pids = 1:16; end
    if nargin < 4 || isempty(runs),     runs = 1; end
    if nargin < 5, params = struct(); end

    addpath(fullfile(fileparts(mfilename('fullpath')), 'cec2026'));

    results = struct([]);
    k = 0;
    for d = dims
        for ins = instances
            for pid = pids
                k = k + 1;
                problem = ProblemMM(pid, ins, d);

                agg = struct('n_found', 0, 'pr', 0, 'precision', 0, ...
                             'f1', 0, 'coverage', 0, 'sr', 0);
                for r = 1:runs
                    p = params;
                    p.seed = (k - 1) * 1000 + r;
                    % 官方 CEC2026 口径：MaxFEs = 20000 * dim（与 S-CARD-CMSA 论文一致，
                    % 原默认 10000*dim 与官方不一致，不可对标）。用户显式传入 max_eval 时优先用用户值。
                    if ~isfield(p, 'max_eval')
                        p.max_eval = 20000 * d;
                    end
                    res = dsr_mo(problem, p);
                    m = dsr_mo_metrics(res.X, res.F, problem);
                    agg.n_found   = agg.n_found   + m.n_found;
                    agg.pr        = agg.pr        + m.pr;
                    agg.precision = agg.precision + m.precision;
                    agg.f1        = agg.f1        + m.f1;
                    agg.coverage  = agg.coverage  + m.coverage;
                    agg.sr        = agg.sr        + m.sr;
                end
                agg.n_found   = agg.n_found / runs;
                agg.pr        = agg.pr / runs;
                agg.precision = agg.precision / runs;
                agg.f1        = agg.f1 / runs;
                agg.coverage  = agg.coverage / runs;
                agg.sr        = agg.sr / runs;

                results(k).pid    = pid;
                results(k).insNo  = ins;
                results(k).dim    = d;
                results(k).n_peaks = problem.numGlobMin;
                results(k).n_found = agg.n_found;
                results(k).pr     = agg.pr;
                results(k).precision = agg.precision;
                results(k).f1     = agg.f1;
                results(k).coverage = agg.coverage;
                results(k).sr     = agg.sr;

                fprintf('P%02d I%02d D%2d peaks=%2d found=%6.1f PR=%5.3f Prec=%5.3f F1=%5.3f SR=%3.1f\n', ...
                    pid, ins, d, problem.numGlobMin, agg.n_found, ...
                    agg.pr, agg.precision, agg.f1, agg.sr);
            end
        end
    end

    % ---- 汇总表 ----
    fprintf('\n==== CEC2026 单目标 niching 汇总 ====\n');
    fprintf('%-6s %-6s %-6s %8s %8s %8s %8s\n', ...
        'Dim', 'nProb', 'nPeak', 'Coverage', 'Precision', 'F1', 'SR');
    for d = dims
        sel = [results.dim] == d;
        if sum(sel) == 0, continue; end
        cov = mean([results(sel).coverage]);
        prec = mean([results(sel).precision]);
        f1  = mean([results(sel).f1]);
        sr  = mean([results(sel).sr]);
        fprintf('%-6d %-6d %-6d %8.4f %8.4f %8.4f %8.4f\n', ...
            d, sum(sel), mean([results(sel).n_peaks]), cov, prec, f1, sr);
    end
end