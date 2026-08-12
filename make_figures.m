function make_figures()
% make_figures --- Generate all DSR-MO publication figures.
%   Priority 1:
%     fig1_bar_compare.png    per-problem Coverage & Precision grouped bars
%     fig2_radar.png          DSR-MO vs S-CARD-CMSA radar (published metrics)
%     fig3_scatter_peaks.png  decision-space peak distribution P02/P14/P05
%   Priority 2:
%     fig4_graft.png          density-filter on/off (P02/P05/P08)
%     fig5_convergence.png    eval budget vs peaks found (P02)
%     fig6_radius_box.png     self-estimated niche radius distribution
%   Output is written to ./figs

    addpath(pwd);
    addpath(fullfile(pwd, 'cec2026'));
    outdir = fullfile(pwd, 'figs');
    if ~exist(outdir, 'dir'), mkdir(outdir); end

    % ------------------------------------------------------------------
    % 1. Load aggregated D=2 results (16 problems x 15 instances x 5 runs)
    % ------------------------------------------------------------------
    S = load(fullfile('e:\算法\算法\MMO_ND_Octave', 'results_final_d2.mat'));
    R = S.res2;
    pids = 1:16;
    covP = zeros(16, 1); precP = zeros(16, 1); f1P = zeros(16, 1); srP = zeros(16, 1);
    for pid = pids
        sel = [R.pid] == pid;
        covP(pid)  = mean([R(sel).coverage]);
        precP(pid) = mean([R(sel).precision]);
        f1P(pid)   = mean([R(sel).f1]);
        srP(pid)   = mean([R(sel).sr]);
    end

    % ---- Figure 1: per-problem Coverage & Precision grouped bar ----
    fig = figure('Visible', 'off', 'Position', [100 100 980 520]);
    x = 1:16;
    b = bar(x, [covP, precP]);
    b(1).FaceColor = [0.20 0.45 0.85];
    b(2).FaceColor = [0.90 0.45 0.20];
    set(gca, 'XTick', x, 'XTickLabel', arrayfun(@(p) sprintf('P%02d', p), x, 'UniformOutput', false));
    xlabel('Problem ID'); ylabel('Score');
    title('DSR-MO: Peak Coverage & Precision per Problem (D=2, 15 instances x 5 runs)');
    legend({'Coverage', 'Precision'}, 'Location', 'best');
    grid on; ylim([0 1.05]);
    exportgraphics(fig, fullfile(outdir, 'fig1_bar_compare.png'), 'Resolution', 200);
    close(fig);

    % ---- Figure 2: radar DSR-MO vs S-CARD-CMSA (published metrics) ----
    % S-CARD-CMSA publishes Coverage(0.805)/Precision(0.915)/F1(0.831)/official(0.818).
    % SR and runtime are NOT published for S-CARD-CMSA, so the radar uses the
    % four metrics that both algorithms report under the same official rule.
    labels = {'Coverage / RPR', 'Precision', 'F1', 'Official score'};
    ours   = [0.899, 0.992, 0.935, 0.917];
    theirs = [0.805, 0.915, 0.831, 0.818];
    nA = numel(labels);
    ang = linspace(0, 2*pi, nA+1); ang(end) = [];
    fig = figure('Visible', 'off', 'Position', [100 100 620 620]);
    ax = polaraxes('Color', [0.97 0.97 0.97]); hold(ax, 'on');
    polarplot(ax, [ang, ang(1)], [ours, ours(1)], 'o-', 'LineWidth', 2.2, 'Color', [0.20 0.45 0.85]);
    polarplot(ax, [ang, ang(1)], [theirs, theirs(1)], 's-', 'LineWidth', 2.2, 'Color', [0.90 0.45 0.20]);
    ax.ThetaTick = round(ang*180/pi);
    ax.ThetaTickLabel = labels;
    ax.RLim = [0 1];
    title('DSR-MO vs S-CARD-CMSA (D=2, official 20000xDim budget)');
    legend({'DSR-MO', 'S-CARD-CMSA'}, 'Location', 'best');
    exportgraphics(fig, fullfile(outdir, 'fig2_radar.png'), 'Resolution', 200);
    close(fig);

    % ---- Figure 3: decision-space peak scatter P02 / P14 / P05 ----
    targets = struct('pid', {2, 14, 5}, 'title', {'P02 (strong)', 'P14 (weak)', 'P05 (perfect)'});
    fig = figure('Visible', 'off', 'Position', [100 100 1500 480]);
    for t = 1:numel(targets)
        pid = targets(t).pid;
        problem = ProblemMM(pid, 1, 2);
        p = struct('seed', 2026, 'pop_size', 300, 'lambda', 20);
        p.max_eval = 20000 * 2;
        res = dsr_mo(problem, p);
        found = res.X; true_pk = problem.globMinX;
        subplot(1, 3, t); hold on;
        plot(true_pk(:,1), true_pk(:,2), 'o', 'MarkerSize', 7, ...
             'MarkerEdgeColor', [0.35 0.35 0.35], 'LineWidth', 1.4);
        if ~isempty(found)
            plot(found(:,1), found(:,2), 'x', 'MarkerSize', 9, ...
                 'MarkerEdgeColor', [0.85 0.15 0.15], 'LineWidth', 1.8);
        end
        axis equal; xlim([problem.lowBound problem.upBound]); ylim([problem.lowBound problem.upBound]);
        title(targets(t).title);
        xlabel('x_1'); ylabel('x_2'); grid on;
        if t == 1
            legend({'True peaks', 'Found peaks'}, 'Location', 'best');
        end
    end
    sgtitle('Decision-space peak distribution (D=2): true peaks (o) vs DSR-MO found (x)');
    exportgraphics(fig, fullfile(outdir, 'fig3_scatter_peaks.png'), 'Resolution', 200);
    close(fig);

    % ---- Figure 4: density-filter on/off on P02 / P05 / P08 ----
    gd = [2 5 8];
    nRun = 5;
    covOff = zeros(1,3); precOff = zeros(1,3);
    covOn  = zeros(1,3); precOn  = zeros(1,3);
    for g = 1:3
        problem = ProblemMM(gd(g), 1, 2);
        co = 0; po = 0; cn = 0; pn = 0;
        for r = 1:nRun
            pbase = struct('seed', 1000+r, 'pop_size', 300, 'lambda', 20, 'report_density', false);
            pbase.max_eval = 20000 * 2;
            ro = dsr_mo(problem, pbase);
            mo = dsr_mo_metrics(ro.X, ro.F, problem);
            co = co + mo.coverage; po = po + mo.precision;

            pon = pbase; pon.report_density = true;
            rn = dsr_mo(problem, pon);
            mn = dsr_mo_metrics(rn.X, rn.F, problem);
            cn = cn + mn.coverage; pn = pn + mn.precision;
        end
        covOff(g) = co/nRun; precOff(g) = po/nRun;
        covOn(g)  = cn/nRun; precOn(g)  = pn/nRun;
    end
    fig = figure('Visible', 'off', 'Position', [100 100 900 480]);
    Xc = categorical({'P02', 'P05', 'P08'}); Xc = reordercats(Xc, {'P02', 'P05', 'P08'});
    data = [covOff(:) covOn(:) precOff(:) precOn(:)];
    b = bar(Xc, data);
    b(1).FaceColor = [0.55 0.55 0.55]; b(2).FaceColor = [0.20 0.45 0.85];
    b(3).FaceColor = [0.85 0.60 0.60]; b(4).FaceColor = [0.90 0.45 0.20];
    ylim([0 1.05]); grid on; ylabel('Score');
    title('Density-filtered reporting: before vs after (D=2, 5 runs avg)');
    legend({'Coverage (no filter)', 'Coverage (filtered)', ...
            'Precision (no filter)', 'Precision (filtered)'}, 'Location', 'northwest');
    exportgraphics(fig, fullfile(outdir, 'fig4_graft.png'), 'Resolution', 200);
    close(fig);

    % ---- Figure 5: convergence curve, P02, increasing budget ----
    problem = ProblemMM(2, 1, 2);
    fullBudget = 20000 * 2;
    frac = [0.05 0.075 0.1 0.15 0.2 0.3 0.4 0.6 0.8 1.0];
    nFound = zeros(size(frac));
    for k = 1:numel(frac)
        p = struct('seed', 2026, 'pop_size', 300, 'lambda', 20);
        p.max_eval = round(fullBudget * frac(k));
        res = dsr_mo(problem, p);
        mm = dsr_mo_metrics(res.X, res.F, problem);
        nFound(k) = mm.n_found;
    end
    fig = figure('Visible', 'off', 'Position', [100 100 760 520]);
    plot(frac*100, nFound, 'o-', 'LineWidth', 2.2, 'Color', [0.20 0.45 0.85], 'MarkerSize', 7);
    xlabel('Evaluation budget (% of 20,000xDim)'); ylabel('Global peaks found');
    ylim([0 problem.numGlobMin + 1]); grid on;
    title(sprintf('Convergence: peaks found vs budget (P02, D=2, %d true peaks)', problem.numGlobMin));
    exportgraphics(fig, fullfile(outdir, 'fig5_convergence.png'), 'Resolution', 200);
    close(fig);

    % ---- Figure 6: self-estimated niche radius distribution ----
    radii = [];
    for pid = 1:16
        for ins = 1:5
            problem = ProblemMM(pid, ins, 2);
            p = struct('seed', pid*100+ins, 'pop_size', 300, 'lambda', 20);
            p.max_eval = 20000 * 2;
            res = dsr_mo(problem, p);
            if ~isempty(res.X)
                r = estimate_radius_pub(res.X);
                radii = [radii; r]; %#ok<AGROW>
            end
        end
    end
    fig = figure('Visible', 'off', 'Position', [100 100 760 520]);
    sv = sort(radii); nR = numel(sv);
    q = [sv(max(1, ceil(0.25*nR))), sv(max(1, ceil(0.50*nR))), sv(min(nR, ceil(0.75*nR)))];
    iqr_ = q(3) - q(1);
    lo = max(min(radii), q(1) - 1.5*iqr_);
    hi = min(max(radii), q(3) + 1.5*iqr_);
    hold on;
    boxc = [0.20 0.45 0.85];
    rectangle('Position', [0.55 q(1) 0.9 (q(3)-q(1))], 'FaceColor', boxc, 'EdgeColor', 'k');
    plot([0.55 1.45], [q(2) q(2)], '-', 'Color', 'k', 'LineWidth', 2);
    plot([1 1], [lo q(1)], '-', 'Color', 'k');
    plot([1 1], [q(3) hi], '-', 'Color', 'k');
    plot([0.775 1.225], [lo lo], '-', 'Color', 'k');
    plot([0.775 1.225], [hi hi], '-', 'Color', 'k');
    xlim([0 2]); ylim([0 max(radii)*1.1]);
    set(gca, 'XTick', [1], 'XTickLabel', {'Self-estimated niche radius'});
    ylabel('Radius'); grid on;
    title(sprintf('Self-estimated niche radius across problems (D=2, n=%d)', numel(radii)));
    text(1.6, q(2), sprintf('median=%.3g', q(2)), 'FontSize', 10);
    exportgraphics(fig, fullfile(outdir, 'fig6_radius_box.png'), 'Resolution', 200);
    close(fig);

    fprintf('All figures written to %s\n', outdir);
end

% Oracle-free radius: median nearest-neighbour distance (mirrors dsr_mo.estimate_radius)
function r = estimate_radius_pub(Xtop)
    N = size(Xtop, 1);
    if N <= 1, r = 0.1; return; end
    D = pdist2(Xtop, Xtop);
    D(1:N+1:end) = inf;
    mn = min(D, [], 2);
    r = median(mn);
    if ~isfinite(r) || r < 1e-9, r = 0.1; end
end