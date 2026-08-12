function [xb, fb] = dsr_mo_nelder_mead(fun, x0, lb, ub, step, tol, maxiter)
    % dsr_mo_nelder_mead --- 无导数 Nelder-Mead 单纯形最小化（Octave 原生实现）。
    % 用于对每个候选峰代表做局部精修，能沿弯曲谷道（如 Rosenbrock）收敛。
    %
    % fun      : 目标函数句柄，fun(x) -> 标量（最小化）
    % x0       : 起始点 (1 x D)
    % lb, ub   : 边界（1 x D），越界顶点被夹回
    % step     : 初始单纯形边长
    % tol      : 收敛步长阈值
    % maxiter  : 最大迭代次数
    %
    % 返回 [xb, fb] 最优解及其函数值。

    D = numel(x0);
    X = zeros(D + 1, D);
    X(1, :) = x0(:)';
    for i = 1:D
        X(i + 1, :) = x0;
        X(i + 1, i) = x0(i) + step;
    end

    Y = zeros(D + 1, 1);
    for i = 1:D + 1
        Y(i) = fun(X(i, :));
    end

    alpha = 1.0; gamma = 2.0; rho = 0.5; sigma = 0.5;

    for it = 1:maxiter
        [Y, idx] = sort(Y);
        X = X(idx, :);

        % 收敛判据：单纯形最大边长
        srange = max(max(abs(X - X(1, :))));
        if srange < tol || abs(Y(1) - Y(end)) < 1e-12
            break;
        end

        bary = sum(X(1:D, :), 1) / D;   % 去掉最差点的质心

        % 反射
        xr = bary + alpha * (bary - X(D + 1, :));
        xr = min(max(xr, lb), ub);
        yr = fun(xr);

        if yr < Y(1)
            % 扩张
            xe = bary + gamma * (xr - bary);
            xe = min(max(xe, lb), ub);
            ye = fun(xe);
            if ye < yr
                X(D + 1, :) = xe; Y(end) = ye;
            else
                X(D + 1, :) = xr; Y(end) = yr;
            end
        elseif yr < Y(D)
            X(D + 1, :) = xr; Y(end) = yr;
        else
            % 收缩
            if yr < Y(D + 1)
                xc = bary + rho * (xr - bary);
                xc = min(max(xc, lb), ub);
                yc = fun(xc);
                if yc <= yr
                    X(D + 1, :) = xc; Y(end) = yc;
                else
                    X(2:D + 1, :) = X(1, :) + sigma * (X(2:D + 1, :) - X(1, :));
                    X(2:D + 1, :) = min(max(X(2:D + 1, :), lb), ub);
                    for i = 2:D + 1, Y(i) = fun(X(i, :)); end
                end
            else
                X(2:D + 1, :) = X(1, :) + sigma * (X(2:D + 1, :) - X(1, :));
                X(2:D + 1, :) = min(max(X(2:D + 1, :), lb), ub);
                for i = 2:D + 1, Y(i) = fun(X(i, :)); end
            end
        end
    end

    [fb, idx] = min(Y);
    xb = X(idx, :);
end