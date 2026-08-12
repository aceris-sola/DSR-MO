classdef BasicFun
    % BasicFun --- CEC2026 的 8 个基础函数（最小化）。
    % 支持批量评估：X 为 N x D 矩阵（每行一个解），返回 N x 1 列向量。
    % 向量化后大幅降低 Octave 的解释开销。
    methods (Static)
        function f = evaluate_batch(X, hGO, funID)
            switch funID
                case 1, f = BasicFun.belliptic(X, hGO);
                case 2, f = BasicFun.bdiffPow(X, hGO);
                case 3, f = BasicFun.bschwefelN02Skew(X, hGO);
                case 4, f = BasicFun.brosenbrock(X, hGO);
                case 5, f = BasicFun.backleySkew(X, hGO);
                case 6, f = BasicFun.brastrigin(X, hGO);
                case 7, f = BasicFun.bweierstrass(X, hGO);
                case 8, f = BasicFun.bschwefelN26(X, hGO);
                otherwise, error('This function is not defined');
            end
        end

        % 单行兼容接口（内部转矩阵）
        function f = evaluate(x, hGO, funID)
            f = BasicFun.evaluate_batch(x(:)', hGO, funID);
        end

        function f = backleySkew(X, hGO)
            y = 5.0 .^ (sign(X) * hGO) .* X;
            term1 = sqrt(mean(y .^ 2, 2));
            term2 = mean(cos(2 * pi * y), 2);
            f = -20 * exp(-0.2 * term1) - exp(term2) + 20 + exp(1);
        end

        function f = bdiffPow(X, hGO)
            D = size(X, 2);
            H = hGO * 4.0;
            p = 2.0 + H * (0:(D - 1)) / (D - 1);   % 1 x D
            f = sqrt(sum(abs(X) .^ p, 2));
        end

        function f = belliptic(X, hGO)
            D = size(X, 2);
            pow0 = (0:(D - 1)) / (D - 1) * hGO * 3.0;
            f = 10000.0 ^ (0.5 - hGO) * sum((10.0 .^ pow0 .* X) .^ 2, 2);
        end

        function f = brosenbrock(X, hGO)
            fsphere = 20 * sum(X .^ 2, 2);
            y = X + 1;
            term1 = 100 * sum((y(:, 2:end) - y(:, 1:end - 1) .^ 2) .^ 2, 2);
            term2 = sum((y(:, 1:end - 1) - 1) .^ 2, 2);
            f0 = term1 + term2;
            f = hGO * f0 + (1 - hGO) * fsphere;
        end

        function f = brastrigin(X, hGO)
            base = 5;
            A = 10.0 * (base ^ hGO - 1) / (base - 1);
            f = sum(X .^ 2 + A * (1 - cos(2 * pi * X)), 2);
        end

        function f = bschwefelN02Skew(X, hGO)
            y = 5.0 .^ (sign(X) * hGO) .* X;
            % f = sum_k (sum_{j<=k} y_j)^2 ，对每行(每个解)计算。用累积和。
            cy = cumsum(y, 2);          % N x D
            f = sum(cy .^ 2, 2);
        end

        function f = bschwefelN26(X, hGO)
            D = size(X, 2);
            base = 5.0;
            H = (base ^ hGO - 1) / (base - 1);
            xstar = 420.96874635998202731184436501869;
            fshift = 418.9828872724337062747864351956;
            y = X + xstar;
            p1 = sum((-300 - y) .* (y < -500), 2);
            p2 = sum((y > 500) .* (y - 420), 2);
            P = p1 + p2;
            g = 1.5 * P + sum(-y .* sin(abs(y) .^ 0.5), 2) + fshift * D;
            f = H * g + (1 - H) * sum(abs(y - xstar), 2);
        end

        function f = bweierstrass(X, hGO)
            base = 5.0;
            H = (base ^ hGO - 1) / (base - 1);
            D = size(X, 2);
            N = size(X, 1);
            a = 0.5;
            b = 3.0;
            k = (0:20)';
            Bcoef = 2 * pi * b .^ k;      % 21 x 1
            Acoef = a .^ k;               % 21 x 1
            % sum_d sum_k a^k cos(2*pi*b^k * (X(i,d)+0.5))
            M = cos((X(:) + 0.5) * Bcoef');   % (N*D) x 21
            h1 = sum(reshape(M * Acoef, N, D), 2);   % N x 1
            h2 = sum(a .^ k .* cos(pi * b .^ k));
            P = sum((abs(X) - 0.5) .* (abs(X) > 0.5), 2);
            f0 = h1 - D * h2 + P;
            f = f0 * H + (1 - H) * sum(abs(X), 2);
        end
    end
end