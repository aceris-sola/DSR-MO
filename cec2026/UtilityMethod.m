classdef UtilityMethod
    % UtilityMethod --- 从官方 Python 代码 (UtilityMethod.py) 移植到 Octave。
    % 提供 CEC2026 多模态基准构造所需的三个工具函数。
    methods (Static)
        % -----------------------------------------------------------------
        function R = gen_rot_mat_pseudo(u0, v0, alp)
            % 给定两个随机向量 uu, vv 生成一个随机旋转矩阵。
            % 见 [1] "On the Rigid Rotation Concept in n-Dimensional Spaces"
            %     Daniele Mortari (2001)
            D = numel(u0);
            u = u0 / norm(u0);                 % 1 x D
            v = v0 - (u * v0') * u;            % Gram-Schmidt 正交化
            v = v / norm(v);
            R = eye(D) + sin(alp) * (v' * u - u' * v) + ...
                (cos(alp) - 1) * (u' * u + v' * v);
        end

        % -----------------------------------------------------------------
        function [Y, minDis] = keep_farthest(X, n)
            % 迭代移除距离最近的点，最终保留 n 个尽可能分散的点。
            % X : N x D，每行一个解，坐标在 [0,1]。
            [N, D] = size(X);
            dis  = pdist2(X, X);
            dis2 = dis + eye(N) * sqrt(D) * 2;     % 自身距离置为极大
            keepInd = [1];
            candidInd = 2:N;
            for k = 2:n
                closest = min(dis(keepInd, candidInd), [], 1);
                [~, idx] = max(closest);           % 到已选集合最远的候选
                keepInd(end + 1) = candidInd(idx); %#ok<AGROW>
                candidInd(idx) = [];
            end
            Y = X(keepInd, :);
            minDis = min(min(dis2(keepInd, keepInd)));
        end

        % -----------------------------------------------------------------
        function [Y, countTry] = redist_glob_min(X, Xref, hardNU, disTol)
            % 按到 Xref 的距离排序，并按 rank 非线性收缩距离，重分布全局最小点，
            % 同时尽量保证任意两解之间距离 >= disTol。
            maxTry = 100;
            tauNU = 0.2;
            [N, D] = size(X);
            Y = X;
            countTry = zeros(N, 1);
            if N > 1
                Vlength = zeros(N, 1);
                V = zeros(N, D);
                for k = 1:N
                    V(k, :) = X(k, :) - Xref;
                    Vlength(k) = norm(V(k, :));
                end
                [~, ind] = sort(Vlength);          % 升序
                rnk = zeros(N, 1);
                rnk(ind) = (1:N)';                 % 按距离排序的名次
                rnk = rnk / N * (1 - tauNU) + tauNU;
                targetCoef = rnk .^ hardNU;        % 到 Xref 的距离按名次收缩
                for kk = 1:numel(ind)
                    j = ind(kk);
                    for tryNo = 0:maxTry
                        term1 = (maxTry - tryNo) / maxTry;
                        finCoef = targetCoef(j) ^ term1;
                        Y(j, :) = Xref + finCoef * V(j, :);
                        tmp = setdiff(1:N, j);
                        dis1 = pdist2(Y(j, :), Y(tmp, :));
                        countTry(j) = tryNo;
                        if kk == 1 || kk == N || min(dis1) >= disTol
                            break;                 % 接受该位置
                        end
                    end
                end
            end
        end
    end
end