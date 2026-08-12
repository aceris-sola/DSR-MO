classdef ProblemMM
    % ProblemMM --- 从官方 Python 代码 (ProblemMM.py) 移植到 Octave。
    % CEC2026 多模态优化基准问题类（单目标最小化，多个全局最小点）。
    properties
        pid
        funID
        insNo
        dim
        maxEval
        lowBound
        upBound
        hardGO
        hardNU
        rotAngleLocalCoef
        rotAngleGlobalCoef
        lambda0
        usedEval
        sigmaW
        dMin
        maxEvalCoef
        numGlobMin
        globMinX
        globMinF
        globMinHard
        globMinRangeCoef
        rotMatGlobal
        rotMat
        nichRad
        crowdBasinInd
        specialNo
    end

    methods
        function obj = ProblemMM(pid, insNo, dim)
            % 构造函数：读取 pid 元数据并调用 form() 生成问题数据。
            dataDir = fileparts(mfilename('fullpath'));
            data = csvread(fullfile(dataDir, 'data', 'pidData.csv'));
            obj.pid = pid;
            obj.funID = data(pid, 1);
            obj.insNo = insNo;
            obj.dim = dim;
            obj.lowBound = -5;
            obj.upBound = 5;
            obj.numGlobMin = data(pid, 2);
            obj.hardGO = data(pid, 3:4);
            obj.hardNU = data(pid, 5);
            obj.rotAngleLocalCoef = pi / 2;
            obj.rotAngleGlobalCoef = pi;
            obj.lambda0 = data(pid, 6);
            obj.dMin = 0.3 * dim ^ 0.5;
            obj.usedEval = 0;
            obj.sigmaW = 0.5;
            obj.maxEvalCoef = 50000;
            obj.globMinRangeCoef = 0.9;
            obj = obj.form();
        end

        function obj = form(obj)
            dataDir = fileparts(mfilename('fullpath'));
            numUniform = csvread(fullfile(dataDir, 'data', 'num-uniform.csv'));
            numNormal  = csvread(fullfile(dataDir, 'data', 'num-normal.csv'));
            sequences  = csvread(fullfile(dataDir, 'data', 'sequences.csv'));
            fstarData  = csvread(fullfile(dataDir, 'data', 'fstarData.csv'));

            indUni = 1; indNorm = 1;
            obj.specialNo = 15 * (obj.pid - 1) + obj.insNo;
            useSeq = sequences(obj.specialNo, :);   % 已是 1-based 索引

            % 选取全局最小点的位置
            temp = numUniform(useSeq(indUni : indUni + 5 * obj.dim * obj.numGlobMin - 1));
            indUni = indUni + 5 * obj.dim * obj.numGlobMin;
            randX = reshape(temp, obj.numGlobMin * 5, obj.dim);

            % 重分布参考点 Xref 的索引
            obj.crowdBasinInd = ceil(numUniform(useSeq(indUni)) * obj.numGlobMin);
            indUni = indUni + 1;

            % 设置全局最小点（先最远点挑选，再按硬度收缩重分布）
            [uniformX, ~] = UtilityMethod.keep_farthest(randX, obj.numGlobMin);
            uniformX = uniformX(1:obj.numGlobMin, :) * obj.globMinRangeCoef + ...
                       (1 - obj.globMinRangeCoef) / 2;
            uniformX = uniformX * (obj.upBound - obj.lowBound) + obj.lowBound;
            [obj.globMinX, ~] = UtilityMethod.redist_glob_min( ...
                uniformX, uniformX(obj.crowdBasinInd, :), obj.hardNU, obj.dMin);

            % 每个全局最小点的硬度
            [~, ind0] = sort(numUniform(useSeq(indUni : indUni + obj.numGlobMin - 1)));
            indUni = indUni + obj.numGlobMin;
            coef = (ind0 - 1) / (obj.numGlobMin - 1);
            obj.globMinHard = coef * (obj.hardGO(2) - obj.hardGO(1)) + obj.hardGO(1);

            % niching 半径
            if obj.numGlobMin == 1
                obj.nichRad = 5 * sqrt(obj.dim);
            else
                tmp = pdist2(obj.globMinX, obj.globMinX);
                tmp = tmp + max(tmp(:)) * eye(obj.numGlobMin);
                obj.nichRad = min(tmp, [], 1) / 2.0;
            end

            obj.maxEval = round(obj.maxEvalCoef * obj.dim);
            obj.globMinF = fstarData(obj.funID);

            % 旋转矩阵
            temp0 = 2 * obj.dim * (obj.numGlobMin + 1);
            tempUV = numNormal(useSeq(indNorm : indNorm + temp0 - 1));
            indNorm = indNorm + temp0;
            allUV = reshape(tempUV, 2 * (obj.numGlobMin + 1), obj.dim);
            allAngleData = numNormal(useSeq(indNorm : indNorm + obj.numGlobMin));
            indNorm = indNorm + obj.numGlobMin + 1;

            if obj.rotAngleGlobalCoef ~= 0 && obj.dim > 1
                rotAngleGlobal = obj.rotAngleGlobalCoef * allAngleData(1);
                u0 = allUV(1, :); v0 = allUV(2, :);
                obj.rotMatGlobal = UtilityMethod.gen_rot_mat_pseudo(u0, v0, rotAngleGlobal);
            else
                obj.rotMatGlobal = eye(obj.dim);
            end

            obj.rotMat = cell(obj.numGlobMin, 1);
            if obj.rotAngleLocalCoef > 0 && obj.dim > 1
                for k = 1:obj.numGlobMin
                    rotAngleLocal = obj.rotAngleLocalCoef * allAngleData(k + 1);
                    u0 = allUV(2 * k + 1, :);
                    v0 = allUV(2 * k + 2, :);
                    R0 = UtilityMethod.gen_rot_mat_pseudo(u0, v0, rotAngleLocal);
                    obj.rotMat{k} = obj.rotMatGlobal * R0;
                end
            else
                for k = 1:obj.numGlobMin
                    obj.rotMat{k} = obj.rotMatGlobal;
                end
            end
        end

        % 目标函数：接受 N x D 矩阵（每行一个解），返回 N x 1 列向量。
        % 全向量化：同时对 N 个解和 K 个全局峰做批量评估，降低解释开销。
        function f = func_eval(obj, x0)
            x = x0;
            if isvector(x0)
                x = x0(:)';
            end
            N = size(x, 1);
            K = obj.numGlobMin;

            % ---- 各峰基础函数值 F (N x K) ----
            F = zeros(N, K);
            for k = 1:K
                shift = obj.globMinX(k, :);
                xr = (x - shift) * obj.rotMat{k};          % N x D
                F(:, k) = BasicFun.evaluate_batch(xr / obj.lambda0, ...
                                                  obj.globMinHard(k), obj.funID);
            end

            % ---- 权重 W (N x K) ----
            dis = zeros(N, K);
            for k = 1:K
                dis(:, k) = sqrt(sum((x - obj.globMinX(k, :)) .^ 2, 2));
            end
            if numel(obj.nichRad) > 1
                nich = obj.nichRad(:)';                    % 1 x K
            else
                nich = repmat(obj.nichRad, 1, K);
            end
            normDis2 = (dis ./ (obj.sigmaW * nich)) .^ 2;
            normDis2min = min(normDis2, [], 2);            % N x 1
            C0 = 1 - normDis2min;
            C0(normDis2min <= 1) = 0;
            W = exp(-normDis2 - C0);
            maxW = max(W, [], 2);                          % N x 1
            term = abs(W - maxW) < 1e-14;
            W = W .* (1 - maxW .^ 10) .* (1 - term) + W .* term;
            W = W ./ sum(W, 2);

            f = sum(F .* W, 2) + obj.globMinF;
            obj.usedEval = obj.usedEval + N;
        end

        % 单点评估兼容接口（内部转批量）
        function f = func_eval_single(obj, x)
            f = obj.func_eval(x);
        end
    end
end