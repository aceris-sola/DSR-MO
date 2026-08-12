function D = pdist2(X, Y)
    % pdist2 --- 计算 X (nx x d) 与 Y (ny x d) 两两之间的欧氏距离，返回 nx x ny。
    % 替代 Octave 的 statistics 包，避免额外依赖。
    nx = size(X, 1);
    ny = size(Y, 1);
    D = zeros(nx, ny);
    for i = 1:nx
        D(i, :) = sqrt(sum((X(i, :) - Y) .^ 2, 2))';
    end
end