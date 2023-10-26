function [Q, R] = GS(A, normalized)
    [m, n] = size(A);
    Q = zeros(m, n);
    R = zeros(n, n);
    
    for j = 1:n
        v = A(:, j);
        for i = 1:j-1
            R(i, j) = Q(:, i)' * A(:, j);
            v = v - R(i, j) * Q(:, i);
        end
        if normalized
            R(j, j) = norm(v);
            Q(:, j) = v / R(j, j);
        else
            Q(:, j) = v;
            R(j, j) = norm(v);
        end
    end
end
