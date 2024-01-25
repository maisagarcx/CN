function r = vetor_r(A, b, x)
    n = length(b);
    r = zeros(n, 1);
    for i = 1:n
        r(i) = b(i);
        for j = 1:n
            r(i) = r(i) - A(i, j) * x(j);
        end
    end
end
