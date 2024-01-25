function matrix_transposta = find_transp(A)
    [m, n] = size(A);
    matrix_transposta = zeros(n, m);
    for i = 1:m
        for j = 1:n
            matrix_transposta(j, i) = A(i, j);
        end
    end
end
