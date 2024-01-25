function A_inversa = inver(A)
    [m, n] = size(A);
    if m ~= n
        error('A matriz não é quadrada. A inversão não é possível.');
    end
    if det(A) == 0
        error('A matriz tem determinante igual a zero, não é possível inverter.');
    end
    n = length(A);
    augmented_matrix = [A eye(n)];
    for col = 1:n
        pivot = augmented_matrix(col, col);
        augmented_matrix(col, :) = augmented_matrix(col, :) / pivot;
        for row = 1:n
            if row ~= col
                factor = augmented_matrix(row, col);
                augmented_matrix(row, :) = augmented_matrix(row, :) - factor * augmented_matrix(col, :);
            end
        end
    end
    A_inversa = augmented_matrix(:, n+1:end);
end
