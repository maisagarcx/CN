function [pivot, matrix, det] = dec_LU_EYE(order, matrix)

    pivot = eye(order);
    det = 1;
    for j = 1:order - 1
        % Choice of pivot element
        p = j;
        A_max = abs(matrix(j, j));
        for k = j + 1:order
            if abs(matrix(k, j)) > A_max
                A_max = abs(matrix(k, j));
                p = k;
            end
        end
        if p ~= j
            % Change rows
            matrix([j, p], :) = matrix([p, j], :);
            pivot([j, p], :) = pivot([p, j], :);
            det = -det;
        end
        det = det * matrix(j, j);
        if (abs(matrix(j, j)) ~= 0)
            % Gauss elimination
            r = 1 / matrix(j, j);
            for i = j + 1:order
                m = matrix(i, j) * r;
                matrix(i, j) = m;
                for k = j + 1:order
                    matrix(i, k) = matrix(i, k) - m * matrix(j, k);
                end
            end
        end
    end
    det = det * matrix(order, order);
end
