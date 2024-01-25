function [canBeDiagonallyDominant, permutationMatrix] = canBeMadeDiagonalDominant(matrix)
    [numRows, numCols] = size(matrix);
    canBeDiagonallyDominant = true;
    permutationMatrix = eye(numRows); % Inicializa a matriz de permutações como a matriz de identidade
    
    for col = 1:numCols
        [~, maxRowIndex] = max(abs(matrix(col:numRows, col)));
        maxRowIndex = maxRowIndex + col - 1;

        if maxRowIndex ~= col
            % Troca as linhas na matriz de permutações
            permutationMatrix([col, maxRowIndex], :) = permutationMatrix([maxRowIndex, col], :);
            % Troca as linhas na matriz original
            matrix([col, maxRowIndex], :) = matrix([maxRowIndex, col], :);
        end

        diagonalElement = abs(matrix(col, col));
        sumOffDiagonal = sum(abs(matrix(col, :))) - diagonalElement;

        if diagonalElement <= sumOffDiagonal
            canBeDiagonallyDominant = false;
            break;
        end
    end
end
