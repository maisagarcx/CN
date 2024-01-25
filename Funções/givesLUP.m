function [L, U, P] = givesLUP(A)
    %P usado na funcao de substituicoes sucessivas pivotal (suc_subst_piv)
    [~, order] = size(A); %encontra a ordem da matriz A
    P = 1:order; %inicializa um vetor de permutação do tamanho de A para controlar as pivotacoes 
    U = A; %inicializa Upper como identica a A
    L = eye(order); %inicializa Lower como uma matriz identidade de ordem n
    for j = 1:order-1 %loop responsavel por fazer a eliminacao de Gauss e pivotacao parcial
        [~,position] = max(abs(U(j:order, j))); %encontra a o indice da coluna onde esta o maior valor absoluto usando submatrizes de U
        position = position + j - 1; %ajusta m para refletir a posicao do pivo na matriz original U
        if position ~= j %se o indice encontrado nao estiver na diagonal principal, fazemos a troca de linhas
            U([j position], :) = U([position j], :); %trocamos a linha position pela linha j, para garantir que o pivo esteja na diagonal de U
            %L([j position], 1:j-1) = L([position j], 1:j-1); SE DER ERRADO DESCOMENTA
            P([j position]) = P([position j]); %troca os elementos position e j
        end
        for i = j+1:order %eliminacao de Gauss
            L(i, j) = U(i, j) / U(j, j); %guarda os multiplicadores na Lower
            U(i, :) = U(i, :) - L(i, j) * U(j, :); %zerando os elementos abaixo do pivo usando os multiplicadores
        end
    end
end
