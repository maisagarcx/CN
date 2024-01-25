% prova 1.1 de cálculo numérico
format default;
fprintf("\n\n");

% *primeira questão*
fprintf("Primeira questão: ***************************\n");

A = [0.2 1.8 8.1 3.1; 2.1 3.7 1.5 9.7; 2.6 6.4 3.0 1.2; 5.1 1.0 0.0 2.3];
b = [2.5; 2.0; 4.7; 6.1];
order = 4;
precisao = 1e-4;

% letra a)
fprintf("letra a: ***************************\n");

% Determine as matrizes L, U e P resultantes da fatoração LU com pivotação 
% parcial e o determinante

[LU, dete, pivot] = dec_LU(order, A);
[L, U] = find_L_U_from_A(LU);
eye_pivot = dec_LU_EYE(order, A);
fprintf("A matriz L é:\n");
disp(L);
fprintf("A matriz U é:\n");
disp(U);
fprintf("A matriz de pivotação (identidade) é:\n");
disp(eye_pivot);
fprintf("O determinante da matriz A é:\n");
disp(dete);

% letra b) (LUx = b)
fprintf("letra b: ***************************\n");

% Mostre as soluções dos sistemas triangulares inferior e superior resultantes
% da solução por fatoração LU com pivotação parcial do sistema acima

y = suc_subst_piv(order, L, b, pivot);
x = ret_subst(order, U, y);  

fprintf("A solução (transposta) do sistema triangular inferior (L) é:\n");
disp(y);
fprintf("A solução (transposta) do sistema triangular superior (U) é:\n");
disp(x);

% letra c)
fprintf("letra c: ***************************\n");

% O sistema acima pode ser resolvido por algum método iterativo? Se não for 
% possível, existe alguma troca de linhas que torna possível a sua solução? 
% Por qual método e quantas iterações são necessárias para uma tolerância de 10-5?

% para ser resolvido por método iterativo, precisamos saber se a matriz A é
% diagonalmente positiva (ou pode ser feita diagonalmente positiva)

[canBeDiagonallyDominant, permutationMatrix] = canBeMadeDiagonalDominant(A);
fprintf("A matriz pode ser diagonalmente dominante (0 para falso, 1 para verdadeiro)? ");
disp(canBeDiagonallyDominant);

% nesse caso, permutationMatrix vai mostrar onde estão os maiores valores
% de cada linha de A, mas eles não são grandes o suficiente para A ser
% resolvida por método iterativo

% *segunda questão*
fprintf("Segunda questão: ***************************\n");

F = [4 1 2; 1 3 1; 2 1 5];
g = [1.1; 5.1; 0.7];
cho_order = 3;

% letra a) 
fprintf("letra a: ***************************\n");

% Encontre a solução por Cholesky, apresente a matriz L, os 
% resultados dos sistemas triangulares e determine o resíduo

[cho_L, cho_dete, error] = cholesky(cho_order, F);
fprintf("A matriz cho_L é: \n");
disp(cho_L);
cho_U = find_transp(cho_L);
fprintf("A matriz cho_U é: \n");
disp(cho_U);

cho_y = suc_subst(cho_order, cho_L, g);
cho_x = ret_subst(cho_order, cho_U, cho_y);

fprintf("A solução (transposta) do sistema triangular inferior (cho_L) é:\n");
disp(cho_y);
fprintf("A solução (transposta) do sistema triangular superior (cho_U) é:\n");
disp(cho_x);

% determinando o vetor residuo
cho_r = vetor_r(F, g, cho_x);

fprintf("O vetor resíduo é:\n");
disp(cho_r);

% *funções*

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


function [L, U] = find_L_U_from_A(matriz_composta)
    ordem_da_matriz = length(matriz_composta);
    L = zeros(ordem_da_matriz);
    U = zeros(ordem_da_matriz);

    for linha = 1:ordem_da_matriz
        for coluna = 1:ordem_da_matriz
            if linha == coluna
                L(linha, coluna) = 1.0;
            end
            if coluna < linha
                L(linha, coluna) = matriz_composta(linha, coluna);
            end
            if coluna >= linha
                U(linha, coluna) = matriz_composta(linha, coluna);
            end
        end
    end
end

function X = suc_subst_piv(order, lower_tri_matrix_uni, ind_vet, pivot)
    X=zeros(1,order);  
    k=pivot(1);
    X(1)=ind_vet(k);
    for i=2:order
        summy=0;
        for j=1:i-1
            summy=summy+lower_tri_matrix_uni(i,j)*X(j);
        end
        k=pivot(i);
        X(i)=ind_vet(k)-summy;
    end
end

function [matrix, det, pivot] = dec_LU(order, matrix)
    pivot = zeros(order,1);
    for i=1:order
        pivot(i)=i;
    end
    det=1;
    for j=1:order-1
        p=j;
        A_max=abs(matrix(j,j));
        for k=j+1:order
            if abs(matrix(k,j))>A_max
                A_max=abs(matrix(k,j));
                p=k;
            end
        end
        if p~=j
            for k=1:order
                t=matrix(j,k);
                matrix(j,k)=matrix(p,k);
                matrix(p,k)=t;
            end
            t=pivot(j);
            pivot(j)=pivot(p);
            pivot(p)=t;
            det=-det;
        end
        det=det*matrix(j,j);
        if (abs(matrix(j,j))~=0)
            r=1/matrix(j,j);
            for i=j+1:order
                m=matrix(i,j)*r;
                matrix(i,j)=m;
                for k=j+1:order
                    matrix(i,k)=matrix(i,k)-m*matrix(j,k);
                end
            end
        end
    end
    det=det*matrix(order,order);
end

function X = ret_subst(order, upper_tri_matrix, vet_ans)
    X=zeros(1,order);
    X(order) = vet_ans(order)/upper_tri_matrix(order,order);
    for i=order-1:-1:1
        summy=0;
        for j=i+1:order
            summy=summy+upper_tri_matrix(i,j)*X(j);
        end
        X(i)=(vet_ans(i)-summy)/upper_tri_matrix(i,i);
    end
end

function [canBeDiagonallyDominant, permutationMatrix] = canBeMadeDiagonalDominant(matrix)
    [numRows, numCols] = size(matrix);
    canBeDiagonallyDominant = true;
    permutationMatrix = eye(numRows);
    
    for col = 1:numCols
        [~, maxRowIndex] = max(abs(matrix(col:numRows, col)));
        maxRowIndex = maxRowIndex + col - 1;

        if maxRowIndex ~= col
            permutationMatrix([col, maxRowIndex], :) = permutationMatrix([maxRowIndex, col], :);
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

function r = vetor_r(A, b, x)
    length_b = length(b);
    r = zeros(length_b, 1);
    for i = 1:length_b
        r(i) = b(i);
        for j = 1:length_b
            r(i) = r(i) - A(i, j) * x(j);
        end
    end
end

function [L, det, error] = cholesky(order, matrix)
    %does the decomposition LLᵗ of a matrix
    det=1;
    L=zeros(order,order);

    for j=1:order
        summy=0;
        for k = 1:j-1
            summy=summy + L(j,k)^2;
        end
        t=matrix(j,j)-summy;
        det=det*t;
        error= t<=0;
        if error
            prompt = "A matriz não é definida positiva.";
            error(prompt);
        else
            L(j,j)=sqrt(t);
            r=1/L(j,j);
        end
        for i=j+1:order
            summy=0;
            for k=1:j-1
                summy=summy+L(i,k)*L(j,k);
            end
            L(i,j)=(matrix(i,j)-summy)*r;
        end
    end
end

function matrix_transposta = find_transp(any)
    [m, n] = size(any);
    matrix_transposta = zeros(n, m);
    for i = 1:m
        for j = 1:n
            matrix_transposta(j, i) = any(i, j);
        end
    end
end
