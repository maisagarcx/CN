clear; clc;

tolerancia = 10^(-6);
maxIter = 30;
%solucao = [0.001158700814233 -0.123895211968722 0.073039656247272];
solucao = [0.10251921584124164724846158739186, -0.11967786050292097079129821563958, 0.075466196327028590414204970851599];
%solucao = [0 0 0];
n = 3;

syms x1 x2 x3

simbolos = [x1, x2, x3];

F1 = f1(simbolos);
F2 = f2(simbolos);
F3 = f3(simbolos);

vet_func = [F1 F2 F3];

[solucao, k] = metNewtonRaphsonN(vet_func, simbolos, n, solucao, tolerancia, maxIter);

function F1 = f1(v)
    F1 = (9 - v(1) + sin(v(1) + v(2)));
end

function F2 = f2(v)
    F2 = cos(v(2) - v(3))^2 - 8*v(2) + 1;
end

function F3 = f3(v)
    F3 = sin(v(3)) - 1 + 12*v(3);
end

function [solucao, k] = metNewtonRaphsonN(vet_func, simbolos, n, solucao, tolerancia, maxIter)
    
    matrixN = zeros(n,n);
    vetN = zeros(1,n);
    k = 0;
    while (true)
        k = k + 1;
        if (k > maxIter)
            fprintf("\n\nO número de iterações excedeu o permitido!");
            break
        end

        for i=1:n
            vetN(i) = subs(vet_func(i), simbolos, solucao);
        end

        for i=1:n
            for j=1:n
                matrixN(i,j) = vpa(subs(diff(vet_func(i),simbolos(j)), simbolos, solucao));
            end
        end

        vet_b = (-1) * vetN;
        %disp(vet_b);

        % [LU, ~, pivos] = calculaPivotacaoParcial(matrixA, vet_b, n, pivos);
        %[matrixA, ~, pivot] = calculaPivotacaoParcial(matrixN, vet_b, n);
        % [L, U, P] = givesLUP(matrixN);
        %disp(matrixN);
        [matrixA, ~, pivot] = dec_LU(n, matrixN);
        %disp(matrixA);
        % matrix = dec_LU(n, matrixA);
        % disp(matrix);
        [L, U] = find_L_U_from_A(matrixA);
        %disp(L);
        %disp(U);
        Y = suc_subst_piv(L, vet_b, pivot);
        incremento = ret_subst(U, Y);
        
        % disp(incremento);
        
        for k=1:n
            solucao(k) = solucao(k) + incremento(k);
        end

        max_n = maxMag(incremento, n);
        %disp(max_n);

        if (max_n <= tolerancia)
            fprintf("Foi encontrada a solução com a precisão desejada\n");
            %fprintf("A última iteração é: %d\n", iteracao);
            fprintf("O último dx é aproximadamente: %.15f\n", round(solucao(1),15));
            fprintf("O último dy é aproximadamente: %.15f\n", round(solucao(2),15));
            fprintf("O último dz é aproximadamente: %.15f\n", round(solucao(3),15));
            break
        end
        k=k+1;
    end
end

function [decomposicaoPorPivotacaoParcial, vetPermutado, Pivots] = calculaPivotacaoParcial(matriz, vetor, tamMat)
    pivots = 1:tamMat;
    for j=1:1:tamMat-1
        p = j;
        maximo = round((abs(matriz(j,j))),15);

        for k=j+1:1:tamMat      % armazena o maior elemento de cada coluna coluna
            if abs(matriz(k,j)) > maximo
                maximo = round((abs(matriz(k,j))),15);
                p=k;
            end
        end

        if p~=j    % permuta as linhas se o elemento pivô não estiver na sua devida posição
            for k=1:1:tamMat
                aux1 = matriz(j,k);
                matriz(j,k) = matriz(p,k);
                matriz(p,k) = aux1;
            end
            
            aux4 = vetor(j,1);
            vetor(j,1) = vetor(p,1);
            vetor(p,1) = aux4;

            aux2 = pivots(j,1);
            pivots(j,1) = pivots(p,1);
            pivots(p,1) = aux2;
    
        end

        if abs(matriz(j,j)) ~= 0
            aux3 = round((1 / matriz(j,j)),15);
            for l=j+1:1:tamMat
                multiplicador = round((matriz(l,j) * aux3),15);
                matriz(l,j) = multiplicador;
                for m=j+1:1:tamMat
                    matriz(l,m) = round((matriz(l,m) - multiplicador*matriz(j,m)),15);
                end
            end
        end
    end

    decomposicaoPorPivotacaoParcial = matriz;
    vetPermutado = vetor;
    Pivots = pivots;

end

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
            L([j position], 1:j-1) = L([position j], 1:j-1); %SE DER ERRADO DESCOMENTA
            P([j position]) = P([position j]); %troca os elementos position e j
        end
        for i = j+1:order %eliminacao de Gauss
            L(i, j) = U(i, j) / U(j, j); %guarda os multiplicadores na Lower
            U(i, :) = U(i, :) - L(i, j) * U(j, :); %zerando os elementos abaixo do pivo usando os multiplicadores
        end
    end
end

function [matrix, det, pivot] = dec_LU(order, matrix)
    pivot = zeros(order,1);
    for i=1:order
        pivot(i)=i;
    end
    det=1;
    for j=1:order-1
        %choice of pivot element
        p=j;
        A_max=abs(matrix(j,j));
        for k=j+1:order
            if abs(matrix(k,j))>A_max
                A_max=abs(matrix(k,j));
                p=k;
            end
        end
        if p~=j
            %change rows
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
            %Gauss elimination
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

function max_n = maxMag(vetor, n)
    for i=1:1:n
        vetor(i) = round((abs(vetor(i))),15);
    end
    max_n = max(abs(vetor));
end

function [matriz_L, matriz_U] = find_L_U_from_A(matriz_composta)
    ordem_da_matriz = length(matriz_composta);
    matriz_L = zeros(ordem_da_matriz);
    matriz_U = zeros(ordem_da_matriz);

    for linha = 1:ordem_da_matriz
        for coluna = 1:ordem_da_matriz
            if linha == coluna
                matriz_L(linha, coluna) = 1.0;
            end
            if coluna < linha
                matriz_L(linha, coluna) = matriz_composta(linha, coluna);
            end
            if coluna >= linha
                matriz_U(linha, coluna) = matriz_composta(linha, coluna);
            end
        end
    end
end

function X = suc_subst_piv(lower_tri_matrix_uni, ind_vet, pivot)
    %pode receber vetor linha ou coluna no ind_vet
    %o P usado deve ser um com indices das linhas trocados e nao as linhas da matriz identidade    
    % solves lower triangular matrix systems using LX=Pc    
    % lower_tri_matrix*X=permutation_matrix*ind_vet
    
    [order,~] = size(lower_tri_matrix_uni);
    %X=zeros(order,1); %para devolver vetor coluna
    X=zeros(1,order); %para devolver vetor linha
    k=pivot(1);
    X(1)=ind_vet(k);
    for i=2:order
        summy=0;
        for j=1:i-1
            summy= round((summy+ round((lower_tri_matrix_uni(i,j)*X(j)),15)),15);
        end
        k=pivot(i);
        X(i)=round((ind_vet(k)-summy), 15);
    end
end

function X = ret_subst(upper_tri_matrix, vet_ans)
    %aceita vetor linha ou coluna como vet_ans
    %solves upper triangular matrix systems using UX=D
    %UX=D, is upper_tri_matrix*X=vet_ans
    
    [order,~] = size(upper_tri_matrix);
    %X=zeros(order,1); %para devolver vetor coluna
    X=zeros(1,order); %para devolver vetor linha
    X(order) = vet_ans(order)/upper_tri_matrix(order,order);
    for i=order-1:-1:1
        summy=0;
        for j=i+1:order
            summy=summy+upper_tri_matrix(i,j)*X(j);
        end
        X(i)=(vet_ans(i)-summy)/upper_tri_matrix(i,i);
    end
end
