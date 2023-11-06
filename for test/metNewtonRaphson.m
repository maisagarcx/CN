% NEWTON-RAPHSON
fprintf("\nAPLICAÇÃO DO MÉTODO DE NEWTON-RAPHSON\n");

% Decomposição LU com pivotação parcial para sistemas de ordem N

syms f1(x1,x2);
f1(x1,x2) = x1^10 + x2^10 - 1024;

syms f2(x1,x2);
f2(x1,x2) = exp(x1) - exp(x2) - 1;

% x1 = -5:0.1:5;
% x2 = -5:0.1:5;
% [X1, X2] = meshgrid(x1, x2);
% F1 = X1 + X2 - 3;
% F2 = X1.^2 + X2.^2 - 9;
% figure;
% surf(X1, X2, F, 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'FaceColor', 'b');
% hold on;
% surf(X1, X2, G, 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'FaceColor', 'r');
% xlabel('X / x1');
% ylabel('Y / x2');
% zlabel('F1(x1, x2) / F2(x1, x2)');
% legend('F1(x1, x2)', 'F2(x1, x2)');

% usando que Ax = b

func = [f1(x1,x2), f2(x1,x2)];
var = [x1 , x2];

% for i=1:n
%     for j=1:n
%         matrixA(i, j) = diff(func(j), var(i));
%     end
% end

% ind_vet = [dx; dy];
vet_b = (-1)*func; 

dDeF1X1(x1,x2) = diff(f1(x1,x2), x1);
dDeF2X1(x1,x2) = diff(f2(x1,x2), x1);

dDeF1X2(x1,x2) = diff(f1(x1,x2), x2);
dDeF2X2(x1,x2) = diff(f2(x1,x2), x2);

n = 2; % número de variaveis
solucao = [0.1 -2];
tolerancia = 10^(-2);
maxIter = 30;
k = 0;

while (true)
    k = k + 1;
    if (k > maxIter)
        fprintf("\n\nO número de iterações excedeu o permitido!");
        break;
    end

    % matrix de derivadas parciais
    dDeF1X1X0 = round(dDeF1X1(solucao(1),solucao(2)),15);
    dDeF1X2Y0 = round(dDeF1X2(solucao(1),solucao(2)),15);
    dDeF2X1X0 = round(dDeF2X1(solucao(1),solucao(2)),15);
    dDeF2X2Y0 = round(dDeF2X2(solucao(1),solucao(2)),15);

    matrixA = [dDeF1X1X0 dDeF1X2Y0; dDeF2X1X0 dDeF2X2Y0];
    % matrixA(1, 1) = dDeF1X1X0;
    % matrixA(1, 2) = dDeF1X2Y0;
    % matrixA(2, 1) = dDeF2X1X0;
    % matrixA(2, 2) = dDeF2X2Y0;

    % vetor b

    F1X0Y0 = round(f1(solucao(1),solucao(2)),15);
    F2X0Y0 = round(f2(solucao(1),solucao(2)),15);

    vet_b = (-1) * [F1X0Y0; F2X0Y0];

    % vetor x
    pivos(1,1) = 1;
    pivos(2,1) = 2;

    % [LU, ~, pivos] = calculaPivotacaoParcial(matrixA, vet_b, n, pivos);
    % disp(LU);

    [~, ~, P] = givesLUP(matrixA);
    P = P';
    % disp(P);
    matrix = dec_LU(n, matrixA);
    % disp(matrix);
    [L, U] = find_L_U_from_A(matrix);
    % disp(L);
    % disp(U);
    Y = suc_subst_piv(L, vet_b, P);
    disp(Y);
    vetorSolucao = ret_subst(U, Y);
    disp(vetorSolucao);


    fprintf("\n\nIteração: %d", k);

    dx = round(vetorSolucao(1),15);
    fprintf("\nO incremento dx é aproximadamente: %.15f\n", round(dx,15));
    solucao(1) = solucao(1) + dx;

    dy = round(vetorSolucao(2),15);
    fprintf("\nO incremento dy é aproximadamente: %.15f\n", round(dy,15));
    solucao(2) = solucao(2) + dy;

    max_n = maxMag(vetorSolucao, n);

    if (max_n < tolerancia)
        fprintf("\n\nFoi encontrada a solução com a precisão desejada\n");
        fprintf("A última iteração é: %d\n", k);
        fprintf("O último dx é aproximadamente: %.15f\n", round(dx,15));
        fprintf("O último dy é aproximadamente: %.15f\n", round(dy,15));
        fprintf("A solução com a precisão desejada é aproximadamente:\n");
        fprintf("[%.15f, %.15f]\n", solucao(1), solucao(2));
        break;
    end
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
            %L([j position], 1:j-1) = L([position j], 1:j-1); SE DER ERRADO DESCOMENTA
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
    max_n = max(vetor);
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
