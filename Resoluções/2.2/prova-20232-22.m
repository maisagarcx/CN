% PARTE 2 DA PROVA 2 DE CALCULO NUMERICO
% ALUNA: MAISA GARCIA NEPOMUCENO CORREA

minimo1 = [1 1 1]; % mínimo local
tolerancia1 = 0.05;
n = 3; % números de variáveis
iterMAX1 = 100;

fprintf("Para encontrar uma solução inicial pelo método de gradiente descendente\n");
fprintf("inicializei o minimo local e chamei a função com uma tolerância de\n");
fprintf("0.05 e número de iterações máximas alto.\n");

[aproximacao,~,~] = metGradienteDescendente(n, minimo1, tolerancia1, iterMAX1);
fprintf("\n\nA solução inicial, ao sair do método gradiente é:\n");
disp(aproximacao);
fprintf("Usando o método de Newton-Raphson para refinar a solução temos:\n")

tolerancia2 = 10^(-6);
iterMAX2 = 100;

syms x1 x2 x3
simbolos = [x1, x2, x3];
F1 = f1(simbolos);
F2 = f2(simbolos);
F3 = f3(simbolos);
vet_func = [F1 F2 F3];

[solucao, ~] = metNewtonRaphsonN(vet_func, simbolos, n, aproximacao, tolerancia2, iterMAX2);

function F1 = f1(v)
    F1 = v(1)^2 + v(2)^2 + v(3)^2 - 1;
end

function F2 = f2(v)
    F2 = sin(v(1)) - v(3) + cos(v(2));
end

function F3 = f3(v)
    F3 = exp(v(1)) - v(2) + v(3)^3;
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
            vetN(i) = vpa(subs(vet_func(i), simbolos, solucao));
        end

        for i=1:n
            for j=1:n
                matrixN(i,j) = vpa(subs(diff(vet_func(i),simbolos(j)), simbolos, solucao));
            end
        end

        vet_b = (-1) * vetN;

        [matrixA, ~, pivot] = dec_LU(n, matrixN);
        [L, U] = find_L_U_from_A(matrixA);
        Y = suc_subst_piv(L, vet_b, pivot);
        incremento = ret_subst(U, Y);
        
        for m=1:n
            solucao(m) = solucao(m) + incremento(m);
        end

        max_n = maxMag(incremento, n);

        % disp(solucao);
        % disp(max_n);

        if (max_n <= tolerancia)
            fprintf("\nFoi encontrada a solução com a precisão desejada\n");
            fprintf("A última iteração é: %d\n", k);
            fprintf("O último dx é aproximadamente: %.15f\n", round(solucao(1),15));
            fprintf("O último dy é aproximadamente: %.15f\n", round(solucao(2),15));
            fprintf("O último dz é aproximadamente: %.15f\n", round(solucao(3),15));
            break
        end
        k=k+1;
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

function [aproximacao,gNew,k] = metGradienteDescendente(~, minimo, tolerancia, iterMAX)
    % GRADIENTE DESCENDENTE
    %fprintf("\nAPLICAÇÃO DO MÉTODO DE GRADIENTE DESCENDENTE\n");

    % a solucão do sistema de equação não linear é o ponto de interseção 
    % entre as N equações

    syms f1(x1,x2,x3);
    f1(x1,x2,x3) = x1^2 + x2^2 + x3^2 - 1;

    syms f2(x1,x2,x3);
    f2(x1,x2,x3) = sin(x1) + cos(x2) - x3;
    %f2(x1,x2,x3) = 8*x2 - 1 + cos(x3 - x2)^2;

    syms f3(x1,x2,x3);
    f3(x1,x2,x3) = exp(x1) - x2 - x3^3;

    syms gs(x1,x2,x3);
    gs(x1,x2,x3) = ((f1(x1,x2,x3))^2) + ((f2(x1,x2,x3))^2) + ((f3(x1,x2,x3))^2);

    % cálculo das derivadas parciais para montar o gradiente

    dDeF1X1(x1,x2,x3) = diff(f1(x1,x2,x3),x1);
    dDeF2X1(x1,x2,x3) = diff(f2(x1,x2,x3),x1);
    dDeF3X1(x1,x2,x3) = diff(f3(x1,x2,x3),x1);

    dDeF1X2(x1,x2,x3) = diff(f1(x1,x2,x3),x2);
    dDeF2X2(x1,x2,x3) = diff(f2(x1,x2,x3),x2);
    dDeF3X2(x1,x2,x3) = diff(f3(x1,x2,x3),x2);

    dDeF1X3(x1,x2,x3) = diff(f1(x1,x2,x3),x3);
    dDeF2X3(x1,x2,x3) = diff(f2(x1,x2,x3),x3);
    dDeF3X3(x1,x2,x3) = diff(f3(x1,x2,x3),x3);

    % declaração do vetor gradiente

    syms gradiente(x1,x2,x3);
    gradiente(x1,x2,x3) = [2*f1(x1,x2,x3)*dDeF1X1(x1,x2,x3) + 2*f2(x1,x2,x3)*dDeF2X1(x1,x2,x3) + 2*f3(x1,x2,x3)*dDeF3X1(x1,x2,x3), ...
                           2*f1(x1,x2,x3)*dDeF1X2(x1,x2,x3) + 2*f2(x1,x2,x3)*dDeF2X2(x1,x2,x3) + 2*f3(x1,x2,x3)*dDeF3X2(x1,x2,x3), ...
                           2*f1(x1,x2,x3)*dDeF1X3(x1,x2,x3) + 2*f2(x1,x2,x3)*dDeF2X3(x1,x2,x3) + 2*f3(x1,x2,x3)*dDeF3X3(x1,x2,x3)];

    % INPUTS
    aproximacao = minimo;

    k = 1;
    while (k <= iterMAX)

        % STEP 3
        g1 = round(gs(aproximacao(1),aproximacao(2),aproximacao(3)),5); % APLICA NO PONTO
        z = round(gradiente(aproximacao(1),aproximacao(2),aproximacao(3)),5); % APLICA NO PONTO
        z0 = norm_p(z, 2); % NORMA EUCLIDIANA

        % STEP 4
        if (z0 == 0)
            fprintf("Gradiente zero!\n"); 
            fprintf("O ponto de aproximação é:\n");
            fprintf("%.15f \n", aproximacao);
            fprintf("\nA função G (soma das funções) aplicada nesse ponto é: %.15f", g1);
            fprintf("\nProcedimento concluído!");   
            return
        end

        % STEP 5
        z = z/z0; % tornando z um vetor unidade
        alpha1 = 0;
        alpha3 = 1;
        aux3 = aproximacao - alpha3*z;
        g3 = round(gs(aux3(1),aux3(2),aux3(3)), 5);

        % STEP 6
        while (g3 >= g1) % DO STEPS 7 AND 8

            % STEP 7
            alpha3 = alpha3/2;
            aux3 = aproximacao - alpha3*z;
            g3 = round(gs(aux3(1),aux3(2),aux3(3)), 5);

            % STEP 8
            if (alpha3 < (tolerancia/2))
                fprintf("\n\nNenhuma melhora!\n");
                fprintf("O ponto de aproximação é:\n");
                fprintf("%.15f ", aproximacao);
                fprintf("A função G (soma das funções) aplicada nesse ponto é: %.15f", g1);
                fprintf("\n\nProcedimento concluído! Talvez tenha um mínimo\n");   
                return
            end
        end

        % STEP 9 (precisa do IF?)
        if (g3 < g1)
            alpha2 = alpha3/2;
            aux2 = aproximacao - alpha2*z;
            g2 = round(gs(aux2(1),aux2(2),aux2(3)), 5);
        end

        % STEP 10
        % a diferença divida de Newton é usada para encontrar a quadrática
        % P(alpha) = g1 + alpha1*h1 + h3*alpha1(alpha1-alpha2)
        % que interpola h(alpha1) em alpha1 = alpha2 = alpha3 = 0
        h1 = (g2 - g1)/alpha2;
        h2 = (g3 - g2)/(alpha3 - alpha2);
        h3 = (h2 - h1) /alpha3;

        % STEP 11
        alpha0 = 0.5*(alpha2 - (h1/h3)); % o ponto crítico de P ocorre em alpha0
        aux0 = aproximacao - (alpha0 * z);
        g0 = round(gs(aux0(1),aux0(2),aux0(3)), 5);

        % STEP 12
        % encontrar alpha E [alpha0, alpha3] para que
        % g = g(aproximacao - alpha*z) = min[g0, g3]
        if(g0 <= g3)    
            if (g0 <= g2) 
                alpha = alpha0; 
            end
        elseif (g3 <= g2)
            alpha = alpha3;
        else 
            alpha = alpha2;
        end

        % STEP 13
        aproximacao = aproximacao - alpha*z; % atualiza aproximação inicial
        fprintf("\nAtualização: %d\n", k);
        fprintf("O ponto de aproximação é:\n");
        fprintf("%.15f ", aproximacao);
        gNew = round(gs(aproximacao(1),aproximacao(2),aproximacao(3)),5);
        fprintf("\nA função G (soma das funções) aplicada nesse ponto é: %.15f", gNew);

        % STEP 14
        if (abs(gNew - g1) < tolerancia)
            fprintf("\n\nProcedimento concluído!\n");
            fprintf("O último ponto de aproximação é:\n");
            fprintf("%.15f ", aproximacao);
            fprintf("\nA função G (soma das funções) aplicada nesse ponto é: %.15f", gNew);
            return;
        end

        % STEP 15
        k = k + 1;
    end
    if (k > maxIter)
        fprintf("Número máximo de iterações excedido!\n");
        fprintf("O procedimento falhou!");
    end
end

function norm_p = norm_p(vet, p) 
     n = length(vet); 
     summy = 0; 
     for i=1:n 
         summy = summy + (abs(vet(i)))^p; 
     end 
     norm_p = summy^(1/p);  
end

