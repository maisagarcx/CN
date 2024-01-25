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
