% GRADIENTE DESCENDENTE
fprintf("\nAPLICAÇÃO DO MÉTODO DE GRADIENTE DESCENDENTE\n");

% x + y = 3
% x² + y² = 9

% a solucão do sistema de equação não linear é o ponto de interseção 
% entre as N equações

syms f1(x1,x2);
f1(x1,x2) = x1 + x2 - 3;

syms f2(x1,x2);
f2(x1,x2) = x1^2 + x2^2 - 9;

syms g(x1,x2,x3);
g(x1,x2,x3) = ((f1(x1,x2)^2) + (f2(x1,x2))^2);

% cálculo das derivadas parciais para montar o gradiente

dDeF1X1 = diff(f1(x1,x2), x1);
dDeF2X1 = diff(f2(x1,x2), x1);

dDeF1X2 = diff(f1(x1,x2), x2);
dDeF2X2 = diff(f2(x1,x2), x2);

% declaração do vetor gradiente

syms gradiente(x1,x2);
gradiente(x1,x2,x3) = [2*f1(x1,x2)*dDeF1X1 + 2*f2(x1,x2)*dDeF2X1, ...
                      2*f1(x1,x2)*dDeF1X2 + 2*f2(x1,x2)*dDeF2X2];

% INPUTS
aproximacao = [0,0]; % mínimo local
% segundo o Joselias, plota o gráfico e espitualiza a reposta
tolerancia = 10^(-2);
n = 2; % números de variáveis
iterMAX = 30;

k = 1;
while (k <= iterMAX)

    % STEP 3
    g1 = round(g(aproximacao(1),aproximacao(2), 5)); % APLICA NO PONTO
    z = round(gradiente(aproximacao(1),aproximacao(2),5)); % APLICA NO PONTO
    z0 = norm_p(z, 2); % NORMA EUCLIDIANA

    % STEP 4
    if (z0 == 0)
        fprintf("Gradiente zero!\n"); 
        fprintf("O ponto de aproximção é:\n");
        fprintf("%.15f ", aproximacao);
        fprintf("A função G (soma das funções) aplicada nesse ponto é: %.15f", g1);
        fprintf("Procedimento concluído!");   
        return
    end

    % STEP 5
    z = z/z0; % tornando z um vetor unidade
    alpha1 = 0;
    alpha3 = 1;
    aux3 = aproximacao - alpha3*z;
    g3 = round(g(aux3(1),aux3(2), 5));

    % STEP 6
    while (g3 >= g1) % DO STEPS 7 AND 8

        % STEP 7
        alpha3 = alpha3/2;
        g3 = round(g(aux3(1),aux3(2), 5));

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
        g2 = round(g(aux2(1),aux2(2), 5));
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
    g0 = round(g(aux0(1),aux0(2), 5));

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
    fprintf("\n\nAtualização: %d\n", k);
    fprintf("O ponto de aproximação é:\n");
    fprintf("%.15f ", aproximacao);
    gNew = round(g(aproximacao(1),aproximacao(2),5));
    fprintf("\nA função G (soma das funções) aplicada nesse ponto é: %.15f", gNew);

    % STEP 14
    if (abs(gNew - g1) < tolerancia)
        fprintf("\nProcedimento concluído!\n");
        fprintf("O último ponto de aproximção é:\n");
        fprintf("%.15f ", aproximacao);
        fprintf("\n\nA função G (soma das funções) aplicada nesse ponto é: %.15f", gNew);
        return;
    end

    % STEP 15
    k = k + 1;
end
fprintf("Número máximo de iterações excedido!\n");
fprintf("O procedimento falhou!");

function norm_p = norm_p(vet, p) 
     n = length(vet); 
     summy = 0; 
     for i=1:n 
         summy = summy + abs(vet(i))^p; 
     end 
     norm_p = summy.^(1/p);  
 end
