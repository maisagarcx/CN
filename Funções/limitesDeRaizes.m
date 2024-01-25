function L = limitesDeRaizes(n, coef)
    % parâmetros de entrada:
    % n = grau do polinômio 
    % coef = coeficientes
    
    % parâmetros de saída:
    % L = limites inferior e superior das raízes positivas e negativas
    
    if coef(1) == 0
        disp('O coef(1) é nulo');
        return;
    end
    
    t = n+1;
    coef(t+1) = 0;
    
    while(1) % se coef(n+1) for nulo, então o polinômio é deflacionado
        if coef(t) ~= 0
            break;    
        end
        t = t-1;
    end
    
    for i=1:4 % cálculo dos quatro limites da raízes reais
        if (i == 2)||(i == 4) % inversão da ordem dos coeficientes
            for j=1:(t/2)  
                aux = coef(j);
                coef(j) = coef(t-j+1);
                coef(t-j+1) = aux;
            end
        elseif i == 3 % reinversão da ordem e troca de sinais dos coeficientes
            for j=1:(t/2)  
                aux = coef(j);
                coef(j) = coef(t-j+1);
                coef(t-j+1) = aux;
            end  
            for j=(t-1):-2:1
                coef(j) = -coef(j); % troca sinal
            end    
        end
        
        % se coef(1) for negativo, então é trocado o sinal de todos os coef
        if coef(1) < 0 
            for j=1:t
                coef(j) = -coef(j);
            end    
        end
        
        % cálculo de k, o maior índice dos coeficientes negativos
        k = 2;
        while(1)
            if (coef(k) < 0)||(k > t)
                break;
            end
            k = k+1;
        end    
        
        % cálculo de B, o maior coeficiente negativo em módulo
        if k <= t
            B = 0;
            for j=2:t
                if (coef(j) < 0) && (abs(coef(j)) > B)
                    B = abs(coef(j));
                end    
            end
            
            % limite das raízes positivas de P(x) = 0 e das equações auxiliares
            L(i) = 1 + (B/coef(1))^(1/(k-1));
        else
            L(i) = 10^100;   % limite infinito
        end      
    end
    
    % limites das raízes positivas e negativas de P(x) = 0
    aux = L(1);
    L(1) = 1/L(2);
    L(2) = aux;
    L(3) = -L(3);
    L(4) = -1/L(4);
end
