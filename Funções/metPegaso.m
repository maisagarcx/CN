function [raiz, iter, erro] = metPegaso(a, b, toler, maxIter)
    % parâmetros de entrada:
    % a = limite inferior do intervalo
    % b = limite superior do intervalo
    % toler = tolerância
    % maxIter = número máximo de iterações
    
    % parâmetros de saída:
    % raiz = raiz da equação
    % iter = número de iterações gastas
    % se erro == 0 = a raiz foi encontrada
    % se erro == 1 = a raiz não foi encontrada com a tolerância e o maxIter fornecidos
    
    x = 0;
    Fx = 2*x^3 - cos(x + 1) - 3; % função analisada
    a = -1;
    b = 2;
    toler = 0.01;
    maxIter = 100;
    
    Fa = 2*a^3 - cos(a + 1) - 3;
    Fb = 2*b^3 - cos(b + 1) - 3; % avaliar a função em a e b
    
    iter = 0;
    x = b;
    Fx = Fb;

    while(1)
        deltaX = -Fx/(Fb - Fa)*(b - a);
        x = x + deltaX;
        Fx = 2*x^3 - cos(x + 1) - 3;  % avaliando a função em x
        fprintf("iterações = %d \n", iter);
        fprintf("intervalo [a,b] = [%f,%f] \n", a,b);
        fprintf("F(a) = %f \n", Fa);
        fprintf("F(b) = %f \n", Fb);
        fprintf("x = %f \n", x);
        fprintf("F(x) = %f \n", Fx);
        fprintf("deltaX = %f \n", deltaX);
        fprintf("\n");
        
        if ((abs(deltaX) <= toler) && (abs(Fx) <= toler)) || iter >= maxIter
            break;
        end
        
        if (Fx * Fb) < 0
            a = b;
            Fa = Fb;
        else
            Fa = (Fa * Fb)/(Fb + Fx);
        end
        
        b = x;
        Fb = Fx;
        iter = iter + 1;
    end
    
    raiz = x;
    % teste de convergência
    if (abs(deltaX) <= toler) && (abs(Fx) <= toler)
        erro = 0;
    else
        erro = 1;
    end
end
