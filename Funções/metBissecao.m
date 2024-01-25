function [raiz, iter, erro] = metBissecao(a, b, toler, maxIter)
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
   
    % syms x;
    % F1 = 2*x^3 - cos(x + 1) - 3;
   
    % a = -1;
    % b = 2;
    % toler = 0.01;
    % maxIter = 100;
    
    Fa = subs(F1, x, a);
    Fb = subs(F1, x, b);
    
    if (Fa * Fb) > 0
        error('A função não muda de sinal nos extremos do intervalo dado');
    end
    
    deltaX = abs(b-a)/2;
    iter = 0;

    while(1)
        u = (a + b)/2;
        Fx = subs(F1, x, u);
        fprintf("iterações = %d \n", iter);
        fprintf("intervalo [a,b] = [%f,%f] \n", a,b);
        fprintf("F(a) = %f \n", Fa);
        fprintf("F(b) = %f \n", Fb);
        fprintf("x = %f \n", u);
        fprintf("F(x) = %f \n", Fx);
        fprintf("deltaX = %f \n", deltaX);
        fprintf("\n");
        
        if ((deltaX <= toler) && (abs(Fx) <= toler)) || iter >= maxIter
            break;
        end
        
        if (Fa * Fx) > 0
            a = u;
            Fa = Fx;
        else
            b = u;
        end
        
        deltaX = deltaX/2;
        iter = iter + 1;
    end
    
    raiz = u;
    
    % teste de convergência
    if (deltaX <= toler) && (abs(Fx) <= toler)
        erro = 0;
    else
        erro = 1;
    end 
    
    % plot(Fx,[a,b]);
    % hold on;
    % plot(raiz,[a,b]);
    % grid on;
end
