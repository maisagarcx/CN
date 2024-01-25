function [raiz, iter, erro] = metSecante(a, b, toler, maxIter)
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
    
    % a = -1;
    % b = 2;
    % toler = 0.01;
    % maxIter = 100;
    
    syms x
    % F1 = 2*x^3 - cos(x + 1) - 3;

    Fa = subs(F1, x, a);
    Fb = subs(F1, x, b);
    
    if abs(Fa) < abs(Fb) % troca valores
        t = a;
        a = b;
        b = t;
        t = Fa;
        Fa = Fb;
        Fb = t;
    end
    
    iter = 0;
    u = b;
    Fx = Fb;

    while(1)
        deltaX = -Fx/(Fb - Fa)*(b - a);
        u = u + deltaX;
        Fx = subs(F1, x, u);  % avaliando a função em u
        fprintf("iterações = %d \n", iter);
        fprintf("intervalo [a,b] = [%f,%f] \n", a,b);
        fprintf("F(a) = %f \n", Fa);
        fprintf("F(b) = %f \n", Fb);
        fprintf("x = %f \n", u);
        fprintf("F(x) = %f \n", Fx);
        fprintf("deltaX = %f \n", deltaX);
        fprintf("\n");
        
        if ((abs(deltaX) <= toler) && (abs(Fx) <= toler)) || iter >= maxIter
            break;
        end
        
        a = b;
        Fa = Fb;
        b = u;
        Fb = Fx;
        iter = iter + 1;
    end
    
    raiz = u;
    % teste de convergência
    if (abs(deltaX) <= toler) && (abs(Fx) <= toler)
        erro = 0;
    else
        erro = 1;
    end
end
