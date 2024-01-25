function [raiz, iter, erro] = metMuller(a, c, toler, maxIter)
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
    c = 2;
    toler = 0.01;
    maxIter = 100;
    
    Fa = 2*a^3 - cos(a + 1) - 3;
    Fc = 2*c^3 - cos(c + 1) - 3; % avaliar a função em a e c
    
    
    b = (a + c)/2;
    Fb = 2*b^3 - cos(b + 1) - 3; % avaliar a função em a, b e c
    
    x = b;
    Fx = Fb;
    deltaX = c - a;
    iter = 0;
    
    while(1)
        h1 = c - b;
        h2 = b - a;
        r = h1/h2;
        t = x;
    
        A = (Fc - (r + 1)*Fb + r*Fa)/(h1*(h1 + h2));
        B = (Fc - Fb)/h1 - A*h1;
        C = Fb;
        z = (-B + sign(B)*sqrt(B^2 - 4*A*C))/(2*A);
        x = b + z;
        deltaX = x - t;
        Fx = 2*x^3 - cos(x + 1) - 3; % avaliar a função em x
        fprintf("iterações = %d \n", iter);
        fprintf("intervalo [a,b] = [%f,%f] \n", a,b);
        fprintf("F(a) = %f \n", Fa);
        fprintf("F(b) = %f \n", Fb);
        fprintf("x = %f \n", x);
        fprintf("F(x) = %f \n", Fx);
        fprintf("deltaX = %f \n", deltaX);
        fprintf("\n");
        
        if (abs(deltaX) <= toler && abs(Fx) <= toler) || iter >= maxIter
            break;
        end
    
        if x > b
            a = b;
            Fa = Fb;
        else
            c = b;
            Fc = Fb;
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
