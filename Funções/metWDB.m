function [raiz, iter, erro] = metWDB(a, b, toler, maxIter)
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
    Fx = x^4 + 2*x^3 - 13*x^2 - 14*x + 24; % função analisada
    a = -5;
    b = -3;
    toler = 1.0e-10
    maxIter = 100;
    
    Fa = a^4 + 2*a^3 - 13*a^2 - 14*a + 24; 
    Fb = b^4 + 2*b^3 - 13*b^2 - 14*b + 24; % avaliar a função em a e b
    
    if (Fa * Fb) > 0;
        error('A função não muda de sinal nos extremos do intervalo dado');
    end
    
    c = b;
    Fc = Fb;
    iter = 0;

    while(1) % altera a, b e c para que b seja a melhor estimativa da raiz
        if (Fb * Fc) > 0
            c = a;
            Fc = Fa;
            d = b - a;
            e = d;
        end
        if abs(Fc) < abs(Fb)
            a = b;
            b = c;
            c = a;
            Fa = Fb;
            Fb = Fc;
            Fc = Fa;
        end
        
        tol = 2*toler*max(abs(b),1);
        z = (c - b)/2;
        fprintf("iterações = %d \n", iter);
        fprintf("a = %f \n", a);
        fprintf("c = %f \n", c);
        fprintf("b = %f \n", b);
        fprintf("F(b) = %f \n", Fb);
        fprintf("z = %f \n", z);
        fprintf("\n");
        
        % teste de convergência 
        if (abs(z) <= toler) || Fb == 0 || iter >= maxIter
            break;
        end
        
        % escolha entre interpolação e bisseção
        if (abs(e) >= tol) && (abs(Fa) > abs(Fb))
            s = Fb/Fa;
            if a == c % interpolação linear
                p = 2*z*s;
                q = 1 - s;
            else % interpolação inversa quadrática
                q = Fa/Fc;
                r = Fb/Fc;
                p = s*(2*z*q*(q-r) - (b-a)*(r-1));
                q = (q-1)*(r-1)*(s-1);
            end
            if p > 0
                q = -q;
            else
                p = -p;
            end
            if 2*p < min(3*z*q-abs(tol*q), abs(e*q)) % aceita interpolação
                e = d;
                d = p/q;
            else % usa bisseção devido à falha na interpolação
                d = z;
                e = z;
            end
        else % bisseção
            d = z;
            e = z;
        end
        
        a = b;
        Fa = Fb;
        
        if abs(d) > tol
            b = b + d;
        else
            b = b + sign(z)*tol;
        end
        
        iter = iter + 1;
        Fb = b^4 + 2*b^3 - 13*b^2 - 14*b + 24; % avaliar a função em b
    end
    
    raiz = b;
    
    if (abs(z) <= toler) || (Fb == 0)
        erro = 0;
    else
        erro = 1;
    end 
end
