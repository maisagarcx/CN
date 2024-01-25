function [a, b, erro] = ondeTrocaSinal(z)
    % parâmetro de entrada:
    % z = ponto a partir do qual o intervalo será gerado
    
    % parâmetros de saída:
    % a = limite inferior do intervalo
    % b = limite superior do intervalo
    % se erro == 0 = existe um número ímpar de raízes no intervalo [a,b]
    % se erro == 1 = ou não tem raízes, ou tem número par de raízes no intervalo [a,b]
    
    %z = 5;
    
    if z == 0
        a = -0.05; 
        b = 0.05;
    else
        a = 0.95 * z;
        b = 1.05 * z;
    end
    
    iter = 0;
    aureo = 2/(sqrt(5) - 1);
    Fa = 2*a^3-cos(a+1)-3;
    Fb = 2*b^3-cos(b+1)-3; % avaliar a função em a e b
    fprintf("iterações = %d \n", iter);
    fprintf("intervalo [a,b] = [%f,%f] \n", a,b);
    fprintf("F(a) = %f \n", Fa);
    fprintf("F(b) = %f \n", Fb);
    fprintf("\n");
    
    while(1)
        if ((Fa * Fb) <= 0) || (iter >= 20)
            break;
        end
        iter = iter+1;
        
        if(abs(Fa) < abs(Fb))
            a = a - aureo*(b-a);
            Fa = 2*a^3-cos(a+1)-3;   % avaliar a função em a
        else
            b = b + aureo*(b-a);
            Fb = 2*b^3-cos(b+1)-3;   % avaliar a função em b
        end
        
        fprintf("iterações = %d \n", iter);
        fprintf("intervalo [a,b] = [%f,%f] \n", a,b);
        fprintf("F(a) = %f \n", Fa);
        fprintf("F(b) = %f \n", Fb);
        fprintf("\n");
    end
    
    if((Fa * Fb) <= 0)
        erro = 0;
    else
        erro = 1;
    end
    %plot(fx,[a,b]);
end
