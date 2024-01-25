function [raiz, iter, erro] = metSchroder(m, x0, toler, maxIter)
    % parâmetros de entrada:
    % x0 = valor inicial
    % toler = tolerância
    % maxIter = número máximo de iterações
    
    % parâmetros de saída:
    % raiz = raiz da equação
    % iter = número de iterações gastas
    % se erro == 0 = a raiz foi encontrada
    % se erro == 1 = a raiz não foi encontrada com a tolerância e o maxIter fornecidos
    
    toler = 1e-5;
    x0 = 1.5;
    maxIter = 100;
    m = 3;
    
    Fx = x0^4 + 2*x0^3 - 12*x0^2 + 14*x0 - 5;
    DFx = 4*x0^3 + 6*x0^2 - 24*x0 + 14; % avaliar a função e sua derivada em x0
    x = x0;
    iter = 0;
    fprintf("iterações = %d \n", iter);
    fprintf("x = %f \n", x);
    fprintf("DFx = %f \n", DFx);
    fprintf("Fx = %f \n", Fx);
    fprintf("\n");
    
    while(1)
        deltaX = -m*Fx/DFx;
        x = x + deltaX;
        Fx = x^4 + 2*x^3 - 12*x^2 + 14*x - 5;
        DFx = 4*x^3 + 6*x^2 - 24*x + 14; % avaliar a função e sua derivada em x
        iter = iter + 1;
        
        disp(['iteração: ', num2str(iter), ' x = ', num2str(x), ' Fx = ', num2str(Fx), ' DeltaX = ', num2str(deltaX)]);
        
        fprintf("iterações = %d \n", iter);
        fprintf("x = %f \n", x);
        fprintf("DFx = %f \n", DFx);
        fprintf("Fx = %f \n", Fx);
        fprintf("deltaX = %f \n", deltaX);
        fprintf("\n");
        
        if ((abs(deltaX) <= toler) && (abs(Fx) <= toler)) || DFx == 0 || iter>= maxIter
            break;
        end
    end
    
    raiz = x;
    % teste de convergência
    if (abs(deltaX) <= toler) && (abs(Fx) <= toler)
        erro = 0;
    else
        erro = 1;
    end
end
