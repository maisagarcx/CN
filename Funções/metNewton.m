function [raiz, iter, erro] = metNewton(x0, toler, maxIter)
    % parâmetros de entrada:
    % x0 = valor inicial
    % toler = tolerância
    % maxIter = número máximo de iterações
    
    % parâmetros de saída:
    % raiz = raiz da equação
    % iter = número de iterações gastas
    % se erro == 0 = a raiz foi encontrada
    % se erro == 1 = a raiz não foi encontrada com a tolerância e o maxIter fornecidos
    
    toler = 0.00001;
    x0 = -2;
    maxIter = 100;
    eul = 2.718281828459045;
    Fx = 12*x0 - eul^x0 + 6;
    DFx = 12 - eul^x0; % avaliar a função e sua derivada em x0
    x = x0;
    iter = 0;
    
    fprintf("iterações = %d \n", iter);
    fprintf("x = %f \n", x);
    fprintf("DFx = %f \n", DFx);
    fprintf("Fx = %f \n", Fx);
    fprintf("\n");
    
    while(1)
        deltaX = -Fx/DFx;
        x = x + deltaX;
        Fx = 12*x - eul^x + 6;
        DFx = 12 - eul^x; % avaliar a função e sua derivada em x
        iter = iter + 1;
        
        disp(['iteração: ', num2str(iter), ' x = ', num2str(x), ' Fx = ', num2str(Fx), ' DeltaX = ', num2str(deltaX)]);
    
        fprintf("iterações = %d \n", iter);
        fprintf("x = %f \n", x);
        fprintf("DFx = %f \n", DFx);
        fprintf("Fx = %f \n", Fx);
        fprintf("deltaX = %f \n", deltaX);
        fprintf("\n");
        
        if ((abs(deltaX) <= toler) && (abs(Fx) <= toler)) || (DFx == 0 || iter>= maxIter)
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
