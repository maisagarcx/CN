function [resultado, delta, erro] = metGaussLegendreIterativo(a, b, toler, maxIter, symK, F)
    a = 0;
    b = 20;
    toler = 1e-10;
    maxIter = 10;
    syms k;
    symK = sym('k'); 
    F = k*sin(15*k);
    
    iter = 1;
    n1 = 5;
    n2 = 8;
    [resultado1, erro] = metGaussLegendre(a, b, n2, symK, F);
    fprintf("iter: %d \nn2: %d \nintegral: %f \n", iter, n2, resultado1);
    
    % Sucessivos cálculos das integrais
    while (1)
        iter = iter + 1;
        n = n1 + n2;
        [resultado, erro] = metGaussLegendre(a, b, n, symK, F);
        
        if (resultado ~= 0)
            delta = abs((resultado - resultado1)/resultado);
        else
            delta = abs(resultado - resultado1);
        end
        fprintf("iter: %d \nn: %d \nintegral: %f \ndelta: %f\n", iter, n, resultado, delta);
        if (delta <= toler) || (iter == maxIter)
            return;
        end
        resultado1 = resultado;
        n1 = n2;
        n2 = n;
    end
    if (delta <= toler)
        erro = 0;
    else
        erro = 1;
    end
end
