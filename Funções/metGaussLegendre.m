function [resultado, erro] = metGaussLegendre(a, b, n, symK, F)
    
    syms k;
    symK = sym('k');
    F = (15^2/0.25)*exp((-2*k)/0.025);
    a = 0; b = 7; n = 1000;
    
    resultado = 0;
    [A, T, condErro] = givesPesoAbsGL(n);
    if (condErro ~= 0)
        erro = 1;
        error('Houve erro na função givesPesoAbsGL.');
    end
    erro = 0;
    
    % Cálculo da integral
    e1 = (b-a)/2; 
    e2 = (a+b)/2;
    
    if(mod(n,2) == 0)
        c1 = 1;
        c2 = 0.5;
    else
        c1 = 0;
        c2 = 1;
    end
    
    for i = 1:n
        k = fix(i - 0.5*(n+1) + sign(i - 0.5*(n + c1))*c2);
        t = sign(k)*T(abs(k));
        x = e1*t + e2;
        y = vpa(subs(F, symK, x));
        c = A(abs(k));
        resultado = resultado + y*c;
        % fprintf('i: %d \nt: %f \nx: %f \ny: %f\nc: %f \n\n', i, t, x, y, c);
    end
    resultado = e1*resultado;
end
