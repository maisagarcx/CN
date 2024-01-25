function [resultado, erro] = metNewtonCotes(a, b, n, m, F)
    
    syms k;
    symK = sym('k');
    
    a = 0; b = 7; n = 1; m = 1000;
    F = (15^2/0.25)*exp((-2*k)/0.025);
    
    d(1)=2;d(2)=6;d(3)=8;d(4)=90;d(5)=288;d(6)=840;d(7)=17280;d(8)=28350;
    c(1)=1;c(2)=1;c(3)=4;c(4)=1;c(5)=3;c(6)=7;c(7)=32;c(8)=12;c(9)=19;c(10)=75;
    c(11)=50;c(12)=41;c(13)=216;c(14)=27;c(15)=272;c(16)=751;c(17)=3577;
    c(18)=1323;c(19)=2989;c(20)=989;c(21)=5888;c(22)=-928;c(23)=10496;c(24)=-4540;

    if ((n<1)||(n>8))
        erro = 1;
        error('Grau do polinômio interpolador fora do permitido.');   
    end

    if((mod(m,n) ~= 0))   
        erro = 2;
        error('Número de subintervalos não é múltiplo do grau do polinômio interpolador.');
    end
    
    p = fix(0.25*(n*(n + 2)+ mod(n,2)));
    h = (b-a)/m; % passo
    resultado = 0;
    
    for i = 0:m
        x = a + i*h;
        y = vpa(subs(F, symK, x));
        j = p + fix(0.5*n - abs(mod(i,n) - 0.5*n));
        k = 1 + fix((n - mod(i,n))/n) - fix((m - mod(i,m))/m);
        resultado = resultado + y*c(j)*k;
        % fprintf('i: %d \nx: %f \ny: %f \n(c(j)*k): %d\n\n', i, x, y, c(j)*k);
    end
    resultado = resultado * n * h/d(n);
end
