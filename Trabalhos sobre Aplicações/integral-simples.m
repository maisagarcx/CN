clear; clc; close all;

% TRABALHO 3.1 DE CÁLCULO NUMÉRICO 2023.2
% ALUNA: MAÍSA GARCIA NEPOMUCENO CORRÊA

syms t;
max = 50; toler = 1e-10;

fprintf("Implementação referente ao Trabalho 3.1 de Cálculo Numérico 2023.2.\n\n");
fprintf("Calcular, usando métodos computacionais de integração, a energia\n");
fprintf("dissipada em um circuto RC.");

t0 = 0; tf = 7; % Limites de integração
C = 0.1; % Capacitância 
R = 1/(8+12) + 1/5; % Resistência
tau = R*C; % Constante de tempo capacitiva
V0 = 15; % Voltagem inicial

fprintf("\n\nTemos como limites de integração o tempo inicial e final, já que a potência\n");
fprintf("é dada em função do tempo t. Logo, t₀ = %i e tf = %i.", t0, tf);
fprintf("\nTemos também que a capacitância C é %.1f F, e que a resistência equivalente R é %.2f ohms.\n", C, R);
fprintf("A DDP inicial é %i V e a constante de tempo capacitiva (R*C) é %.3f.", V0, tau);

F = (V0^2/R)*exp((-2*t)/tau); % Função da potência

plotFF(t0, 0.1, 500, t, F); % Plotando P(t) x t

fprintf("Para o método de Newton-Cotes, usamos os graus 1,2 e 3 com 400-300 pontos.\n");
fprintf("Para o método de Gauss-Legendre não iterativo, também foi utilizado 400 pontos.\n");
fprintf("Já para o método iterativo de Gauss-Legendre, o máximo de iterações foi %i\n", max);
fprintf("e a tolerância foi de %.1i.", toler);

fprintf("\n\nResultados:");
fprintf("\nNewton-Cotes com polinômio interpolador de grau 1: %.5f", metNewtonCotes(t0, tf, 1, 400, F));
fprintf("\nNewton-Cotes com polinômio interpolador de grau 2: %.5f", metNewtonCotes(t0, tf, 2, 400, F));
fprintf("\nNewton-Cotes com polinômio interpolador de grau 3: %.5f", metNewtonCotes(t0, tf, 3, 300, F));
fprintf("\nGauss-Legendre não-iterativo: %.5f", metGaussLegendre(t0, tf, 400, F));
fprintf("\nGauss-Legendre iterativo: %.5f", metGaussLegendreIterativo(t0, tf, toler, max, F));

function plotFF(a, b, n, t, y)
    x = linspace(a, b, n);
    Y = vpa(subs(y,t,x));
    figure;
    plot(x, Y, 'LineWidth', 2);
    hold on;
    fill([x,fliplr(x)], [Y,zeros(size(Y))], 'b', 'FaceAlpha', 0.3);
    xlabel('t');
    ylabel('P(t)');
    legend('Função da Potência', 'Energia Dissipada');
    hold off;
end

function [resultado, erro] = metNewtonCotes(a, b, n, m, F)
    
    syms t;
    t = sym('t');
    
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
        y = vpa(subs(F, t, x));
        j = p + fix(0.5*n - abs(mod(i,n) - 0.5*n));
        k = 1 + fix((n - mod(i,n))/n) - fix((m - mod(i,m))/m);
        resultado = resultado + y*c(j)*k;
        % fprintf('i: %d \nx: %f \ny: %f \n(c(j)*k): %d\n\n', i, x, y, c(j)*k);
    end
    resultado = resultado * n * h/d(n);
end

function [A, T, erro] = givesPesoAbsGL(n)
    
    if(n < 1)
        erro = 1;
        error('Número de pontos menor que 1.');
    end
    
    m = fix(0.5*(n+1));
    erro = 0;

    for i = 1:m    
        z = cos(pi*(i-0.25)/(n+0.5));
        while (1)
            p1 = 1;     
            p2 = 0;      
            for j = 1:n % Polinômio de Legendre no ponto z          
                p3 = p2;          
                p2 = p1;            
                p1 = ((2*j-1)*z*p2 - (j-1)*p3)/j;      
            end      
            % Derivada do polinômio de Legendre no ponto z       
            pp = n*(z*p1-p2)/((z^2)-1);       
            z1 = z;      
            % Método de Newton para calcular os zeros do polinômio    
            z = z1 - (p1/pp);     
            if abs(z-z1) < 10e-15        
                break;
            end  
        end  
        T(m+1-i) = z; % Abscissa
        A(m+1-i) = 2/((1-z^2)*(pp^2)); % Peso
    end
end

function [resultado, erro] = metGaussLegendre(a, b, n, F)
    
    syms t;
    t = sym('t');
    
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
        ti = sign(k)*T(abs(k));
        x = e1*ti + e2;
        y = vpa(subs(F, t, x));
        c = A(abs(k));
        resultado = resultado + y*c;
        % fprintf('i: %d \nt: %f \nx: %f \ny: %f\nc: %f \n\n', i, ti, x, y, c);
    end
    resultado = e1*resultado;
end

function [resultado, delta, erro] = metGaussLegendreIterativo(a, b, toler, maxIter, F)
    
    syms t;
    t = sym('t'); 
    
    iter = 1;
    n1 = 5;
    n2 = 8;
    [resultado1, erro] = metGaussLegendre(a, b, n2, F);
    % fprintf("iter: %d \nn2: %d \nintegral: %f \n", iter, n2, resultado1);
    
    % Sucessivos cálculos das integrais
    while (1)
        iter = iter + 1;
        n = n1 + n2;
        [resultado, erro] = metGaussLegendre(a, b, n, F);
        
        if (resultado ~= 0)
            delta = abs((resultado - resultado1)/resultado);
        else
            delta = abs(resultado - resultado1);
        end
        % fprintf("iter: %d \nn: %d \nintegral: %f \ndelta: %f\n", iter, n, resultado, delta);
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
