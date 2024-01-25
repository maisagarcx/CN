% TRABALHO 3.4 DE CÁLCULO NUMÉRICO 2023.2
% ALUNA: MAÍSA GARCIA NEPOMUCENO CORRÊA

fprintf("Implementação referente ao Trabalho 3.4 de Cálculo Numérico 2023.2.\n\n");

syms z k w;
R = 200; L = 0.1; C = 10e-7;

fprintf("Um circuito RLC série, com R = %d ohms, L = %i H e C = %d F,\n", R, L, C);
fprintf("tem como valor inicial I(0) = 0 e I'(0) = 1. Encontre o gráfico da corrente");

a = 0; b = 6e-3; m = 50; y10 = 0; y20 = 1;
FA = exp(-1000*z)*(sin(3000*z)/3000); % Solução exata, obtida analiticamente
F1 = w;
F2 = -R/L * w - k/(C*L);

fprintf("\n\nResultados:");
fprintf("\nUsando o método de Runge-Kutta 4 para sistema de ordem 2, obtemos:");
fprintf("\nPodemos também, comparar com o resultado exato (Última coluna):\n\n")
[X, Y1, Y2] = metRK42(a, b, m, y10, y20, F1, F2, FA);

plotFF(a, b, 200, sym('z'), FA, X, Y1);

function plotFF(a, b, num, x, F, X, Y1)
    hold on;
    u = linspace(a, b, num);
    Fx = vpa(subs(F, x, u));
    plot(u, Fx);
    grid on;
    hold on;
    scatter(X,Y1);
end

function [X, Y1, Y2] = metRK42(a, b, m, y10, y20, F1, F2, FA)
    syms z k w; 
    symZ = sym('z');
    symK = sym('k');
    symW = sym('w');
    vet = [symZ, symK, symW];

    X = zeros(1,m+1);
    Y1 = zeros(1,m+1);
    Y2 = zeros(1,m+1);

    h = (b-a)/m; 
    xt = a; y1t = y10; y2t = y20;
    X(1) = xt; Y1(1) = y1t; Y2(1) = y2t;
    fprintf("i\t x \t\t y1 \t\t y2 \t\t yA");
    fprintf('\n 0 \t %f \t %f \t %f', xt, y1t, y2t);
    for i=1:m
        x = xt; y1 = y1t; y2 = y2t;
        k11 = vpa(subs(F1, vet, [x,y1, y2]));
        k12 = vpa(subs(F2, vet, [x,y1, y2])); 
        x = xt + h/2; y1 = y1t + h/2*k11; y2 = y2t + h/2*k12;
        k21 = vpa(subs(F1, vet, [x,y1, y2]));
        k22 = vpa(subs(F2, vet, [x,y1, y2])); 
        y1 = y1t + h/2*k21; y2 = y2t + h/2*k22;
        k31 = vpa(subs(F1, vet, [x,y1, y2]));
        k32 = vpa(subs(F2, vet, [x,y1, y2])); 
        x = xt + h; y1 = y1t + h*k31; y2 = y2t + h*k32;
        k41 = vpa(subs(F1, vet, [x,y1, y2]));
        k42 = vpa(subs(F2, vet, [x,y1, y2])); 
        xt = a + i*h;
        y1t = y1t + h/6*(k11 + 2*(k21 + k31) + k41);
        y2t = y2t + h/6*(k12 + 2*(k22 + k32) + k42);
        Fa = vpa(subs(FA, vet, [x,y1, y2]));
        fprintf('\n %i \t %f \t %f \t %f \t %f', i, xt, y1t, y2t, Fa);
        X(i+1) = xt; Y1(i+1) = y1t; Y2(i+1) = y2t;
    end
end
