% TRABALHO 3.3 DE CÁLCULO NUMÉRICO 2023.2
% ALUNA: MAÍSA GARCIA NEPOMUCENO CORRÊA

fprintf("Implementação referente ao Trabalho 3.3 de Cálculo Numérico 2023.2.\n\n");

syms k w;
R = 12; L = 4; E = 60;

fprintf("Suponha que em um circuito simples, a resistência seja %d ohms e\n", R);
fprintf("a indutância seja %d H. Uma pilha fornece uma voltagem constante\n", L);
fprintf("de %d V e o interruptor é fechado quanto t = 0, \n", E);
fprintf("a corrente começa com I(0) = 0. Encontre a corrente depois de 1s.");

a = 0; b = 1; m = 50; y0 = 0;
F = (E-R*w)/L; % EDO de primeira ordem
FA = -(15*exp(-3*k))/3 + 5; % Solução exata, obtida analiticamente

fprintf("\n\nResultados:");
fprintf("\nUsando o método de Euler e comparando com a resposta exata (4º coluna), obtemos:\n\n");
[X, Y] = metEuler(a, b, 100, y0, F, FA);

plotF(a, b, 200, sym('k'), FA, X, Y);

function plotF(a, b, num, x, F, X, Y)
    hold on;
    u = linspace(a, b, num);
    Fx = vpa(subs(F, x, u));
    plot(u, Fx);
    grid on;
    hold on;
    scatter(X,Y);
end

function [X, Y] = metEuler(a, b, m, y0, F, FA)
    syms k w;
    symK = sym('k');
    symW = sym('w');
    vet = [symK, symW];
    
    X = zeros(1,m+1);
    Y = zeros(1,m+1);

    h = (b-a)/m;
    x = a; 
    y = y0;
    Fa = vpa(subs(FA, vet, [x,y]));
    Fxy = vpa(subs(F, vet, [x,y]));
    fprintf("i\t x \t\t y \t\t FA(x,y) \t F(x,y)\n")
    fprintf('0\t %f\t %f \t %f \t %f\t', x, y, Fa, Fxy);
    fprintf("\n");
    X(1) = x;
    Y(1) = y;

    for i=1:m
        x = a + i*h;
        Fa = vpa(subs(FA, vet, [x,y]));
        y = y + h*Fxy;
        Fxy = vpa(subs(F, vet, [x,y]));
        fprintf('%d\t %f \t %f \t %f \t %f\n', i, x, y, Fa, Fxy);
        X(i+1) = x;
        Y(i+1) = y;
    end
end
