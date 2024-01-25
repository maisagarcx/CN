function [X, Y] = metEuler(a, b, m, y0, F)
    syms k w;
    symK = sym('k');
    symW = sym('w');
    vet = [symK, symW];
    
    X = zeros(1,m+1);
    Y = zeros(1,m+1);

    h = (b-a)/m;
    x = a; 
    y = y0;
    Fxy = vpa(subs(F, vet, [x,y]));
    fprintf("i\t x\t \t y \t\t F(x,y)\n")
    fprintf('0\t %f\t %f \t %f \t', x, y, Fxy);
    fprintf("\n");
    X(1) = x;
    Y(1) = y;

    for i=1:m
        x = a + i*h;
        y = y + h*Fxy;
        Fxy = vpa(subs(F, vet, [x,y]));
        % fprintf("\n");
        fprintf('%d\t %f \t %f \t %f\n', i, x, y, Fxy);
        X(i+1) = x;
        Y(i+1) = y;
    end
end
