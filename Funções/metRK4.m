function [X, Y] = metRK4(a, b, m, y0, F)
    syms k; syms w;
    symK = sym('k');
    symW = sym('w');
    vet = [symK, symW];

    a = 0; b = 6e-3; m = 100; y0 = 0; 
    V = 10; R = 1; L = 0.001;
    F = 1/L * (V - R*L*k);

    X = zeros(1,m+1);
    Y = zeros(1,m+1);

    h = (b-a)/m;
    xt = a;
    yt = y0;
    X(1) = xt;
    Y(1) = yt;

    % fprintf(' i \t\t x \t\t\t y\n');
    % fprintf(' 0 \t %f \t %f\n', xt, yt);

    for i=1:m
        x = xt;
        y = yt;
        k1 = vpa(subs(F, vet, [x,y]));
        x = xt + h/2;
        y = yt + h/2*k1;
        k2 = vpa(subs(F, vet, [x,y]));
        y = yt + h/2*k2;
        k3 = vpa(subs(F, vet, [x,y]));
        x = xt + h;
        y = yt + h*k3;
        k4 = vpa(subs(F, vet, [x,y]));
        xt = a + i*h;
        yt = yt + h/6*(k1 + 2*(k2+k3) + k4);
        % fprintf(' %d   %f \t %f\n', i, xt, yt);
        X(i+1) = xt;
        Y(i+1) = yt;
    end
end

