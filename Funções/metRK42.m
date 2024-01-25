function [X, Y1, Y2] = metRK42(a, b, m, y10, y20, F1, F2)
    syms z; syms k; syms w; 
    symZ = sym('z');
    symK = sym('k');
    symW = sym('w');
    vet = [symZ, symK, symW];

    a = 0; b = 1; m = 10; y10 = 4; y20 = -20;
    F1 = w; % w = v
    F2 = - W*R/L;

    X = zeros(1,m+1);
    Y1 = zeros(1,m+1);
    Y2 = zeros(1,m+1);

    h = (b-a)/m; 
    xt = a; y1t = y10; y2t = y20;
    X(1) = xt; Y1(1) = y1t; Y2(1) = y2t;
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
        fprintf('\n %i \t %f \t %f \t %f', i, xt, y1t, y2t);
        X(i+1) = xt; Y1(i+1) = y1t; Y2(i+1) = y2t;
    end
end
