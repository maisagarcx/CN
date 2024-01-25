function [X, Y, erro] = metDormandPrince(a, b, m, y0, F)
    syms k w;
    symK = sym('k');
    symW = sym('w');
    vet = [symK, symW];

    X = zeros(1,m+1);
    Y = zeros(1,m+1);
    erro = zeros(1,m+1);

    a21 = 1/5; a31 = 3/40; a32 = 9/40;
    a41 = 44/45;  a42 = -56/15; a43 = 32/9;
    a51 = 19372/6561; a52 = -25360/2187; a53 = 64448/6561; a54 = -212/729;
    a61 = 9017/3168; a62 = -355/33; a63 = 46732/5247; a64 = 49/176;
    a65 = -5103/18656; a71 = 35/384; a73 = 500/1113; a74 = 125/192;
    a75 = -2187/6784; a76 = 11/84;
    
    c2 = 1/5; c3 = 3/10; c4 = 4/5; c5 = 8/9; c6 = 1; c7 = 1;

    e1 = 71/57600; e3 = -71/16695; e4 = 71/1920; e5 = -17253/339200;
    e6 = 22/525; e7 = -1/40;

    h = (b-a)/m;
    xt = a;  
    yt = y0;
    X(1) = xt;   
    Y(1) = yt;  
    erro(1) = 0;
    fprintf("i \tx \ty \terro");
    fprintf("0 \t%f \t%f", xt, yt);
    for i=1:m
        x = xt; y = yt;
        k1 = h*vpa(subs(F, vet, [x,y]));
        x = xt + c2*h; y = yt + a21*k1;
        k2 = h*(subs(F, vet, [x,y]));     
        x = xt + c3*h; y = yt + a31*k1 + a32*k2;      
        k3 = h*(subs(F, vet, [x,y]));
        x = xt + c4*h; y = yt + a41*k1 + a42*k2 + a43*k3;
        k4 = h*(subs(F, vet, [x,y]));
        x = xt + c5*h; y = yt + a51*k1 + a52*k2 + a53*k3 + a54*k4;       
        k5 = h*(subs(F, vet, [x,y]));    
        x = xt + c6*h; y = yt + a61*k1 + a62*k2 + a63*k3 + a64*k4 + a65*k5;
        k6 = h*(subs(F, vet, [x,y])); 
        x = xt + c7*h; y = yt + a71*k1 + a73*k3 + a74*k4 + a75*k5 + a76*k6;        
        k7 = h*(subs(F, vet, [x,y]));       
        xt = a + i*h;     
        yt = yt + a71*k1 + a73*k3 + a74*k4 + a75*k5 + a76*k6;
        erroGlobal = e1*k1 + e3*k3 + e4*k4 + e5*k5 + e6*k6 + e7*k7;
        X(i + 1) = xt;
        Y(i + 1) = yt; 
        erro(i + 1) = erroGlobal;
        fprintf("%i \t%f \t%f \t%f", i, xt, yt, erro);
    end
end
