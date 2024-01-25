function [X, Y, erro] = metAdams4(a, b, m, y0, F)
    syms k; syms w;
    symK = sym('k');
    symW = sym('w');
    vet = [symK, symW];

    % a = 0; b = 1; m = 10; y0 = 0; 
    % F = 15 - 3*k;
    
    h = (b - a)/m;
    [X, Y, erro] = metDormandPrince(a, (a + 3*h), 3, y0, F);

    for i=1:4
        fprintf('\n %d \t %f \t %f \t %f', i-1, X(i), Y(i), erro(i));    
    end

    for i=4:m
        x = X(i-3); y = Y(i-3);
        f0 = vpa(subs(F, vet, [x,y]));   
        x = X(i-2); y = Y(i-2); 
        f1 = vpa(subs(F, vet, [x,y]));
        x = X(i-1); y = Y(i-1); 
        f2 = vpa(subs(F, vet, [x,y]));  
        x = X(i); y = Y(i); 
        f3 = vpa(subs(F, vet, [x,y]));
        Ypre = h*(55*f3 - 59*f2 + 37*f1 - 9*f0)/24 + Y(i);
        Y(i+1) = Ypre; X(i+1) = a + i*h; x = X(i+1);
        
        for j=1:2     
            y = Y(i+1); 
            f4 = vpa(subs(F, vet, [x,y]));
            Ycor = h*(9*f4 + 19*f3 - 5*f2 + f1)/24 + Y(i);
            Y(i+1) = Ycor;
        end

        erro = abs(Ycor - Ypre)*19/270;  
        fprintf('\n %d \t %f \t %f \t %f', i, X(i+1), Y(i+1), erro); 
    end
end
