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
