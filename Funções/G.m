function [L, Ge, Gl, dia] = G(numbe, X, Y, v)
    X = [0.1, 0.6, 0.8];
    Y = [1.221, 3.320, 4.953];
    numbe = 3;
    v = 0.2;
    for i = 1:numbe
        for j = 1:numbe
            if i == j
                Ge(i, j) = v - X(i);
            else
                Ge(i, j) = X(i) - X(j);
            end
        end
    end
    
    % encontrando a diagonal
    dia = prod(diag(Ge));
    
    % multiplicando cada linha
    Gl = ones(1,numbe);
    for i=1:numbe
        for j=1:numbe
            Gl(i)=Gl(i)*Ge(i,j);
        end
    end
    % encontrando o valor
    L=1;
    for i=1:numbe
        L = L + Y(i)/Gl(i);
    end
    L=L*dia;
end
