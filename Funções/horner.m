function Pa = horner(n, coef, a)
    % n = grau do polinômio
    % coef = coeficientes do polinômio
    % a = ponto de avaliação 
    % P(x) = c(1)x^(n)+ c(2)x^(n-1) + ... + c(n)x + c(n+1)  
    Pa = coef(1);
    for i=2:(n+1)
        Pa = Pa*a + coef(i); 
    end  
end
