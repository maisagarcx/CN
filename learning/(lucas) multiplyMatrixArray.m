% begin
l= 3;
c= 3;
% vetor entrada
A= [1, 2, 3; 4, 5, 6; 7, 8, 9];
v= [2, 2, 4];
% vetor saída
r= [0, 0, 0];

for i=1:1:l
    soma=0;
    for j=1:1:c
        soma= soma+(A(i,j)*v(j));
    end
    r(i)=soma;
end

% Saídas
disp(r)
