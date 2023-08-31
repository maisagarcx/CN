% begin
i=0;
soma=0;
soma2=0;

% vetor entrada
n=5;
vet= [10 8 7 9 6];

for i=1:1:n
    soma=soma+vet(i);
    soma2=soma2+vet(i)^(2);
end

med=soma/n;
dev=sqrt(((soma2-soma^(2)/n)/(n-1)));

% Saídas
disp(med)
disp(dev)
