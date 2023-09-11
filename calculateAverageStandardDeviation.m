% calcule a média e desvio

% vetor entrada

prompt = "Input size of array: ";
n = input(prompt);
vet = zeros(1, n);

prompt = "Input element: ";
for i=1:1:n
    vet(i) = input(prompt)
end

X = ['Your array is: ', newline, num2str(vet)];
disp(X);

soma = sum(vet);
soma2 = sum(vet.^2);

med=soma/n;
dev=sqrt(((soma2-soma^(2)/n)/(n-1)));

Y = ['Your avarage is: ', newline, num2str(med), newline, 'Your deviation is: ', newline, num2str(dev)];
disp(Y);


