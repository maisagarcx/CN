% norm-p

prompt = "Input the size of array: ";
n = input(prompt);
vet = zeros(1, n);

prompt = "Input the norm P: ";
p = input(prompt);

for i=1:1:n
    prompt = sprintf("Input element %d: ", i);
    vet(i) = input(prompt);
end

summy=0;
for i=1:1:n
    summy = summy + abs(vet(i)).^p;
end

pnorm = summy .^ (1/p);

fprintf('The norm-p of your array is %f', pnorm);

