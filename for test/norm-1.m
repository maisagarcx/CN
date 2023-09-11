% norm-1

prompt = "Input the size of array: ";
n = input(prompt);
vet = zeros(1, n);

for i=1:1:n
    prompt = sprintf("Input element %d: ", i);
    vet(i) = input(prompt);
end

summy=0;
for i=1:1:n
    summy = summy + abs(vet(i));
end

fprintf('The norm-1 of your array is %f', summy);
