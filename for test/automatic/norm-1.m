% automatic norm-1

n = randi([2 10]);
vet = randi([-10 10], 1, n);
summy=0;

fprintf("The size of your array is %d \n", n);
sprintf("And your array is: ", newline);

for(i=1:1:n)
    fprintf("%d ", vet(i));
end

for(i=1:1:n)
    summy=summy+abs(vet(i));
end

fprintf("\nThe norm-1 of your array is %d", summy);
