function norm = automatic_norm_2(vet)
    n = length(vet);
    summy = 0;

    fprintf("The size of your array is %d \n", n);
    disp("And your array is: ");

    for i = 1:n % the same as i=1:1:1
        fprintf("%d ", vet(i));
    end

    for i = 1:n
        summy = summy + abs(vet(i)).^2;
    end
    
     norm = summy.^(1/2);

    %fprintf("\nThe norm-2 of your array is %d", norm);
end
