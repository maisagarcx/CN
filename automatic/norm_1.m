function summy = norm_1(vet)
    n = length(vet); % gets the size of the array
    summy = 0;

    fprintf("The size of your array is %d \n", n); % use fprintf if you wanna display variables
    disp("And your array is: "); % use disp if you gonna display just a string

    for i = 1:n % the same as i=1:1:1
        fprintf("%d ", vet(i));
    end

    for i = 1:n
        summy = summy + abs(vet(i)); % abs gives the absolute value of each element
    end

    fprintf("\nThe norm-1 of your array is %d", summy);

end
