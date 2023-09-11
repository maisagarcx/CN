function summy = automatic_norm_1(vet)
    % automatic norm-1
    n = length(vet);
    summy = 0;

    fprintf("The size of your array is %d \n", n);
    disp("And your array is: ");

    for i = 1:n
        fprintf("%d ", vet(i));
    end

    for i = 1:n
        summy = summy + abs(vet(i));
    end

    fprintf("\nThe norm-1 of your array is %d", summy);

end
