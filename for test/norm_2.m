function norm = norm_2()
    prompt = "Input the size of array: ";
    n = input(prompt);
    vet = zeros(1, n);

    for i = 1:n
        prompt = sprintf("Input element %d: ", i);
        vet(i) = input(prompt);
    end

    summy=0;
    for i = 1:n
        summy = summy + abs(vet(i)).^2;
    end

    norm = summy .^ (1/2);

    %fprintf('The norm-2 of your array is %f', norm);
end
    
