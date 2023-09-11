function summy = norm_1()
    prompt = "Please input the size of your array: ";
    n = input(prompt);
    vet = zeros(1, n);
    
    for i = 1:n
        prompt = sprintf("Input element %d: ", i);
        vet(i) = input(prompt);
    end
    
    summy=0;
    for i = 1:n
        summy = summy + abs(vet(i));
    end

    %fprintf('The norm-1 of your array is %f\n', summy);
end

    
    
