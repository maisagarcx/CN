function norm = norm_p() 
    prompt = "Input the size of array: ";
    n = input(prompt);
    vet = zeros(1, n);
    prompt = "Please input the norm p: ";
    
    p = -1;
    while p<1
        p = input(prompt);
        if p<1
            disp("The norm p must be equal or higher then 1");
        end
    end      
    
    for i = 1:n
        prompt = sprintf("Input element %d: ", i);
        vet(i) = input(prompt);
    end

    summy=0;
    for i = 1:n
        summy = summy + abs(vet(i)).^p;
    end

    norm = summy.^(1/p);

    %fprintf('The norm-p of your array is %f', norm);
end
