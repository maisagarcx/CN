function norm = automatic_norm_p(vet, p)
    if p<1
        disp("The norm p must be equal or higher then 1");
        return
    end
    
    n = length(vet);
    summy = 0;

    fprintf("The size of your array is %d \n", n);
    disp("And your array is: ");

    for i = 1:n
        fprintf("%d ", vet(i));
    end

    for i = 1:n
        summy = summy + abs(vet(i)).^p;
    end
    
    norm = summy.^(1/p);

    %fprintf('The norm-p of your array is %f', norm);
end
