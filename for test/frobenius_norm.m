function norm = frobenius_norm()
    prompt = "Input the number of rows: ";
    m = input(prompt);
    prompt = "Input the number of columns: ";
    n = input(prompt);
    summy=0;
    matrix = zeros(m,n);

    for i = 1:m
        for j = 1:n
            prompt = sprintf("Input element [%d][%d]: ", i, j);
            matrix(i,j) = input(prompt);
        end
    end
    
    disp("Your matrix is: ");
      
    for i = 1:m
        for j = 1:n
            fprintf("%d ", matrix(i,j));
        end
        fprintf("\n");
    end
    
    for i = 1:m
        for j = 1:n
            summy = summy + abs(matrix(i,j))^2;
        end
    end
    
    norm = summy^(1/2);

    %fprintf("\nThe Frobenius-norm of your matrix is %d", norm);

end
