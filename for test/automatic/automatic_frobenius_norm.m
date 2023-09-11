function norm = automatic_frobenius_norm(matrix)
    summy=0;
    [m,n] = size(matrix);
    
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
