function frobenius_norm = frobenius_norm(A)
    summy=0;
    [m,n]=size(A);

    for i = 1:m
        for j = 1:n
            summy=summy+abs(A(i,j))^2;
        end
    end
    
    frobenius_norm=summy^(1/2);
end
