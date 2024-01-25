function sumax_cols = sumax_cols(matrix)
    [m,n] = size(matrix);
    vet_summy = zeros(1,n);
    
    for i=1:n
        summy=0;
        for j=1:m
            summy = summy + abs(matrix(i,j));
        end
        vet_summy(i)=summy;
    end
    
    sumax_cols = vet_summy(1,1);
     for i = 1:n
         if vet_summy(n)>=sumax_cols
             sumax_cols=vet_summy(n);
         end
    end
end
