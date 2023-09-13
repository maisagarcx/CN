function sumax_rows = sumax_rows(matrix)
    [m,n] = size(matrix);
    vet_summy = zeros(1,m);
    
    for i=1:m
        summy=0;
        for j=1:n
            summy = summy + abs(matrix(i,j));
        end
        vet_summy(i)=summy;
    end
    
    sumax_rows = vet_summy(1,1);
     for i = 1:m
         if vet_summy(m)>=sumax_rows
             sumax_rows=vet_summy(m);
         end
    end
end
