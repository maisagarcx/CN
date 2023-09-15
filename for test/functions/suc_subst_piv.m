function X = suc_subst_piv(order, lower_tri_matrix_uni, ind_vet, pivot)
    %solves lower triangular matrix systems using LX=Pc
    %lower_tri_matrix*X=permutation_matrix*ind_vet
    if order<1 
         prompt = "The order must be equal or higher then 1, try again."; 
         error(prompt);
    end
    
    X=zeros(1,order);
    
    k=pivot(1);
    X(1)=ind_vet(k);
    for i=2:order
        summy=0;
        for j=1:i-1
            summy=summy+lower_tri_matrix_uni(i,j)*X(j);
        end
        k=pivot(i);
        X(i)=ind_vet(k)-summy;
    end
end
