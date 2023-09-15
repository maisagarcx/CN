function X = suc_subst(order, lower_tri_matrix, vet_ans)
    %solves lower triangular matrix systems using LX=C
    %LX=C, is lower_tri_matrix*X=vet_ans
    if order<1 
         prompt = "The order must be equal or higher then 1, try again."; 
         error(prompt);
    end
    
    X=zeros(1,order); %just of speed
    
    X(1) = vet_ans(1)/lower_tri_matrix(1,1);
    for i=2:order
        summy=0;
        for j=1:i-1
            summy=summy+lower_tri_matrix(i,j)*X(j);
        end
        X(i)=(vet_ans(i)-summy)/lower_tri_matrix(i,i);
    end
end
