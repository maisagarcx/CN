function X = ret_subst(order, upper_tri_matrix, vet_ans)
    %solves upper triangular matrix systems using UX=D
    %UX=D, is upper_tri_matrix*X=vet_ans
    if order<1 
         prompt = "The order must be equal or higher then 1, try again."; 
         error(prompt);
    end
    
    X=zeros(1,order); %just of speed
    
    X(order) = vet_ans(order)/upper_tri_matrix(order,order);
    for i=order-1:-1:1
        summy=0;
        for j=i+1:n
            summy=summy+upper_tri_matrix(i,j)*X(j);
        end
        X(i)=(vet_ans(i)-summy)/upper_tri_matrix(i,i);
    end
end
