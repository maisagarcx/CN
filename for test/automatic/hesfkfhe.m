function [X, iter, error] = jacobi(order, matrix, ind_vet, toler, max_iter)
    %construction of matrices for iterations
    for i=1:order
        r=1/matrix(i,i);
        for j=1:order
            if i~=j
                matrix(i,j)=matrix(i,j)*r;
            end
        end
        ind_vet(i)=ind_vet(i)*r;
        X(i)=ind_vet(i);
    end
    iter=0;
    %Jacobi iterations
    while ~DifMax<toler||iter>=max_iter
        iter=iter+1;
        for i=1:order
            summy=0;
            for j=1:order
                if i~=j
                    summy=summy+matrix(i,j)*X(j);
                end
            end
            v(i)=ind_vet(i)-summy;
        end
        norm_1=0;
        norm_2=0;
        for 
        
    end
end
