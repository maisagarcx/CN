function X = ret_subst(upper_tri_matrix, vet_ans)
    %aceita vetor linha ou coluna como vet_ans
    %solves upper triangular matrix systems using UX=D
    %UX=D, is upper_tri_matrix*X=vet_ans
    [order,~] = size(upper_tri_matrix);
    %X=zeros(order,1); %para devolver vetor coluna
    X=zeros(1,order); %para devolver vetor linha
    X(order) = vet_ans(order)/upper_tri_matrix(order,order);
    for i=order-1:-1:1
        summy=0;
        for j=i+1:order
            summy=summy+upper_tri_matrix(i,j)*X(j);
        end
        X(i)=(vet_ans(i)-summy)/upper_tri_matrix(i,i);
    end
end
