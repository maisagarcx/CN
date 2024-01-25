function X = suc_subst_piv(lower_tri_matrix_uni, ind_vet, pivot)
    %pode receber vetor linha ou coluna no ind_vet
    %o P usado deve ser um com indices das linhas trocados e nao as linhas da matriz identidade    
    % solves lower triangular matrix systems using LX=Pc    
    % lower_tri_matrix*X=permutation_matrix*ind_vet
    
    [order,~] = size(lower_tri_matrix_uni);
    %X=zeros(order,1); %para devolver vetor coluna
    X=zeros(1,order); %para devolver vetor linha
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
