function [X, iter, error] = gauss_seidel(order, matrix, ind_vet, toler, iterMax)
    %solve systems matrix*X=ind_vet using Gauss-Siedel
    
    X=zeros(1,order);
    v=zeros(1,order);
    
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
    
    %Gauss-Seidel iterations
    difMax=toler+1;
    while ~(difMax<toler||iter>=iterMax)
        iter=iter+1;
        for i=1:order
            summy=0;
            for j=1:order
                if i~=j
                    summy=summy+matrix(i,j)*X(j);
                end
            end
            v(i)=X(i);
            X(i)=ind_vet(i)-summy;
        end
        norm_1=0;
        norm_2=0;
        for i=1:order
            if abs(X(i)-v(i))>norm_1
                norm_1=abs(X(i)-v(i));
            end
            if abs(X(i))>norm_2
                norm_2=abs(X(i));
            end
        end
        difMax=norm_1/norm_2;
        fprintf("O número de iterações é %d\nA diferença máxima é %d", iter, difMax);
        disp(X);
        %teste de convergencia
    end
    error=difMax>=toler;
end
