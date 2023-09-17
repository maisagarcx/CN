function [fator, det, error] = cholesky(order, matrix)
    %does the decomposition LLᵗ of a matrix
    det=1;
    fator=zeros(1,order);

    for j=1:n
        summy=0;
        for k=1:j-1
            summy=summy+fator(j,k)^2;
        end
        t=matrix(j,j)-summy;
        det=det*t;
        error=t<=0;
        if error
            prompt = "A matriz não é definida positiva.";
            error(prompt);
        else
            fator(j,j)=sqrt(t);
            r=1/fator(j,j);
        end
        for i=j+1:order
            summy=0;
            for k=1:j-1
                summy=summy+fator(i,k)*fator(j,k);
            end
            fator(i,j)=(matrix(i,j)-summy)*r;
        end
    end
end