function [L, det, error] = cholesky(order, matrix)
    %order = 3;
    %matrix = [4 -2 2; -2 10 -7; 2 -7 30];
    %does the decomposition LLᵗ of a matrix
    det=1;
    L=zeros(order,order);

    for j=1:order
        summy=0;
        for k = 1:j-1
            summy=summy + L(j,k)^2;
        end
        t=matrix(j,j)-summy;
        det=det*t;
        error = t<=0;
        if error==1
            return;
            %prompt = "A matriz não é definida positiva";
            %error(prompt);
        else
            L(j,j)=sqrt(t);
            r=1/L(j,j);
        end
        for i=j+1:order
            summy=0;
            for k=1:j-1
                summy=summy+L(i,k)*L(j,k);
            end
            L(i,j)=(matrix(i,j)-summy)*r;
        end
    end
end
