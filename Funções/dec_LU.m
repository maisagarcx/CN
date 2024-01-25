function [matrix, det, pivot] = dec_LU(order, matrix)
    pivot = zeros(order,1);
    for i=1:order
        pivot(i)=i;
    end
    det=1;
    for j=1:order-1
        %choice of pivot element
        p=j;
        A_max=abs(matrix(j,j));
        for k=j+1:order
            if abs(matrix(k,j))>A_max
                A_max=abs(matrix(k,j));
                p=k;
            end
        end
        if p~=j
            %change rows
            for k=1:order
                t=matrix(j,k);
                matrix(j,k)=matrix(p,k);
                matrix(p,k)=t;
            end
            t=pivot(j);
            pivot(j)=pivot(p);
            pivot(p)=t;
            det=-det;
        end
        det=det*matrix(j,j);
        if (abs(matrix(j,j))~=0)
            %Gauss elimination
            r=1/matrix(j,j);
            for i=j+1:order
                m=matrix(i,j)*r;
                matrix(i,j)=m;
                for k=j+1:order
                    matrix(i,k)=matrix(i,k)-m*matrix(j,k);
                end
            end
        end
    end
    det=det*matrix(order,order);
end
