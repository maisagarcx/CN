function [X, iter, error] = jacobi(order, matrix, ind_vet, toler, max_iter) 
     %solve systems matrix*X=ind_vet using Jacobi
     
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
     
     %Jacobi iterations 
     while true 
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
         for  i=1:order
             if abs(v(i)-X(i))>norm_1
                 norm_1=abs(v(i)-X(i));
             end
             if abs(v(i))>norm_2
                 norm_2=abs(v(i));
             end
             X(i)=v(i);
         end
         difMax=norm_1/norm_2;
         fprintf('Iteração: %d\n', iter);
         fprintf('x: %s\n', mat2str(X));
         fprintf('DifMax: %f\n', difMax);
         disp(X);
         % if difMax<toler || iter>=max_iter
         %     break;
         % end
         if (difMax <= toler) || (iter >= max_iter)
            break;  % interrompa
         end
     end
 error=difMax>=toler;
end
