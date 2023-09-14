function [X, iter, error] = jacobi(order, matrix, ind_vet, toler, max_iter) 
     %construction of matrices for iterations 
     
     X=zeros(1,order);
     v=zeros(1,order);
     
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
     while ~difMax<toler|iter>=max_iter 
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
         fprintf("O número de iterações é %d\nA diferença máxima é %d", iter, difMax);
         disp(X);
     end 
     error=difMax>=toler;
 end
