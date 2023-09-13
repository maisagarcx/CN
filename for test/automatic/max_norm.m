function max_norm(vet) 
     n = length(vet); 
     
     %fprintf("The size of your array is %d \n", n); 
  
     max_norm = vet(1);
     for i = 1:n 
         if vet(i)>=max_norm
             max_norm=vet(i)
         end
     end 
  
     %fprintf('The max norm of your array is %f', norm); 
 end
