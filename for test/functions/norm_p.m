function norm_p = norm_p(vet, p) 
     if p<1 
         prompt = "The norm-p must be equal or higher then 1, try again."; 
         error(prompt);
     end 
  
     n = length(vet); 
     summy = 0; 

     for i=1:n 
         summy = summy + abs(vet(i)).^p; 
     end 
  
     norm_p = summy^(1/p);  
 end
