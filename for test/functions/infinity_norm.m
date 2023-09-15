function infinity_norm = infinity_norm(A) 
     [m,n] = size(A);

     infinity_norm = 0;
     for i=1:m
         summy=0;
         for j=1:n
            summy=summy+abs(A(i,j));
         end
         if summy>=infinity_norm
                infinity_norm=summy;
         end
     end  
end
