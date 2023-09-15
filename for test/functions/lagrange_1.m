function interpolated = lagrange_1(number_of_points, X, Y, value_to_interpolate)
    %X, Y = arrays of points i have
    interpolated=0;
    for i=1:number_of_points
       P=Y(i);
       for j=1:number_of_points
           if i~=j
               P = P*((value_to_interpolate-X(j))/X(i)-X(j));
           end
       end
       interpolated=interpolated+P;        
    end
end
