function interpolated = intLagrange2(number_of_points, X, Y, value_to_interpolate)
    %X, Y = arrays of points i have
    interpolated=0;
    for i=1:number_of_points
       numerator=1;
       denominator=1;
       for j=1:number_of_points
           if i~=j
               numerator = numerator*(value_to_interpolate-X(j));
               denominator = denominator*(X(i)-X(j));
           end
       end
       interpolated=interpolated+Y(i)*(numerator/denominator);    
    end
end
