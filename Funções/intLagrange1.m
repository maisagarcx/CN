function P = intLagrange1(number_of_points, X, Y, value_to_interpolate)
    %X,Y sao vetores com os pontos das abscissas e ordenadas
    interpolated=0;
    for i=1:number_of_points
       %P=Y(i);
       for j=2:number_of_points
           if i~=j
               P = Y(1)+(Y(j)-Y(i))/(X(j)-X(i))*(value_to_interpolate-X(i));
           end
       end
       %interpolated=P;        
    end
end
