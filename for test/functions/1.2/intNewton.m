function interpolated = intNewton(number_of_points, X, Y, value_to_interpolate)
    for i=1:number_of_points
        Dely(i)=Y(i);
    end
    % construindo as diferencas divididas
    for k=1:number_of_points-1
        for i=number_of_points:-1:k+1
            Dely(i)=(Dely(i)-Dely(i-1))/(X(i)-X(i-k));
        end
    end
    % avaliacao pelo processo de Horner
    interpolated=Dely(number_of_points);
    for i=number_of_points-1:-1:1
        interpolated=interpolated*(value_to_interpolate-X(i))+Dely(i);
    end
end
