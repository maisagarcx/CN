function interpolated = intGregoryNewton(number_of_points, X, Y, value_to_interpolate)
    for i=1:number_of_points
        Dely(i)=Y(i);
    end
    % construindo as diferencas finitas
    for k=1:number_of_points-1
        for i=number_of_points:-1:k+1
            Dely(i)=Dely(i)-Dely(i-1);
        end
    end
    % avaliacao pelo processo de Horner
    u=(value_to_interpolate-X(1))/(X(2)-X(1));
    interpolated=Dely(number_of_points);
    for i=number_of_points-1:-1:1
        interpolated=interpolated*(u-i+1)/i+Dely(i);
    end
end
