% begin
l=3;
c=3;
% vetor entrada
A=[1, 2, 3; 4, 5, 6; 7, 8, 9];
% vetor saída
bigger=[0 ,0, 0];

for i=1:1:l
    bigger(i)=A(i, 1);
    for j=2:1:c
        if A(i,j)>bigger(i)
            bigger(i)=A(i,j);
        end
    end
    % Saídas
    disp(i)
    disp(bigger(i))
end

