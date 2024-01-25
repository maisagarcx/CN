% RESULUCAO DA PROVA DE CALCULO NUMERICO 1 PARTE 2 DE 2020.2
HORAS = [2.0 4.0 6.0 6.4 6.5 7.00];
CARGAS = [0.8 1.6 5.4 4.1 5.9 3.2];
% R2 é a determinacao e SIGMA é a variancia residual

% scatter(HORAS, CARGAS, 'filled');
% xlabel('HORAS');
% ylabel('CARGAS');
% grid on;

% Os dados da Tabela 1 são para a demanda de energia (Y) de um imóvel 
% residencial em função da hora do dia (X). 

% Encontre o polinômio regressor mais adequado ao conjunto
% de pontos da Tabela 1, entre os graus de 2 (manual), 3 e 4 (computador).

R2 = zeros(3,1);
sigma = zeros(3,1);
param = [2; 3; 4];

HORAS1 = [1 1 1 1 1 1; 2.0 4.0 6.0 6.4 6.5 7.0];

[~, R2(1,1), sigma(1,1)] = regressaoLinearMultipla(6, 1, 2, HORAS1', CARGAS);
[~, R2(2,1), sigma(2,1)] = regressaoLinearMultipla(6, 1, 3, HORAS1', CARGAS);
[PARAM1, R2(3,1), sigma(3,1)] = regressaoLinearMultipla(6, 1, 4, HORAS1', CARGAS);
% TEM QUE ESCOLHER OS PARAMETROS?????? Acho que é melhor usar os 4 ultimos.

% Mostre os parâmetros, a matriz das equações normais do modelo mais 
% adequado e o motivo da escolha do polinômio [30%].

C = cat(2,param,R2,sigma);
% disp(C);
% disp(PARAM1);
% Usaremos o grau 4 pois tem o menor SIGMA e R2 esta mais proximo de 1
% EQUACAO DE QUADRADOS MININOS = 19.1908 - 16.6854X1 + 4.4199X2 - 0.3372X3

% Determine os polinômios interpoladores P3(5.44) de Lagrange 

% P³(5.44) precisamos dos 4 pontos mais pertos de 5.44
% usarei os pontos 2, 3, 4 e 5
HORAL = [4.0 6.0 6.4 6.5]; % X
CARGAL = [1.6 5.4 4.1 5.9]; % Y
INTLAGRANGE = intLagrange2(4, HORAL, CARGAL, 5.44);
% INTLAGRANGE = 19.8913

% e P²(3.67) Gregory-Newton para os pontos da Tabela 1 
% (manual e computador) [30%].

% P²(3.67) precisamos dos 3 pontos mais pertos de 3.67
% usarei os pontos 1, 2 e 3
HORAGN = [2.0 4.0 6.0];
CARGAGN = [0.8 1.6 5.4];
INTGREGNEWTON = intGregoryNewton(3, HORAGN, CARGAGN, 3.67);
% INTGREGNEWTON = 1.2613

% Determine os parâmetros do modelo y = ae^(bx) (computador) para os 
% pontos acima [40%].
% Y = ae^(bx) => ln(y) = ln(a) + bx
% isso é a coisa mais errada que eu ja fiz na vida
[PARAM2, R22, sigma22] = regressaoLinearMultipla(6, 1, 3, HORAS1', CARGAS);
% PARAM2 =  [-2.7414 1.8022 -0.1065]
% R22 = 0.6574
% sigma22 = 2.3664

function [b, R2, sigma2] = regressaoLinearMultipla(n, v, p, x, y)
    
    
    if (v>1)&&((v+1)~=p)
        prompt = "Modelo Invalido";
        error(prompt);
        return;
    end
    %for i=1:n % inclusao de uma coluna de 1s relativa a b⁰
        %for j=v+1:-1:2
            %x(i,j)=x(i,j-1);
        %end
        %x(i,1)=1;
    %end
    % disp(x);
    if (v==1)&&(p>2) % ou seja, se regressao polinomial, entao gera potenciais de x
        for j=2:p-1
            for i=1:n
                x(i,j+1)=x(i,2)^j; % ^j
            end
        end
    end
    % U = zeros(p, p);
    for i=1:p % equacoes normais
        for j=1:p
            summy=0;
            for k=1:n
                summy=summy+x(k,i)*x(k,j);
            end
            sxx(i,j)=summy; % matriz dos coeficientes
        end
        summy=0;
        for k=1:n
            summy=summy+x(k,i)*y(k);
        end
        sxy(i)=summy; % vetor dos termos independentes
    end
    L = cholesky(p, sxx); % p = ordem, sxx = matriz
    t = suc_subst(p, L, sxy); % L = lower
    for i=1:p
        for j=1:i
            U(j,i)=L(i,j); % U = L transposta
        end
    end
    b = ret_subst(U, t);
    D=0;
    sy2 = 0;
    for i=1:n
        summy=0;
        for j=1:p
            summy=summy+b(j)*x(i,j);
        end
        u(i)=summy;
        d(i)=y(i)-u(i);
        D=D+d(i)^2;
        sy2=sy2+y(i)^2;
    end
    R2=1-D/(sy2-sxy(1)^2/n); % coeficiente de determinacao
    sigma2= D/(n-p); % variancia residual
end

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
