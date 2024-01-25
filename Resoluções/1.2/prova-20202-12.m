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

% RESULUCAO DA PROVA DE CALCULO NUMERICO 1 PARTE 2 DE 2022.2

CARGAS = [1004 1026 989 971 993 1069 1270 1474 1554 1565];
PERDAS = [36 34 36 35 33 31 27 32 35 42];

% scatter(CARGAS, PERDAS, 'filled');
% xlabel('Cargas');
% ylabel('Perdas');
% grid on;

% Um Engenheiro Eletricista precisa determinar um modelo para estimar as perdas de 
% potência total (y) em função da carga total (x) de um sistema de potência.
% Para isso, o Engenheiro dispõe dos dados da Tabela 1

% Determine os parâmetros do melhor modelo polinomial (computador) para estimar 
% a perda em função da carga. Em seguida, mostre como encontrar a matriz das 
% equações normais (manual) do melhor modelo [40%].
CARGAS1 = [1 1 1 1 1 1 1 1 1 1; 1004 1026 989 971 993 1069 1270 1474 1554 1565];

R2 = zeros(6,1);
sigma = zeros(6,1);
param = [1; 2; 3; 4; 5; 6];

for i=1:6
    [~, R2(i,1), sigma(i,1)] = regressaoLinearMultipla(10, 1, i, CARGAS1', PERDAS);
end

% [~, R12, sigma12] = regressaoLinearMultipla(10, 1, 1, CARGAS1', PERDAS);
% [~, R22, sigma22] = regressaoLinearMultipla(10, 1, 2, CARGAS1', PERDAS);
% [~, R32, sigma32] = regressaoLinearMultipla(10, 1, 3, CARGAS1', PERDAS);
% [~, R42, sigma42] = regressaoLinearMultipla(10, 1, 4, CARGAS1', PERDAS);
% [~, R52, sigma52] = regressaoLinearMultipla(10, 1, 5, CARGAS1', PERDAS);
% [~, R62, sigma62] = regressaoLinearMultipla(10, 1, 6, CARGAS1', PERDAS);

% construcao da tabela
% MATRIZR2 = [-2.6645e-15; 0.0421; 0.7796; 0.8115; 0.8123; 0.8719];
% SIGMA2 = [15.2111; 16.3925; 4.3107; 4.3001; 5.1384; 4.3854];
C = cat(2,param,R2,sigma); % concatena horizontalmente, use 1 para concatenar verticalmente
% grau mais adequado: 3; pois tem o SIGMA baixo e o R2 tambem

[b, ~, ~] = regressaoLinearMultipla(10, 1, 3, CARGAS1', PERDAS);
% b = [215.3286 -0.3011 0.0001];
% estimat1 = b1(1) + x*b1(2) + x^2*b1(3);

% Determine os parâmetros do modelo linear múltiplo para determinação da perda
% (computador). Estime a perda para 3.7 h e carga de 900 MW. [30%]
X0 = [1 1 1 1 1 1 1 1 1 1; 1 2 3 4 5 6 7 8 9 10; 1004 1026 989 971 993 1069 1270 1474 1554 1554];
X1 = X0';
[b2, R2, sigma2] = regressaoLinearMultipla(10, 2, 3, X1, PERDAS);
estimat2 = b2(1) + 3.7*b2(2) + 900*b2(3);
% b = [25.6470 -0.6449 0.0101]
% R2 = 0.0909
% sigma2 = 17.7789
% estimat2 = 32.3334

% Determine os polinômios interpoladores horários para 𝑃𝑒𝑟𝑑𝑎3(3,7) de Lagrange
% e 𝑃𝑒𝑟𝑑𝑎4(7,3) de Gregory-Newton, para os pontos da Tabela 1 [30%].

% P³(3.7) precisamos dos 4 pontos mais pertos de 3.7
% usarei os pontos (horas) 2, 3, 4 e 5
CARGASL = [2 3 4 5]; % X
PERDASL = [34 36 35 33]; % Y 
INTLAGRANGE = intLagrange2(4, CARGASL, PERDASL, 3.7);
% INTLAGRANGE = 35.4960

% P⁴(7.3) precisaremos dos 5 pontos mais pertos de 7.3
% usarei os pontos (horas) 5, 6, 7, 8, 9
CARGASGN = [5 6 7 8 9]; % X
PERDASGN = [33 31 27 32 35]; % Y
INTGREGNEWTON = intGregoryNewton(5, CARGASGN, PERDASGN, 7.3);
% INTGREGNEWTON = 27.6301

function [L, det, error] = cholesky(order, matrix)
    %order = 3;
    %matrix = [4 -2 2; -2 10 -7; 2 -7 30];
    %does the decomposition LLᵗ of a matrix
    det=1;
    L=zeros(order,order);

    for j=1:order
        summy=0;
        for k = 1:j-1
            summy=summy + L(j,k)^2;
        end
        t=matrix(j,j)-summy;
        det=det*t;
        error = t<=0;
        if error==1
            return;
            %prompt = "A matriz não é definida positiva";
            %error(prompt);
        else
            L(j,j)=sqrt(t);
            r=1/L(j,j);
        end
        for i=j+1:order
            summy=0;
            for k=1:j-1
                summy=summy+L(i,k)*L(j,k);
            end
            L(i,j)=(matrix(i,j)-summy)*r;
        end
    end
end

function Y = suc_subst(order, lower_tri_matrix, vet_ans)
    %solves lower triangular matrix systems using LX=C
    %LX=C, is lower_tri_matrix*X=vet_ans
    if order<1 
         prompt = "The order must be equal or higher then 1, try again."; 
         error(prompt);
    end
    
    Y=zeros(1,order); %just of speed
    
    Y(1) = vet_ans(1)/lower_tri_matrix(1,1);
    for i=2:order
        summy=0;
        for j=1:i-1
            summy=summy+lower_tri_matrix(i,j)*Y(j);
        end
        Y(i)=(vet_ans(i)-summy)/lower_tri_matrix(i,i);
    end
end

function X = ret_subst(upper_tri_matrix, vet_ans)
    %aceita vetor linha ou coluna como vet_ans
    %solves upper triangular matrix systems using UX=D
    %UX=D, is upper_tri_matrix*X=vet_ans
    
    [order,~] = size(upper_tri_matrix);
    %X=zeros(order,1); %para devolver vetor coluna
    X=zeros(1,order); %para devolver vetor linha
    X(order) = vet_ans(order)/upper_tri_matrix(order,order);
    for i=order-1:-1:1
        summy=0;
        for j=i+1:order
            summy=summy+upper_tri_matrix(i,j)*X(j);
        end
        X(i)=(vet_ans(i)-summy)/upper_tri_matrix(i,i);
    end
end

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
