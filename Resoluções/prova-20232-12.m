% PARTE 2 DA PROVA 1 DE CALCULO NUMERICO
% ALUNA: MAISA GARCIA NEPOMUCENO CORREA

t = [9.0 10.5 12.8 14.0 16.3 18.5 21.0 23.0 25.0 27.0]; % X
Rt = [1.00 0.88 0.72 0.53 0.47 0.35 0.26 0.21 0.16 0.07]; % Y

% PRIMEIRA QUESTAO
fprintf("Primeira Questão:\n\n");

fprintf("Letra A:\n");
% LETRA A: Determine os parâmetros do modelo linear, quadrático e 
% exponencial da Equação 1
t1 = [1 1 1 1 1 1 1 1 1 1; 9.0 10.5 12.8 14.0 16.3 18.5 21.0 23.0 25.0 27.0];

tab(1,1) = "R2";
tab(1,2) = "SIGMA2";

% MODELO LINEAR:
% numero de parametros = 2
[b2, tab(2,1), tab(2,2)] = regressaoLinearMultipla(10, 1, 2, t1', Rt);

% MODELO QUADRATICO:
% numero de parametros = 3
[b3, tab(3,1), tab(3,2)] = regressaoLinearMultipla(10, 1, 3, t1', Rt);

% MODELO EXPONENCIAL:
% numero de parametros = 4
[b4, tab(4,1), tab(4,2)] = regressaoLinearMultipla(10, 1, 4, t1', Rt);

fprintf("Tabela com R2 e SIGMA2 dos modelos testados (cada linha um grau):\n");
disp(tab);
fprintf("Parâmetros do modelo linear:\n");
disp(b2);
fprintf("Parâmetros do modelo quadrático:\n");
disp(b3);
fprintf("Parâmetros do modelo exponencial:\n");
disp(b4);

fprintf("Vemos então, que o melhor modelo é o exponencial pois tem o SIGMA2 de\n");
fprintf("apenas 0.0011605 e o R2 de 0.99224. Que é o mais próximo de 1.\n");

fprintf("Para estimar a confiabilidade para t = 13.7 meses usando o melhor\n");
fprintf("modelo encontrado, substituirei no polinomio construido com os\n");
fprintf("parametros encontrados (pelo melhor modelo).\n\n")
tempo = 13.7;
POL = b4(1) + b4(2)*tempo + b4(3)*tempo^2 + b4(4)*tempo^3;

fprintf("Substituindo em 'b4(1) + tempo*b4(2) + b4(3)*tp^2 + b4(4)*tempo^3'\n");
fprintf("um tempo = 13.7, temos que a confiabilidade é: ");
disp(POL);

fprintf("Letra B:\n");
% Quantos meses são necessários para a confiabilidade atingir 80%?
% 80% em decimal é: 0.8, portanto R(t) = 0.8

RtLOG = [log(1.00) log(0.88) log(0.72) log(0.53) log(0.47) log(0.35) log(0.26) log(0.21) log(0.16) log(0.07)];
% Equacao 1: Rt = ae^(-kt) => log(a) - kt = log(0.8)
[b1, determinacao, variancia] = regressaoLinearMultipla(10, 1, 2, t1', RtLOG);
temp = (log(0.8) - b1(1))/b1(2);
fprintf("Usando que R(t) = ae^(-kt) => log(0.8) = log(a) - kt, temos que\n");
fprintf("tempo necessario para que R(t) = 0.8, é:");
disp(temp);
fprintf("Temos tambem, que o coeficiente de determinação é (R2): ");
disp(determinacao);
fprintf("E que a variancia residual é: ");
disp(variancia);

% SEGUNDA QUESTAO
fprintf("Segunda Questão:\n\n");
% Determine R³(13.7) usando Lagrange

% Para encontrar os parametros dos polinomios interpoladores de grau 3
% precisamos de 4 pontos (os mais proximos do ponto interpolado)
tLA = [10.5 12.8 14.0 16.3];
RtLA = [0.88 0.72 0.53 0.47];
INTLAGRANGE = intLagrange2(4, tLA, RtLA, 13.7);

fprintf("O resultado da interpolação de R³(13.7), por Lagrange é: ");
disp(INTLAGRANGE);

% Determine R²(24.7) usando Gregory-Newton

% Para encontrar os parametros dos polinomios interpoladores de grau 2
% precisamos de 3 pontos (os mais proximos do ponto interpolado)
tGN = [23.0 25.0 27.0];
RtGN = [0.21 0.16 0.07];
INTGREGNEWTON = intGregoryNewton(3, tGN, RtGN, 24.7);

fprintf("O resultado da interpolação de R²(24.7), por Gregory-Newton é: ");
disp(INTGREGNEWTON);

% FUNCOES

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

