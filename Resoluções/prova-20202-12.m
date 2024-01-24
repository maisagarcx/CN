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
