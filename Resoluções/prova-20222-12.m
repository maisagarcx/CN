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
