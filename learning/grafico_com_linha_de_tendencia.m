% Dados
distancias = [0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0];
tempo = [0.415 0.551 0.668 0.826 0.885 0.984 1.128 1.169];

% Ajuste linear
coeffs = polyfit(tempo, distancias, 1);
ajuste_linear = polyval(coeffs, tempo);

% Criar o gráfico
figure;
plot(tempo, distancias, 'o', tempo, ajuste_linear, '-');
xlabel('Tempo (s)');
ylabel('Distância (M)');
title('Experimento II');
legend('Dados', 'Ajuste Linear', 'Location', 'Northwest');
grid on;

% Exibir os coeficientes da reta de ajuste linear
disp('Coeficientes da reta de ajuste linear:');
disp(['Inclinação (a): ' num2str(coeffs(1))]);
disp(['Interceptação (b): ' num2str(coeffs(2))]);
