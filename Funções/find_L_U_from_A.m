function [matriz_L, matriz_U] = find_L_U_from_A(matriz_composta)
    ordem_da_matriz = length(matriz_composta);
    matriz_L = zeros(ordem_da_matriz);
    matriz_U = zeros(ordem_da_matriz);

    for linha = 1:ordem_da_matriz
        for coluna = 1:ordem_da_matriz
            if linha == coluna
                matriz_L(linha, coluna) = 1.0;
            end
            if coluna < linha
                matriz_L(linha, coluna) = matriz_composta(linha, coluna);
            end
            if coluna >= linha
                matriz_U(linha, coluna) = matriz_composta(linha, coluna);
            end
        end
    end
end
