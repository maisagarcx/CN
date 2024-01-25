function [solucao, k] = metNewtonRaphsonN(vet_func, simbolos, n, solucao, tolerancia, maxIter)
    
    matrixN = zeros(n,n);
    vetN = zeros(1,n);
    k = 0;
    while (true)
        k = k + 1;
        if (k > maxIter)
            fprintf("\n\nO número de iterações excedeu o permitido!");
            break
        end

        for i=1:n
            vetN(i) = subs(vet_func(i), simbolos, solucao);
        end

        for i=1:n
            for j=1:n
                matrixN(i,j) = vpa(subs(diff(vet_func(i),simbolos(j)), simbolos, solucao));
            end
        end

        vet_b = (-1) * vetN;
        %disp(vet_b);

        % [LU, ~, pivos] = calculaPivotacaoParcial(matrixA, vet_b, n, pivos);
        %[matrixA, ~, pivot] = calculaPivotacaoParcial(matrixN, vet_b, n);
        % [L, U, P] = givesLUP(matrixN);
        %disp(matrixN);
        [matrixA, ~, pivot] = dec_LU(n, matrixN);
        %disp(matrixA);
        % matrix = dec_LU(n, matrixA);
        % disp(matrix);
        [L, U] = find_L_U_from_A(matrixA);
        %disp(L);
        %disp(U);
        Y = suc_subst_piv(L, vet_b, pivot);
        incremento = ret_subst(U, Y);
        
        % disp(incremento);
        
        for k=1:n
            solucao(k) = solucao(k) + incremento(k);
        end

        max_n = maxMag(incremento, n);
        %disp(max_n);

        if (max_n <= tolerancia)
            fprintf("Foi encontrada a solução com a precisão desejada\n");
            %fprintf("A última iteração é: %d\n", iteracao);
            fprintf("O último dx é aproximadamente: %.15f\n", round(solucao(1),15));
            fprintf("O último dy é aproximadamente: %.15f\n", round(solucao(2),15));
            fprintf("O último dz é aproximadamente: %.15f\n", round(solucao(3),15));
            break
        end
        k=k+1;
    end
end
