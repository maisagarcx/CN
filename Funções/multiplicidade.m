function [m, erro] = multiplicidade(raiz, x, F1)
    DFx = diff(F1, x);
    m = 1;
    max = 50;
    while (vpa(subs(DFx, x, raiz)) == 0)
        DFx = diff(DFx, x);
        m = m + 1;
        if m == max
            erro = 1;
            break;
        end
    end
    erro = 0;
end
