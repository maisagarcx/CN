function [resultado, erro] = metGaussLegendreDupla(ax, bx, nx, ay, by, ny, F)
    syms k;
    ax = 2; bx = 6; nx = 6; ay = 1; by = 3; ny = 4;
    symK = sym('k');
    
    syms w;
    symW = sym('w');
    
    F = 
    vet = [symK, symW];
    
    % Cálculo dos pesos e abscissas
    [A, T, condErro] = givesPesoAbsGL(nx);
    resultado = 0;
    if condErro ~= 0
        return;
    end
    if ny == nx
        for j=1:fix(0.5*(nx + 1))
            B(j) = A(j);
            U(j) = T(j);
        end
    else
        [B, U, condErro] = givesPesoAbsGL(ny);
        if condErro ~= 0
            return;
        end
    end
    
    % Cálculo da integral dupla
    ex1 = (bx-ax)/2;
    ex2 = (ax+bx)/2;
    ey1 = (by-ay)/2;
    ey2 = (ay+by)/2;
    
    if mod(nx,2) == 0
        cx1 = 1;
        cx2 = 0.5;
    else
        cx1 = 0;
        cx2 = 1;
    end
    
    if mod(ny,2) == 0
        cy1 = 1;
        cy2 = 0.5;
    else
        cy1 = 0;
        cy2 = 1;
    end

    for i = 1:nx
        kx = fix(i-0.5*(nx+1) + sign(i - 0.5*(nx + cx1))*cx2);
        tx = sign(kx)*T(abs(kx));
        Axx = A(abs(kx));
        x = ex1*tx + ex2;
        soma = 0;
        for j = 1:ny
            ky = fix(j-0.5*(ny+1) + sign(j-0.5*(ny+cy1))*cy2);
            ty = sign(ky)*U(abs(ky));
            Ayy = B(abs(ky));
            y = ey1*ty + ey2;
            fxy = vpa(subs(F, vet, [x,y]));
            soma = soma + Ayy*fxy;
            if (j == 1)
                % fprintf("i: %d\ntx: %f\nx: %d\nAxx: %d\nj: %f\nty: %d\ny: %f\nAyy: %f\nfxy: %f\n\n", i, tx, x, Axx, j, ty, y, Ayy, fxy);
            else
                % fprintf("j: %d\nty: %f\ny: %d\nAyy: %f\nfxy: %f\n\n", j, ty, y, Ayy, fxy);
            end
        end
        resultado = resultado + Axx*soma;
    end
    resultado = ex1*ey1*resultado;
end
