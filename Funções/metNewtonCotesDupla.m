function [resultado, erro] = metNewtonCotesDupla(ax, bx, nx, mx, ay, by, ny, my, F)
    
    syms k;
    symK = sym('k');
    
    syms w;
    symW = sym('w');
    
    vet = [symK, symW];
    
    d(1)=2;d(2)=6;d(3)=8;d(4)=90;d(5)=288;d(6)=840;d(7)=17280;d(8)=28350;
    c(1)=1;c(2)=1;c(3)=4;c(4)=1;c(5)=3;c(6)=7;c(7)=32;c(8)=12;c(9)=19;c(10)=75;
    c(11)=50;c(12)=41;c(13)=216;c(14)=27;c(15)=272;c(16)=751;c(17)=3577;
    c(18)=1323;c(19)=2989;c(20)=989;c(21)=5888;c(22)=-928;c(23)=10496;
    c(24)=-4540;
    
    % Consistência dos parâmetros
    if (((nx<1) || (nx>8)) || ((ny<1) || (ny>8)))
        erro = 1;
        error('Grau do polinômio interpolador fora do permitido.');   
    end

    if (mod(mx,nx) ~= 0) || (mod(my,ny ~= 0))  
        erro = 2;
        error('Número de subintervalos não é múltiplo do grau do polinômio interpolador.');
    end
    
    % Cálculo da integral
    px = fix(0.25*(nx*(nx + 2)+ mod(nx,2)));
    py = fix(0.25*(ny*(ny + 2)+ mod(ny,2)));
    hx = (bx-ax)/mx;
    hy = (by-ay)/my;
    resultado = 0;
    
    for i = 0:mx
        x = ax + i*hx;
        jx = px + fix(0.5*nx - abs(mod(i,nx) - 0.5*nx));
        kx = 1 + fix((nx - mod(i,nx))/nx) - fix((mx - mod(i,mx))/mx);
        for j = 0:my
            y = ay + j*hy;
            jy = py + fix(0.5*ny - abs(mod(j,ny) - 0.5*ny));
            ky = 1 + fix((ny - mod(j,ny))/ny) - fix((my - mod(j,my))/my);
            fxy = vpa(subs(F, vet, [x,y]));
            resultado = resultado + fxy*c(jx)*kx*c(jy)*ky;
            if (j == 0)
                fprintf("i: %d\nx: %f\nc(i): %d\nj: %d\ny: %f\nc(j): %d\nfxy: %f\n\n", i, x, c(jx)*kx, j, y, c(jy)*ky, fxy);
            else
                fprintf("j: %d\ny: %f\nc(j): %d\nfxy: %f\n\n", j, y, c(jy)*ky, fxy);
            end
        end
    end
    resultado = resultado*nx*ny*hx*hy/(d(nx)*d(ny));
end
