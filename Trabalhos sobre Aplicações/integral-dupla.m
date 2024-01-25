
% TRABALHO 3.2 DE CÁLCULO NUMÉRICO 2023.2
% ALUNA: MAÍSA GARCIA NEPOMUCENO CORRÊA

fprintf("Implementação referente ao Trabalho 3.2 de Cálculo Numérico 2023.2.\n\n");
fprintf("Calcular, usando métodos computacionais de integração, a densidade de carga\n");
fprintf("em um retângulo.");

syms k;
syms w;

F = 2*k + 4*w; % Função da densidade de carga
a = 0; b = 5; c = 2; d = 5; % Limites do retângulo
plotR(a, b, c, d); % Plotando o retângulo

fprintf("Para o método de Newton-Cotes, usamos os graus 1,2 e 3 com 9-10 pontos.\n");
fprintf("Para o método de Gauss-Legendre, foi utilizado 20 pontos.\n");

fprintf("\n\nResultados:");
fprintf("\nNewton-Cotes Dupla com polinômio interpolador de grau 1: %.5f", metNewtonCotesDupla(a, b, 1, 10, c, d, 1, 10, F));
fprintf("\nNewton-Cotes Dupla com polinômio interpolador de grau 2: %.5f", metNewtonCotesDupla(a, b, 2, 10, c, d, 2, 10, F));
fprintf("\nNewton-Cotes Dupla com polinômio interpolador de grau 3: %.5f", metNewtonCotesDupla(a, b, 3, 9, c, d, 3, 9, F));
fprintf("\nGauss-Legendre Dupla: %.5f", metGaussLegendreDupla(a, b, 20, c, d, 20, F));

function plotR(a, b, c, d)
    x = linspace(0, 10, 100);
    y = linspace(0, 10, 100);
    x1 = a * ones(size(y));
    x2 = b * ones(size(y));
    y1 = c * ones(size(x));
    y2 = d * ones(size(x));
    figure;
    plot(x, y1, 'r-', 'LineWidth', 2);
    hold on;
    plot(x, y2, 'r-', 'LineWidth', 2);
    hold on;
    plot(x1, y, 'b-', 'LineWidth', 2);
    hold on;
    plot(x2, y, 'b-', 'LineWidth', 2);
    xlim([-1,10]);
    ylim([0,10]);
    grid on;
end

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
                % fprintf("i: %d\nx: %f\nc(i): %d\nj: %d\ny: %f\nc(j): %d\nfxy: %f\n\n", i, x, c(jx)*kx, j, y, c(jy)*ky, fxy);
            else
                % fprintf("j: %d\ny: %f\nc(j): %d\nfxy: %f\n\n", j, y, c(jy)*ky, fxy);
            end
        end
    end
    resultado = resultado*nx*ny*hx*hy/(d(nx)*d(ny));
end

function [A, T, erro] = givesPesoAbsGL(n)
    if(n < 1)
        erro = 1;
        error('Número de pontos menor que 1.');
    end
    m = fix(0.5*(n+1));
    erro = 0;

    for i = 1:m    
        z = cos(pi*(i-0.25)/(n+0.5));
        while (1)
            p1 = 1;     
            p2 = 0;      
            for j = 1:n % Polinômio de Legendre no ponto z          
                p3 = p2;          
                p2 = p1;            
                p1 = ((2*j-1)*z*p2 - (j-1)*p3)/j;      
            end      
            % Derivada do polinômio de Legendre no ponto z       
            pp = n*(z*p1-p2)/((z^2)-1);       
            z1 = z;      
            % Método de Newton para calcular os zeros do polinômio    
            z = z1 - (p1/pp);     
            if abs(z-z1) < 10e-15        
                break;
            end  
        end  
        T(m+1-i) = z; % Abscissa
        A(m+1-i) = 2/((1-z^2)*(pp^2)); % Peso
    end
end

function [resultado, erro] = metGaussLegendreDupla(ax, bx, nx, ay, by, ny, F)
    syms k;
    symK = sym('k');
    
    syms w;
    symW = sym('w');
   
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
