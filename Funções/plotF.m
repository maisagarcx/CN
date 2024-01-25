function plotF(a, b, num, x, F)
    u = linspace(a, b, num);
    Fx = vpa(subs(F, x, u));
    figure;
    plot(u, Fx);
    grid on;
    hold on;
    plot(u, zeros(size(u)), 'r--');
end
