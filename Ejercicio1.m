syms x sigx;
%pi = sym(pi);
phi(x,sigx)=(1/((2 * pi())^(1/4) * sigx^0.5)) * exp(-0.25 * (x/sigx)^2);

x = -5:0.1:5;
sigx = -2:0.5:2;

for i = 1:length(sigx)
    if sigx(i) ~= 0
        plot(x,phi(x,sigx(i)))
        title('Gausiana para muchos valores de sig_x')
        xlabel('Valor de x')
        ylabel('Valor de Gaussiana')
        hold on
    end
end

%syms y sig pi;
psi = @(y,sig) (1/((2 * 3.14159)^(1/4) * (sig).^0.5))^2 * exp(-0.5 * (y./sig).^2);

%ahora integro la función 
I = integral(@(y)psi(y,0.5),-inf, inf)