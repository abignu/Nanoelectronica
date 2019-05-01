function F=fermiD(E, mu_D, T)
    k = 8.617e-5; %eV/K
    F = 1 ./ (1 + exp((E-mu_D)./(k*T)));
end