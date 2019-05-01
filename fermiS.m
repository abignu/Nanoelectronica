function F=fermiS(E, mu_S, T)
    k = 8.617e-5; %eV/K
    F = 1 ./ (1 + exp((E-mu_S)./(k*T)));
end