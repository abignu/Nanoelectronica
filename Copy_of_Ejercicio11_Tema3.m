%constantes
h = 4.135e-15; %eV * s
eta = 0.5; %s
Ef = -5.0; %eV
HOMO = -5.5; %eV
gammaS = 0.1; %eV
ts = h / gammaS;
gammaD = 0.1; %eV
td = h / gammaD; %s
T = 298; %K
q = 1.6e-19;
gamma = gammaS + gammaD;
%mu_S = 0;
%mu_D = 0;
%damos unos valores a Vds
Vds = -4:0.1:4;
I1 = zeros(1,length(Vds));
U = 0;
g = @(E) (1/3.14159) .* (gamma) ./ ((E - HOMO).^2 + (gamma/2).^2); %Lorentziana

No = @(E) ((fermi(E+HOMO,mu_S,T))).*g(E-U);
%apartado a - no broadening; Calculamos la I-V y G-V
N1o = integral(@(E)No(E), -inf , inf);
for i = 1:length(Vds)
    %calculo mu_S y mu_d
    mu_S = Ef + eta * Vds(i);
    mu_D = Ef +(eta - 1) * Vds(i);
    
    g2 = @(E) (gamma/3.14159) ./ ((E - HOMO).^2 + (gamma/2).^2); %Lorentziana
    N = @(E) ((td.*fermiS(E,mu_S,T) + ts.*fermiD(E,mu_D,T))./(ts+td)).*g2(E-U);
    
    
    N1 = integral(@(E)N(E), HOMO+U-15*gamma, HOMO+U+15*gamma);
   
    Ucalc = N1 - N1o;
    Uold = 0.1;
    while abs(Ucalc-Uold) > 1e-2
        U=Uold + 1e-2 * (Ucalc - Uold);
        N = @(E) ((td.*fermiS(E,mu_S,T) + ts.*fermiD(E,mu_D,T))./(ts+td)).*g2(E-U);
        
        N1 = integral(@(E)N(E), HOMO+U-15*gamma, HOMO+U+15*gamma);
        Ucalc = N1 - N1o;
        Uold = U;
    end
    I = @(E) ((fermiS(E,mu_S,T) - fermiD(E,mu_D,T)) ./ (ts + td)).*g2(E-U);
    I1(i) = q * integral(@(E)I(E), HOMO+U-15*gamma, HOMO+U+15*gamma); 
   % I1(i)
end


length(Vds)
length(I1)
%figure(1)
%plot(Vds,I1)

h = 0.00001;
G1 = diff(I1)/h;
G1(length(Vds)) = 0;
%figure(2)
%plot(Vds,G1)

%HOMO = linspace(-1,1,length(Vds)); %valor continuo de la energía
I2 = zeros(1,length(Vds));
%apartado a - no broadening; Calculamos la I-V y G-V
N2o = fermi(HOMO,Ef,T);
for i = 1:length(Vds)
    %calculo mu_S y mu_d
    mu_S = Ef + eta * Vds(i);
    mu_D = Ef +(eta - 1) * Vds(i);
    
    N2 = ((td*fermiS(HOMO+U,mu_S,T) + ts*fermiD(HOMO+U,mu_D,T))./(ts+td));
    Ucalc = N2 - N2o;
    Uold = 0.1;
    while abs(Ucalc-Uold) > 1e-3
        U=Uold + 1e-2 * (Ucalc - Uold);
        N2 = ((td*fermiS(HOMO+U,mu_S,T) + ts*fermiD(HOMO+U,mu_D,T))./(ts+td));
        Ucalc = N2 - N2o;
        Uold = U;
    end
    I2(i) = (q * (fermiS(HOMO+U,mu_S,T) - fermiD(HOMO+U,mu_D,T)) ./ (ts + td)); 
end
length(Vds)
length(I2)
%figure(3)
%plot(Vds,I2)
%saco conductancia
h = 0.00001;
G2 = diff(I2)/h;
G2(length(Vds)) = 0;
%figure(4)
%plot(Vds,G2)


figure(6)
title('Ids-Vds characteristic')
plot(Vds,I1)
hold on
plot(Vds,I2)
hold on 
legend('I1','I2')

figure(7)
title('G-Vds characteristic')
plot(Vds,G1)
hold on
plot(Vds,G2)
hold on 
legend('G1','G2')