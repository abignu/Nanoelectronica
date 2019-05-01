%constantes
h = 4.135e-15; %eV * s
eta = 0.5; %s
Ef = [-2.5 -3.5 -5.0]; %eV
HOMO = -5.5; %eV
LUMO = -1.5; %eV
gammaS = 0.1; %eV
ts = h / gammaS;
gammaD = 0.1; %eV
td = h / gammaD; %s
T = 298; %K
U = 0; %ya que suponemos Ces --> infinity; valor inicial
q = 1e-19;
%Ces = inf;
%damos unos valores a Vds
Vds = -6:0.1:6;
I = zeros(1,length(Vds));
I1 = zeros(1,length(Vds));
I2 = zeros(1,length(Vds));
I3 = zeros(1,length(Vds));
%apartado a

%apartado b; Calculamos la I-V y G-V
j = 1;
while j < 4 %le doy 3 vueltas, una por cada apartado
    No = 2 * fermi(HOMO,Ef(j),T) + 2 * fermi(LUMO,Ef(j),T);
    for i = 1:length(Vds)
        %calculo mu_S y mu_d
        mu_S = Ef(j) + eta * Vds(i);
        mu_D = Ef(j) +(eta - 1) * Vds(i);

        N = 2 * (td*fermiS(HOMO+U,mu_S,T) + ts*fermiD(HOMO+U,mu_D,T))/(ts+td) + 2 * (td*fermiS(LUMO+U,mu_S,T) + ts*fermiD(LUMO+U,mu_D,T))/(ts+td);
        Ucalc = N - No;
        Uold = 0.1;
        while abs(Ucalc-Uold) > 1e-3
            U=Uold + 1e-2 * (Ucalc - Uold);
            N = 2 * (td*fermiS(HOMO+U,mu_S,T) + ts*fermiD(HOMO+U,mu_D,T))/(ts+td) + 2 * (td*fermiS(LUMO+U,mu_S,T) + ts*fermiD(LUMO+U,mu_D,T))/(ts+td);
            Ucalc = N - No;
            Uold = U;
        end
        I(i) = 2 * q * (fermiS(HOMO+U,mu_S,T) - fermiD(HOMO+U,mu_D,T)) / (ts + td) + 2 * q * (fermiS(LUMO+U,mu_S,T) - fermiD(LUMO+U,mu_D,T)) / (ts + td); 
    end
    
    %chequeo para hacer las gráficas
    if j == 1
        I1 = I;
    end
    
    if j == 2
        I2 = I;
    end
    
    if j == 3
        I3 = I;
    end    
     
    j = j + 1; %continua el bucle
end


%saco conductancia
h = 0.00001;

G1 = diff(I1)/h;
G1(length(Vds)) = 0;

G2 = diff(I2)/h;
G2(length(Vds)) = 0;

G3 = diff(I3)/h;
G3(length(Vds)) = 0;

figure(1)
plot(Vds,I1)
hold on
plot(Vds,I2)
hold on
plot(Vds,I3)
hold on 
title('Ids-Vds characteristic')
xlabel('Vds (V)')
ylabel('I (A)')
legend('I1','I2','I3')

figure(2)
plot(Vds,G1)
hold on
plot(Vds,G2)
hold on 
plot(Vds,G3)
hold on
title('G-Vds characteristic')
xlabel('Vds (V)')
ylabel('G (s)')
legend('G1','G2','G3')