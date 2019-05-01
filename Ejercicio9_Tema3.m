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
U = 0; %ya que suponemos Ces --> infinity; valor inicial

Ces = inf;
%damos unos valores a Vds
Vds = -4:0.1:4;
I1 = zeros(1,length(Vds));
%N = @(E) DOS(E,U, HOMO)* (td * fermiS(E,mu_S,T) + ts * fermiD(E,mu_D,T)) / (ts + td);
%I = @(E) q * DOS(E,U, HOMO)* (fermiS(E,mu_S,T) - fermiD(E,mu_D,T)) / (ts + td);

%apartado a (acá U = 0)
for i = 1:length(Vds)
    %calculo mu_S y mu_d
    mu_S = Ef + eta .* Vds(i);
    mu_D = Ef +(eta - 1) .* Vds(i);
   
    I1(i) = 2 * q * (fermiS(HOMO,mu_S,T) - fermiD(HOMO,mu_D,T)) / (ts + td);
end
%ploteo
figure(1)
plot(Vds,I1)
%saco conductancia
h = 0.00001;
G1 = diff(I1)/h;
G1(length(Vds)) = 0;
figure(2)
plot(Vds,G1)

q = 1e-19;
Von = (Ef-HOMO) / (q * eta); %puesto en julios y no en eV; Son 0.5eV
Idsmax = (2*q)/(ts+td);
Vdsmax = (Idsmax * td) / (Ces * eta) + Von; %Ces = inf

figure(3)
plot([0 Von],[0 0], 'r') %primer tramo Ids = 0; Vds = Von
title('Ids-Vds characteristic')
xlabel('Vds (V)')
ylabel('Ids (A)')
hold on
plot([Von Vdsmax],[0 Idsmax], 'r') %segundo tramo
hold on
plot([Vdsmax Vdsmax+1],[Idsmax Idsmax], 'r') %tercer tramo

I2 = zeros(1,length(Vds));
%apartado b; Calculamos la I-V y G-V
No = 2 * fermi(HOMO,Ef,T);
for i = 1:length(Vds)
    %calculo mu_S y mu_d
    mu_S = Ef + eta * Vds(i);
    mu_D = Ef +(eta - 1) * Vds(i);
    
    N = 2 * (td*fermiS(HOMO+U,mu_S,T) + ts*fermiD(HOMO+U,mu_D,T))/(ts+td);
    Ucalc = N - No;
    Uold = 0.1;
    while abs(Ucalc-Uold) > 1e-3
        U=Uold + 1e-2 * (Ucalc - Uold);
        N = 2 * (td*fermiS(HOMO+U,mu_S,T) + ts*fermiD(HOMO+U,mu_D,T))/(ts+td);
        Ucalc = N - No;
        Uold = U;
    end
    I2(i) = 2 * q .* (fermiS(HOMO+U,mu_S,T) - fermiD(HOMO+U,mu_D,T)) / (ts + td); 
end
length(Vds)
length(I2)
figure(4)
plot(Vds,I2)


%saco conductancia
h = 0.00001;
G2 = diff(I2)/h;
G2(length(Vds)) = 0;
figure(5)
plot(Vds,G2)

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