%Ejercicio 6 tema 3
%constantes
Epsilon = 8.85e-12;
ts = 1e-12;
td = ts;
r = 2e-9;
Cd = 1e-19;
Cs = 1e-19;
T = 0;

%Falta apartado a
gammaS = 4.14e-15 / ts; %en eV
gammaD = 4.14e-15 / td; %en eV

%apartado b
pi = sym(pi);
Ces = 4 * pi * r * Epsilon;
q = 1.602e-19;
Uc = q^2 / Ces; %Charging energy per electron
%me falta sacarla utilizando Cd y CS pero no se hacerlo

%apartado c
eta = 0.5;
Von = (8.0109e-20) / (q * eta); %puesto en julios y no en eV; Son 0.5eV
Idsmax = (2*q)/(ts+td);
Vdsmax = (Idsmax * td) / (Ces * eta) + Von;

figure(1)
plot([0 Von],[0 0], 'r') %primer tramo Ids = 0; Vds = Von
title('Ids-Vds characteristic')
xlabel('Vds (V)')
ylabel('Ids (A)')
hold on
plot([Von Vdsmax],[0 Idsmax], 'r') %segundo tramo
hold on
plot([Vdsmax Vdsmax+1],[Idsmax Idsmax], 'r') %tercer tramo

%apratado d
%q^2/Ces = 1eV --> ¿cómo altera a Von e Idsmax?
Uc = 1.60218e-19; %1eV/e-
Nmax = 1;
Imax = q * Nmax / td;
Imax/Idsmax %Como vemos Imax no altera su valor
%Por su parte Von tampoco altera su valor

%apartado e
%para obtener el número de e en el HOMO, N, tenemos en cuenta que la
%carga que entra y sale de la molécula es la misma: Id = Is
%Is = q/ts * (Ns-N); Id = q/td * (N-Nd)
%igualo y despejo N obteniendo N = (td*Ns+ts*Nd)/(ts+td)
syms N x;
N = 2 * dirac(x) * (td*1+td*0)/(td+ts); %las funciones de distribución pasan a ser la funcion escalon
%el source está ocupado por eso vale 1 y el drain vacio por eso vale 0
N = int(N,x,-inf,inf)

%apartado f
%Cuando q^2/Ces = 1eV/e 
Ucmax = 1.60218e-19 * Nmax; %da como resultad 1eV

%apartado g
Vds = 2; %en voltios
HOMO = -8.0109e-20; %-0.5 eV
Ef = 0; %0 eV
mu_S = (Cd * q * Vds) / (Cd + Cs);
mu_D = -(Cs * q * Vds) / (Cd + Cs);
figure(2)
plot([HOMO HOMO], 'y')
title('Niveles energéticos')
ylabel('Energía (Julios)')
hold on
%plot([Ef Ef], 'r')
hold on
plot([mu_S mu_S], 'b')
hold on
plot([mu_D mu_D], 'g')
hold on
legend('HOMO','mu_S','mu_D')

%apartado h

figure(3)
plot([0 Von],[0 0], 'r') %primer tramo Ids = 0; Vds = Von
title('Ids-Vds characteristic')
xlabel('Vds (V)')
ylabel('Ids (A)')
hold on
plot([Von Vdsmax],[0 Idsmax], 'r') %segundo tramo
hold on
plot([Vdsmax Vdsmax+1],[Idsmax Idsmax], 'r') %tercer tramo
hold on 
%ahora la otra gráfica
plot([0 Von],[0 0], 'b') %primer tramo Ids = 0; Vds = Von
hold on
Vdsmax = (Idsmax * td * q)/(eta) + Von;
plot([Von Vdsmax],[0 Idsmax], 'b') %segundo tramo
hold on
plot([Vdsmax Vdsmax+1],[Idsmax Idsmax], 'b') %tercer tramo
hold on
legend('q^2/Ces = 0.72 eV','q^2/Ces = 1 eV')
