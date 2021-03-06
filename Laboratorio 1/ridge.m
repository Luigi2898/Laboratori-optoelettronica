tic
clc
clear all
close all

n_core=3.475;
n_cladding=1.44;
neff = linspace(n_cladding,n_core,1e4);
d=220e-9;
%per modo TE
n_1=n_core;
n_2=n_cladding;
n_3=n_2;
lambda=1.55e-6;
k_0=(2*pi/lambda);
delta = k_0.*sqrt(neff.^2-n_3^2);
k_x = k_0.*sqrt(n_1^2-neff.^2);
gamma = k_0*sqrt(neff.^2 - n_2^2);
for m = 0:1:5
    d_TE = (atan(gamma./k_x) + atan(delta./k_x) + m*pi)./(k_x);
    figure(1)
    grid on
    plot(d_TE,neff);
    hold on
    xlabel Spessore[m]
    ylabel neff
    xlim ([0 5e-6])
    ylim ([n_cladding n_core])
    title 'Plot di neff in funzione dello spessore per i modi TE'
end
d_TE = (atan(gamma./k_x) + atan(delta./k_x))./(k_x);

for m = 0:1:5
    omega_TE = (3e8*(atan(delta./k_x)+atan(gamma./k_x) + m*pi))./(d*sqrt(n_1^2-neff.^2));
    figure(2)
    plot(omega_TE,neff)
    hold on
    title 'Plot di neff in funzione di omega per i modi TE'
    xlabel omega[rad/s]
    ylabel neff
    grid on
    xlim ([0 2e16])
    ylim ([n_cladding n_core])
end
omega_TE = (3e8*(atan(delta./k_x)+atan(gamma./k_x)))./(d*sqrt(n_1^2-neff.^2));
for m = 0:1:5
    d_TM =  (atan(gamma.*n_1^2./(k_x*n_2^2)) + atan(delta.*n_1^2./(k_x*n_2^2)) + m*pi)./(k_x);
    figure(3)
    plot(d_TM,neff);
    hold on
    grid on
    xlabel Thickness[m]
    ylabel neff
    xlim ([0 5e-6])
    ylim ([n_cladding n_core])
    title 'Plot di neff in funzione dello spessore per il modo TM'
end
d_TM =  (atan(gamma.*n_1^2./(k_x*n_2^2)) + atan(delta.*n_1^2./(k_x*n_2^2)))./(k_x);
%plot(d_TE_1,neff,'r');
%hold on
%plot(d_TE_2,neff,'y');
for m = 0:1:5
    omega_TM =  (3e8 * (atan(gamma.*n_1^2./(k_x*n_2^2)) + atan(delta.*n_1^2./(k_x*n_2^2)) + m*pi))./(d*sqrt(n_1^2-neff.^2));
    figure(4)
    plot(omega_TM,neff);
    hold on
    grid on
    title 'Plot di neff in funzione di omega per i modi TM'
    xlabel omega[rad/s]
    ylabel neff
    xlim ([0 2e16])
    ylim ([n_cladding n_core])
end
omega_TM =  (3e8 * (atan(gamma.*n_1^2./(k_x*n_2^2)) + atan(delta.*n_1^2./(k_x*n_2^2))))./(d*sqrt(n_1^2-neff.^2));

 %CALCOLO DEL BETA
beta_TE0=neff.*omega_TE./3e8;
beta_TEco = n_core.*omega_TE./3e8;
beta_TEcl = n_cladding.*omega_TE./3e8;
figure(5)
plot(omega_TE,beta_TE0)
hold on
plot(omega_TE,beta_TEco)
plot(omega_TE,beta_TEcl)
grid on
xlim ([0 10e15])
title 'Plot di \beta(\omega) in funzione di \omega per i modi TE'
xlabel \omega[rad/s]
ylabel \beta(\omega)

beta_TM0=neff.*omega_TM./3e8;
beta_TMco = n_core.*omega_TM./3e8;
beta_TMcl = n_cladding.*omega_TM./3e8;
figure(6)
plot(omega_TM,beta_TM0)
hold on
plot(omega_TM,beta_TMco)
plot(omega_TM,beta_TMcl)
grid on
xlim ([0 10e15])
title 'Plot di \beta(\omega) in funzione di \omega per i modi TM'
xlabel \omega[rad/s]
ylabel \beta(\omega)

%velocit� di gruppo
group_velocity_TE=1./(diff(beta_TE0)./(diff(omega_TE)));
figure(7)
plot(omega_TE(1:9999),group_velocity_TE)
hold on
grid on
xlim ([0 4e16])
title 'Plot della velocità di gruppo in funzione di \omega TE'
xlabel \omega[rad/s]
ylabel 'indice di gruppo'

%indice di gruppo
group_index_TE=(3e8)./(group_velocity_TE);
figure(8)
plot(omega_TE(1:9999),group_index_TE)
hold on
grid on
xlim ([0 4e16])
title 'Plot dell indice di gruppo in funzione di \omega TE'
xlabel \omega[rad/s]
ylabel 'indice di gruppo'

indiceTE220 = 6912; %220nm
%profilo di campo TE 220nm

x_TE=linspace(-550e-9,380e-9,1e4);
x_TE1 = x_TE(find(x_TE < -d));
x_TE2 = x_TE(find(x_TE<=0 & x_TE>=-d));
x_TE3 = x_TE(find(x_TE > 0));

n_eff_fix = neff(indiceTE220);
delta_fix = delta(indiceTE220);
gamma_fix = gamma(indiceTE220);
kx_fix = k_x(indiceTE220);
B=-(delta_fix)/(kx_fix);

Ey_TE = [(cos(kx_fix.*d)-(B.*sin(kx_fix.*d))).*exp(gamma_fix.*(x_TE1+d)) cos(kx_fix.*x_TE2)+(B.*sin(kx_fix.*x_TE2)) exp(-delta_fix.*x_TE3)];

zz = zeros(1,1000);
x = linspace(0,2.7,1000);
figure(9)
int_norm_Ey_TE = abs(Ey_TE).^2./max(abs(Ey_TE).^2);
plot(int_norm_Ey_TE, x_TE)
hold on
plot(x,zz)
plot(x,zz-d)
title 'Plot del profilo del modo TE fondamentale 220nm'
xlim([0 1])
grid on

%profilo di campo TM 220nm
indiceTM220 = 2994;
x_TM=linspace(-1100e-9,1100e-9,1e4);
x_TM1 = x_TM(find(x_TM < -d));
x_TM2 = x_TM(find(x_TM <= 0 & x_TM >= -d));
x_TM3 = x_TM(find(x_TM > 0));
n_eff_fix = neff(indiceTM220);
delta_fix = delta(indiceTM220);
gamma_fix = gamma(indiceTM220);
kx_fix = k_x(indiceTM220);
B=-(delta_fix*n_1^2)/(n_2^2*kx_fix);%porco il cazzo
Hy_TM = [(cos(kx_fix*d)-(B.*sin(kx_fix*d)))*exp(gamma_fix*(x_TM1+d)) cos(kx_fix*x_TM2)+B*sin(kx_fix*x_TM2) exp(-delta_fix*x_TM3)];

int_norm_Hy_TM = abs(Hy_TM).^2./max(abs(Hy_TM).^2);

figure(10)
zz = zeros(1,1000);
x = linspace(0,11,1000);
plot(int_norm_Hy_TM, x_TM)
hold on
plot(x,zz)
plot(x,zz-d)
title 'Plot del profilo del modo TM fondamentale 220nm'
xlim([0 1])
grid on

%profilo di campo TE 150nm
indiceTE150 = 5395; %150nm
d = 150e-9;
x_TE=linspace(-575e-9,425e-9,1e4);
x_TE1 = x_TE(find(x_TE < -d));
x_TE2 = x_TE(find(x_TE<=0 & x_TE>=-d));
x_TE3 = x_TE(find(x_TE > 0));
n_eff_fix = neff(indiceTE150);
delta_fix = delta(indiceTE150);
gamma_fix = gamma(indiceTE150);
kx_fix = k_x(indiceTE150);
B=-(delta_fix)/(kx_fix);
Ey_TE_150 = [(cos(kx_fix.*d)-(B.*sin(kx_fix.*d))).*exp(gamma_fix.*(x_TE1+d)), cos(kx_fix.*x_TE2)+(B.*sin(kx_fix.*x_TE2)), exp(-delta_fix.*x_TE3)];
zz = zeros(1,1000);
x = linspace(0,2,1000);

int_norm_Ey_TE_150 = abs(Ey_TE_150).^2./max(abs(Ey_TE_150).^2);

figure(11)
plot(int_norm_Ey_TE_150, x_TE)
hold on
plot(x,zz)
plot(x,zz-d)
title 'Plot del profilo del modo TE fondamentale 150nm'
xlim([0 1])
grid on

%profilo di campo TM 150nm
indiceTM150 = 834;
x_TM=linspace(-575e-9,425e-9,1e4);
x_TM1 = x_TM(find(x_TM < -d));
x_TM2 = x_TM(find(x_TM <= 0 & x_TM >= -d));
x_TM3 = x_TM(find(x_TM > 0));
n_eff_fix = neff(indiceTM150);
delta_fix = delta(indiceTM150);
gamma_fix = gamma(indiceTM150);
kx_fix = k_x(indiceTM150);
B=-(delta_fix*n_1^2)/(n_2^2*kx_fix);%porco il cazzo
Hy_TM_150 = [(cos(kx_fix*d)-(B.*sin(kx_fix*d)))*exp(gamma_fix*(x_TM1+d)) cos(kx_fix*x_TM2)+B*sin(kx_fix*x_TM2) exp(-delta_fix*x_TM3)];

int_norm_Hy_TM_150 = abs(Hy_TM_150).^2./max(abs(Hy_TM_150).^2);

figure(12)
zz = zeros(1,1000);
x = linspace(0,3,1000);
plot(int_norm_Hy_TM_150, x_TM)
hold on
plot(x,zz)
plot(x,zz-d)
title 'Plot del profilo del modo TM fondamentale 150nm'
xlim([0 1])
grid on

%Metodo ERI per analisi ridge
n_coreERI = neff(indiceTE220);
n_claddingERI = neff(indiceTE150);
n_1ERI = n_coreERI;
n_2ERI = n_claddingERI;
n_3ERI = n_2ERI;
neffERI = linspace(n_claddingERI, n_coreERI, 1e4);
dERI=500e-9;
k_0ERI = (2*pi/lambda);
deltaERI = k_0ERI.*sqrt(neffERI.^2-n_3ERI^2);
k_xERI = k_0ERI.*sqrt(n_1ERI^2-neffERI.^2);
gammaERI = k_0ERI*sqrt(neffERI.^2 - n_2ERI^2);
d_TMERI =  (atan(gammaERI.*n_1ERI^2./(k_xERI*n_2ERI^2)) + atan(deltaERI.*n_1ERI^2./(k_xERI*n_2ERI^2)))./(k_xERI);
omega_TMERI =  (3e8 * (atan(gammaERI.*n_1ERI^2./(k_xERI*n_2ERI^2)) + atan(deltaERI.*n_1ERI^2./(k_xERI*n_2ERI^2))))./(dERI*sqrt(n_1ERI^2-neffERI.^2));
x_TMERI=linspace(-1000e-9,500e-9,15000);
indiceTMERI = 5274; %500nm
n_eff_fixERI = neffERI(indiceTMERI);
delta_fixERI = deltaERI(indiceTMERI);
gamma_fixERI = gammaERI(indiceTMERI);
kx_fixERI = k_xERI(indiceTMERI);
BERI=-(delta_fixERI*n_1ERI^2)/(n_2ERI^2*kx_fixERI);%porco il cazzo
%Calcolo Hy per rami
Hy_TMERI(find(x_TMERI > 0)) = exp(-delta_fixERI*x_TMERI(find(x_TMERI > 0)));
Hy_TMERI(find(x_TMERI<=0 & x_TMERI>=-dERI))= cos(kx_fixERI.*x_TMERI(find(x_TMERI<=0 & x_TMERI>=-dERI)))+(BERI.*sin(kx_fixERI.*x_TMERI(find(x_TMERI<=0 & x_TMERI>=-dERI))));
Hy_TMERI(find(x_TMERI < -dERI))= (cos(kx_fixERI.*dERI)-(BERI.*sin(kx_fixERI.*dERI))).*exp(gamma_fixERI.*(x_TMERI(find(x_TMERI < -dERI))+dERI));

int_norm_Hy_TMERI = abs(Hy_TMERI).^2./max(abs(Hy_TMERI).^2);

figure(13)
plot(int_norm_Hy_TMERI,x_TMERI)
hold on
zzERI = zeros(1,1000);
xERI = linspace(0,3.5,1000);
plot(xERI,zzERI)
plot(xERI,zzERI-dERI)
plot(xERI,zzERI-dERI/2)
title 'Plot del profilo del modo TM 500nm ERI Hy'
grid on
xlim([0 1])

omega_TEERI = 3e8*2*pi/lambda;
epsi0 = 8.9e-12;
%Ex_TMERI = (kx_fixERI.*Hy_TMERI)./(omega_TEERI*epsi0*n_claddingERI^2);
Ex_TMERI(find(x_TMERI > 0))=(kx_fixERI.*Hy_TMERI(find(x_TMERI > 0)))./(omega_TEERI*epsi0*n_claddingERI^2);
Ex_TMERI(find(x_TMERI<=0 & x_TMERI>=-dERI))= (kx_fixERI.*Hy_TMERI(find(x_TMERI<=0 & x_TMERI>=-dERI)))./(omega_TEERI*epsi0*n_coreERI^2);
Ex_TMERI(find(x_TMERI < -dERI))= (kx_fixERI.*Hy_TMERI(find(x_TMERI < -dERI)))./(omega_TEERI*epsi0*n_claddingERI^2);

int_norm_Ex_TMERI = abs(Ex_TMERI).^2./max(abs(Ex_TMERI).^2);

figure(14)
plot(x_TMERI,int_norm_Ex_TMERI)
hold on
zzERI = zeros(1,1000);
xERI = linspace(0,5000,1000);
plot(zzERI,xERI)
plot(zzERI-dERI,xERI)
title 'Plot del profilo del modo TM 500nm ERI Ex'
grid on
ylim([0 1])

prodot = Ey_TE'.*Ex_TMERI;
prodot_norm = int_norm_Ey_TE'.*int_norm_Ex_TMERI;
% guide = zeros(size(prodot));
% prodot(5950:6050,:) = ones(101,numel(prodot(1,:)));
% prodot(5950-1500:6050-1500,1:5000) = ones(101,5000);
% prodot(5950-2200:6050-2200,5001:10000) = ones(101,5000);
% prodot(5950-1500:6050-1500,10001:15000) = ones(101,5000);
% prodot(3750:4550,4950:5050) = ones(801,101);
% prodot(3750:4550,4950+5000:5050+5000) = ones(801,101);

prodot_norm(5950:6050,:) = ones(101,numel(prodot(1,:))).*0.5;
 prodot_norm(5950-1500:6050-1500,1:5000) = ones(101,5000).*0.5;
 prodot_norm(5950-2200:6050-2200,5001:10000) = ones(101,5000).*0.5;
 prodot_norm(5950-1500:6050-1500,10001:15000) = ones(101,5000).*0.5;
 prodot_norm(3750:4550,4950:5050) = ones(801,101).*0.5;
 prodot_norm(3750:4550,4950+5000:5050+5000) = ones(801,101).*0.5;
%FIXME Aggiungere contorno guida (matrice con tutti 0 e 100 (da vedere in base al colore) dove voglio vedere il contorno
%FIXME Sistemare in modo da vederla simmetrica
figure (15)
imagesc(prodot)
figure(16)
imagesc(prodot_norm)

%figure
%mesh(x_TMERI, x_TE, prodot_norm)

toc
