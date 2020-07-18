clear all
close all
clc
load('lab4W.mat')
L = 5e-3; %m
lambda = 1550e-9; %m
open('meas_phase.fig')
D = get(gca, 'Children');
XData = get(D, 'XData');
YData = get(D, 'YData');
V = cell2mat(XData(1));
Dphi = cell2mat(YData(1));
%FIT CON IL TOOLBOX PER OTTENERE LA CARATTERISTICA Dphi(V)
Vfit = linspace(0,6,1000);
Dphifit = Dphi_V(Vfit);
Dneff_V = lambda .* Dphifit ./ (L * 2 * pi);
Dneff_exp = lambda .* Dphi ./ (2 * pi * L);
close

figure
plot(Vfit, Dneff_V);
hold on
plot(V, Dneff_exp, 'o');
grid on
title '\Deltan_{eff} vs V'

open('TdB_meas.fig')
D1 = get(gca, 'Children');
Vfig = get(D1, 'XData');
TLfig = get(D1, 'YData');
%FIT CON IL TOOLBOX PER OTTENERE LA CARATTERISTICA T(V)
Tfit = T_Vfit(Vfit);
alpha_fit = Tfit./(4.34*L);
alpha_exp = TLfig./(4.34*L);
close

figure
plot(Vfit, alpha_fit./100);
hold on
plot(Vfig, alpha_exp./100, 'o');
grid on
title 'ALPHA vs V'
%ALPHA DA m^-1 A cm^-1
alpha_fit = alpha_fit./100;


%MZI modulator: transmission versus bias voltage
%tutto convertito in cm
neff = 3.75; %Dall'articolo
DL = 1e-4*100; %m
reverse_bias = [0:0.1:2]; %V
alpha1 = -log(10^(-8.7/20))/(L*100);
Ton = 10^(-8.7/20);
lambdavect = linspace(1545e-9, 1565e-9, 1000) .* 100;
lambda = lambda * 100;
L = L * 100;
DneffV = lambda .* Dphi_V(reverse_bias) ./ (L * 2 * pi);
Dphi_rev = Dphi_V(reverse_bias);

for j = 1:1:numel(reverse_bias)
  Dneff1 = DneffV(j);
  T_lambda = (1/4).*Ton.*abs(exp(-(2i*pi).*Dneff1.*L./(lambdavect))+exp(-alpha1*DL)*exp((-2i*pi).*(neff*DL)./(lambdavect))).^2;
  %T_lambda = Ton.*(1/2).*(1+cos(2*pi*Dneff1*L*reverse_bias(j)./(lambdavect)))
  %T_lambda = (1/2).*(1+cos(Dphi_rev))
  figure(3)
  plot(lambdavect*10^7, 20*log10(T_lambda))
  hold on
  grid on
end
figure(3)
ylabel("Transmission spectrum [dB]")
xlabel("\lambda [nm]")

FSR = 1559 - 1553; %6nm Dal grafico

Dneff1 = Dneff_V(1);
T_lambda_0 = (1/4).*Ton.*abs(exp(-(2i*pi).*Dneff1.*L./(lambdavect))+exp(-alpha1*DL)*exp((-2i*pi).*(neff*DL)./(lambdavect))).^2;
figure
plot(lambdavect, 20*log10(T_lambda_0))
%lambda_T_max = lambdavect(find(T_lambda == max(T_lambda)));
lambda_T_max = 0.0001556;
Vrange = [0:0.1:6];
DneffV = lambda_T_max .* Dphi_V(Vrange) ./ (L * 2 * pi);
T_lambda_max = (1/4).*Ton.*abs(exp(-(2i*pi).*DneffV.*L./(lambda_T_max))+exp(-alpha1*DL)*exp((-2i*pi).*(neff*DL)./(lambda_T_max))).^2;
figure
plot(Vrange, 20*log10(T_lambda_max))
ylabel("Transmission coefficient [dB]")
xlabel("V_{bias} [V]")

V_pi = 3.3; %V

%ER
int_V_ER = [0.5:0.01:5.5];
VH = int_V_ER + 0.5;
VL = int_V_ER - 0.5;
%DneffV_H = lambda_T_max .* Dphi_V(VH) ./ (L * 2 * pi);
%DneffV_L = lambda_T_max .* Dphi_V(VL) ./ (L * 2 * pi);
%T_H = (1/4).*Ton.*abs(exp(-(2i*pi).*DneffV_H.*L./(lambda_T_max))+exp(-alpha1*DL)*exp((-2i*pi).*(neff*DL)./(lambda_T_max))).^2;
%T_L = (1/4).*Ton.*abs(exp(-(2i*pi).*DneffV_L.*L./(lambda_T_max))+exp(-alpha1*DL)*exp((-2i*pi).*(neff*DL)./(lambda_T_max))).^2;
%ER = abs(T_H./T_L);
for i = 1:numel(int_V_ER)
	DneffV_H = lambda_T_max .* Dphi_V(int_V_ER(i)+0.5) ./ (L * 2 * pi);
	DneffV_L = lambda_T_max .* Dphi_V(int_V_ER(i)-0.5) ./ (L * 2 * pi);
	T_H(i) = (1/4).*Ton.*abs(exp(-(2i*pi).*DneffV_H.*L./(lambda_T_max))+exp(-alpha1*DL)*exp((-2i*pi).*(neff*DL)./(lambda_T_max))).^2;
	T_L(i) = (1/4).*Ton.*abs(exp(-(2i*pi).*DneffV_L.*L./(lambda_T_max))+exp(-alpha1*DL)*exp((-2i*pi).*(neff*DL)./(lambda_T_max))).^2;
	ER(i) = abs(T_H(i)./T_L(i));
end
figure
plot(int_V_ER, abs(10*log10(ER)))
figure
plot(int_V_ER, 10*log10(T_H))
hold on
plot(Vrange, 10*log10(T_lambda_max))
plot(int_V_ER, 10*log10(T_L))
grid on

%Modulation bandwidth
f_vect = linspace(10e7,10e12,100000);
omega = f_vect*2*pi;
alpha_m = 0;
n_m = 3.75;
n_o = 3.86;
Z0 = 37;
c = 3e8;
v_m = c/n_m;
Z_L = 37;
R_G = 50;
R_L = 37;
Z_G = 50;
up = alpha_m*L + 1i*omega*(n_m-n_o)*L./c;
um = -alpha_m*L + 1i*omega*(-n_m-n_o)*L./c;
gamma_m = alpha_m + 1i*omega/v_m;
Z_in = Z0*(Z_L + Z0*tanh(gamma_m*L))/(Z0 + Z_L*tanh(gamma_m*L));
F_up = (1-exp(up))./up;
F_um = (1-exp(um))./um;
m_omega = ((R_L + R_G) / R_L) * abs(Z_in / (Z_in + Z_G)) * abs((((Z_L+Z0) .* F_up) + ((Z_L - Z0) .* F_um)) ./ ((Z_L + Z0) * exp(gamma_m * L) + (Z_L - Z0) * exp(-gamma_m * L)));
velMismatch = (n_m - n_o) * 100 / n_o


velMismatch = 1;
up = alpha_m*L + 1i*omega*velMismatch*n_o*L./(c*100);
n_m = n_o + velMismatch*n_o/100;
um = -alpha_m*L + 1i*omega*(-n_m-n_o)*L./c;
gamma_m = alpha_m + 1i*omega/v_m;
Z_in = Z0*(Z_L + Z0*tanh(gamma_m*L))/(Z0 + Z_L*tanh(gamma_m*L));
F_up = (1-exp(up))./up;
F_um = (1-exp(um))./um;
m_omega_0 = ((R_L + R_G) / R_L) * abs(Z_in / (Z_in + Z_G)) * abs((((Z_L+Z0) .* F_up) + ((Z_L - Z0) .* F_um)) ./ ((Z_L + Z0) * exp(gamma_m * L) + (Z_L - Z0) * exp(-gamma_m * L)));

velMismatch = 10;
up = alpha_m*L + 1i*omega*velMismatch*n_o*L./(c*100);
n_m = n_o + velMismatch*n_o/100;
um = -alpha_m*L + 1i*omega*(-n_m-n_o)*L./c;
gamma_m = alpha_m + 1i*omega/v_m;
Z_in = Z0*(Z_L + Z0*tanh(gamma_m*L))/(Z0 + Z_L*tanh(gamma_m*L));
F_up = (1-exp(up))./up;
F_um = (1-exp(um))./um;
m_omega_10 = ((R_L + R_G) / R_L) * abs(Z_in / (Z_in + Z_G)) * abs((((Z_L+Z0) .* F_up) + ((Z_L - Z0) .* F_um)) ./ ((Z_L + Z0) * exp(gamma_m * L) + (Z_L - Z0) * exp(-gamma_m * L)));

figure
semilogx(f_vect/1e9, 10*log10(m_omega))
hold on
semilogx(f_vect/1e9, 10*log10(m_omega_0))
semilogx(f_vect/1e9, 10*log10(m_omega_10))
grid on
xlabel("f [GHz]")
ylabel("Modulation index [dB]")
title("\alpha_m = 0 - 37 \Omega")
ylim([-20 2])


f_vect = linspace(10e7,10e12,100000);
omega = f_vect*2*pi;
a = 12.02/sqrt(10e9);
alpha_m = a.*sqrt(omega/(2*pi));
n_m = 3.75;
n_o = 3.86;
Z0 = 37;
c = 3e8;
v_m = c/n_m;
Z_L = 37;
R_G = 50;
R_L = 37;
Z_G = 50;
up = alpha_m.*L + 1i*omega*(n_m-n_o)*L./c;
um = -alpha_m.*L + 1i*omega*(-n_m-n_o)*L./c;
gamma_m = alpha_m + 1i*omega/v_m;
Z_in = Z0*(Z_L + Z0*tanh(gamma_m*L))/(Z0 + Z_L*tanh(gamma_m*L));
F_up = (1-exp(up))./up;
F_um = (1-exp(um))./um;
m_omega_alpha = ((R_L + R_G) / R_L) * abs(Z_in / (Z_in + Z_G)) * abs((((Z_L+Z0) .* F_up) + ((Z_L - Z0) .* F_um)) ./ ((Z_L + Z0) * exp(gamma_m * L) + (Z_L - Z0) * exp(-gamma_m * L)));
velMismatch = (n_m - n_o) * 100 / n_o


figure
semilogx(f_vect/1e9, 10*log10(m_omega_alpha))
grid on
xlabel("f [GHz]")
ylabel("Modulation index [dB]")
ylim([-20 2])
title("\alpha_m \div 0 - 37 \Omega")

f_vect = linspace(10e5,10e10,1000000);
omega = f_vect*2*pi;
alpha_m = 0;
n_m = 3.75;
n_o = 3.86;
Z0 = 37;
c = 3e8;
v_m = c/n_m;
%Z_L = 0;
R_G = 50;
%R_L = 0;
Z_G = 50;
up = alpha_m*L + 1i*omega*(n_m-n_o)*L./c;
um = -alpha_m*L + 1i*omega*(-n_m-n_o)*L./c;
gamma_m = alpha_m + 1i*omega/v_m;
Z_in = Z0./tanh(gamma_m*L);
F_up = (1-exp(up))./up;
F_um = (1-exp(um))./um;
m_omega = abs((F_up + F_um) ./ (exp(gamma_m * L) + exp(-gamma_m * L)));
velMismatch = (n_m - n_o) * 100 / n_o

figure
semilogx(f_vect/1e9, 10*log10(m_omega))
grid on
xlabel("f [GHz]")
ylabel("Modulation index [dB]")
title("\alpha_m = 0 - open")
%ylim([-20 2])


f_vect = linspace(10e5,10e10,100000);
omega = f_vect*2*pi;
a = 12.02/sqrt(10e9);
alpha_m = a.*sqrt(omega/(2*pi));
n_m = 3.75;
n_o = 3.86;
Z0 = 37;
c = 3e8;
v_m = c/n_m;
%Z_L = 0;
R_G = 50;
%R_L = 0;
Z_G = 50;
up = alpha_m.*L + 1i*omega*(n_m-n_o)*L./c;
um = -alpha_m.*L + 1i*omega*(-n_m-n_o)*L./c;
gamma_m = alpha_m + 1i*omega/v_m;
Z_in = Z0./tanh(gamma_m*L);
F_up = (1-exp(up))./up;
F_um = (1-exp(um))./um;
m_omega_alpha = abs((F_up + F_um) ./ (exp(gamma_m * L) + exp(-gamma_m * L)));
velMismatch = (n_m - n_o) * 100 / n_o

figure
semilogx(f_vect/1e9, 10*log10(m_omega_alpha))
grid on
xlabel("f [GHz]")
ylabel("Modulation index [dB]")
title("\alpha_m \div 0 - open")
%ylim([-20 2])

f_vect = linspace(10e7,10e12,1000000);
omega = f_vect*2*pi;
alpha_m = 0;
n_m = 3.75;
n_o = 3.86;
Z0 = 37;
c = 3e8;
v_m = c/n_m;
Z_L = 50;
R_G = 50;
R_L = 50;
Z_G = 50;
up = alpha_m*L + 1i*omega*(n_m-n_o)*L./c;
um = -alpha_m*L + 1i*omega*(-n_m-n_o)*L./c;
gamma_m = alpha_m + 1i*omega/v_m;
Z_in = Z0*(Z_L + Z0*tanh(gamma_m*L))/(Z0 + Z_L*tanh(gamma_m*L));
F_up = (1-exp(up))./up;
F_um = (1-exp(um))./um;
m_omega = ((R_L + R_G) / R_L) * abs(Z_in / (Z_in + Z_G)) * abs((((Z_L+Z0) .* F_up) + ((Z_L - Z0) .* F_um)) ./ ((Z_L + Z0) * exp(gamma_m * L) + (Z_L - Z0) * exp(-gamma_m * L)));
velMismatch = (n_m - n_o) * 100 / n_o


velMismatch = 1;
up = alpha_m*L + 1i*omega*velMismatch*n_o*L./(c*100);
n_m = n_o + velMismatch*n_o/100;
um = -alpha_m*L + 1i*omega*(-n_m-n_o)*L./c;
gamma_m = alpha_m + 1i*omega/v_m;
Z_in = Z0*(Z_L + Z0*tanh(gamma_m*L))/(Z0 + Z_L*tanh(gamma_m*L));
F_up = (1-exp(up))./up;
F_um = (1-exp(um))./um;
m_omega_0 = ((R_L + R_G) / R_L) * abs(Z_in / (Z_in + Z_G)) * abs((((Z_L+Z0) .* F_up) + ((Z_L - Z0) .* F_um)) ./ ((Z_L + Z0) * exp(gamma_m * L) + (Z_L - Z0) * exp(-gamma_m * L)));

velMismatch = 10;
up = alpha_m*L + 1i*omega*velMismatch*n_o*L./(c*100);
n_m = n_o + velMismatch*n_o/100;
um = -alpha_m*L + 1i*omega*(-n_m-n_o)*L./c;
gamma_m = alpha_m + 1i*omega/v_m;
Z_in = Z0*(Z_L + Z0*tanh(gamma_m*L))/(Z0 + Z_L*tanh(gamma_m*L));
F_up = (1-exp(up))./up;
F_um = (1-exp(um))./um;
m_omega_10 = ((R_L + R_G) / R_L) * abs(Z_in / (Z_in + Z_G)) * abs((((Z_L+Z0) .* F_up) + ((Z_L - Z0) .* F_um)) ./ ((Z_L + Z0) * exp(gamma_m * L) + (Z_L - Z0) * exp(-gamma_m * L)));

figure
semilogx(f_vect/1e9, 10*log10(m_omega))
hold on
semilogx(f_vect/1e9, 10*log10(m_omega_0))
semilogx(f_vect/1e9, 10*log10(m_omega_10))
grid on
title("\alpha_m = 0 - 50 \Omega")
xlabel("f [GHz]")
ylabel("Modulation index [dB]")
ylim([-20 2])


f_vect = linspace(10e7,10e12,100000);
omega = f_vect*2*pi;
a = 12.02/sqrt(10e9);
alpha_m = a.*sqrt(omega/(2*pi));
n_m = 3.75;
n_o = 3.86;
Z0 = 37;
c = 3e8;
v_m = c/n_m;
Z_L = 50;
R_G = 50;
R_L = 50;
Z_G = 50;
up = alpha_m.*L + 1i*omega*(n_m-n_o)*L./c;
um = -alpha_m.*L + 1i*omega*(-n_m-n_o)*L./c;
gamma_m = alpha_m + 1i*omega/v_m;
Z_in = Z0*(Z_L + Z0*tanh(gamma_m*L))/(Z0 + Z_L*tanh(gamma_m*L));
F_up = (1-exp(up))./up;
F_um = (1-exp(um))./um;
m_omega_alpha = ((R_L + R_G) / R_L) * abs(Z_in / (Z_in + Z_G)) * abs((((Z_L+Z0) .* F_up) + ((Z_L - Z0) .* F_um)) ./ ((Z_L + Z0) * exp(gamma_m * L) + (Z_L - Z0) * exp(-gamma_m * L)));
velMismatch = (n_m - n_o) * 100 / n_o


figure
title("\alpha_m \div 0 - 50 \Omega")
semilogx(f_vect/1e9, 10*log10(m_omega_alpha))
grid on
xlabel("f [GHz]")
ylabel("Modulation index [dB]")
ylim([-20 2])