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
Dneff_V = Dphifit .* lambda./(2*pi);

figure(2)
plot(Vfit, Dneff_V);
title '\Deltan_{eff} vs V'

open('TdB_meas.fig')
D1 = get(gca, 'Children');
Vfig = get(D1, 'XData');
TLfig = get(D1, 'YData');
%FIT CON IL TOOLBOX PER OTTENERE LA CARATTERISTICA T(V)
Tfit = T_Vfit(Vfit);
alpha = Tfit./(4.34*L);

figure(4)
plot(Vfit, alpha);
title 'ALPHA vs V'
%ALPHA DA m^-1 A cm^-1
alpha = alpha./100;



%MZI modulator: transmission versus bias voltage
neff = 3.86;
DL = 1e-4;
reverse_bias = [0:0.1:2];
alpha = T_Vfit(reverse_bias)./(100 * 4.34 * L);
Dneff = Dphi_V(reverse_bias).* lambda./(2*pi);
%SPETTRO DEL COEFF DI TRASMISSIONE >> FARE CICLO SU V

for i = 1:1:numel(alpha)
  alpha1 = alpha(i);
  Dneff1 = Dneff_V(i);
  lambdavect = linspace(1545e-9, 1565e-9, 1000);
  T_lambda = (1/4).*exp(-alpha1*L).*abs(exp(-(1i*2*pi)./(lambdavect).*Dneff1.*L)+exp(-alpha1*L)*exp((-1i*2*pi)./(lambdavect).*(neff*DL)));

  figure(5)
  plot(lambdavect, 10*log10(T_lambda))
  hold on
  grid on
end
alpha1 = alpha(1);
Dneff1 = Dneff_V(1);
T_lambda = (1/4).*exp(-alpha1*L).*abs(exp(-(1i*2*pi)./(lambdavect).*Dneff1.*L)+exp(-alpha1*L)*exp((-1i*2*pi)./(lambdavect).*(neff*DL)));



lambdain = lambdavect(find(T_lambda == max(T_lambda)));
T_lambda = (1/4).*exp(-alpha1.*L).*abs(exp(-(1i*2*pi)./(lambdain).*Dneff_V.*L)+exp(-alpha.*L)*exp((-1i*2*pi)./(lambdain).*(neff*DL)));

figure(6)
plot(Vfit,10*log10(T_lambda));
