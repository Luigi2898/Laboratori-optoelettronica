clear all
close all
clc

%TODO plot diagramma a bande

m_e = 4.385358e-32; %kg
m_h = 3.9351705e-31; %kg
h_t = 6.582199e-16; %eV*s
U1e = 0; %eV
U2e = 0.162; %eV
Ee = linspace(U1e, U2e, 1000000);
gammae = (sqrt(2*m_e)./h_t).*sqrt(U2e - Ee);
deltae = (sqrt(2*m_e)./h_t).*sqrt(U2e - Ee);
k_xe = (sqrt(2*m_e)./h_t).*sqrt(Ee - U1e);
d0e = (1./k_xe).*(atan(gammae./k_xe)+ atan(deltae./k_xe)); %nm
d_TE = (atan(gammae./k_xe) + atan(deltae./k_xe) + pi)./(k_xe);
d1e = (1./k_xe).*(atan(gammae./k_xe)+ atan(deltae./k_xe)+pi); %nm


U1h = 0.108; %eV
U2h = 0; %eV
Eh = linspace(U2h, U1h, 1000000);
gammah = (sqrt(2*m_h)./h_t).*sqrt(Eh - U2h);
deltah = (sqrt(2*m_h)./h_t).*sqrt(Eh - U2h);
k_xh = (sqrt(2*m_h)./h_t).*sqrt(U1h - Eh);
d0h = (1./k_xh).*(atan(gammah./k_xh)+ atan(deltah./k_xh)); %nm
d1h = (1./k_xh).*(atan(gammah./k_xh)+ atan(deltah./k_xh)+pi); %nm

x = 20 * ones(1, 50);
y = linspace(0, 0.18, 50);

figure
plot(d0e,Ee,'r',d1e,Ee,'b')
title('E vs d - ELETTRONI')
grid on
xlim([0 100])

figure
plot(d0h,Eh,'r',d1h,Eh,'b')
title('E vs d - LACUNE')
grid on
xlim([0 100])

figure
plot(d0h,Eh,'r',d1h,Eh,'b', d0e,Ee, 'g', d1e, Ee, 'm', x, y, '.-')
title('E vs d - ELETTRONI e LACUNE')
grid on
xlim([0 100])

E0e = 0.048556; %eV
E1e = 0.155420; %eV
E1h = 0.072684; %eV
E0h = 0.099008; %eV

lambda_ric = h_t*2*pi*3e8/(E0e+E0h+0.775437);
