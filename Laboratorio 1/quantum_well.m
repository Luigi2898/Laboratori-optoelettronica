clear all
close all
clc

m_e = 4.385358e-32;
m_h = 3.9351705e-31;
h_t = 6.582199e-16*1.60217e-19;
U1e = 0;
U2e = 0.162*1.60217e-19;
Ee = linspace(U1e, U2e, 10000);
gammae = (sqrt(2*m_e)./h_t).*sqrt(U2e - Ee);
deltae = (sqrt(2*m_e)./h_t).*sqrt(U2e - Ee);
k_xe = (sqrt(2*m_e)./h_t).*sqrt(Ee - U1e);
d0e = (1./k_xe).*(atan(gammae./k_xe)+ atan(deltae./k_xe));
d1e = (1./k_xe).*(atan(gammae./k_xe)+ atan(deltae./k_xe)+pi);


U1h = (-0.7-0.108)*1.60217e-19;
U2h = -0.7 * 1.60217e-19;
Eh = linspace(U1h, U2h, 10000);
gammah = (sqrt(2*m_h)./h_t).*sqrt(U2h - Eh);
deltah = (sqrt(2*m_h)./h_t).*sqrt(U2h - Eh);
k_xh = (sqrt(2*m_h)./h_t).*sqrt(Eh - U1h);
d0h = (1./k_xh).*(atan(gammah./k_xh)+ atan(deltah./k_xh));
d1h = (1./k_xh).*(atan(gammah./k_xh)+ atan(deltah./k_xh)+pi);

figure(1)
plot(d0e,Ee,'r',d1e,Ee,'b')
title('E vs d - ELETTRONI')
grid on

figure(2)
plot(d0h,Eh,'r',d1h,Eh,'b')
title('E vs d - LACUNE')
grid on

d = 10e-9;
E1e = 2.0875e-20;
E0e = 5.8665e-21;
E1h = -1.2553e-19;
E0h = -1.2847e-19;

lambda_ric = h_t*2*pi*3e8/(E0e-E0h);