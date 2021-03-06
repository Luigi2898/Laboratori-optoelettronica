clc
close all
clear variables

%esercizio 1

r_1=0.3;
r_2=0.9;
La=900e-4;
alpha_m_fix = (1/La)*(log(1/(sqrt(r_1*r_2)))); % rircordati che � cmetri^-1
alpha_i_fix = 6; %unit� di misura cmetri^-1
gamma_per_g_threshold_fix = alpha_m_fix + alpha_i_fix;
N_tr = 1.8e18;%elettroni per centrimetro cubo
a_0 =  5.34e-16; %sono cm^2
guadagno_modale = a_0 * gamma_per_g_threshold_fix;
gamma = 0.06;
N_th_fix=N_tr + ((alpha_i_fix + alpha_m_fix) / (gamma * a_0));%elettroni per centrimetro cubo
V_fix = (900e-6)*(100e-10)*(2.5e-6); %metri al cubo
V_fix = V_fix*1e6;%cm al cubo
q = 1.60217e-19;
eta_i = 0.8; % � adimensionato
A = 3.57e8;%secondi alla meno 1
B = 0.8e-10;%cm^3/secondo
C = 3.5e-30; % cm^6/secondo
n_g = 4.2;
c = 3e8;
v_g = c/n_g;
I_th = ((A*N_th_fix + B*N_th_fix^2 + C*N_th_fix^3))*((q*V_fix)/(eta_i));

n_g = 4.2;
c = 3e8;
v_g = c/n_g;
I = 0.010;
beta_sp = 1e-4;
g_th = gamma_per_g_threshold_fix/gamma;
Np = eta_i*(I-I_th)/(q*V_fix*v_g*g_th);
Rsp_p = gamma*beta_sp*B*N_th_fix^2;
alpha_h = 3;
linewidth = gamma*Rsp_p*(1+alpha_h)^2/(4*pi*Np)

tau_p = 1/(v_g * (alpha_i_fix * 10^2 + alpha_m_fix * 10^2));
beta_sp = 1e-4;
I = 6e-3;
Np = eta_i*(I - I_th)/(q * V_fix * v_g * gamma_per_g_threshold_fix/gamma);
delta_nu = ((gamma^2*beta_sp)/(4*pi*Np));

%ci manca la linewidth

%esercizio 2

La=linspace(300e-4,1500e-4,10000);
alpha_m = (1./La).*(log(1/(sqrt(r_1*r_2)))); % rircordati che � cmetri^-1
alpha_i = 3; %unit� di misura cmetri^-1
gamma_per_g_threshold=alpha_m+alpha_i;
N_tr = 1.8e18;%elettroni per centrimetro cubo
a_0 =  1e-16; %sono cm^2
guadagno_modale = a_0 .* gamma_per_g_threshold;
gamma = 0.06;
N_th=N_tr + ((alpha_i+alpha_m)./(gamma.* a_0));%elettroni per centrimetro cubo
V = (La).*(100e-8)*(2.5e-4); %Cmetri al cubo
%V = V*1e6;%cm al cubo
q = 1.60217e-19;
eta_i = 0.8; % � adimensionato
A = 3.57e8;%secondi alla meno 1
B = 0.8e-10;%cm^3/secondo
C = 3.5e-30; % cm^6/secondo

I_th = ((A.*N_th + B.*N_th.^2 + C.*N_th.^3)).*((q*V)./(eta_i));
figure(1)
plot(La,I_th,'r')
grid on
title('Plot della corrente di threshold in funzione della lunghezza')
xlabel('La [cm]')
ylabel('I_{th} [A]')

eta_d = eta_i .* alpha_m ./(alpha_i + alpha_m );

figure(2)
plot(La,eta_d,'r')
grid on
title('Plot efficienza differenziale in funzione della lunghezza')
xlabel('La [cm]')
ylabel('\eta_{d}')

%Esercizio 3
%riportare sullo stesso grafico le diverse risposte

h1_05 = openfig('esercizio_3_punto_6__1_05_IM.fig');
t1_05 = findobj(h1_05,'type', 'line');

x1_05 = get(t1_05,'XData');
y1_05 = get(t1_05,'YData');

h2 = openfig('esercizio_3_punto_6__2_IM.fig');
t2 = findobj(h2,'type', 'line');

x2 = get(t2,'XData');
y2 = get(t2,'YData');

h6 = openfig('esercizio_3_punto_6__6_IM.fig');
t6 = findobj(h6,'type', 'line');

x6 = get(t6,'XData');
y6 = get(t6,'YData');

h10 = openfig('esercizio_3_punto_6__10_IM.fig');
t10 = findobj(h10,'type', 'line');

x10 = get(t10,'XData');
y10 = get(t10,'YData');
close all
figure(3)
semilogx(x1_05, y1_05, 'b')
grid on
hold on
semilogx(x2, y2, 'r')
semilogx(x6, y6, 'g')
semilogx(x10, y10, 'y')
legend('I_{bias}/I_{th} = 1.05','I_{bias}/I_{th} = 2','I_{bias}/I_{th} = 6','I_{bias}/I_{th} = 10')
title('Confronto tra le risposte IM per diverse correnti di bias')
ylabel('Risposta IM [dB]')
xlabel('Frequenza [GHz]')

%banda come funzione della corrente di bias
f1_05_3dB = x1_05(20);
f2_3dB = x2(89);
f6_3dB = x6(186);
f10_3dB = x10(232);
I_th_fix = 9.677; %mA
Ibias = [I_th_fix*1.05 I_th_fix*2 I_th_fix*6 I_th_fix*10];
f_3dB = [f1_05_3dB f2_3dB f6_3dB f10_3dB];

%Banda di modulazione massima
n_g = 4.2;
c = 3e8;
v_g = c/n_g;
tau_p = 1/(v_g * (alpha_i_fix * 10^2 + alpha_m_fix * 10^2));
banda_max_ideale = sqrt(2)/(tau_p * 2 * pi);

one_ov_tau_Dn = 2 * B * N_th_fix + A + 3 * C * N_th_fix;
Np = eta_i*(Ibias - I_th_fix)./(q * V_fix * v_g * gamma_per_g_threshold_fix/gamma); %FIXME Calcolare la densità di fotoni come funzione della corrente di bias
damping = v_g * a_0 * Np + one_ov_tau_Dn;
wr_square = (v_g*a_0*Np)./tau_p;
wp_square = wr_square*(1-1/2*(damping/sqrt(wr_square)));
w3dB = sqrt(wp_square + sqrt(wp_square.^2 + wr_square.^2));

figure(4)
plot(Ibias, f_3dB, '-o')
grid on
%hold on
%plot(Ibias, w3dB./(2*pi)./1e9)

%esercizio 3.4 (ricalcolo con epsilon = 1.5e-14)
epsilon = 1.5e-14;
h1_05 = openfig('esercizio_3_punto_4__1_05_IM.fig');
t1_05 = findobj(h1_05,'type', 'line');

x1_05 = get(t1_05,'XData');
y1_05 = get(t1_05,'YData');

h2 = openfig('esercizio_3_punto_4__2_IM.fig');
t2 = findobj(h2,'type', 'line');

x2 = get(t2,'XData');
y2 = get(t2,'YData');

h6 = openfig('esercizio_3_punto_4__6_IM.fig');
t6 = findobj(h6,'type', 'line');

x6 = get(t6,'XData');
y6 = get(t6,'YData');

h10 = openfig('esercizio_3_punto_4__10_IM.fig');
t10 = findobj(h10,'type', 'line');

x10 = get(t10,'XData');
y10 = get(t10,'YData');

figure(5)
semilogx(x1_05, y1_05, 'b')
grid on
hold on
semilogx(x2, y2, 'r')
semilogx(x6, y6, 'g')
semilogx(x10, y10, 'y')
legend('I_{bias}/I_{th} = 1.05','I_{bias}/I_{th} = 2','I_{bias}/I_{th} = 6','I_{bias}/I_{th} = 10')
title('Confronto tra le risposte IM per diverse correnti di bias (\epsilon = 1.5e-14)')
ylabel('Risposta IM [dB]')
xlabel('Frequenza [GHz]')

%banda come funzione della corrente di bias
f1_05_3dB_eps = x1_05(6);
f2_3dB_eps = x2(8);
f6_3dB_eps = x6(14);
f10_3dB_eps = x10(20);
I_th_fix_eps = 9.6672; %mA
Ibias_eps = [I_th_fix_eps*1.05 I_th_fix_eps*2 I_th_fix_eps*6 I_th_fix_eps*10];
f_3dB_eps = [f1_05_3dB_eps f2_3dB_eps f6_3dB_eps f10_3dB_eps];

figure(6)
plot(Ibias_eps, f_3dB_eps, '-o')
grid on

%Banda di modulazione massima
%FIXME da rifare con epsilon = 1.5e-14
n_g = 4.2;
c = 3e8;
v_g = c/n_g;
tau_p = 1/(v_g * (alpha_i_fix * 10^2 + alpha_m_fix * 10^2));
banda_max_eps = sqrt(2)*(1/(1+(gamma * epsilon/a_0)))/(tau_p * 2 * pi);

%esercizio 5
P_bias = 5; %mW
I_bias = 22.2; %mA 4.04 * i_th
P_0 = 1; %mW
P_1 = 10; %mW
I_0 = 8.7; %mA 1.58 * I_th
I_1 = 36.6; %mA 6.65 * I_th
I_th_NRZ = 5.5001; %mA
f_3dB = 4.0981; %GHz
