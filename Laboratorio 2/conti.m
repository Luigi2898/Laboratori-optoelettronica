clc
close all
clear variables

%esercizio 1
%%
r_1=0.3;
r_2=0.9;
La=900e-4;
alpha_m = (1/La)*(log(1/(sqrt(r_1*r_2)))); % rircordati che � cmetri^-1
alpha_i = 6; %unit� di misura cmetri^-1
gamma_per_g_threshold=alpha_m+alpha_i;
N_tr=1.8e18;%elettroni per centrimetro cubo
a_0 =  5.34e-16; %sono cm^2
guadagno_modale = a_0 * gamma_per_g_threshold;
gamma = 0.06;
N_th=N_tr + ((alpha_i+alpha_m)/(gamma * a_0));%elettroni per centrimetro cubo
V = (900e-6)*(100e-10)*(2.5e-6); %metri al cubo
V = V*1e6;%cm al cubo
q = 1.60217e-19;
eta_i = 0.8; % � adimensionato
A = 3.57e8;%secondi alla meno 1
B = 0.8e-10;%cm^3/secondo
C = 3.5e-30; % cm^6/secondo

I_th = ((A*N_th + B*N_th^2 + C*N_th^3))*((q*V)/(eta_i));

%ci manca la linewidth

%esercizio 2

La=linspace(300e-4,1500e-4,10000);
alpha_m = (1./La).*(log(1/(sqrt(r_1*r_2)))); % rircordati che � cmetri^-1
alpha_i = 3; %unit� di misura cmetri^-1
gamma_per_g_threshold=alpha_m+alpha_i;
N_tr=1.8e18;%elettroni per centrimetro cubo
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
ylabel('I_{th} [A]')
