%PER T A 20 GRADI
I_bias_laser_diode_controller_T20=[2 4 6 8 10 12 14 16 18 20];%mA prendiamo da 2mA a 20 mA con passo 2mA
P_analizzatore_di_spettro_T20=[0.1872e-6 0.6495e-6 1.5457e-6 3.090e-6 33.83e-6 73.59e-6 115.89e-6 159.91e-6 0.2055e-3 0.2524e-3];% W
I_TE_T20=[0.062 0.063 0.064 0.065 0.066 0.067 0.068 0.070 0.0.072 0.073];%I_TE in Ampere

figure('Name','P-I_bias')
plot(P_analizzatore_di_spettro_T20,I_bias_laser_diode_controller_T20);

figure('Name','I_TE-I_bias')
plot(I_TE_T20,I_bias_laser_diode_controller_T20);

%PER T A 40 GRADI
I_bias_laser_diode_controller_T40=[20 18 16 14 12 10 8 6 4 2];%mA prendiamo da 2mA a 20 mA con passo 2mA
P_analizzatore_di_spettro_T40=[68.14e-6 42.89e-6 19,72e-6 3.60e-6 2.56e-6 1.120e-6 0.643e-6 0.314e-6 109.28e-9];% W
I_TE_T40=[-0.116 -0.116 -0.116 -0.117 -0.118 -0.120 -0.121 -0.122 -0.123];%I_TE in Ampere

%Partiamo a T=20 e poin rifaccimo le misure senza il controllo sulla
%temperatura senza controllo
I_bias_laser_diode_controller_free=[2 4 6 8 10 12 14 16 18 20];%mA prendiamo da 2mA a 20 mA con passo 2mA
P_analizzatore_di_spettro_free=[0.1560e-6 0.5024e-6 1.1325e-6 2.102e-6 3.804e-6 33.75e-6 67.89e-6 104.22e-6 142.14e-6 180.50e-6];% W
T=[27.3 27.2 27.1 27.1 27.2 27.2 27.3 27.4 27.5 27.6]%gradi centigradi


