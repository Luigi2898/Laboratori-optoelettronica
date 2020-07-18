%% SCRIPT per il calcolo di SMSR @20°C
%Ratio between the power of the lasing  mode and the power of the highest power longitudinal mode below threshold
%%% SMSR=Plaserante/P(primo modo più alto sotto soglia)
I=[1:40];
SMSR=[0 65.1-64.4 61.69-60.62 58.63-57 55.5-53.27 52.51-49.36 ...
      49.49-44.96 46.13-39.26 42.87-31.44 41.6-18.34 42.17-13.83 ...
      42.58-11.88 42.92-10.45 43.21-9.37 43.45-8.59 43.66-7.94 ...
      43.86-7.31 44.03-6.76 44.16-6.28 44.29-5.9 44.42-5.52 ...
      44.53-5.15 44.62-4.84 44.7-4.58 44.8-4.33 44.86-4.01 ...
      44.92-3.75 44.99-3.49 45.05-3.32 45.12-3.08 45.16-2.88 ...
      45.21-2.75 45.27-2.53 45.3-2.34 45.35-2.18 45.36-2.07 ...
      45.44-1.86 45.44-1.71 45.49-1.62 45.55-1.45];

figure(1)
plot(I, SMSR)
hold on;
grid on;
xlabel('I [mA]');
ylabel('SMSR [dBm]');
%title('SMSR function of variable stimulus current at T=20°C');

%% Variazione spettri a diverse temperature (sopra soglia)

Istart=20; %mA
Istop=30; %mA
Istep=1; %mA
Icurrent=[Istart:Istep:Istop];
Mmeas_20=load('OSA_Temp20.000000.txt');
Mmeas_25=load('OSA_Temp25.000000.txt');
Mmeas_30=load('OSA_Temp30.000000.txt');
legendInfo=cell(3,1);

legendInfo{1}=sprintf('T=20°C');
legendInfo{2}=sprintf('T=25°C');
legendInfo{3}=sprintf('T=30°C');


col1=[0.54, 0.780, 1];
col2=[0.16, 0.56, 0.94];
col3=[0.02, 0.31, 0.58];

figure(2)
xlabel('\lambda [nm]');
ylabel('Potenza [dBm]');
grid on;
hold on
num_punti=800; % numero di lunghezze d'onda a cui è stato misurato lo spettro ottico ad una corrente fissata
for icurr=Istart:(Istart+length(Icurrent))
    Spettro_20=plot(Mmeas_20((icurr-1)*num_punti+1:icurr*num_punti,1),Mmeas_20((icurr-1)*num_punti+1:icurr*num_punti,2), 'Color', col1);
end
hold on
for icurr=Istart:(Istart+length(Icurrent))
    Spettro_25=plot(Mmeas_25((icurr-1)*num_punti+1:icurr*num_punti,1),Mmeas_25((icurr-1)*num_punti+1:icurr*num_punti,2), 'Color', col2);
end
hold on
for icurr=Istart:(Istart+length(Icurrent))
    Spettro_30=plot(Mmeas_30((icurr-1)*num_punti+1:icurr*num_punti,1),Mmeas_30((icurr-1)*num_punti+1:icurr*num_punti,2), 'Color', col3);
end

legend([Spettro_20,Spettro_25,Spettro_30],legendInfo);