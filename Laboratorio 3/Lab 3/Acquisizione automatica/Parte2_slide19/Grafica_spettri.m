% grafica spettri
% imposta le correnti a cui sono stati misurati gli spettri e che si
% vogliono riportare sul grafico
Istart=1% mA
Istop=40% mA
Istep=1; %mA
Icurrent=[Istart:Istep:Istop];
Mmeas=load('OSA_Temp20.000000.txt'); % carica il file contenente le misure di spettri ottici

%%

figure(1)
hold on
num_punti=800; % numero di lunghezze d'onda a cui è stato misurato lo spettro ottico ad una corrente fissata
legendInfo=cell(40,1); % creo cella per salvare i valori di corrente per la legenda
for icurr=1:length(Icurrent)
    plot(Mmeas((icurr-1)*num_punti+1:icurr*num_punti,1),Mmeas((icurr-1)*num_punti+1:icurr*num_punti,2))
    %pause  
    legendInfo{icurr}=sprintf('I=%d mA',icurr);
end
legend(legendInfo);
    