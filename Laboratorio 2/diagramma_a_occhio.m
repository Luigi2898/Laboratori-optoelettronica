%La frequenza di taglio a 3dB che ho trovato dalla simulazione è: f_3dB = 4.0981 GHz

%Scommentare la riga della simulazione che si vuole analizzare
%obj = openfig('esercizio_4_NRZ_1_1024');
%obj = openfig('esercizio_4_NRZ_2_5_1024');
obj = openfig('esercizio_4_NRZ_5_1024');
%obj = openfig('esercizio_4_NRZ_10_1024');
%obj = openfig('esercizio_4_NRZ_20_1024');
data = findobj(obj, 'type', 'line');

tempo = get(data(1),'XData'); %ns
potenza = get(data(1),'YData'); %mW

%Cambiare in base alla simulazione scelta
bitrate = 5e9; %bit/s
numb = 1024;
Tbit = 1e9/bitrate; %ns
%Calcolo la durata dello stream di simboli
durata = 1e9 * numb / bitrate; %ns
%Trovo il punto di partenza dello stream nel vettore di tempi
ini = tempo(numel(tempo)) - durata; %ns
%Trovo l'indice del punto di partenza nel vettore di tempi
%(quello che si avvicina di più, non posso trovare l'istante esatto)
[M I] = min((tempo - ini).^2);
%Seleziono la parte del vettore di tempi e delle potenze che contengono lo stream
tempo = tempo(I : numel(tempo));
potenza = potenza(I : numel(potenza));
%Scalo il vettore di tempi così da portare lo zero del tempo al punto di partenza della stream
tempo = tempo - tempo(1);

%Plot dello stream nel tempo
figure(2)
plot(tempo,potenza)
grid on

%Imposto la durata dell'asse delle ascisse del diagramma ad occhio (una volta e mezzo il tempo di bit così da poter visualizzare tutta la transizione, non solo quella iniziale)
Tdisp = Tbit * 1.5;
%Il for si muove sull'asse del tempo dividendolo in multipli del tempo di bit
for i = 0 : Tbit : tempo(end)
  %Trovo come prima l'indice del punto di partenza del simbolo di interesse
  [M Istart] = min((tempo - i).^2);
  %Trovo come prima l'indice un tempo di bit e mezzo dopo l'inizio della trasmissione del simbolo di interesse
  [M Istop] = min((tempo - (Tdisp + i)).^2);
  figure(4)
  %Faccio il plot del diagramma simbolo per simbolo usando come riferimento sulle ascisse il tempo in cui è iniziata la trasmissione del simbolo
  %Per la potenza uso gli stessi indici (sono in parallelo)
  plot(tempo(Istart : Istop) - tempo(Istart), potenza(Istart : Istop))
  hold on
  grid on
end
