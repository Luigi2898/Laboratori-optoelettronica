clear all;
close all;
clc;

obj = openfig('esercizio_4_NRZ_1');
data = findobj(obj, 'type', 'line');

tempo = get(data(1),'XData'); %ns
potenza = get(data(1),'YData'); %mW
bitrate = 1e9; %bit/s
numb = 128;
durata = 10^9 * numb / bitrate; %ns
ini = tempo(numel(tempo)) - durata; %ns
[M I] = min((tempo - ini).^2);
tempo = tempo(I : numel(tempo));
potenza = potenza(I : numel(potenza));
tempo = tempo - tempo(1);

figure(988)
plot(tempo,potenza)

eyediagram(potenza, tempo, numb, bitrate, 11)
