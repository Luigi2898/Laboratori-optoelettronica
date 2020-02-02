%clear all;
%close all;
clc;

obj = openfig('esercizio_4_NRZ_10_1024');
data = findobj(obj, 'type', 'line');

tempo = get(data(1),'XData'); %ns
potenza = get(data(1),'YData'); %mW
bitrate =  10e9; %bit/s
numb = 1024;
durata = 1e9 * numb / bitrate; %ns
ini = tempo(numel(tempo)) - durata; %ns
[M I] = min((tempo - ini).^2);
tempo = tempo(I : numel(tempo));
potenza = potenza(I : numel(potenza));
tempo = tempo - tempo(1);

figure(988)
plot(tempo,potenza)

eyediagram(potenza, tempo, numb, bitrate, f)
f = f + 1;
