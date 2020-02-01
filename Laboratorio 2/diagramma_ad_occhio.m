%clear all;
%close all;
clc;

obj = openfig('esercizio_4_NRZ_1');
data = findobj(obj, 'type', 'line');

tempo = get(data(1),'XData');
potenza = get(data(1),'YData');
tempo = tempo(221 : numel(tempo));
potenza = potenza(221 : numel(potenza));

figure(988)
plot(tempo,potenza)

eyediagram(potenza, 128, 11)
