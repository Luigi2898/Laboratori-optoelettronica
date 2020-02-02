%clear all;
%close all;
clc;

obj = openfig('esercizio_4_NRZ_5');
data = findobj(obj, 'type', 'line');

tempo = get(data(1),'XData'); %ns
potenza = get(data(1),'YData'); %mW
bitrate =  5e9; %bit/s
numb = 13;
durata = 1e9 * numb / bitrate; %ns
ini = tempo(numel(tempo)) - durata; %ns
[M I] = min((tempo - ini).^2);
tempo = tempo(I : numel(tempo));
potenza = potenza(I : numel(potenza));
tempo = tempo - tempo(1);

figure(x + 2)
plot(tempo,potenza)

Tbit = 1e9/bitrate; %ns
%inter = ceil(Tbit*1.5/(tempo(2) - tempo(1))); %ns
x = x + 1;
Tdisp = Tbit * 1.5;
n = 0;
for i = 0 : Tbit : tempo(end)
  [M Istart] = min((tempo - i).^2);
  [M Istop] = min((tempo - (Tdisp + i)).^2);
  figure(x)
  plot(tempo(Istart : Istop) - tempo(Istart), potenza(Istart : Istop))
  hold on
  n = n + 1;
end
%for j = 1 : inter : numel(potenza)
%  if (inter + j < numel(potenza))
%    figure(x)
%    plot(tempo(j : inter + j - 1) - tempo(j), potenza(j : inter + j - 1))
%    hold on
%  end
%end
figure(x)
grid on
