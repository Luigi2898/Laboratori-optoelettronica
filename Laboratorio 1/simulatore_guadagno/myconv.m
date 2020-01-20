function [Econv,gconv,lambda,gconvl]=myconv(htw,g,x,L,Egp,Egpmin,step)

% Funzione che implementa la convoluzione tra il guadagno in materiali quantum well e un'altra funzione (Lorenziana o altro)
% 
% Parametri ingresso:
% 
% htw = ascisse corrispondenti al guadagno [eV]
% g = guadagno [cm^-1]
% x = ascisse corrispondenti alla funzione su cui fare la convoluzione [eV]
% L = funzione da usare per la convoluzione [eV^-1]
% Egp = somma dell'energy gap + minimo dello stato confinato nella banda di conduzione + massimo dello stato confinato nella banda di valenza [eV]
% Egpmin = Eg + levelsc(1) + levelsh(1) [eV]
% step = passo di campionamento;
% 
% Parametri uscita
% 
% Econv = energie - ascisse di gconv
% gconv = guadagno convoluto rispetto all'energia
% lambda = lunghezza d'onda [um]
% gconvl = guadagno convoluto rispetto a lambda

lhtw=length(htw);
lx=length(x);

%Innanzitutto mi accerto che prima che il guadagno assuma valori diversi da
%zero ci siano un numero di zeri pari al numero dei campioni della
%lorenziana
gz=find(htw>Egpmin);
gaux=[zeros(1,lx),g(gz(1):end)];

%Mantenendo lo stesso passo di campionamento, l'energia minima si trova
%partendo da Egp andando indietro di deltaE*length(L)
Emin=Egp-step*lx;

%Ciclo esterno: da 1 alla lunghezza del guadagno meno la lunghezza della
%lorenziana
for k=1:length(gaux)-lx-1
    %ciclo interno: da 1 alla lunghezza della lorenziana
    temp(1:lx)=L(1:lx).*gaux(k:k+lx-1);
    gconv(k)=sum(temp)*step;
    Econv(k)=Emin+step*(k-1+lx/2);
    lambda(k)=1.24/Econv(k);
end

lambda=fliplr(lambda);
gconvl=fliplr(gconv);