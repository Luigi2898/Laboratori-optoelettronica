function nw=nquantumwellsing(Econf,Efc,k,T,meff,h)

%Stessa funzione di nquantumwell ma calcola la concentrazione di portatori
%[cm^-2] per ogni singolo livello


sumliv=0;
nb=0;

sumliv=log(1+exp((Efc-Econf)/(k*T)));

% sumliv=sumliv+log(1+exp((Efc-levelsc(1))/(k*T)));

nw=1.602e15 * k * T * meff/(pi * h^2) * sumliv; %cm^-2
