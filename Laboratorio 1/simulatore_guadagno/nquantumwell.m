function nw=nquantumwell(levelsc,Efc,k,T,meff,h)



% 
% Funzione che restituisce la concentrazione dei portatori in ogni quantum well.
% 
% Input:
% 
% levelsc: stati confinati
% Efc: [eV]
% k: costante di Boltzmann
% T: temperatura [K]
% meff: massa efficace
% h: costante di Plank
% 
% Output
% 
% nw: concentrazione portatori in ogni well [cm^-2]

sumliv=0;
nb=0;

for(i=1:length(levelsc))
    sumliv=sumliv+log(1+exp((Efc-levelsc(i))/(k*T)));
end


nw=1.602e15 * k * T * meff/(pi * h^2) * sumliv; %cm^-2
