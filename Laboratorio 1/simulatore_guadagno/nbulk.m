function Nb=nbulk(meff,Efc,k,T,h,Einit)

% [int]=approxintegr(Estart,Eend,meff,Efc,k,T)
% 
% int = restituisce il della concentrazione dei portatori nel bulk [cm^-3]
% meff = massa efficace
% Efc = valore che deve essere calcolato dalla fzero
% k = costante di Boltzmann
% T = temperatura
% h = h tagliato
% Einit = valore energia iniziale

% int=0;
% E=Einit*1000; %valore dell'energia in meV
% Evett(1)=0;
% ind=2;
% k=k*1000;  %per avere Boltzmann in meV/K
% maxim=0;
% incr=1;
% NE(1)=0;
% 
% const=0.29e18*(meff/h^2)^(3/2);
% 
% %Calcolo l'integrale finché non ottengo un valore 1000 volte più piccolo
% %del massimo trovato
% while(NE(ind-1)>=(maxim/1000) || (incr==1))
%     NE(ind)=0.1*sqrt(E-Einit*1000)/(1+exp((E-Efc*1000)/(k*T)));
%     int=int+NE(ind);
%     if(NE(ind)>maxim)
%         maxim=NE(ind);
%     end
%     if((NE(ind)<NE(ind-1)) && incr==1)
%        incr=0;
%     end
%     
%     Evett(ind)=E;
%     E=E+0.1;
%     ind=ind+1;
% end
% 
% int=const*int;

global material
global const_phys



% Portatori nella barriera calcolati con approx. di Boltzman: Coldren pag.
% 415

% Nc=2.51e19*(meff/const_phys.m0.data)^(3/2)*(material.T/300)^(3/2); % in cm^-3
% Fapprox=sqrt(pi)/2*exp((Efc-Einit)/(const_phys.K.data*material.T));
% N=Nc*2/sqrt(pi)*Fapprox;

Nb=0;
