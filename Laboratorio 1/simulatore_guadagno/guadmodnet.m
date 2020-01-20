function [gmn]=guadmodnet(N0,gammaqw,gmat,alfai,gammabarr,kn,n,kp,p)

% [gmn]=guadmodnet(N0,gammaqw,gmat,alfai,gammabarr,kn,n,kp,p)
% Funzione che calcola il guadagno modale netto.
% 
% Parametri ingresso:
% 
% N0 = numero quantum well
% gammaqw = fattore di confinamento per ogni singolo quantum well
% gmat = guadagno calcolato nel materiale [cm^-1]
% alfai = 7 cm^-1
% gammabar = fattore di confinamento nella barriera
% kn = coefficiente IVBA per elettroni [cm^2]
% n = elettroni per unità di volume nella barriera [cm^-3]
% kp = coefficiente IVBA per lacune [cm^2]
% p = lacune per unità di volume nella barriera [cm^-3]
% 
% Parametri uscita:
% 
% gmn = guadagno modale netto

gammaqwtot=N0*gammaqw;

alfamatbarr=kn*n+kp*p;

gmn=gammaqwtot.*gmat-alfai-gammabarr*alfamatbarr;