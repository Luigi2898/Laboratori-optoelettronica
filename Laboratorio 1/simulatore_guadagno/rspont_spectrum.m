function   rsp=rspont_spectrum(htw,gmat,Efc,Efv);
% spettro di emissione spontanea in cm^-3*eV^-1*s^-1
% htw: energie in EV
% gmat: gain materiale in cm^-1

% calcolo dello spettro di emissione spontanea nell' ipotesi di quasi eq.


global material
global const_phys



rsp1=4*material.nr^2/(const_phys.hteV.data)*gmat*1e24;
sep_Ef=Efv+Efc+material.Eg;
e_term=1-exp((htw-sep_Ef)/(const_phys.K.data*material.T));
rsp=rsp1./e_term;