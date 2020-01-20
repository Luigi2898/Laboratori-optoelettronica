function  dN=carrier_RE(N,I,S_satF,S_satR,dz_slice,L,htw,wSOA)

% RE dei portatori di cui si cerca la soluzione stazionaria

% S_satF: spettro fotoni progressivi
% S_satR: spettro fotoni regressivi
% Rsp: spettro di emissione spontanea
% I in mA
% dz_slice: dimensione slice i-esima in um
% L: lunghezza SOA in um
% htw: energie di campionamento dello spettro;

global material
global const_phys

% dato N viene calcolato lo spettro di guadagno ed emissione spontanea


out_gain=gain_spectrum(N,htw); % entrambi in cm^-1
rsp=rspont_spectrum(htw,out_gain.gmat,out_gain.Efc,out_gain.Efv); % spettro di emissione spontanea in cm^-3*eV^-1*s^-1

Rsp=sum(rsp)*(htw(2)-htw(1))*(dz_slice*1e-4*material.W*material.N0*1e-7*wSOA*1e-4); % rate di emissione spontanea totale nella slice dz_slice
S_tot=S_satF+S_satR; % numero totale di fotoni nel volumetto dz_slice 
Rst=const_phys.c.data*1e10/material.nr*material.gammaqw*sum(out_gain.gmat.*S_tot);


dN=I/const_phys.e.data*1e16*dz_slice/L-Rsp-N/(material.tau_nr*1e-9)-Rst;
