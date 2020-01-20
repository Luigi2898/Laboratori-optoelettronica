function  dN_rel=carrier_RE_STAZ(N,I,S_noise,S_sgn,dz_slice,L,htw,wSOA,matrix_L,isgn)

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
global guess_Fermi

% dato N viene calcolato lo spettro di guadagno ed emissione spontanea


Ndens=N*1e6/(dz_slice*1e-4*wSOA*1e-4); % portatori in cm^-2

out_gain=gain_spectrum(Ndens,htw,matrix_L); % entrambi in cm^-1




% figure(1)
% plot(htw,out_gain.gnet)
% figure(2)
% plot(htw,out_gain.rsp)
% pause

Rsp=sum(out_gain.rsp)*(htw(2)-htw(1))*(dz_slice*1e-4*material.W*material.N0*1e-7*wSOA*1e-4); % rate di emissione spontanea totale nella slice dz_slice
Rst_noise=const_phys.c.data*1e10/material.nr*material.gammaqw*sum(out_gain.gmat.*(S_noise.'));
Rst_sgn=const_phys.c.data*1e10/material.nr*material.gammaqw*(out_gain.gmat(isgn)*S_sgn);





% Rsp=Rsp/10
% N
% I/const_phys.e.data*1e16*dz_slice/L
% Rsp
% N/(material.tau_nr*1e-9)
% pause

% I/const_phys.e.data*1e16*dz_slice/L
% Rsp
% N/(material.tau_nr*1e-9)
% Rst_sgn
% pause

dN_rel=(I/const_phys.e.data*1e16*dz_slice/L-Rsp-N/(material.tau_nr*1e-9)-Rst_noise-Rst_sgn)*1e-9/(N*1e6);

