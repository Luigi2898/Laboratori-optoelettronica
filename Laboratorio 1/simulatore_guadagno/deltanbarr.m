function [dn]=deltanbarr(lambda,nbarr,nr,meff)

% [dn]=deltanbarr(lambda,nbarr,nr,meff)
% 
% Input: 
% 
% lambda = lunghezza d'onda [um]
% nbarr = numero di portatori nella barriera [cm^-3]
% nr = coefficiente di rifrazione del mezzo
% meff = massa efficace
% 
% output
% 
% dn = deltan nella barriera

define_constants
for i=1:length(lambda)
    dn(i)=(-1e-18*lambda(i)^2*const_phys.e.data^2*nbarr)/(8*pi*const_phys.eps0.data*const_phys.c.data^2*nr*meff);
end