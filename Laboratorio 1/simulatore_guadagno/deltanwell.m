function [lambda,nconvl]=deltanwell(htw,g,tauin,Egp,Egpmin,step)

% [lambda,nconvl]=deltanwell(htw,g,tauin,Egp,Egpmin,step)
% 
% Input: 
% 
% htw: vettore energie [eV]
% g: guadagno materiale non convoluto [cm^-1]
% tauin: tempo di rilassamento [sec]
% Egp = somma dell'energy gap + minimo dello stato confinato nella banda di conduzione + massimo dello stato confinato nella banda di valenza [eV]
% Egpmin = somma dell'energy gap + minimo del primo stato confinato nella banda di conduzione + massimo del primo stato confinato nella banda di valenza [eV]
% step = passo di campionamento
% 
% Output
% 
% lambda = lunghezza d'onda [um]
% nconvl = deltan nei well


define_constants
for ihtw=1:length(htw)
    w(ihtw)=htw(ihtw)/(const_phys.hteV.data*1e-16);
    dnmat(ihtw)=1e10*g(ihtw)*const_phys.c.data/(2*w(ihtw));
end

x=htw(ihtw)-10*const_phys.hteV.data*1e-16/tauin:step:htw(ihtw)+const_phys.hteV.data*1e-16/tauin*10;
for count=1:length(x)
    L(count)=((htw(ihtw)-x(count))/pi)/((x(count)-htw(ihtw))^2+(const_phys.hteV.data*1e-16/tauin)^2);
end


[Econv,nconv,lambda,nconvl]=myconv(htw,dnmat,x,L,Egp,Egpmin,step);