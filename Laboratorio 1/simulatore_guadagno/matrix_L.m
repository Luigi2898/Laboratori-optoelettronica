function   ML=matrix_L(htw,tauin)

%ML: matrice delle Lorenziane; su ogni colonna ho una Lorentziana centrata
%in htw

global const_phys

for ie=1:length(htw)
   ML(:,ie)=((const_phys.hteV.data*1e-16/(tauin*pi))./((htw-htw(ie)).^2+(const_phys.hteV.data*1e-16/tauin)^2)).';% in eV^-1
end