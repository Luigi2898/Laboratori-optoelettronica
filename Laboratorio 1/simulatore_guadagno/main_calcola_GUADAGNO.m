% calcola il guadagno di una guida attiva MQW

%Dati per il calcolo delle Efc, Efv e guadagno nel caso del quantum well

% Materiali: InGaAsP lattice matched on InP.  QW: InGaAs lattice matached on InP; barrier: InGaAsP $y=0.39
% Parametri materiale da Libro Chuang, p. 711

%% In base alla buca QW scelta impostare i seguenti parametri di ingresso

material.mstar = 0.041; %mc/m0 nel well: massa efficace elettrone in BC nel well (NB: normalizzata rispetto a massa elettrone libero)
material.m2star = 0.065; %mc/m0 nella barriera:  massa efficace elettrone in BC nella barriera (NB: normalizzata rispetto a massa elettrone libero)

material.mstarhh = 0.465; %mhh/m0 nel well per le lacune
material.m2starhh = 0.46; %mhh/m0 nella barriera per le lacune

material.Eg = 0.75; %eV - energy gap del QW 

material.W = 10; %nm - spessore d del QW

material.T = 300; %K - temperatura

material.N0 = 7; %Numero totale di wells

material.Ebulk = 0.130; %ev - altezza della buca nella banda di conduzione

material.Ebulkh = 0.195; %ev - altezza buca nella banda di valenza