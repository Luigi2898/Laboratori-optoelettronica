function error=errmqw_mg(Efc,N)


global material
global const_phys

% Funzione per determinare gli Efc e Efv tali per cui la concentrazione di portatori N è uguale alla somma dei portatori
% nei quantum well più quelli nei bulk
% Efc: [eV]
% N: [cm^-2]
% meff: massa efficace
% k: costante di Boltzmann
% T: temperatura [K]
% h: costante di Plank h tagliata
% levelsc: stati confinati
% W: larghezza della buca
% N0: numero di well
% Ebulk: profondità buca [eV]
% totlen: lunghezza totale struttura [nm]


%Numero di portatori in ogni quantum well
nw=nquantumwell_mg(Efc); %cm^-2

% %Altro modo di calcolarlo
% nw2=0;
% for(iconf=1:length(levelsc))
%     nw2(iconf,:)=nquantumwellsing(levelsc(iconf),Efc,k,T,meff,h);
% end

%Numero di portatori in tutti i quantum well
ntotwell=nw*material.N0; %cm^-2

%Numero di portatori nei bulk
nb=(nbulk(material.meff,Efc,material.k,material.T,material.h,material.Ebulk))*material.totlen*1e-7; %cm^-2


error=N-ntotwell-nb;


