function [gconv,rsp_conv]=mqw_guad(N,htw,mstar,mstarhh,T,levelsc,levelsh,Efc,Efh,W,N0,Ebulk,Ebulkh,totlen,Eg,ilev,cp,delta,tauin,nr,flag,M_Ltz)

global const_phys
global material

% Funzione che calcola il guadagno in strutture multi quantum well.
% 
% 
% Parametri d'ingresso:
% 
% N = numero dei portatori [cm^-2]
% htw = vettore con le energie su cui calcolare il guadagno
% mstar = mc/m0 nel well
% mstarhh = mhh/m0 nel well
% T = temperatura [K]
% levelsc = stati confinati banda di conduzione
% levelsh = stati confinati banda di valenza
% Efc = quasi livello di Fermi nella banda di conduzione
% Efh = quasi livello di Fermi nella banda di valenza
% W = larghezza del well [nm]
% N0 = numero dei well
% Ebulk = altezza della buca nella banda di conduzione [eV]
% Ebulkh = altezza buca nella banda di valenza [eV]
% totlen = lunghezza totale della barriera [nm]
% Eg = energy gap [eV]
% Egp = somma dell'energy gap + minimo dello stato confinato nella banda di conduzione + massimo dello stato confinato nella banda di valenza [eV]
% cp = costante di polarizzazione
% delta = [eV]
% tauin = tempo di rilassamento [sec]
% nr = indice di rifrazione del mezzo
% flag = quando = 1 la convoluzione viene effettuata con la lorenziana, quando = 2 si usa la secante iperbolica 
% 
% Parametri uscita:
% 
% gnc = valore del guadagno senza la convoluzione funzione dell'energia
% gconv = valore del guadagno con la convoluzione funzione dell'energia [cm^-1]


pl=0; %flag per decidere se plottare la curva con cui fare la convoluzione (=1) o no (=0)
define_constants;
mrstar=(mstar*mstarhh)/(mstar+mstarhh);

Egpmin = Eg;% + levelsc(ilev) + levelsh(ilev);

meffel=mstar*const_phys.m0.data;
meffh=mstarhh*const_phys.m0.data;
mr=((mstar*mstarhh)/(mstar+mstarhh))*const_phys.m0.data;



% NB: i livelli di Fermi per il calcolo delle funzioni di Fermi sono
% riferiti rispetto al fondo della banda di conduzione ed al fondo della
% banda di valenza

Ec = (mrstar/mstar) * (htw - Egpmin)+levelsc(ilev);
Eh = (mrstar/mstarhh) * (htw - Egpmin)+levelsh(ilev);  




Cg=1/nr*(1e11/1.602)*(1./htw)*(pi*const_phys.e.data^2*const_phys.ht.data/(const_phys.eps0.data*const_phys.c.data*const_phys.m0.data^2)); %m^2/kg
    
Mtsq=2.1e-48/10; %kg*J - da tabella II pag.49 libro quantum well lasers, GaAS trasformato in kg*J


igp=find(htw>=Egpmin);
rho(1:igp(1))=0; % per energie minori di Egp
rho(igp(1)+1:length(htw))=1e47*mr/(pi*W*const_phys.ht.data^2); %s^2/(kg*m^5)
fc=1./(1+exp((Ec-Efc)/(const_phys.K.data*T)));
fvb=1./(1+exp((Eh-Efh)/(const_phys.K.data*T)));


gnc=1e-2*Cg.*Mtsq.*rho.*(fc+fvb-1); %cm^-1
ll=1.24./htw; % lunghezza d'onda in um

rsp1=4*material.nr^2/(const_phys.hteV.data)*gnc*1e24./(ll.^2);
sep_Ef=Efh+Efc+material.Eg;
e_term=1-exp((htw-sep_Ef)/(const_phys.K.data*material.T));
rspnc=rsp1./e_term;   % in 
  
% figure(1)
% plot(htw,gnc,'r')
% figure(2)
% plot(htw,rspnc,'r')
% pause


%Convoluzione con la Lorentziana


gconv=gnc(1:end-1)*M_Ltz(1:end-1,1:end-1)*(htw(2)-htw(1)); % in cm^-1
rsp_conv=rspnc(1:end-1)*M_Ltz(1:end-1,1:end-1)*(htw(2)-htw(1)); 
gsgn=interp1(htw(1:end-1),gconv,htw(end),'linear');
rsp_sgn=interp1(htw(1:end-1),rsp_conv,htw(end),'linear');
if isfinite(gsgn)
    gconv=[gconv,gsgn];
    rsp_conv=[rsp_conv,rsp_sgn];
end

