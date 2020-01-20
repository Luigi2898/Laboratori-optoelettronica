function   out_gain=gain_spectrum(N,htw,matrix_L) % guadagno modale netto in cm^{-1}

global material
global const_phys
global guess_Fermi


% out_gain.gmat: guadagno materiale in cm^-1
% out_gain.gnet: gudagno modale netto in cm^-1
% out_gain.rsp: emissione spontanea in in cm^-3*eV^-1*s^-1




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     Calcolo del guadagno
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%Calcolo della Efc e Efh
    opt=optimset('fzero');
    
    out_gain.Efc=fzero('errmqw',guess_Fermi.Efc,opt,N,material.meffel,const_phys.K.data,material.T,const_phys.ht.data,material.levelsc,material.W,material.N0,...
                 material.Ebulk,material.totlen);

   
    out_gain.Efh=fzero('errmqw',guess_Fermi.Efh,opt,N,material.mefflac,const_phys.K.data,material.T,const_phys.ht.data,material.levelsh,...
                material.W,material.N0,material.Ebulkh,material.totlen);

    guess_Fermi.Efc=out_gain.Efc;
    guess_Fermi.Efh=out_gain.Efh;

    
    %Concentrazione portatori nei quantum well [cm^-3]
    Nqw=(nquantumwell(material.levelsc,out_gain.Efc,const_phys.K.data,material.T,material.meffel,const_phys.ht.data))/(material.W*1e-7);
    Hqw=(nquantumwell(material.levelsh,out_gain.Efh,const_phys.K.data,material.T,material.mefflac,const_phys.ht.data))/(material.W*1e-7);   

    %Calcolo del guadagno per ogni stato confinato
    for ilev=1:length(material.levelsc)
        Egp = material.Eg + material.levelsc(ilev) + material.levelsh(ilev);

        
        %Rinormalizzazione dell'energy gap
        %deltaEg=-material.xi*((Nqw+Hqw)/2)^(1/3);
        deltaEg=0;; %=-material.xi*((Nqw+Hqw)/2)^(1/3);
        Egpn=Egp+deltaEg;
        
        %Calcolo del guadagno
        [gconv_lev,rspconv_lev]=mqw_guad(N,htw,material.mstar,...
material.mstarhh,material.T,material.levelsc,material.levelsh,out_gain.Efc,out_gain.Efh,material.W,material.N0,...
material.Ebulk,material.Ebulkh,material.totlen,...
Egp,ilev,material.cp,material.deltaSO,material.tauin,material.nr,1,matrix_L);

        gconv(ilev,:)=gconv_lev;
        rspconv(ilev,:)=rspconv_lev; % in cm^-3*eV^-1*s^-1

end


    
    gconvt=sum(gconv); %Guadagno convoluto funzione dell'energia
    rspt=sum(rspconv);
    out_gain.gmat=gconvt; %Guadagno materiale in cm^-1 in funzione  dell'energia
    out_gain.rsp=rspt;
    
    ebarr=nbulk(material.meffel,out_gain.Efc,const_phys.K.data,material.T,const_phys.ht.data,material.Ebulk); %numero di elettroni nel bulk [cm^-3]
    hbarr=nbulk(material.mefflac,out_gain.Efh,const_phys.K.data,material.T,const_phys.ht.data,material.Ebulkh); %numero di lacune nel bulk [cm^-3]
    
    out_gain.gnet=guadmodnet(material.N0,material.gammaqw,out_gain.gmat,material.alfai,material.gammabarr,material.kn,ebarr,material.kp,hbarr); %Guadagno modale netto [cm^-1]
