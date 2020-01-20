
beta_sponMAT=repmat(material.beta_sp,ns,1).';

% set condizioni di partenza

%inclusione nelle RE del segnale saturante
Nini=Nstz_norm*1e6;
S_sgn_ini=dz_slice*Pin_var(1)/(const_phys.c.data/material.nr)/htw_sgn/const_phys.e.data*100;
[t,Nt_sat]=ode45('carrier_RE_TIME',[0 5e-9],Nini,opt_ode,I,0,S_sgn_ini,dz_slice,L,htw,wSOA,M_Lorentzian,isgn);
gain_ini=gain_spectrum(Nt_sat(end)/(dz_slice*1e-4*wSOA*1e-4),htw,M_Lorentzian); % portatori in cm^-2,htw,M_Lorentzian); % entrambi in cm^-1
%carrier_RE_time(t,N,opt,I,S_noise,S_sgn,dz_slice,L,htw,wSOA,matrix_L,isgn)



for is=1:ns
    slice.gnet(:,is)=gain_ini.gnet;
    slice.rsp(:,is)=gain_ini.rsp; % in  cm^-3*eV^-1*s^-1
    slice.N(is)=Nt_sat(end);
    slice.Efc(is)=gain_ini.Efc;
    slice.Efh(is)=gain_ini.Efh;
end



for  iP=1: length(Pin_var)
    Pin=Pin_var(iP);
    % Propagazione dei fotoni e del segnale
    
    it=1;
    conv_cond=0; % condizione di convergenza su emissione spontanea amplificata


    %condizione al contorno


    SF_sgn_prec(1:ns)=dz_slice*Pin/(const_phys.c.data/material.nr)/htw_sgn/const_phys.e.data*100;
    SF_sgn=SF_sgn_prec;

    % condizioni iniziali
    if iP==1
        SF_noise_prec(1:length(htw),1:ns)=0;
        SR_noise_prec(1:length(htw),1:ns)=0;
    end
       

while (conv_cond==0) 
    Pin 
    it
    % propagazione del segnale
     SF_sgn(2:ns)=SF_sgn_prec(1:ns-1).*exp(slice.gnet(isgn,1:ns-1)*dz_slice*1e-4);
    
     %propagazione emissione spontanea:
     if incl_spon==1
         spon_z(:,1:ns)=0.5*beta_sponMAT.*slice.rsp*wSOA*1e-4*material.N0*material.W*1e-7*(dz_slice*1e-4)^2/(const_phys.c.data/material.nr*1e10)*(htw(2)-htw(1));
         SF_noise(:,1)=spon_z(:,1);
         SF_noise(:,2:ns)=SF_noise_prec(:,1:ns-1).*exp(slice.gnet(:,1:ns-1)*dz_slice*1e-4)+spon_z(:,2:ns);
         SR_noise(:,ns)=spon_z(:,ns);
         SR_noise(:,1:ns-1)=SR_noise_prec(:,2:ns).*exp(slice.gnet(:,2:ns)*dz_slice*1e-4)+spon_z(:,1:ns-1);   
         Stot_noise=SR_noise+SF_noise;
    else
         SR_noise=zeros(1,ns);
         SF_noise=zeros(1,ns);
         Stot_noise=zeros(1,ns);
    end
     
     for iz=1:ns
        guess_Fermi.Efc=slice.Efc(iz);
        guess_Fermi.Efh=slice.Efh(iz);
        newN(iz)=fsolve('carrier_RE_STAZ',slice.N(iz)*1e-6,opt_fsolve,I,Stot_noise(:,iz),SF_sgn(iz),dz_slice,L,htw,wSOA,M_Lorentzian,isgn);
        %carrier_RE_STAZ(N,I,S_noise,S_sgn,dz_slice,L,htw,wSOA,matrix_L,isgn)

        new_gain=gain_spectrum(newN(iz)*1e6/(dz_slice*1e-4*wSOA*1e-4),htw,M_Lorentzian);
        slice.gnet(:,iz)=new_gain.gnet;
        slice.rsp(:,iz)=new_gain.rsp;
        slice.Efc(iz)=new_gain.Efc;
        slice.Efh(iz)=new_gain.Efh;
        slice.N(iz)=newN(iz)*1e6;
        
    end



    SF_sgn_prec=SF_sgn;
    SR_noise_prec=SR_noise;
    SF_noise_prec=SF_noise;
    Pout_noise(it)=sum(SF_noise(:,ns).*(htw.'))/(dz_slice*1e-4)*(const_phys.c.data/material.nr*1e10)*const_phys.e.data*1e-19*1e3; % in mW

    % controllo della condizione di convergenza
    if (it>=20) & (incl_spon==1)
        dP_rel=gradient(Pout_noise(it-10:it))./Pout_noise(it-10:it);
        inc=find(dP_rel>1e-4);
        if isempty(inc)
            conv_cond=1;
        end
    end

    if (it>(ns+1)) & (incl_spon==0)
        conv_cond=1;
    end
            
    it=it+1;
end


    % calcolo della figura di rumore 
    Gz=exp(slice.gnet(isgn,:)*dz_slice*1e-4); % guadagno per ogni singola sotto-sezsione
    nspz=1./(1-exp((htw(isgn)-(slice.Efc+slice.Efh+material.Eg))/(const_phys.K.data*material.T))); % spontaneous emission inversion factor; dipende da z
    Fz=2*nspz.*(Gz-1)./Gz+1./Gz; % formula valida sempre, anche quando non vale Gz>>1!!

   
    Gtot=prod(Gz);
    
    Gz_spectrum=exp(slice.gnet(:,1:ns)*dz_slice*1e-4);
    Gtot_spectrum=prod(Gz_spectrum,2);
    
    if (incl_spon==1) % calcolo della figura di rumore
        Ftot=Fz(1);
        for iz=2:ns
            Ftot=Ftot+(Fz(iz)-1)/(prod(Gz(1:iz-1)));
        end
    else
        Ftot=1;
    end
    out_SOA(iP).slice=slice;
    out_SOA(iP).Gtot=Gtot;
    out_SOA(iP).Gtot_spectrum=Gtot_spectrum;
    out_SOA(iP).F=Ftot;
    out_SOA(iP).SF_noise=SF_noise(:,ns);
    out_SOA(iP).SF_sgn=SF_sgn(ns);

end



