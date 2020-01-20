% soluzione RE in assenza di saturazione: risolve Eq. diff portatoti nel
% dominio del tempo con ODE45

I=[0.1:0.2:20];

L=500;
wSOA=4;


% Vettore energia htw
step=2e-3;
htw=material.Egp-0.01:step:material.Egp2+0.2;

ns=10;
dz_slice=L/ns;


J=I/(wSOA*1e-4*L*1e-4)*1e-6; % densità di corrente in kA/cm^2

% Matrice LORENTZIANe per la convoluzione

M_Lorentzian=matrix_L(htw,material.tauin);



Nini=1e15*dz_slice*1e-4*wSOA*1e-4*material.N0*material.W*1e-7;

opt_ode=odeset('Reltol',1e-5);
guess_Fermi.Efc=-material.Eg/2;
guess_Fermi.Efh=-material.Eg/2;

ll=1.24./(htw(1:end-1));
material.beta_sp=1./(2*pi^2*material.nr^2).*(ll*1e-3).^2*material.gammaqw*material.N0/(wSOA*material.N0*material.W*1e-3);
[t,Nt]=ode45('carrier_RE_TIME',[0 5e-9],Nini,opt_ode,I(1),0,0,dz_slice,L,htw,wSOA,M_Lorentzian,1);
Nstz(1)=Nt(end);



%calcolo dello spettro di guadagno ed emissione spontanea

out_gain(1)=gain_spectrum(Nstz/(dz_slice*1e-4*wSOA*1e-4),htw,M_Lorentzian); % portatori in cm^-2,htw,M_Lorentzian); % entrambi in cm^-1


opt_fsolve=optimset('fsolve');
%opt_fzero=optimset('fzero');
guess_Fermi.Efc=out_gain(1).Efc;
guess_Fermi.Efh=out_gain(1).Efh;

% ciclo per calcolo del guadagno al variare della corrente
for icurr=1:length(I)
        if (icurr==1)
            Nini=Nstz(1);
        else
            Nini=Nstz(icurr-1);
        end
        Nstz_norm=fsolve('carrier_RE_STAZ',Nini*1e-6,opt_fsolve,I(icurr),0,0,dz_slice,L,htw,wSOA,M_Lorentzian,1);   
        Nstz(icurr)=Nstz_norm*1e6;
        out_gain(icurr)=gain_spectrum(Nstz(icurr)/(dz_slice*1e-4*wSOA*1e-4),htw,M_Lorentzian); % portatori in cm^-2,htw,M_Lorentzian); % entrambi in cm^-1


        guess_Fermi.Efc=out_gain(icurr).Efc;
        guess_Fermi.Efh=out_gain(icurr).Efh;


        ampl_term=(exp(out_gain(icurr).gnet*L*1e-4)-1)./(out_gain(icurr).gnet*L*1e-4);
        ASE(icurr).ase_out=material.beta_sp.*out_gain(icurr).rsp.*ampl_term;  % in eV^-1*cm^-3*s^-1

        
        figure(2)
        hold on
        plot(1.24./htw(1:end-1)*1000,out_gain(icurr).gmat)
        
        figure(3)
        hold on
        plot(1.24./htw(1:end-1)*1000,out_gain(icurr).gnet)
        


        figure(4)
        hold on
        plot(1.24./htw(1:end-1)*1000,out_gain(icurr).rsp)

        figure(5)
        hold on
        plot(1.24./htw(1:end-1)*1000,ASE(icurr).ase_out)
        
        Efh_vs_curr(icurr)=out_gain(icurr).Efh;
        
        Efe_vs_curr(icurr)=out_gain(icurr).Efc;
end

figure(1)
plot(J,(Efe_vs_curr-Efh_vs_curr)+material.Eg)
xlabel('Current density in kA/cm^{2} ')
ylabel('Separation between quasi-Fermi levels (E_{F,e}-E_{F,h}) [eV]')




figure(2)
axis([1.24./(material.Egp2+0.1)*1000  1.24/material.Eg*1000  -3000  max(out_gain(length(I)).gmat)*1.1])
xlabel('Wavelength  (nm)')
ylabel('Material gain cm^{-1}')
title(['Current density in kA/cm^{2} ', num2str(J)])


figure(3)
axis([1.24./(material.Egp2+0.1)*1000  1.24/material.Eg*1000  -100  max(out_gain(length(I)).gnet)*1.1])
xlabel('Wavelength  (nm)')
ylabel('Net modal gain cm^{-1}')
title(['Current density in kA/cm^{2} ', num2str(J)])


figure(4)
axis([1.24./(material.Egp2+0.1)*1000  1.24/material.Eg*1000  0  max(out_gain(length(I)).rsp)*1.1])
xlabel('Wavelength  (nm)')
ylabel('Spontaneous emission rate   r_{sp}(\lambda)  (cm^-3*eV^-1*s^-1) ')
title(['Current density in kA/cm^{2} ', num2str(J)])


figure(5)
axis([1.24./(material.Egp2+0.1)*1000  1.24/material.Eg*1000  0  max(ASE(length(I)).ase_out)*1.1])
xlabel('Wavelength  (nm)')
ylabel('Amplified spontaneous emission rate   ASE(\lambda)  (cm^-3*eV^-1*s^-1) ')
title(['Current density in kA/cm^{2} ', num2str(J)])











