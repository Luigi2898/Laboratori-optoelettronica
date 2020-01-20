
M_Lorentzian=matrix_L(htw,material.tauin);



Nini=1e15*dz_slice*1e-4*wSOA*1e-4*material.N0*material.W*1e-7;

opt_ode=odeset('Reltol',1e-5);

opt_fsolve=optimset('fsolve');

guess_Fermi.Efc=-material.Eg/2;
guess_Fermi.Efh=-material.Eg/2;

ll=1.24./(htw);
material.beta_sp=1./(2*pi^2*material.nr^2).*(ll).^2*material.gammaqw*material.N0/(wSOA*material.N0*material.W*1e-3);
[t,Nt]=ode45('carrier_RE_TIME',[0 5e-9],Nini,opt_ode,I,0,0,dz_slice,L,htw,wSOA,M_Lorentzian,1);
Nstz=Nt(end);


% Refine (per essere sicuri di aver scelto la soluzione stazionaria)

Nstz_norm=fsolve('carrier_RE_STAZ',Nstz*1e-6,opt_fsolve,I,0,0,dz_slice,L,htw,wSOA,M_Lorentzian,1); 
out_gain=gain_spectrum(Nstz_norm*1e6/(dz_slice*1e-4*wSOA*1e-4),htw,M_Lorentzian); % portatori in cm^-2,htw,M_Lorentzian); % entrambi in cm^-1
guess_Fermi.Efc=out_gain.Efc;
guess_Fermi.Efh=out_gain.Efh;

figure
plot(1.24./htw(1:end-1)*1000,out_gain.gnet(1:end-1))
axis([1.24./(material.Egp2+0.1)*1000  1.24/material.Eg*1000  -100  max(out_gain.gnet)*1.1])
xlabel('Wavelength  (nm)')
ylabel('Net modal gain cm^{-1}')
hold on
plot(1.24./htw(end)*1000,out_gain.gnet(end),'r*')
text(1.24./htw(end)*1000,out_gain.gnet(end)*1.01,'Signal')

