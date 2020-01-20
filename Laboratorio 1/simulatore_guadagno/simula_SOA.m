% Simulazione comportamento statico SOA in MQW
% M. Gioannini - Maggio 2007



global material
global const_phys
global guess_Fermi


define_constants

% -------------------------------------------------------------------------
% scelta del materiale del SOA

% caricamento del file con parametri del materiale SOA

% load_data: caricamento dei dati di ingresso 

[filegeo,pathgeo]=uigetfile('*.in','Choose file with material parameters');
DATA=textread([pathgeo filegeo],'%s','commentstyle','matlab','whitespace','\n');
for idata=1:length(DATA)
   eval(DATA{idata});
end

% -------------------------------------------------------------------------

QW_properties
% calcolo di proprietà QW




% Vettore energia htw
step=1e-3;
htw=material.Egp-0.1:step:material.Egp2+0.2;

ns=50;






% inserimento geometria del SOA

L=input('SOA length  (um)? ');
wSOA=input('SOA waveguide width (um)? '); % um fissata ad esempio a  4 um 
dz_slice=L/ns;
I=input('Injected current  (mA)? ');






lin=input('Wavelength of input signal (nm)? ');




disp('Choose simulation type: ')
disp('1: Calculates SOA gain versus input power (no noise included)')
disp('2: Calculates SOA gain versus input power INCLUDING NOISE')
%disp('3: Calculates SOA output spectrum including noise')
sim_type=input('Insert your choice ','s');




switch  sim_type
    case '1'  
        Prange=input('Power variation: [Pmin  Pmax] in dBm? ');
        DP=input('Power variation step in dBm? ');
        Pin_vardBm=Prange(1):DP:Prange(2);
        Pin_var=10.^(Pin_vardBm/10);
        htw_sgn=1.24./(lin)*1000;
        htw=[htw,htw_sgn]; % l'ultimo elemento è il segnale
        material_Ifixed  % Calculation of material properties for fixed injected current
        isgn=length(htw);
        incl_spon=0; % non viene incluso il rumore per emissione spontanea amplificata
        SOA_Pin_Iin
        % grafici
        Gplot=cat(1,out_SOA.Gtot);
        Pout_fine=linspace(Pin_vardBm(1)+10*log10(Gplot(1).'),Pin_vardBm(end)+10*log10(Gplot(end).'),length(Pin_var)*10);
        GSOA_fine=interp1(Pin_vardBm+10*log10(Gplot.'),10*log10(Gplot),Pout_fine,'pchip');
        figure
        plot(Pin_vardBm+10*log10(Gplot.'),10*log10(Gplot),'*');
        hold on
        plot(Pout_fine,GSOA_fine)
        xlabel('Output power (dBm)')
        ylabel('Signal Gain  G (dB)')
    case '2'
        Prange=input('Power variation: [Pmin  Pmax] in dBm? ');
        DP=input('Power variation step in dBm? ');
        Pin_vardBm=Prange(1):DP:Prange(2);
        Pin_var=10.^(Pin_vardBm/10);
        htw_sgn=1.24./(lin)*1000;
        htw=[htw,htw_sgn]; % l'ultimo elemento è il segnale
        material_Ifixed  % Calculation of material properties for fixed injected current
        isgn=length(htw);
        incl_spon=1; %  viene incluso il rumore per emissione spontanea amplificata
        SOA_Pin_Iin
        % grafici
        Gplot=cat(1,out_SOA.Gtot);
        Pout_fine=linspace(Pin_vardBm(1)+10*log10(Gplot(1).'),Pin_vardBm(end)+10*log10(Gplot(end).'),length(Pin_var)*10);
        GSOA_fine=interp1(Pin_vardBm+10*log10(Gplot.'),10*log10(Gplot),Pout_fine,'pchip');
        figure
        plot(Pin_vardBm+10*log10(Gplot.'),10*log10(Gplot),'*');
        hold on
        plot(Pout_fine,GSOA_fine)
        xlabel('Output power (dBm)')
        ylabel('Signal Gain  G (dB)')


%         Fplot=cat(1,out_SOA.F);
%         figure
%         plot(Pin_vardBm+10*log10(Gplot.'),10*log10(Fplot),'-*');
%         hold on
%         xlabel('Output power (dBm)')
%         ylabel('Noise Figure  G (dB)')

        sgn_vett=zeros(1,length(htw));
        sgn_vett(isgn)=1;

        % spettri in uscita
        figure(100)
        
        for iP=1:length(Pin_var)
            PN(:,iP)=(out_SOA(iP).SF_noise.*(htw.'))/(dz_slice*1e-4)*(const_phys.c.data/material.nr*1e10)*const_phys.e.data*1e-19*1e3; % spettro della potenza di rumore in uscita
            Pout_sgn(iP)=(out_SOA(iP).SF_sgn.*(htw(isgn).'))/(dz_slice*1e-4)*(const_phys.c.data/material.nr*1e10)*const_phys.e.data*1e-19*1e3; % potenza di segnale in uscita
            Pout_SOA_spectrum(:,iP)=PN(:,iP)+Pout_sgn(iP)*sgn_vett.';
            [htw_sort,isort]=sort(htw);
            
            figure(100)
            %h(iP)=plot(ll*1e3,10*log10(Pout_SOA_spectrum(:,iP)),'Linewidth',[2]);
            h(iP)=plot(1.24./htw_sort*1e3,10*log10(Pout_SOA_spectrum(isort,iP)),'Linewidth',[2]);
            axis([1.24./(material.Egp2+0.1)*1000  1.24/material.Eg*1000  max(10*log10(Pout_SOA_spectrum(:,iP)))-80   max(10*log10(Pout_SOA_spectrum(:,iP)))+5])
            outP=num2str(Pin_vardBm(iP)+10*log10(Gplot(iP).'));
            outG=num2str(10*log10(Gplot(iP)))
            ht=title(['Input power ',num2str(Pin_vardBm(iP)),'dBm; ','Output power ',outP, 'dBm; ', 'Gain ', outG,'dB' ]);
            set(ht,'Fontweight','bold')
            pause
            hold on
            set(h(iP),'Linewidth',[0.5])
            xlabel('Wavelength  (nm)')
            ylabel('Output power (dBm)')
        end
   maxGG=0;   
  % spettri di guadagno del chip
        figure(101)
        for iP=1:length(Pin_var)
            %plot(ll*1e3,10*log10(out_SOA(iP).Gtot_spectrum));
            plot(1.24./htw_sort*1e3,10*log10(out_SOA(iP).Gtot_spectrum(isort)))
            hold on
            maxGG=max(maxGG,max(10*log10(out_SOA(iP).Gtot_spectrum)));
            text(1550,max(10*log10(out_SOA(iP).Gtot_spectrum)),['Pin=',num2str(Pin_vardBm(iP)),'dBm'])
        end
           xlabel('Wavelength  (nm)')
           ylabel('Chip Gain G (dB)')
          axis([1.24./(material.Egp2+0.1)*1000  1.24/material.Eg*1000  -15   maxGG+5])


%     case '3'
%         
%         Pin_var=input('Input power  (mW)? ')
%         htw_sgn=1.24./(lin)*1000;
%         in=find(htw>=htw_sgn);
%         isgn=in(1);
%         incl_spon=1; % non viene incluso il rumore per emissione spontanea amplificata
%         SOA_Pin_Iin
%         % grafici
% 
%         PN=(SF_noise(:,ns).*(htw.'))/(dz_slice*1e-4)*(const_phys.c.data/material.nr*1e10)*const_phys.e.data*1e-19*1e3; % spettro della potenza di rumore in uscita
%         sgn_vett=zeros(1,length(htw));
%         sgn_vett(isgn)=1;
%         Pout_sgn=(SF_sgn(ns).*(htw(isgn).'))/(dz_slice*1e-4)*(const_phys.c.data/material.nr*1e10)*const_phys.e.data*1e-19*1e3; % potenza di segnale in uscita
%         Pout_SOA_spectrum=PN+Pout_sgn*sgn_vett.';
%         figure
%         plot(ll*1e3,10*log10(Pout_SOA_spectrum));
%         xlabel('Wavelength  (nm)')
%         ylabel('Output power (dBm)')
end

        
        
 





