%Script per calcolare guadagno nel caso di multi quantum well.

clear all;
close all;
clc;

%%Caricamento delle costanti
define_constants;



plotta=1; 
ogniquanto=1; %due flag per determinare quante curve del guadagno plottare

lambdafix=0.84 %[um] - lambda fissa su cui calcolare il guadagno differenziale

%Calcolo degli stati confinati
levelsc=solv_disp_eq1D(mstar,m2star,Ebulk,10*W,0,0);
levelsh=solv_disp_eq1D(mstarhh,m2starhh,Ebulkh,10*W,0,0);

Ec0=levelsc(1);
Eh0=levelsh(1);

Egp = Eg + levelsc(1) + levelsh(1);
Egp2 = Eg + levelsc(2) + levelsh(2);

mrstar=(mstar*mstarhh)/(mstar+mstarhh);
meffel=mstar*const_phys.m0.data;
mefflac=mstarhh*const_phys.m0.data;

%Vettore concentrazione portatori N
%stepN=1e12;

%N = 1e12:stepN:1e13;
N=1e13;

% Vettore energia htw
step=0.001;
htw=Eg:step:Egp+0.3;

figure(1);
figure(2);
h = waitbar(0,'Please wait...');
for iN=1:length(N)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     Calcolo del guadagno
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Calcolo della Efc e Efh
    guess=0; 
    opt=optimset('fzero');
    meffel=mstar*const_phys.m0.data;
    Efc=fzero('errmqw',guess,opt,N(iN),meffel,const_phys.K.data,T,const_phys.ht.data,levelsc,W,N0,Ebulk,totlen);

    meffh=mstarhh*const_phys.m0.data;
    Efh=fzero('errmqw',guess,opt,N(iN),meffh,const_phys.K.data,T,const_phys.ht.data,levelsh,W,N0,Ebulkh,totlen);
    
    %Concentrazione portatori nei quantum well [cm^-3]
    Nqw(iN)=(nquantumwell(levelsc,Efc,const_phys.K.data,T,meffel,const_phys.ht.data))/(W*1e-7);
    Hqw(iN)=(nquantumwell(levelsh,Efh,const_phys.K.data,T,meffh,const_phys.ht.data))/(W*1e-7);
    
    %Inizializzazione variabili per il guadagno    
    gnc(1,:)=zeros(1,length(htw));
    gconv(1,:)=zeros(1,299);
    gconvl(1,:)=zeros(1,299);
    gnc(2,:)=zeros(1,length(htw));
    gconv(2,:)=zeros(1,299);
    gconvl(2,:)=zeros(1,299);
    
    %Calcolo del guadagno per ogni stato confinato
    for ilev=1:length(levelsc)
        Egp = Eg + levelsc(ilev) + levelsh(ilev);
        
        %Rinormalizzazione dell'energy gap
        deltaEg(iN)=-xi*((Nqw(iN)+Hqw(iN))/2)^(1/3);
        Egpn(iN)=Egp+deltaEg(iN);
        
        %Calcolo del guadagno
        [gnc(ilev,:),Econvp,gconv(ilev,:),lambdap,gconvl(ilev,:),valmaxp,lambdavalmaxp]=mqw_guad(N(iN),htw,mstar,mstarhh,T,levelsc,levelsh,Efc,Efh,W,N0,Ebulk,Ebulkh,totlen,Eg,Egpn(iN),cp,delta,tauin,nr,1);
        if(ilev==1)
            Econv=Econvp;
            lambda=lambdap;
            Egpmin=Egp;
        end
    end
    
    gnct=sum(gnc); %Guadagno non convoluto
    gconvt=sum(gconv); %Guadagno convoluto funzione dell'energia
    gconvlt=sum(gconvl); %Guadagno convoluto funzione della lunghezza d'onda
    gmat=gconvlt; %Guadagno materiale
    
    ebarr=nbulk(meffel,Efc,const_phys.K.data,T,const_phys.ht.data,Ebulk); %numero di elettroni nel bulk [cm^-3]
    hbarr=nbulk(mefflac,Efh,const_phys.K.data,T,const_phys.ht.data,Ebulkh); %numero di lacune nel bulk [cm^-3]
    
    [gmn]=guadmodnet(N0,gammaqw,gmat,alfai,gammabarr,kn,ebarr,kp,hbarr); %Guadagno modale netto [cm^-1]
    
    %Per determinare il guadagno a una lunghezza d'onda fissata:
    indlamb=find(lambda>=lambdafix);
    glf(iN)=gmn(indlamb(1));
    lstr=num2str(lambda(indlamb(1)));
    rl(iN)=lambda(indlamb(1));
    
        
    %Per un determinato N plotto il guadagno con e senza le perdite (ma con
    %il fattore di confinamento)
    if(N(iN)==8e12)
        [gmnnoperdbarr]=guadmodnet(N0,gammaqw,gmat,alfai,gammabarr,0,ebarr,0,hbarr);
        [gmnnoperd]=guadmodnet(N0,gammaqw,gmat,0,gammabarr,0,ebarr,0,hbarr);
        Nstr=num2str(N(iN));
        figure;
        plot(lambda,gmnnoperd);
        hold on;
        plot(lambda,gmnnoperdbarr,'r');
        hold on;
        plot(lambda,gmn,'g');        
        legend('Gain without losses','Gain without losses in the barrier','Guadagno with losses');
        xlabel('lambda [um]');
        ylabel('G [cm^-^1]');
        title(['Gain vs lambda con N =' Nstr 'cm^-^2']);
        grid on;
        
        %Inoltre, sempre per questo determinato N, plotto la differenza tra
        %il guadagno convoluto con la lorenziana e quello convoluto con la
        %secante iperbolica
        for ilev=1:length(levelsc)
            Egp = Eg + levelsc(ilev) + levelsh(ilev);
        
            %Rinormalizzazione dell'energy gap
            deltaEg(iN)=-xi*((Nqw(iN)+Hqw(iN))/2)^(1/3);
            Egpn(iN)=Egp+deltaEg(iN);
        
            %Calcolo del guadagno
            [gncs(ilev,:),Econvps,gconvs(ilev,:),lambdaps,gconvls(ilev,:),valmaxps,lambdavalmaxps]=mqw_guad(N(iN),htw,mstar,mstarhh,T,levelsc,levelsh,Efc,Efh,W,N0,Ebulk,Ebulkh,totlen,Eg,Egpn(iN),cp,delta,tauin,nr,2);
            if(ilev==1)
                Econvs=Econvps;
                lambdas=lambdaps;
                valmaxs(iN)=valmaxps;
                lambdavalmaxs(iN)=lambdavalmaxps;
                Egpmins=Egp;
            end
        end
    
        gncts=sum(gncs);
        gconvts=sum(gconvs);
        gconvlts=sum(gconvls);
        gmats=gconvlts;
        
        [gmns]=guadmodnet(N0,gammaqw,gmats,alfai,gammabarr,kn,ebarr,kp,hbarr);
        
        figure
        plot(lambda,gmn);
        hold on;
        plot(lambda,gmns,'r');
        legend('Convoluzione con lorenziana','Convoluzione con secante iperbolica');
        xlabel('lambda [um]');
        ylabel('G [cm^-^1]');
        title(['Gain vs lambda con N =' Nstr 'cm^-^2']);
        grid on;
        
    end
        
    %Calcolo dell'efficienza (Considero il caso specifico di 2 livelli
    %confinati nella banda di conduzione e 3 nella banda di valenza
    ne1(iN)=nquantumwellsing(levelsc(1),Efc,const_phys.K.data,T,meffel,const_phys.ht.data);
    ne2(iN)=nquantumwellsing(levelsc(2),Efc,const_phys.K.data,T,meffel,const_phys.ht.data);
    
    nh1(iN)=nquantumwellsing(levelsh(1),Efh,const_phys.K.data,T,mefflac,const_phys.ht.data);
    nh2(iN)=nquantumwellsing(levelsh(2),Efh,const_phys.K.data,T,mefflac,const_phys.ht.data);
    nh3(iN)=nquantumwellsing(levelsh(3),Efh,const_phys.K.data,T,mefflac,const_phys.ht.data);
    
    effe1(iN)=N0*ne1(iN)/N(iN);
    effe2(iN)=N0*ne2(iN)/N(iN);
    effh1(iN)=N0*nh1(iN)/N(iN);
    effh2(iN)=N0*nh2(iN)/N(iN);
    effh3(iN)=N0*nh3(iN)/N(iN);
    
    [valmax(iN),ascmax]=max(gmn);
    lambdavalmax(iN)=lambda(ascmax);
    
    %Plot del guadagno
    if(plotta==1)
        Nstr=num2str(N(iN));
        figure(1);
        hold on;
		plot(lambda,gmn);
		text(lambdavalmax(iN),valmax(iN),['N = ' Nstr 'cm^-^2']);  
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     Variazione di Neff
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Rinormalizzazione dell'energy gap
    Egp=Eg+levelsc(1)+levelsh(1);
    deltaEg(iN)=-xi*((Nqw(iN)+Hqw(iN))/2)^(1/3);
    Egpn(iN)=Egp+deltaEg(iN);
    [dnb]=deltanbarr(lambda,ebarr,nr,meffel);
    [lambda,nconvl]=deltanwell(htw,gnc(1,:),tauin,Egpn(iN),Egpmin,step);
    
    %Calcolo di deltaneff
    for il=1:length(lambda)
        termone(il)=gammaqw*N0*nconvl(il);
        termtwo(il)=gammabarr*dnb(il);
        dneff(il)=termone(il)+termtwo(il);
    end
    
    %Deltaneff per un determinato lambda
    indlamb=find(lambda>=lambdafix);
    dnefflf(iN)=dneff(indlamb(1));
    lstr2=num2str(lambda(indlamb(1)));
    
 
    [valmaxn(iN),ascmaxn]=min(dneff);
    lambdavalmaxn(iN)=lambda(ascmaxn);
    Nstr=num2str(N(iN));
    if(plotta==1)
        figure(2);
        hold on;
		plot(lambda,dneff,'r');
		text(lambdavalmaxn(iN),valmaxn(iN),['N = ' Nstr 'cm^-^2']);  
    end

    plotta=plotta+1;
    if(plotta>ogniquanto)
        plotta=1;
    end
    
    waitbar(iN/length(N),h)
end

close(h);

%Sistemazione figure guadagno e variazione neff
figure(1);
grid on;
xlabel('lambda [um]');
ylabel('G [cm^-^1]');
title('Gain vs lambda with differen values of N');
prism;

figure(2);
grid on;
xlabel('lambda [um]');
ylabel('deltaneff');
title('Variazione di neff');
prism;

%Incrementiamo il numero di punti con la spline
NN=1e12:1e10:1e13;
glfspl = spline(N,glf,NN);

%Figure guadagno differenziale e alfaH
figure;
plot(N,glf,'*',NN,glfspl);
grid on;
xlabel('N - [cm^-^2]');
ylabel('Gain [cm^-^1]');
title(['Gain corresponding to lambda = ' lstr 'um']);

dnefflfspl = spline(N,dnefflf,NN);
figure;
plot(N,dnefflf,'*',NN,dnefflfspl);
grid on;
xlabel('N - [cm^-^2]');
ylabel('deltaneff');
title(['Deltaneff corresponding to lambda = ' lstr 'um']);

%Guadagno differenziale:
dg=diff(glfspl)/stepN;
figure;
plot(NN(2:end),dg);
title(['Differential gain - lambda = ' lstr 'um']);
xlabel('Ntot [cm^-2]');
ylabel('dg/dN [cm]');
grid on;

%d(deltaneff)/dN
dneff=diff(dnefflfspl)/stepN;
figure;
plot(NN(2:end),dneff);
title(['Differential delta eff - lambda = ' lstr 'um']);
xlabel('Ntot [cm^-2]');
ylabel('d(deff)/dN [cm^2]');
grid on;

alfah=-1e4*4*pi/lambda(indlamb(1))*dneff./dg;
figure;
plot(NN(2:end),alfah);
title(['alfah(N) - lambda = ' lstr 'um']);
xlabel('N [cm^-2]');
ylabel('alfah');
grid on;

%Altri grafici
figure;
plot(N,valmax,'-*');
grid on;
xlabel('N - [cm^-^2]');
ylabel('Highest gain [cm^-^1]');
title('Highest values of gain function of total number of N');

figure;
plot(Nqw,valmax,'-*');
grid on;
xlabel('N - [cm^-^3]');
ylabel('Highest gain [cm^-^1]');
title('Highest values of gain function of N in quantum wells');

figure;
plot(Nqw,lambdavalmax,'-*');
grid on;
xlabel('N - [cm^-^3]');
ylabel('Lambda [um]');
title('Lambda corresponding to the maximum gain with different values of N');

figure;
plot(N,Egpn);
grid on;
xlabel('N - [cm^-^2]');
ylabel('Egp [eV]');
title('Rinormalization of the Energy Gap');


%Plot dell'efficienza
figure;
plot(N,effe1);
hold on;
plot(N,effe2,'r');
legend('First level','Second level');
grid on;
title('Efficiency - conduction band');

figure;
plot(N,effh1);
hold on;
plot(N,effh2,'r');
hold on;
plot(N,effh3,'g');
legend('First level','Second level','Third level');
grid on;
title('Efficiency - valence band');