close all
clear all
clc

%% Script per lo studio delle variazioni della caratteristica PI al variare di T

SourceFiles=["PI_Temp20,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "PI_Temp21,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "PI_Temp22,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "PI_Temp23,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "PI_Temp24,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "PI_Temp25,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "PI_Temp26,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "PI_Temp27,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "PI_Temp28,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "PI_Temp29,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "PI_Temp30,00_IStart0,00_IStop40,00_IStep0,50.txt"];
   
legendInfoPI=cell(length(SourceFiles),1);
figure(1)
for i=1:length(SourceFiles)
    DataIN=importdata(SourceFiles(i));
    Power=DataIN(:,2);
    Current=DataIN(:,1);
    Power=Power.*10; %Tengo così conto delle perdite dovute ai connettori
    
    plot(Current, Power.*1e3);
    xlabel('Current [mA]');
    ylabel('Power [mW]');
    grid on;
    legendInfoPI{i}=sprintf('@T= %d°C',(19+i));
    hold on;
end

legend(legendInfoPI);

% %% Script per lo studio delle variazioni della caratteristica VI al variare di T

SourceFiles=["VI_Temp20,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "VI_Temp21,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "VI_Temp22,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "VI_Temp23,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "VI_Temp24,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "VI_Temp25,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "VI_Temp26,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "VI_Temp27,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "VI_Temp28,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "VI_Temp29,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "VI_Temp30,00_IStart0,00_IStop40,00_IStep0,50.txt"];
   
legendInfoVI=cell(length(SourceFiles),1);
figure(2)
for i=1:length(SourceFiles)
    DataIN=importdata(SourceFiles(i));
    Voltage=DataIN(:,2);
    Current=DataIN(:,1);
    
    grid on;
    plot(Current, Voltage);
    xlabel('Current [mA]');
    ylabel('Voltage [V]');
    legendInfoVI{i}=sprintf('@T= %d°C',(19+i));
    hold on;
end

legend(legendInfoVI);

% %% Script per lo studio delle variazioni della caratteristica IpdI al variare di T

SourceFiles=["IpdI_Temp20,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "IpdI_Temp21,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "IpdI_Temp22,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "IpdI_Temp23,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "IpdI_Temp24,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "IpdI_Temp25,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "IpdI_Temp26,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "IpdI_Temp27,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "IpdI_Temp28,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "IpdI_Temp29,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "IpdI_Temp30,00_IStart0,00_IStop40,00_IStep0,50.txt"];
   
legendInfoIpdI=cell(length(SourceFiles),1);
figure(3)
for i=1:length(SourceFiles)
    DataIN=importdata(SourceFiles(i));
    Ipd=DataIN(:,2);
    Current=DataIN(:,1);

    grid on;
    plot(Current, Ipd);
    xlabel('Current [mA]');
    ylabel('Ipd [uA]');
    legendInfoIpdI{i}=sprintf('@T= %d°C',(19+i));
    hold on;
end

legend(legendInfoIpdI);

% %% Script per lo studio delle variazioni della Wall Plug Efficiency
%    
SourceFiles=["PI_Temp20,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "PI_Temp21,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "PI_Temp22,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "PI_Temp23,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "PI_Temp24,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "PI_Temp25,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "PI_Temp26,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "PI_Temp27,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "PI_Temp28,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "PI_Temp29,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "PI_Temp30,00_IStart0,00_IStop40,00_IStep0,50.txt"];
   
SourceFilesVI=["VI_Temp20,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "VI_Temp21,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "VI_Temp22,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "VI_Temp23,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "VI_Temp24,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "VI_Temp25,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "VI_Temp26,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "VI_Temp27,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "VI_Temp28,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "VI_Temp29,00_IStart0,00_IStop40,00_IStep0,50.txt", ...
             "VI_Temp30,00_IStart0,00_IStop40,00_IStep0,50.txt"];
   
figure(4)
for i=1:length(SourceFiles)
    DataIN=importdata(SourceFiles(i));
    Power=DataIN(:,2);
    Power=Power.*10; %Tengo così conto delle perdite dovute ai connettori

    DataIN=importdata(SourceFilesVI(i));
    Voltage=DataIN(:,2);
    Current=DataIN(:,1);
    
    WPE=(Power.*1e3./(Voltage.*Current)).*0.1*100; %Tengo conto dell'attenuazione di 10 dB e dell'unità di misura della WPE 
    
    grid on
    plot(Current, WPE);
    xlabel('Current [mA]');
    ylabel('WPE %');
    hold on;
    xlim([3 40])
end

legendInfoWPE=cell(11,1);
figure(4)
for i=1:11
   
    hold on;
    legendInfoWPE{i}=sprintf('@T= %d°C',(19+i));
    
end

legend(legendInfoWPE);

%% Script per lo studio dell'andamento delle Ith

% Fissata Poh=0.13mW si ottiene:
Ith=[12.278 12.656 13.083 13.537 14.029 14.517 15.02 15.543 16.044 16.555 17.061];
T=[20 21 22 23 24 25 26 27 28 29 30];

figure(5)
hold on;
grid on;
plot(T, Ith);
xlabel('T [°C]');
ylabel('Ith [mA]');

% %% Script per la valutazione di Rs
% 
% % A partire dal grafico VI è agevole ricavare la resistenza a partire dalla
% % pendenza della caratteristica. In particolare:
Current_ON=Current(find(Current==10):end); %V
Voltage_ON=Voltage(find(Current==10):end); %mA

figure(6)
hold on;
grid on;
plot(Current_ON, Voltage_ON*1e3);
ylabel('V [mV]');
xlabel('I [mA]');
% 

