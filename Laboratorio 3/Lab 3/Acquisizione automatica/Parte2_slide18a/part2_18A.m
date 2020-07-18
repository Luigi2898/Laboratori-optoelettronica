%% catatterisitica PI

DataIN=importdata('PI_Temp20,00_IStart0,10_IStop12,00_IStep0,10.txt');
Power=DataIN(:,2);
Current=DataIN(:,1);
Power=Power.*10; %Tengo così conto delle perdite dovute ai connettori

figure(1)
hold on
grid on;
plot(Current, Power.*1e3);
xlabel('Current [mA]');
ylabel('Power [mW]');
axis([0 12 0 0.9]);

%% Script per plotting catatterisitica VI
DataIN=importdata('VI_Temp20,00_IStart0,10_IStop12,00_IStep0,10.txt');
Voltage=DataIN(:,2);
Current=DataIN(:,1);

figure(2)
hold on
grid on;
plot(Current, Voltage);
xlabel('Current [mA]');
ylabel('Voltage [V]');

% %% Script per plotting catatterisitica IpdI
% DataIN=importdata('IpdI_Temp20,00_IStart0,10_IStop12,00_IStep0,10.txt');
% Ipd=DataIN(:,2);
% Current=DataIN(:,1);
% 
% figure(3)
% plot(Current, Ipd, '.');
% xlabel('Current [mA]');
% ylabel('Ipd [uA]');


% %% Valutazione della WPE
% WPE=(Power./(Current.*Voltage)).*1e3;
% 
% figure(4)
% plot(Current, WPE)
% xlabel('Current [mA]');
% ylabel('WPE ');