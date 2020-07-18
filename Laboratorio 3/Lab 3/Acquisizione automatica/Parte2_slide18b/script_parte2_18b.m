%% catatterisitica PI

DataIN=importdata('PI_IStart0,00_IStop50,00_IStep1,00.txt');
Power=DataIN(:,2);
Current=DataIN(:,1);
%Power=Power.*10; %Tengo così conto delle perdite dovute ai connettori

figure(1)
hold on
grid on;
plot(Current, Power.*1e3);
xlabel('Current [mA]');
ylabel('Power [mW]');


% %% Script per plotting catatterisitica VI
DataIN=importdata('VI_IStart0,00_IStop50,00_IStep1,00.txt');
Voltage=DataIN(:,2);
Current=DataIN(:,1);

figure(2)
hold on
grid on;
plot(Current, Voltage);
xlabel('Current [mA]');
ylabel('Voltage [V]');

% %% Script per lo studio della caratteristica IpdI in assenza di controllo di temperatura
%   
DataIN=importdata('IpdI_IStart0,00_IStop50,00_IStep1,00.txt');
Ipd=DataIN(:,2);
Current=DataIN(:,1);
 
figure(3)
hold on;
grid on;
plot(Current, Ipd);
xlabel('Current [mA]');
ylabel('Ipd [uA]');
% 
% %% Script per ottenere i risultati di confronto con le caratteristiche con il controllo della temperatura
% 
% fig1=hgload('IpdI_no_temp_ctrl.fig');
% fig2=hgload('IpdI_temp_ctrl.fig');
% fig3=hgload('PI_no_temp_ctrl.fig');
% fig4=hgload('PI_temp_ctrl.fig');
% fig5=hgload('VI_no_temp_ctrl.fig');
% fig6=hgload('VI_temp_ctrl.fig');
% 
% L2 = findobj(fig2,'type','line');
% L4 = findobj(fig4,'type','line');
% L6 = findobj(fig6,'type','line');
% copyobj(L2,findobj(fig1,'type','axes'));
% copyobj(L4,findobj(fig3,'type','axes'));
% copyobj(L6,findobj(fig5,'type','axes'));
