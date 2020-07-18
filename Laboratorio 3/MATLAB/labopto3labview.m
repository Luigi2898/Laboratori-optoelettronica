% close all
% clear all
% clc% %%%% Grafico 1 con con temperatura 20 gradi con Istart=0.10mA e Istop=12mA
% % filename = 'Iph-T20.txt';
% % A = importdata(filename);
% % I=A(:,1); %corrente in mA
% % I_t=I';
% % Iph=A(:,2);% corrente del fotodiodo del controllore di temperatura in uA
% % Iph_t=Iph';
% % figure('Name','Iph vs I')
% % hold on
% % plot(I_t,Iph_t)
% % %%%% Grafico 2 con con temperatura 31.22 gradi con Istart=0mA e Istop=50mA
% % filename = 'Iph-T31.txt';
% % A1 = importdata(filename);
% % I1=A1(:,1); %corrente in mA
% % I_t1=I1';
% % Iph1=A(:,2);% corrente del fotodiodo del controllore di temperatura in uA
% % Iph_t1=Iph1';
% % 
% % %%%% Grafico 3 con con temperatura 20 gradi con Istart=0mA e Istop=40mA
% % filename = 'Iph-T20_1.txt';
% % A2 = importdata(filename);
% % I2=A2(:,1); %corrente in mA
% % I_t2=I2';
% % Iph2=A2(:,2);% corrente del fotodiodo del controllore di temperatura in uA
% % Iph_t2=Iph2';
% % %%%% Grafico 4 con con temperatura 21 gradi con Istart=0mA e Istop=40mA
% % filename = 'Iph-T21.txt';
% % A3 = importdata(filename);
% % I3=A3(:,1); %corrente in mA
% % I_t3=I3';
% % Iph3=A3(:,2);% corrente del fotodiodo del controllore di temperatura in uA
% % Iph_t3=Iph3';
% % %%%% Grafico 5 con con temperatura 20 gradi con Istart=0mA e Istop=40mA
% % filename = 'Iph-T22.txt';
% % A4 = importdata(filename);
% % I4=A4(:,1); %corrente in mA
% % I_t4=I4';
% % Iph4=A4(:,2);% corrente del fotodiodo del controllore di temperatura in uA
% % Iph_t4=Iph4';
% % %%%% Grafico 6 con con temperatura 20 gradi con Istart=0mA e Istop=40mA
% % filename = 'Iph-T23.txt';
% % A5 = importdata(filename);
% % I5=A5(:,1); %corrente in mA
% % I_t5=I5';
% % Iph5=A5(:,2);% corrente del fotodiodo del controllore di temperatura in uA
% % Iph_t5=Iph5';
% % %%%% Grafico 7 con con temperatura 20 gradi con Istart=0mA e Istop=40mA
% % filename = 'Iph-T24.txt';
% % A6 = importdata(filename);
% % I6=A6(:,1); %corrente in mA
% % I_t6=I6';
% % Iph6=A6(:,2);% corrente del fotodiodo del controllore di temperatura in uA
% % Iph_t6=Iph6';
% % %%%% Grafico 8 con con temperatura 20 gradi con Istart=0mA e Istop=40mA
% % filename = 'Iph-T25.txt';
% % A7 = importdata(filename);
% % I7=A7(:,1); %corrente in mA
% % I_t7=I7';
% % Iph7=A7(:,2);% corrente del fotodiodo del controllore di temperatura in uA
% % Iph_t7=Iph7';
% % %%%% Grafico 8 con con temperatura 20 gradi con Istart=0mA e Istop=40mA
% % filename = 'Iph-T26.txt';
% % A8 = importdata(filename);
% % I8=A8(:,1); %corrente in mA
% % I_t8=I8';
% % Iph8=A8(:,2);% corrente del fotodiodo del controllore di temperatura in uA
% % Iph_t8=Iph8';
% % %%%% Grafico 8 con con temperatura 20 gradi con Istart=0mA e Istop=40mA
% % filename = 'Iph-T27.txt';
% % A9 = importdata(filename);
% % I9=A9(:,1); %corrente in mA
% % I_t9=I9';
% % Iph9=A9(:,2);% corrente del fotodiodo del controllore di temperatura in uA
% % Iph_t9=Iph9';
% % %%%% Grafico 8 con con temperatura 20 gradi con Istart=0mA e Istop=40mA
% % filename = 'Iph-T28.txt';
% % A10 = importdata(filename);
% % I10=A10(:,1); %corrente in mA
% % I_t10=I10';
% % Iph10=A10(:,2);% corrente del fotodiodo del controllore di temperatura in uA
% % Iph_t10=Iph10';
% % %%%% Grafico 8 con con temperatura 20 gradi con Istart=0mA e Istop=40mA
% % filename = 'Iph-T29.txt';
% % A11 = importdata(filename);
% % I11=A11(:,1); %corrente in mA
% % I_t11=I11';
% % Iph11=A11(:,2);% corrente del fotodiodo del controllore di temperatura in uA
% % Iph_t11=Iph11';
% % %%%% Grafico 8 con con temperatura 20 gradi con Istart=0mA e Istop=40mA
% % filename = 'Iph-T30.txt';
% % A12 = importdata(filename);
% % I12=A12(:,1); %corrente in mA
% % I_t12=I12';
% % Iph12=A12(:,2);% corrente del fotodiodo del controllore di temperatura in uA
% % Iph_t12=Iph12';
%%%%%%%Grafico P-I
%%% Grafico 1 con con temperatura 20 gradi con Istart=0,1mA e Istop=12mA
filename = 'P_I_T20.txt';
B = importdata(filename);
I_1=B(:,1); %corrente in mA
I_1t=I_1';
P=B(:,2);% potenza in W che moltiplico per 10 per considerare delle perdite dei connettori
P_t=10.*P';
figure('Name','P vs I a T=20°C con passo di 0.1mA da 0-12mA')
plot(I_1t,P_t); 
xlabel('I[mA]');
ylabel('P[W]');
grid on
% %%%% Grafico 2 con con temperatura 31 gradi con Istart=0mA e Istop=50mA
filename = 'P_I_T31.txt';
B1= importdata(filename);
I_2=B1(:,1); %corrente in mA
I_2t=I_2';
P1=B1(:,2);% potenza in W che moltiplico per 10 per considerare le perdite dei connettori
P1_t=10.*P1';
figure ('Name','P vs I con temperatura a 31 gradi con passo di 1mA tra 0-50mA')
plot(I_2t,P1_t)
xlabel('I[mA]');
ylabel('P[W]');
grid on
%%%% Grafico 3 con con temperatura 20 gradi con Istart=0mA e Istop=40mA
filename = 'P_I_T20_1.txt';
B2= importdata(filename);
I_3=B2(:,1); %corrente in mA
I_3t=I_3';
P2=B2(:,2);% potenza in W che moltiplico per 10 per considerare le perdite dei connettori
P2_t=10.*P2';
figure ('Name','P vs I al variare della T')
plot(I_3t,P2_t);
hold on
xlabel('I[mA]');
ylabel('P[W]');
grid on
%%%% Grafico 4 con con temperatura 21 gradi con Istart=0mA e Istop=50mA
filename = 'P_I_T21.txt';
B3= importdata(filename);
I_4=B3(:,1); %corrente in mA
I_4t=I_4';
P3=B3(:,2);% potenza in W che moltiplico per 10 per considerare le perdite dei connettori
%figure ('Name','P vs I con T=21°C')
P3_t=10.*P3';
plot(I_4t,P3_t);
grid on
%%%% Grafico 5 con con temperatura 22 gradi con Istart=0mA e Istop=40mA
filename = 'P_I_T22.txt';
B4= importdata(filename);
I_5=B4(:,1); %corrente in mA
I_5t=I_5';
P4=B4(:,2);% potenza in W che moltiplico per 10 per considerare le perdite dei connettori
P4_t=10.*P4';
%figure ('Name','P vs I con T=22°C')
plot(I_5t,P4_t);
grid on
%%%% Grafico 6 con con temperatura 23 gradi con Istart=0mA e Istop=40mA
filename = 'P_I_T23.txt';
B5= importdata(filename);
I_6=B5(:,1); %corrente in mA
I_6t=I_6';
P5=B5(:,2);% potenza in W che moltiplico per 10 per considerare le perdite dei connettori
P5_t=10.*P5';
%figure ('Name','P vs I con T=23°C')
plot(I_6t,P5_t);
grid on
%%%% Grafico 7 con con temperatura 20 gradi con Istart=0mA e Istop=50mA
filename = 'P_I_T24.txt';
B6= importdata(filename);
I_7=B6(:,1); %corrente in mA
I_7t=I_7';
P6=B6(:,2);% potenza in W che moltiplico per 10 per considerare le perdite dei connettori
%figure ('Name','P vs I con T=24°C')
P6_t=10.*P6';
plot(I_7t,P6_t);
grid on
%%%% Grafico 8 con con temperatura 25 gradi con Istart=0mA e Istop=40mA
filename = 'P_I_T25.txt';
B7= importdata(filename);
I_8=B7(:,1); %corrente in mA
I_8t=I_8';
P7=B7(:,2);% potenza in W che moltiplico per 10 per considerare le perdite dei connettori
P7_t=10.*P7';
%figure ('Name','P vs I con T=25°C')
plot(I_8t,P7_t);
grid on
%%%% Grafico 9 con con temperatura 26 gradi con Istart=0mA e Istop=40mA
filename = 'P_I_T26.txt';
B8= importdata(filename);
I_9=B8(:,1); %corrente in mA
I_9t=I_9';
P8=B8(:,2);% potenza in W che moltiplico per 10 per considerare le perdite dei connettori
P8_t=10.*P8';
%figure ('Name','P vs I con T=26°C')
plot(I_9t,P8_t);
grid on
%%%% Grafico 10 con con temperatura 27 gradi con Istart=0mA e Istop=40mA
filename = 'P_I_T27.txt';
B9= importdata(filename);
I_10=B9(:,1); %corrente in mA
I_10t=I_10';
P9=B9(:,2);% potenza in W che moltiplico per 10 per considerare le perdite dei connettori
P9_t=10.*P9';
%figure ('Name','P vs I con T=27°C')
plot(I_10t,P9_t);
%%%% Grafico 11 con con temperatura 28 gradi con Istart=0mA e Istop=40mA
filename = 'P_I_T28.txt';
B10= importdata(filename);
I_11=B10(:,1); %corrente in mA
I_11t=I_11';
P10=B10(:,2);% potenza in W che moltiplico per 10 per considerare le perdite dei connettori
P10_t=10.*P10';
%figure ('Name','P vs I con T=28°C')
plot(I_11t,P10_t);
grid on
%%%% Grafico 12 con con temperatura 29 gradi con Istart=0mA e Istop=40mA
filename = 'P_I_T29.txt';
B11= importdata(filename);
I_12=B11(:,1); %corrente in mA
I_12t=I_12';
P11=B11(:,2);% potenza in W che moltiplico per 10 per considerare le perdite dei connettori
P11_t=10.*P11';
%figure ('Name','P vs I con T=29°C')
plot(I_12t,P11_t);
grid on
%%%% Grafico 13 con temperatura 30 gradi con Istart=0mA e Istop=40mA
filename = 'P_I_T30.txt';
B12= importdata(filename);
I_13=B12(:,1); %corrente in mA
I_13t=I_13';
P12=B12(:,2);% potenza in W che moltiplico per 10 per considerare le perdite dei connettori
P12_t=10.*P12';
%figure ('Name','P vs I con T=30°C')
plot(I_13t,P12_t);
legend('T=20°C','T=21°C','T=22°C','T=23°C','T=24°C','T=25°C','T=26°C','T=27°C','T=28°C','T=29°C','T=30°C','Location','northwest','NumColumns',2);
grid on
hold off
%%%% Calcolo corrente di soglia rispetto alla temperatura
I_th=[12 12.5 13 13.5 14 14.5 15 15.5 16 16.5 17];
T=[20 21 22 23 24 25 26 27 28 29 30];
figure('Name','Ith vs T')
plot(T,I_th);
xlabel('I_{th}[mA]');
ylabel('T[°C]');

% % % % % Caratteristica V-I
%%%%%Grafico 1 con temperatura 20 gradi con Istart=0,1mA e Istop=12mA
filename = 'V_I_T20.txt';
C = importdata(filename);
I_V0=C(:,1); %corrente in mA
I__t=I_V0';
V=C(:,2);% V sul diodo del controllore di temperatura in V
V_t=V';
figure ('Name','V-I con T=20°C con passo di 0.1 mA da 0-12mA')
plot(I__t,V_t);
xlabel('I[mA]');
ylabel('V[V]');
I__t1_1=I__t.*0.001;
pendenza=polyfit(I__t1_1,V_t,1);
grid on
Rs2=121.285;
V_2th=1.689;
%WPE_1=(P_t.*10)/(((I_1t*1000).^2).*Rs2+((I_1t./1000).*V_2th)).*100;
%plot(I_1t,WPE_1);
xlabel('I[mA]');
ylabel('WPE[%]');
grid on
%%%%%Grafico 2 con temperatura 31 gradi con Istart=0mA e Istop=50mA
filename = 'V_I_T31.txt';
C1= importdata(filename);
I_V=C1(:,1); %corrente in mA
I__t1=I_V';
V1=C1(:,2);% V sul diodo del controllore di temperatura in V
V_1t=V1';
figure ('Name','V-I con T=31 °C con passo di 1ma da 0-50mA')
plot(I__t1,V_1t);
xlabel('I[mA]');
ylabel('V[V]');
I__t1_2=I__t1.*0.001;
pendenza1=polyfit(I__t1_2,V_1t,1);
grid on
Rs1=117.1003;
V_1th=1.264;
WPE=(P1_t.*1e4)./(((I_2t).^2).*Rs1+((I_2t./1000).*V_1th)).*100;
plot(I_2t,WPE);
xlabel('I[mA]');
ylabel('WPE[%]');
grid on
%%%%Grafico 3 con temperatura 31 gradi con Istart=0mA e Istop=40mA
filename = 'V_I_T20_1.txt';
C2= importdata(filename);
I_V1=C2(:,1); %corrente in mA
I__t2=I_V1';
V2=C2(:,2);% V sul diodo del controllore di temperatura in V
V_2t=V2';
figure ('Name','V-I al variare della T')
plot(I__t2,V_2t);
grid on
xlabel('I[mA]');
ylabel('V[V]');
hold on
I__t2_3=I__t2.*0.001;
pendenza2=polyfit(I__t2_3,V_2t,1);
%%%%Grafico 1 con temperatura 21 gradi con Istart=0mA e Istop=50mA
filename = 'V_I_T21.txt';
C12= importdata(filename);
I_V11=C12(:,1); %corrente in mA
I__t12=I_V11';
V12=C12(:,2);% V sul diodo del controllore di temperatura in V
V_12t=V12';
%figure ('Name','V-I a 21 gradi')
plot(I__t12,V_12t);
grid on
I__t12_13=I__t12.*0.001;
pendenza12=polyfit(I__t12_13,V_12t,1);
%%%%Grafico 4 con temperatura 31 gradi con Istart=0mA e Istop=40mA
filename = 'V_I_T22.txt';
C3= importdata(filename);
I_V2=C3(:,1); %corrente in mA
I__t3=I_V2';
V3=C3(:,2);% V sul diodo del controllore di temperatura in V
V_3t=V3';
%figure ('Name','V-I a 22 gradi')
plot(I__t3,V_3t);
grid on
I__t3_4=I__t3.*0.001;
pendenza3=polyfit(I__t3_4,V_3t,1);
%%%%Grafico 5 con temperatura 31 gradi con Istart=0mA e Istop=40mA
filename = 'V_I_T23.txt';
C4= importdata(filename);
I_V3=C4(:,1); %corrente in mA
I__t4=I_V3';
V4=C4(:,2);% V sul diodo del controllore di temperatura in V
V_4t=V4';
%figure ('Name','V-I a 23 gradi')
plot(I__t4,V_4t);
grid on
I__t4_5=I__t4.*0.001;
pendenza4=polyfit(I__t4_5,V_4t,1);
%%%%Grafico 1 con temperatura 31 gradi con Istart=0mA e Istop=50mA
filename = 'V_I_T24.txt';
C5= importdata(filename);
I_V4=C5(:,1); %corrente in mA
I__t5=I_V4';
V5=C5(:,2);% V sul diodo del controllore di temperatura in V
V_5t=V5';
%figure ('Name','V-I a 24 gradi')
plot(I__t5,V_5t);
grid on
I__t5_6=I__t5.*0.001;
pendenza5=polyfit(I__t5_6,V_5t,1);
%%%%Grafico 1 con temperatura 31 gradi con Istart=0mA e Istop=50mA
filename = 'V_I_T25.txt';
C6= importdata(filename);
I_V5=C6(:,1); %corrente in mA
I__t6=I_V5';
V6=C6(:,2);% V sul diodo del controllore di temperatura in V
V_6t=V6';
%figure ('Name','V-I a 25 gradi')
plot(I__t6,V_6t);
grid on
I__t6_7=I__t6.*0.001;
pendenza6=polyfit(I__t6_7,V_6t,1);
%%%%Grafico 1 con temperatura 31 gradi con Istart=0mA e Istop=50mA
filename = 'V_I_T26.txt';
C7= importdata(filename);
I_V6=C7(:,1); %corrente in mA
I__t7=I_V6';
V7=C7(:,2);% V sul diodo del controllore di temperatura in V
V_7t=V7';
%figure ('Name','V-I a 26 gradi')
plot(I__t7,V_7t);
grid on
I__t7_8=I__t7.*0.001;
pendenza7=polyfit(I__t7_8,V_7t,1);
%%%%Grafico 1 con temperatura 31 gradi con Istart=0mA e Istop=50mA
filename = 'V_I_T27.txt';
C8= importdata(filename);
I_V7=C8(:,1); %corrente in mA
I__t8=I_V7';
V8=C8(:,2);% V sul diodo del controllore di temperatura in V
V_8t=V8';
%figure ('Name','V-I a 27 gradi')
plot(I__t8,V_8t);
grid on
I__t8_9=I__t8.*0.001;
pendenza8=polyfit(I__t8_9,V_8t,1);
%%%%Grafico 1 con temperatura 31 gradi con Istart=0mA e Istop=50mA
filename = 'V_I_T28.txt';
C9= importdata(filename);
I_V8=C9(:,1); %corrente in mA
I__t9=I_V8';
V9=C9(:,2);% V sul diodo del controllore di temperatura in V
V_9t=V9';
%figure ('Name','V-I a 28 gradi')
plot(I__t9,V_9t);
grid on
I__t9_10=I__t9.*0.001;
pendenza9=polyfit(I__t9_10,V_9t,1);
%%%Grafico 1 con temperatura 31 gradi con Istart=0mA e Istop=50mA
filename = 'nonfunge.txt';
C10= importdata(filename);
I_V9=C10(:,1); %corrente in mA
I__t10=I_V9';
V10=C10(:,2);% V sul diodo del controllore di temperatura in V
V_10t=V10';
%figure ('Name','V-I a 29 gradi')
plot(I__t10,V_10t);
grid on
I__t10_11=I__t10.*0.001;
pendenza10=polyfit(I__t10_11,V_10t,1);
%%%%Grafico 1 con temperatura 31 gradi con Istart=0mA e Istop=50mA
filename = 'V_I_T30txt.txt';
C11= importdata(filename);
I_V10=C11(:,1); %corrente in mA
I__t11=I_V10';
V11=C11(:,2);% V sul diodo del controllore di temperatura in V
V_11t=V11';
%figure ('Name','V-I a 30 gradi')
plot(I__t11,V_11t);
grid on
I__t11_12=I__t11.*0.001;
pendenza11=polyfit(I__t11_12,V_11t,1);
hold off
legend('T=20°C','T=21°C','T=22°C','T=23°C','T=24°C','T=25°C','T=26°C','T=27°C','T=28°C','T=29°C','T=30°C','Location','northwest','NumColumns',2);
%%WPE-WallPlugEffieciency
Rs=[120.359 120.318 120.321 120.317 120.312 120.309 120.329 120.328 120.309 120.137 120.310];
V_th=[1.288 1.296 1.293 1.29 1.29 1.289 1.286 1.288 1.29 1.288 1.289];
WPE1=(P2_t.*1e4)/(((I_3t).^2).*Rs(1)+((I_3t.*0.001)*V_th(1)))*100;
figure('Name','WPE vs I')
plot(I_3t,WPE1);
hold on
xlabel('I[mA]');
ylabel('WPE[%]');
grid on
WPE2=(P3_t.*1e4)./(((I_4t).^2).*Rs(2)+((I_4t)*V_th(2)))*100;
plot(I_4t,WPE2);
WPE3=(P4_t.*1e4)./(((I_5t).^2).*Rs(3)+((I_5t)*V_th(3)))*100;
plot(I_5t,WPE3);
WPE4=(P5_t.*1e4)./(((I_6t).^2).*Rs(4)+((I_6t)*V_th(4)))*100;
plot(I_6t,WPE4);
WPE5=(P6_t.*1e4)./(((I_7t).^2).*Rs(5)+((I_7t)*V_th(5)))*100;
plot(I_7t,WPE5);
WPE6=(P7_t.*1e4)./(((I_8t).^2).*Rs(6)+((I_8t)*V_th(6)))*100;
plot(I_8t,WPE6);
WPE7=(P8_t.*1e4)./(((I_9t).^2).*Rs(7)+((I_9t)*V_th(7)))*100;
plot(I_9t,WPE7);
WPE8=(P9_t.*1e4)./(((I_10t).^2).*Rs(8)+((I_10t)*V_th(8)))*100;
plot(I_10t,WPE8);
WPE9=(P10_t.*1e4)./(((I_11t).^2).*Rs(9)+((I_11t)*V_th(9)))*100;
plot(I_11t,WPE9);
WPE10=(P11_t.*1e4)./(((I_12t).^2).*Rs(10)+((I_12t)*V_th(10)))*100;
plot(I_12t,WPE10);
WPE11=(P12_t.*1e4)./(((I_13t).^2).*Rs(11)+((I_13t)*V_th(11)))*100;
plot(I_13t,WPE11);
hold off
legend('T=20°C','T=21°C','T=22°C','T=23°C','T=24°C','T=25°C','T=26°C','T=27°C','T=28°C','T=29°C','T=30°C','Location','southeast','NumColumns',2);
