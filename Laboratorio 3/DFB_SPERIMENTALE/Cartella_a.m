%%GRAFICI a

%Primo file
text1 = fopen('IpdI_Temp20,00_IStart0,10_IStop12,00_IStep0,10_new.txt');

C1 = textscan(text1, '%f %f', 'Delimiter', ' ', 'EmptyValue', -Inf);
I1 = C1{1};
I_photod= C1{2};

figure(1)
plot(I1, I_photod);
hold;
grid;

xlabel('Current [mA]');
ylabel('I_{pd} [\muA]');

%Secondo file
text2 = fopen('PI_Temp20,00_IStart0,10_IStop12,00_IStep0,10_new.txt');

C2 = textscan(text2, '%f %f', 'Delimiter', ' ', 'EmptyValue', -Inf);
I2 = C2{1};
P = C2{2};
P2 = P .* 10;

figure(2)
plot(I2, P2);
hold;
grid;

xlabel('Current [mA]');
ylabel('Potenza [W]');

%Terzo file
text3 = fopen('VI_Temp20,00_IStart0,10_IStop12,00_IStep0,10_new.txt');

C3 = textscan(text3, '%f %f', 'Delimiter', ' ', 'EmptyValue', -Inf);
I3 = C3{1};
V3 = C3{2};

figure(3)
plot(I3, V3);
hold;
grid;

xlabel('Current [mA]');
ylabel('Tensione [V]');

fclose('all');