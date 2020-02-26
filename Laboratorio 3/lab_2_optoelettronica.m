%%
clc
clear all
close all

%caratterizzazione della potenza e della I_TEC del laser a T = 20°C lambda
%1531.1 nm

I_laser = 5e-3:1e-3:35e-3;
I_TEC =  [0.066,0.066,0.067,0.067,0.068,0.068,0.069,0.069,0.069,0.070,0.070,0.071,0.071,0.072,0.073,0.073,0.074,0.075,0.076,0.077,0.078,0.078,0.079,0.080,0.081,0.082,0.083,0.083,0.085,0.086,0.087]; %ampere
power =[0.72e-6,1.07e-6,1.52e-6,2.05e-6,2.70e-6,3.49e-6,8.08e-6,23.22e-6,40.98e-6,60.10e-6,78.53e-6,98.26e-6,116.42e-6,136.08e-6,155.10e-6,174.85e-6,0.20e-3,0.2144e-3,0.2348e-3,0.2549e-3,0.2754e-3,0.2957e-3,0.3162e-3,0.3371e-3,0.3573e-3,0.3778e-3,0.3986e-3,0.4190e-3,0.4411e-3,0.4620e-3,0.4820e-3];%WATT


figure(1)
plot(I_laser,power)
grid on
figure(2)
plot(I_laser,I_TEC)
%OK

%caratterizzazione della potenza e della I_TEC del laser a T = 35.1°C lambda
%1532.85 nm ricalcolato ricollegando il laser all'OSA

I_TEC_35 = [-0.077,-0.077,-0.077,-0.077,-0.077,-0.077,-0.077,-0.076,-0.076,-0.076,-0.075,-0.075,-0.074,-0.074,-0.073,-0.073,-0.072,-0.071,-0.071,-0.070,-0.070,-0.069,-0.068,-0.068,-0.067,-0.066,-0.066,-0.065,-0.064,-0.063,-0.062,-0.057,-0.054,-0.049,-0.045];
power_35 = [0.410e-6,0.578e-6,0.780e-6,1.003e-6,1.290e-6,1.5981e-6,1.950e-6,2.332e-6,2.753e-6,3.220e-6,3.710e-6,4.415e-6,10.760e-6,21.83e-6,32.53e-6,43.85e-6,56.00e-6,67.93e-6,80.84e-6,93.88e-6,107.11e-6,120.91e-6,134.92e-6,149.50e-6,163.38e-6,178.13e-6,192.70e-6,0.2075e-3,0.2230e-3,0.2378e-3, 0.2526e-3,0.3297e-3,0.4059e-3,0.4840e-3,0.5580e-3];
I_laser_step = 40e-3:5e-3:55e-3; % max corrente prima di self protection
I_plot = horzcat(I_laser,I_laser_step);
% 
% figure(3)
% plot(I_laser,power_35)
% grid on
% figure(4)
% plot(I_laser,I_TEC_35)

figure(5)
plot(I_plot,power_35)
grid on

figure(6)
plot(I_plot,I_TEC_35)
grid on 

%stavolta con step di 5mA
%spegniamo il controllo della temperatura e vediamo T VS I e P VS I
%impostiamo sul power meter una lambda di 1532 nm data dalla media delle
%due precedentemente valutate
I_tec_OFF = 5e-3:5e-3:55e-3;
T =[27.4,27.5,27.7,27.9,28.2,28.6,28.8,29.2,29.7,30.1,31.1];
P =[0.5464e-6,2.341e-6,28.15e-6,100.40e-6,181.01e-6,0.2648e-3,0.3506e-3,0.4340e-3,0.5170e-3,0.5905e-3,0.6616e-3]; %c'è rumore di ampiezza nelle letture

figure(7)
plot(I_tec_OFF,T)
grid on
figure(8)
plot(I_tec_OFF,P)
grid on
