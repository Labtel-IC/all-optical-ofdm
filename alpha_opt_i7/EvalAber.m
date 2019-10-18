close all;clc;clear;
load('Test_4PAM_car_32_offt_0.mat')
figure;
hold on;
grid on;
plot(mean(VetNumCar2),mean(VetOpenO1)./mean(VetOpenO1(:,1)),...
     mean(VetNumCar2),mean(VetOpenO2)./mean(VetOpenO2(:,1)),...
     mean(VetNumCar2),mean(VetOpenO3)./mean(VetOpenO3(:,1)),'LineWidth',2);
axis([0 max(mean(VetNumCar2)) 0 1+1.2*max([max(max(VetOpenO1)) max(max(VetOpenO2)) max(max(VetOpenO3))])]);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
BERt4PAM

load('Test_4PAM_car_32_offt_1.mat')
figure;
hold on;
grid on;
plot(mean(VetNumCar2),mean(VetOpenO1)./mean(VetOpenO1(:,1)),...
     mean(VetNumCar2),mean(VetOpenO2)./mean(VetOpenO2(:,1)),...
     mean(VetNumCar2),mean(VetOpenO3)./mean(VetOpenO3(:,1)),'LineWidth',2);
axis([0 max(mean(VetNumCar2)) 0 1+1.2*max([max(max(VetOpenO1)) max(max(VetOpenO2)) max(max(VetOpenO3))])]);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
BERt4PAM

load('Test_OOK_car_32_offt_0.mat')
figure;
hold on;
grid on;
plot(mean(VetNumCar2),mean(VetOpen2)./mean(VetOpen2(:,1)),'LineWidth',2);
axis([0 max(mean(VetNumCar2)) 0 1+1.2*max(max(VetOpen2))]);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
BERtOOK

load('Test_OOK_car_32_offt_1.mat')
figure;
hold on;
grid on;
plot(mean(VetNumCar2),mean(VetOpen2)./mean(VetOpen2(:,1)),'LineWidth',2);
axis([0 max(mean(VetNumCar2)) 0 1+1.2*max(max(VetOpen2))]);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
BERtOOK