 clear,close all,clc;
%%

NumBit = 2^20;

Dados = randi([0 1],1,NumBit);

tb = 1/(12.5e1);
ta = tb/1;
OvSa = tb/ta;
Data = rectpulse(Dados,OvSa);
tf = tb*NumBit;
t=linspace(0,tf,OvSa*NumBit);
f=time2freq(t);
[ SelecFilt ] = FiltroGaussiano(f,1/tb,0,1);
% [ SelecFilt ] = Filtro_Retangular( 34*fc,16*fc,f);
SelecFilt = fftshift(SelecFilt);
% Data = ifft(fft(Data).*SelecFilt);
figure;
plot(t,Data);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);

figure;
hold all;
plot(f,10*log10(abs(fftshift(fft(Data)./length(Data)))));
plot(f,10*log10(abs(fftshift(SelecFilt))));
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
a=1;