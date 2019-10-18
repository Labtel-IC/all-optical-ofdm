close all;clear;clc;
%        parametros
fc           = 12.5e9;                                                     %center frequency of Electical signal ^10 kilo ^20 mega ^30 giga
Rb           = fc;                                                         %bit rate
NPPB         = 2^6;                                                        %Number of Point Per Bit
fsample      = NPPB*Rb;                                                    %Sampling frequency
T            = 1/Rb;                                                       %period of the filter
Ta           = 1/fsample;                                                  %period of the filter
% Tb           = 1/Rb;                                                     %time of one bit NRZ
NumberOf_T   = 2^13;                                                       %number of periods on the simulation 2^18 is the m�ximum value for numberEyes = 2
FinalTime    = NumberOf_T*T;                                               %Final time of simulation
Nb           = floor(FinalTime/T);                                         %number of bits
% TotalSamples = FinalTime*fsample;                                        %Total Number of Samples
TotalSamples = NumberOf_T*(T/Ta);                                          %Total Number of Samples

t            = linspace(0,FinalTime,round(TotalSamples));                         %Time vector for simulation
f            = time2freq(t); 

CoupConst = 3*pi/4;    %Couples Constant of the DI
ExpCoef   = exp(-1j*pi/4);%Exponential component from Ramaswami equation
MutCoef   = -1;         %Coeficient to control polarity

% E0 = zeros(1,length(t));
E0 = sin(2*pi*fc*t);
E1 = cos(2*pi*fc*t);
figure;plot(t,E0,t,E1);
legend('E0','E1');
axis([-0.1*T 6*T -1.1 1.1]);

Eout11 = ExpCoef.*(cos(CoupConst).*E1 + MutCoef*1j*sin(CoupConst).*E0);        
Eout12 = ExpCoef.*(cos(CoupConst).*E0 + MutCoef*1j*sin(CoupConst).*E1);

figure;plot(t,E0,t,E1,t,abs(Eout11)./max(abs(Eout11)),t,abs(Eout12)./max(abs(Eout12)),t,real(Eout11)./max(real(Eout11)),t,real(Eout12)./max(real(Eout12)),t,imag(Eout11)./max(imag(Eout11)),t,imag(Eout12)./max(imag(Eout12)));
legend('E0','E1','absEout11','absEout12','realEout11','realEout12','imagEout11','imagEout12');
axis([-0.1*T 6*T -1.1 1.1]);

EoutIDt = ifft(fft(Eout11).*exp(-1j*2*pi*f*T));
EoutIDf = exp(1j*pi/4).*Eout12;

figure;plot(t,E0,t,E1,t,abs(Eout11)./max(abs(Eout11)),t,abs(Eout12)./max(abs(Eout12)),t,abs(EoutIDt)./max(abs(EoutIDt)),t,abs(EoutIDf)./max(abs(EoutIDf)),t,real(Eout11)./max(real(Eout11)),t,real(Eout12)./max(real(Eout12)),t,real(EoutIDt)./max(real(EoutIDt)),t,real(EoutIDf)./max(real(EoutIDf)),t,imag(Eout11)./max(imag(Eout11)),t,imag(Eout12)./max(imag(Eout12)),t,imag(EoutIDt)./max(imag(EoutIDt)),t,imag(EoutIDf)./max(imag(EoutIDf)));
legend('E0','E1','absEout11','absEout12','absEoutIDt','absEoutIDf','realEout11','realEout12','realEoutIDt','realEoutIDf','imagEout11','imagEout12','imagEoutIDt','imagEoutIDf');
axis([-0.1*T 6*T -1.1 1.1]);

Eout1 = ExpCoef.*(cos(CoupConst).*EoutIDt + MutCoef*1j*sin(CoupConst).*EoutIDf);  
Eout2 = ExpCoef.*(cos(CoupConst).*EoutIDf + MutCoef*1j*sin(CoupConst).*EoutIDt); 

figure;plot(t,E0,t,E1,t,abs(Eout11)./max(abs(Eout11)),t,abs(Eout12)./max(abs(Eout12)),t,abs(EoutIDt)./max(abs(EoutIDt)),t,abs(EoutIDf)./max(abs(EoutIDf)),t,abs(Eout1)./max(abs(Eout1)),t,abs(Eout2)./max(abs(Eout2)),t,real(Eout11)./max(real(Eout11)),t,real(Eout12)./max(real(Eout12)),t,real(EoutIDt)./max(real(EoutIDt)),t,real(EoutIDf)./max(real(EoutIDf)),t,real(Eout1)./max(real(Eout1)),t,real(Eout2)./max(real(Eout2)),t,imag(Eout11)./max(imag(Eout11)),t,imag(Eout12)./max(imag(Eout12)),t,imag(EoutIDt)./max(imag(EoutIDt)),t,imag(EoutIDf)./max(imag(EoutIDf)),t,imag(Eout1)./max(imag(Eout1)),t,real(Eout2)./max(real(Eout2)));
legend('E0','E1','absEout11','absEout12','absEoutIDt','absEoutIDf','absEout1','absEout2','realEout11','realEout12','realEoutIDt','realEoutIDf','realEout1','realEout2','imagEout11','imagEout12','imagEoutIDt','imagEoutIDf','imagEout1','imagEout2');
axis([-0.1*T 6*T -1.1 1.1]);

a=1;
