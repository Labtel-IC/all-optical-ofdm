close all;clear;clc;
CoupConst = 3*pi/4;    %Couples Constant of the DI
ExpCoef   = exp(-1j*pi/4);%Exponential component from Ramaswami equation
MutCoef   = -1;         %Coeficient to control polarity

%%        parametros
fc           = 12.5e9;                                                     %center frequency of Electical signal ^10 kilo ^20 mega ^30 giga
Rb           = fc;                                                         %bit rate
NPPB         = 2^6;                                                        %Number of Point Per Bit
fsample      = NPPB*Rb;                                                    %Sampling frequency
T            = 1/Rb;                                                       %period of the filter
Ta           = 1/fsample;                                                  %period of the filter
% Tb           = 1/Rb;                                                     %time of one bit NRZ
NumberOf_T   = 2^13;                                                       %number of periods on the simulation 2^18 is the mï¿½ximum value for numberEyes = 2
FinalTime    = NumberOf_T*T;                                               %Final time of simulation
Nb           = floor(FinalTime/T);                                         %number of bits
% TotalSamples = FinalTime*fsample;                                        %Total Number of Samples
TotalSamples = NumberOf_T*(T/Ta);                                          %Total Number of Samples

t            = linspace(0,FinalTime,round(TotalSamples));                         %Time vector for simulation
f            = time2freq(t); 
E0 = zeros(1,length(t));
% CW= cos(2*pi*fc*t) + 1j*sin(2*pi*fc*t);
% CW= exp(1j*2*pi*2*fc*t);
CW= cos(2*pi*2*fc*t);
figure;plot(f,20.*log10(abs(fftshift(fft(CW)./length(CW)))));
a=1;
%%
BitFilBanWid  = 4*12.5e9; %Badnwidth for the bit conformation
BitFiltCenFre = 0;    %Center frequency for the bit conformation
BitFilOrd     = 1;    %The order for the filter which will conform data
[BitFilt,~] = FiltroGaussiano(f,BitFilBanWid,BitFiltCenFre,BitFilOrd);                                 %Creating filter for conformation of the input information
BitFilt = fftshift(BitFilt);
%%             dados

% TxData = randi(2,1,2*Nb)-1;
% PreBits = [0 0];
% TxPoss = linspace(1,length(TxData),length(TxData));

TxDataSeg = [0 0;%01
             0 0;%02
             0 1;%03
             1 0;%04
             1 1;%05
             1 0;%06
             0 0;%07
             1 1;%08
             0 1;%09
             1 1;%10
             1 0;%11
             0 0;%12
             1 0;%13
             0 1;%14
             1 1;%15
            0 1];%16
TxData = [];
PreBits = [0 0];
countseg = 1;
for kk=1:Nb
    TxData = [TxData TxDataSeg(countseg,:)];
    countseg = countseg + 1;
    if countseg > 16
        countseg = 1;
%         TxData = [TxData ones(1,2*Nb-32)];
%         break;
    end
end
TxPoss = linspace(1,length(TxData),length(TxData));
[DataI,DataQ] = DqpskEncodEq(TxData,PreBits);

u = TxData(mod(TxPoss,2)==1);
v = TxData(~mod(TxPoss,2));

[I,Q] = DQPSK_encoder(u,v,PreBits(1),PreBits(2));

TxI = rectpulse(DataI,NPPB);                                
TxQ = rectpulse(DataQ,NPPB);

I = rectpulse(I,NPPB);                                 
Q = rectpulse(Q,NPPB);                                 

TxI(TxI==0) = -0.5;                                     
TxI(TxI==1) = 0.5;                                    
TxQ(TxQ==0) = -0.5;                                    
TxQ(TxQ==1) = 0.5;                                   
% figure;plot(t,I,t,TxI);
% figure;plot(t,Q,t,TxQ);

Vbias         = 3.8;
Vpi2          = 3.8;

Isig1 = ifft(fft(TxI).*BitFilt);
Qsig1 = ifft(fft(TxQ).*BitFilt);
Isig2 = ifft(fft(I).*BitFilt);
Qsig2 = ifft(fft(Q).*BitFilt);

figure;plot(t,Isig1,t,Qsig1);
axis([0 16*T -1.1 1.1]);
%%                Modulação
[Eout,~] = IqModTest (t,CW,Isig1,Qsig1,Vpi2,Vbias);

% Eout = Eout.*conj(Eout);
[ Eout1,Eout2 ] = DelayInterf( t,T,1*pi/4,Eout,E0);
[ Eout3,Eout4 ] = DelayInterf( t,T,-1*pi/4,Eout,E0);
[EoutA,EoutB,~]=mz_interf4(Eout,fftshift(f),Rb,1*pi/4);
[EoutC,EoutD,~]=mz_interf4(Eout,fftshift(f),Rb,-1*pi/4);
A = Eout2.*conj(Eout2);
B = Eout1.*conj(Eout1);

BitFilBanWid  = 1*12.5e9; %Badnwidth for the bit conformation
BitFiltCenFre = 0;    %Center frequency for the bit conformation
BitFilOrd     = 1;    %The order for the filter which will conform data
[BitFilt,~] = FiltroGaussiano(f,BitFilBanWid,BitFiltCenFre,BitFilOrd);                                 %Creating filter for conformation of the input information
BitFilt = fftshift(BitFilt);
A = ifft(fft(A).*BitFilt);
B = ifft(fft(B).*BitFilt);
EoutI1 = A - B;
EoutQ1 = Eout4.*conj(Eout4) - Eout3.*conj(Eout3);
EoutI = EoutB.*conj(EoutB) - EoutA.*conj(EoutA);
EoutQ = EoutD.*conj(EoutD) - EoutC.*conj(EoutC);
figure;
plot(t,Isig1./max(Isig1),...
     t,Qsig1./max(Qsig1),...
     t,abs(EoutI1)./max(abs(EoutI1)),...
     t,real(EoutI1)./max(real(EoutI1)),...
     t,imag(EoutI1)./max(imag(EoutI1)));axis([-0.1*T 16*T -1.1 1.1]);
figure;
plot(t,Isig1./max(Isig1),...
     t,Qsig1./max(Qsig1),...
     t,real(EoutI)./max(real(EoutI)),...
     t,real(EoutQ)./max(real(EoutQ)));





[E,Ei,Eq]=qi_MZM(CW, Isig2, Qsig2,Vbias,Vpi2);
figure;
plot(t,abs(E)./max(abs(E)),...
     t,abs(Eout)./max(abs(Eout)),...
     t,real(E)./max(real(E)),...
     t,real(Eout)./max(real(Eout)),...
     t,imag(E)./max(imag(E)),...
     t,imag(Eout)./max(imag(Eout)));
% Olho(abs(E)./max(abs(E)),T,NPPB,1);
% Olho(abs(Eout)./max(abs(Eout)),T,NPPB,1);
%%           Recepção

E_rec1 = Eout;
E_rec2 = E;
E_rec3 = Eout./max(abs(Eout));
%%
Ui = real(exp(1j*pi/4).*E_rec3.*conj(ifft(fft(E_rec3).*exp(-1j*2*pi*f*T))));
Uq = imag(exp(1j*pi/4).*E_rec3.*conj(ifft(fft(E_rec3).*exp(-1j*2*pi*f*T))));
figure;plot(t,Isig1./max(Isig1),t,Qsig1./max(Qsig1),t,real(Ui)./max(real(Ui)),t,real(Uq)./max(real(Uq)));axis([0 8e-10 -1.1 1.1]);
%%
Eaux1 = exp(-1j*pi/4).*ifft(fft(Eout).*exp(-1j*2*pi*f*T));
componentX = E_rec3.*Eaux1;
componentX = ifft(fft(componentX).*BitFilt);
figure;plot(t,Isig1./max(Isig1),t,Qsig1./max(Qsig1),t,abs(componentX)./max(abs(componentX)));axis([0 8e-10 -1.1 1.1]);
%####################ComponenteI##########################################
Eout11 = ExpCoef.*(cos(CoupConst).*E_rec1 + MutCoef*1j*sin(CoupConst).*E0);        
Eout12 = ExpCoef.*(cos(CoupConst).*E0 + MutCoef*1j*sin(CoupConst).*E_rec1);

EoutIDt = ifft(fft(Eout11).*exp(-1j*2*pi*f*T));
EoutIDf = exp(1j*pi/4).*Eout12;

EoutDemodI1 = ExpCoef.*(cos(CoupConst).*EoutIDt + MutCoef*1j*sin(CoupConst).*EoutIDf);  
EoutDemodI2 = ExpCoef.*(cos(CoupConst).*EoutIDf + MutCoef*1j*sin(CoupConst).*EoutIDt);       

EoutDemodI = EoutIDf.*conj(EoutDemodI2) - EoutIDt.*conj(EoutDemodI1);
EoutDemodI = EoutDemodI./max(EoutDemodI);

[uu1,uu2,~]             = mz_interf4(E_rec2.*(10^(-3.0/20)),fftshift(f),Rb,pi/4);
comp_I1 = uu1.*conj(uu1);
comp_I2 = uu2.*conj(uu2);
comp_I = comp_I2-comp_I1;
comp_I = comp_I./max(comp_I);
% comp_I = ifft(fft(comp_I).*BitFilt);
%###################ComponenteQ###########################################
Eout21 = ExpCoef.*(cos(CoupConst).*E_rec1 + MutCoef*1j*sin(CoupConst).*E0);        
Eout22 = ExpCoef.*(cos(CoupConst).*E0 + MutCoef*1j*sin(CoupConst).*E_rec1);

EoutQDt = ifft(fft(Eout21).*exp(-1j*2*pi*f*T));
EoutQDf = exp(-1j*pi/4).*Eout22;

EoutDemodQ1 = ExpCoef.*(cos(CoupConst).*EoutQDt + MutCoef*1j*sin(CoupConst).*EoutQDf);  
EoutDemodQ2 = ExpCoef.*(cos(CoupConst).*EoutQDf + MutCoef*1j*sin(CoupConst).*EoutQDt); 

EoutDemodQ = EoutQDf.*conj(EoutDemodQ2) - EoutQDt.*conj(EoutDemodQ1);
EoutDemodQ = EoutDemodQ./max(EoutDemodQ);

[vv1,vv2,~]             = mz_interf4(E_rec2.*(10^(-3.0/20)),fftshift(f),Rb,-1*pi/4);
comp_Q1 = vv1.*conj(vv1);
comp_Q2 = vv2.*conj(vv2);
comp_Q = comp_Q2-comp_Q1;
comp_Q = comp_Q./max(comp_Q);
% comp_Q = ifft(fft(comp_Q).*BitFilt);

figure;plot(t,Isig1./max(Isig1),t,Qsig1./max(Qsig1),t,real(EoutDemodI)./max(real(EoutDemodI)),t,real(EoutDemodQ)./max(real(EoutDemodQ)));axis([0 8e-10 -1.1 1.1]);

figure;plot(t,Isig1./max(Isig1),t,Qsig1./max(Qsig1),t,EoutDemodI,t,EoutDemodQ);axis([0 8e-10 -1.1 1.1]);
figure;plot(t,Isig1,t,EoutDemodI,t,comp_I,t,Qsig1,t,EoutDemodQ,t,comp_Q);

Olho(comp_I,T,NPPB,1);
Olho(EoutDemodI,T,NPPB,1);
Olho(comp_Q,T,NPPB,1);
Olho(EoutDemodQ,T,NPPB,1);

a=1;