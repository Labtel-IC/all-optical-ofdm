close all;clear;clc;
CoupConst = 3*pi/4;    %Couples Constant of the DI
ExpCoef   = exp(-1j*1);%Exponential component from Ramaswami equation
MutCoef   = 1;         %Coeficient to control polarity
E0 = 0;
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

% CW= cos(2*pi*fc*t) + 1j*sin(2*pi*fc*t);
CW= exp(-1j*2*pi*fc*t);
% CW= cos(2*pi*fc*t);

%%             dados
U_pi1 = 3.8;
Polirized=1;
TxData = zeros(1,2*Nb);
TxDataPos = linspace(1,2*Nb,2*Nb);
TxData(~mod(TxDataPos,2))=1;
TxData(end/2 - 1) = 1;
TxData = (randi(2,1,2*Nb)-1);
% TxData(1:8)= [0  0];
% TxDataMat(kk,:) = rectpulse(TxData,NPPB);
SigPam = Maping4Pam(TxData,NPPB,Polirized,U_pi1/8);
figure;plot(t,SigPam);

BWD     = 1.5*12.5e9; %Badnwidth for the bit conformation
CenFeq  = 0;    %Center frequency for the bit conformation
FiltOrd = 1;    %The order for the filter which will conform data

[BitFilt,~] = FiltroGaussiano(f,BWD,CenFeq,FiltOrd);                                 %Creating filter for conformation of the input information
BitFilt = fftshift(BitFilt);
TxSig = ifft(fft(SigPam).*BitFilt);
Olho(TxSig,T,NPPB,1);
% TxSigMat(kk,:) = TxSig;
%%
MZ_Input_File = 2;
U.U1t = TxSig;
U.U2t = exp(-1j*pi).*TxSig;
[EoutModAux,~]=Mach_Zehnder_Modulator_simplificado(t,CW,U,MZ_Input_File);   %Modulating all carriers at the same time for simplicity
figure;plot(f,20.*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))))
EoutDemod = EoutModAux.*conj(EoutModAux);
figure;plot(f,20.*log10(abs(fftshift(fft(EoutDemod)./length(EoutDemod)))))
Olho((EoutDemod),T,NPPB,1);
% EoutModAux = ifft(fft(EoutModAux).*BitFilt);
% EoutModAux = EoutModAux+abs(min(real(EoutModAux)))+U_pi1;
% figure;plot(f,20.*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))))
% EoutRec = EoutModAux.*conj(EoutModAux);
% figure;plot(f,20.*log10(abs(fftshift(fft(EoutRec)./length(EoutRec)))))
% EoutRec = EoutRec - mean(EoutRec);
% figure;plot(f,20.*log10(abs(fftshift(fft(EoutRec)./length(EoutRec)))))
% EoutRec = EoutRec.*cos(2*pi*Rb*t);%exp(1j*2*pi*Rb*t);%
% figure;plot(f,20.*log10(abs(fftshift(fft(EoutRec)./length(EoutRec)))))
% EoutRec = ifft(fft(EoutRec).*BitFilt);
% figure;plot(f,20.*log10(abs(fftshift(fft(EoutRec)./length(EoutRec)))))
% EoutRec = 1.*EoutRec./max(EoutRec);
% Olho((EoutRec),T,NPPB,1);
% 
% Olho(real(EoutModAux),T,NPPB,1);
% Olho(abs([(real(EoutModAux)+max(real(EoutModAux)))+1j.*(imag(EoutModAux)+max(imag(EoutModAux)))]),T,NPPB,1);
% Olho(real(EoutModAux),T,NPPB,1);
% Olho(imag(EoutModAux),T,NPPB,1);
% A=1;
% Olho(abs(EoutRec),T,NPPB,1);
% %%  Encontrando Níveis de Decisão
% ErecWindowed = [];
% for kk=1:NPPB:length(EoutRec)
%     ErecWindowed = [ErecWindowed EoutRec(kk+NPPB/4:(kk-1)+NPPB-(NPPB/4))];
% end
Interval = linspace(min(EoutDemod),max(EoutDemod),NPPB);
EyeMax = hist(EoutDemod,Interval);
figure;bar(Interval,EyeMax)
Intervalaux = [0 Interval 0];
EyeMaxaux = [0 EyeMax 0];
% [~,EyeLoc] = findpeaks(EyeMaxaux,'MinPeakHeight',0.5*max(EyeMaxaux));
[~,EyeLoc] = findpeaks(EyeMaxaux,'MinPeakDistance',NPPB/8,'SortStr','descend','NPeaks',4);
figure;hold all;grid on;plot(EyeMax);plot(EyeLoc-1,EyeMax(EyeLoc-1),'xr');
Levels = [Interval(EyeLoc(1)-1) Interval(EyeLoc(2)-1) Interval(EyeLoc(3)-1) Interval(EyeLoc(4)-1)];
Levels = sort(Levels);
LevelDec1 = Levels(1) + (Levels(2)-Levels(1))/2 ;
LevelDec2 = Levels(2) + (Levels(3)-Levels(2))/2 ;
LevelDec3 = Levels(3) + (Levels(4)-Levels(3))/2 ;
%%   Recebendo os bits
DataRec = [];
figure;hold on;grid on;
plot(t,EoutDemod,t,TxSig./max(TxSig));
for kk=1:NPPB:length(EoutDemod)
%     plot(t(kk+NPPB/8:(kk-1)+NPPB-(NPPB/8)),EoutDemod(kk+NPPB/8:(kk-1)+NPPB-(NPPB/8)));
%     plot(t(kk+NPPB/2),EoutDemod(kk+NPPB/2),'xr');
    plot(t(kk:(kk-1)+NPPB),EoutDemod(kk:(kk-1)+NPPB),'xr');
%     aux1 = EoutDemod(kk+NPPB/8:(kk-1)+NPPB-(NPPB/8));
%     aux1 = EoutDemod(kk+NPPB/2);
    aux1 = EoutDemod(kk:(kk-1)+NPPB);
    MeanRec = mean(aux1);
    if MeanRec <= LevelDec1
        DataRec = [DataRec 0 1];
    elseif (MeanRec <= LevelDec2)&&(MeanRec > LevelDec1)
        DataRec = [DataRec 0 0];
    elseif (MeanRec <= LevelDec3)&&(MeanRec > LevelDec2)
        DataRec = [DataRec 1 0];
    elseif MeanRec > LevelDec3
        DataRec = [DataRec 1 1];
    else
        DataRec = [DataRec 0 0];
    end
end

figure;plot(linspace(0,FinalTime,2*Nb),TxData,linspace(0,FinalTime,2*Nb),DataRec);
BitErr = sum(xor(TxData,DataRec));
Ber = BitErr/(2*Nb)
% DataRecCom = rectpulse(DataRec,NPPB);
% LevelDec4 = Levels(1) + (Levels(2)-Levels(1))/2 ;
a=1;