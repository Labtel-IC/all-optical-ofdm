close all;clear;clc;
CoupConst = 3*pi/4;    %Couples Constant of the DI
ExpCoef   = exp(-1j*1);%Exponential component from Ramaswami equation
MutCoef   = 1;         %Coeficient to control polarity
E0 = 0;
NrzMin = -1;
DatMin = 0;
DatMax = 1;
NrzMax = 1;
IQAmp = 0.95;
MZ_Input_File = 1;
p = pwd;                    %Geting the current path of the main program
Local = [p '\input_files\'];%variable that shows where the input files for
                                                     %the mzm will be saved

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
tdqpsk       = linspace(0,FinalTime,round(2*NumberOf_T*(T/Ta)));                         %Time vector for simulation
f            = time2freq(t); 
%%     Variable parameters for evaluation
Varm1_steps = 1;%10+1;%Number of different amplitudes on first arm of the MZM 
Varm2_steps = 1;%10+1;%Number of different amplitudes on second arm of the MZM 
Gain_steps  = 1;%10+2;%Number of different gains used in the optical ring
phase_steps = 1;%20+2;%Number of different phase delays that will be used
V1bias_steps = 1;%10+1;%Number of different Bias voltage that will be used
V2bias_steps = 1;%10+1;%Number of different Bias voltage that will be used
phi_0_steps = 1;
%
%IMPORTANT: The RF amplitude Vpp must respect the device limits!
Varm1         = 1;%linspace(-3.8/2,3.8/2,Varm1_steps);                        %Variation of Amplitude for Arm 1
Varm2         = 1;%linspace(-3.8/2,3.8/2,Varm2_steps);                        %Variation of Amplitude for Arm 2
phase         = pi/2;%linspace(0,pi,phase_steps);                                %Variation of phase betwee arms. [pi/12 pi/4 pi/2 2*pi/3];%Phase Vector;
Vbias         = 3.8;
Vpi2          = 3.8;
V1bias         = 3.1/2;%linspace(-7,7,Vbias_steps);                                %Variation of bias voltage on the simulation%[0.8 1.8 2.8 3.8];%Vbias Vector
V2bias         = 3.25/2;%linspace(-7,7,Vbias_steps);                                %Variation of bias voltage on the simulation%[0.8 1.8 2.8 3.8];%Vbias Vector
phi_0_vet      = 1.85;
Gain_vet      = 0;%linspace(1,2,Gain_steps);                                  %Variation of Gain inside the optical ring

p = pwd;                    %Geting the current path of the main program
Local = [p '\input_files\'];%variable that shows where the input files for
                                                     %the mzm will be saved
Local_Dp = Local;

%%      Parameters of the Continous Wave Lenght
Eo_CW              = 1; %Amplitude of the CW signal
% CW = ones(1,length(t));
CW= cos(2*pi*fc*t) + 1j*sin(2*pi*fc*t);            %Continous Wave (laser)

MZ_Input_Data_DP;

Make_MZ_Input_Files_DP;

TxData  = (randi(2,1,2*Nb)-1);
PreBits = [0 0];
% TxData = zeros(1,2*Nb);
TxPoss = linspace(1,length(TxData),length(TxData));
% TxData(~mod(TxPoss,2))=1;
% PreBits = [1 0];
% TxDataSeg = [0 0;%01
%              0 0;%02
%              0 1;%03
%              1 0;%04
%              1 1;%05
%              1 0;%06
%              0 0;%07
%              1 1;%08
%              0 1;%09
%              1 1;%10
%              1 0;%11
%              0 0;%12
%              1 0;%13
%              0 1;%14
%              1 1;%15
%             0 1];%16
% TxData = [];
% PreBits = [0 0];
% countseg = 1;
% for kk=1:Nb
%     TxData = [TxData TxDataSeg(countseg,:)];
%     countseg = countseg + 1;
%     if countseg > 16
% %         countseg = 1;
%         TxData = [TxData ones(1,2*Nb-32)];
%         break;
%     end
% end
% TxData(1)=1;
% TxData = ones(1,2*Nb);
% TxData(1:2)=[0 1];
% PreBits = [0 0];
% TxData(end/2)=0;
% TxData(end/2 + 1)=0;
TxDataComp = rectpulse(TxData,NPPB);
[DatumI,DatumQ] = DqpskEncodEq(TxData,PreBits);
u = TxData(mod(TxPoss,2)==1);
v = TxData(~mod(TxPoss,2));
[I,Q] = DQPSK_encoder(u,v,PreBits(1),PreBits(2));
% Justifing the information acordingly
TxDataIres = rectpulse(DatumI,NPPB);                                 %Changing the length of the Data acordingly with the time vector
TxDataQres = rectpulse(DatumQ,NPPB);                                 %Changing the length of the Data acordingly with the time vector
I = rectpulse(I,NPPB);                                 %Changing the length of the Data acordingly with the time vector
Q = rectpulse(Q,NPPB);                                 %Changing the length of the Data acordingly with the time vector
%Creating a NRZ Pulses
TxDataIres(TxDataIres==DatMin) = NrzMin;                                     %Where the bits are zero it will be changed to -1
TxDataIres(TxDataIres==DatMax) = NrzMax;                                     %Where the bits are one it will be changed to 1
TxDataQres(TxDataQres==DatMin) = NrzMin;                                     %Where the bits are zero it will be changed to -1
TxDataQres(TxDataQres==DatMax) = NrzMax;                                     %Where the bits are one it will be changed to 1
figure;plot(t,I,t,TxDataIres);
figure;plot(t,Q,t,TxDataQres);
Isig = TxDataIres.*IQAmp;
Qsig = TxDataQres.*IQAmp;

% Isig     = Varm1(1)*sin(2*pi*fc*t);               %The first RF signal at the Uper arm.
% Qsig     = Varm2(1)*sin(2*pi*fc*t + phase(1));
% Isig     = zeros(1,length(t));               %The first RF signal at the Uper arm.
% Qsig     = zeros(1,length(t));
U.U1t = Isig;
U.U2t = Qsig;
[E,Ei,Eq]=qi_MZM(CW, Isig, Qsig,Vbias,Vpi2);
[Eout,~] = IqModTest (t,CW,Isig,Qsig,Vpi2,Vbias);
E = E./max(abs(E));
Eout = Eout./max(abs(Eout));
figure
plot(t,abs(Eout),t,abs(E),t,real(Eout),t,real(E),t,imag(Eout),t,imag(E));
Ui = real(exp(1j*pi/4).*Eout.*ifft(fft(Eout).*exp(-1j*2*pi*f*T)));
Qi = imag(exp(1j*pi/4).*Eout.*ifft(fft(Eout).*exp(-1j*2*pi*f*T)));
figure;hold all;grid on;
plot(t,TxDataIres,t,abs(Ui),t,TxDataQres,t,abs(Qi));
BitFilBanWid  = 10e9; %Badnwidth for the bit conformation
BitFiltCenFre = 0;    %Center frequency for the bit conformation
BitFilOrd     = 1;    %The order for the filter which will conform data
[BitFilt,~] = FiltroGaussiano(f,BitFilBanWid,BitFiltCenFre,BitFilOrd);                                 %Creating filter for conformation of the input information
BitFilt = fftshift(BitFilt);

E_rec = Eout;
% E_rec = E;
[uu1,uu2,~]             = mz_interf4(E_rec.*(10^(-3.0/20)),fftshift(f),Rb,pi/4);
comp_I1 = uu1.*conj(uu1);
comp_I2 = uu2.*conj(uu2);
comp_I = comp_I2-comp_I1;
% comp_I = ifft(fft(comp_I).*BitFilt);
[vv1,vv2,~]             = mz_interf4(E_rec.*(10^(-3.0/20)),fftshift(f),Rb,-1*pi/4);
comp_Q1 = vv1.*conj(vv1);
comp_Q2 = vv2.*conj(vv2);
comp_Q = comp_Q2-comp_Q1;
% comp_Q = ifft(fft(comp_Q).*BitFilt);
figure;plot(t,TxDataIres,t,comp_I)
Olho(comp_I,Ta,NPPB,1);
figure;plot(t,TxDataQres,t,comp_Q)
Olho(comp_Q,Ta,NPPB,1);
% [Eout,~] = Mach_Zehnder_Modulator_DP(t,CW,U,MZ_Input_File,Local);
figure;
hold all;
grid on;
plot(t,abs(Eout),t,real(Eout),t,imag(Eout));

figure; hold all; grid on;
axisaux = 0.5.*ones(1,length(DatumI));
plot(DatumI,axisaux,'xb');
plot(axisaux,DatumQ,'xb');

% Eout = E;
Eout11 = ExpCoef.*(cos(CoupConst).*Eout + MutCoef*1j*sin(CoupConst).*E0);        
Eout12 = ExpCoef.*(cos(CoupConst).*E0 + MutCoef*1j*sin(CoupConst).*Eout);  

    
Eoutaux1 = Eout./2;
Eoutaux2 = Eout./2;
EoutRec = Eout.*conj(Eout);
EoutRec = ifft(fft(EoutRec).*BitFilt);
% EoutRec = EoutRec./max(EoutRec);
Olho(EoutRec,Ta,NPPB,1);

% EoutDt = ifft(fft(Eoutaux1).*exp(-1j*2*pi*f*T));
% EoutDf = exp(1j*pi/4).*Eoutaux2;
EoutDt = ifft(fft(Eout11).*exp(1j*2*pi*f*T));
EoutDf = exp(-1j*pi/4).*Eout12;
EoutDf2 = exp(1j*pi/4).*Eout12;


E11 = ExpCoef.*(cos(CoupConst).*EoutDt + MutCoef*1j*sin(CoupConst).*EoutDf);  
E112 = ExpCoef.*(cos(CoupConst).*EoutDt + MutCoef*1j*sin(CoupConst).*EoutDf2);        
E12 = ExpCoef.*(cos(CoupConst).*EoutDf + MutCoef*1j*sin(CoupConst).*EoutDt);       
E122 = ExpCoef.*(cos(CoupConst).*EoutDf2 + MutCoef*1j*sin(CoupConst).*EoutDt);  
EoutD = EoutDt + EoutDf;  
EoutD2 = EoutDt + EoutDf2;

% EoutDrec = EoutD.*conj(EoutD);
EoutDrec1 = E11.*conj(E11);
EoutDrec12 = E112.*conj(E112);
EoutDrec2 = E12.*conj(E12);
EoutDrec22 = E122.*conj(E122);
EoutDrec = EoutDrec2-EoutDrec1;
EoutDrec2 = EoutDrec22-EoutDrec12;
% EoutDrec = ifft(fft(EoutDrec).*BitFilt);
% EoutDrec2 = ifft(fft(EoutDrec2).*BitFilt);
EoutDrec = EoutDrec./max(EoutDrec);
EoutDrec2 = EoutDrec2./max(EoutDrec2);
figure;plot(t,TxDataIres,t,EoutDrec,t,TxDataIres,t,EoutDrec2)
Olho(EoutDrec,Ta,NPPB,1);
Olho(EoutDrec2,Ta,NPPB,1);
figure;plot(t,TxDataIres,t,EoutDrec,t,comp_I./max(comp_I),t,TxDataQres,t,EoutDrec2,t,comp_Q./max(comp_Q));
figure;plot(t,TxDataQres,t,EoutDrec2,t,comp_Q./max(comp_Q));
for kk=1:NPPB:length(EoutDrec)%length(Irec)%
    MeanI=mean(EoutDrec(kk:(kk-1)+NPPB));
    MeanQ=mean(EoutDrec2(kk:(kk-1)+NPPB));
%     MeanQ=mean(Qrec(kk:(kk-1)+NPPB));
    figure(55);hold all;grid on;
    plot(t(kk:(kk-1)+NPPB),EoutDrec(kk:(kk-1)+NPPB),t(kk:(kk-1)+NPPB),TxDataIres(kk:(kk-1)+NPPB));
    plot(t(kk:(kk-1)+NPPB),rectpulse(MeanI,NPPB));
    figure(44);hold all;grid on;
    plot(t(kk:(kk-1)+NPPB),EoutDrec2(kk:(kk-1)+NPPB),t(kk:(kk-1)+NPPB),TxDataQres(kk:(kk-1)+NPPB));
    plot(t(kk:(kk-1)+NPPB),rectpulse(MeanQ,NPPB));
    a=0;
end
% plot(t,abs(Eout1),t,real(Eout1),t,imag(Eout1));
% plot(t,abs(Eout2),t,real(Eout2),t,imag(Eout2));
% figure;plot(real(Eout1),imag(Eout1))
%%
thisone = 1;
if thisone
    %%
    TimDel = T;
    PhaDel = 1*pi/4;
    [ EI1,EI2 ] = DelayInterf( t,TimDel,PhaDel,Eout);
    %%
    TimDel = T;
    PhaDel = -1*pi/4;
    [ EQ1,EQ2 ] = DelayInterf( t,TimDel,PhaDel,Eout);
    %%
    BitFilBanWid  = 10e9; %Badnwidth for the bit conformation
    BitFiltCenFre = 0;    %Center frequency for the bit conformation
    BitFilOrd     = 1;    %The order for the filter which will conform data
    [BitFilt,~] = FiltroGaussiano(f,BitFilBanWid,BitFiltCenFre,BitFilOrd);                                 %Creating filter for conformation of the input information
    BitFilt = fftshift(BitFilt);
    %%
    I1x=EI1.*conj(EI1);
    I2x=EI2.*conj(EI2);
    I1x = ifft(fft(I1x).*BitFilt);
    I2x = ifft(fft(I2x).*BitFilt);
    I1x = I1x - min(I1x);
    I2x = I2x - min(I2x);
    %%
    Q1x=EQ1.*conj(EQ1);
    Q2x=EQ2.*conj(EQ2);
    Q1x = ifft(fft(Q1x).*BitFilt);
    Q2x = ifft(fft(Q2x).*BitFilt);
    Q1x = Q1x - min(Q1x);%min(abs(Q1x));%
    Q2x = Q2x - min(Q2x);%min(abs(Q2x));%
    %%
    A = I1x;
    B = I2x;
    C = Q1x;
    D = Q2x;
    Irec = B - A;
    Qrec = D - C;
    Irec = Irec./max(abs(Irec));
    Qrec = Qrec./max(abs(Qrec));
    Imean = [];
    Qmean = [];
    for kk=1:NPPB:length(Irec)
        auxI=Irec(kk:(kk-1)+NPPB);
        auxQ=Qrec(kk:(kk-1)+NPPB);
        Imean = [Imean mean(auxI)];
        Qmean = [Qmean mean(auxQ)];
    end
    figure;hold all;grid on;
    plot(Imean);
    plot(Qmean);
    %%
%     figure;plot(t,abs(Irec),t,real(Irec),t,imag(Irec));
    Olho(Irec,T,NPPB,1);
%     figure;plot(t,abs(Qrec),t,real(Qrec),t,imag(Qrec));
    Olho(Qrec,T,NPPB,1);
    BitAnt = PreBits;
    DataRec = [];
    
    for kk=1:NPPB:length(Qrec)%length(Irec)%
        if BitAnt == [0 0]
            ThisTag = 10;
        elseif BitAnt == [0 1]
            ThisTag = 40;
        elseif BitAnt == [1 1]
            ThisTag = 30;
        elseif BitAnt == [1 0]
            ThisTag = 20;
        else
            ThisTag = 00;
        end
        MeanI=mean(Irec(kk:(kk-1)+NPPB));
        MeanQ=mean(Qrec(kk:(kk-1)+NPPB));
        figure(66);hold all;grid on;
        plot(t(kk:(kk-1)+NPPB),Irec(kk:(kk-1)+NPPB));
        plot(t(kk:(kk-1)+NPPB),rectpulse(MeanI,NPPB));
        figure(77);hold all;grid on;
        plot(t(kk:(kk-1)+NPPB),rectpulse(MeanQ,NPPB));
        plot(t(kk:(kk-1)+NPPB),Qrec(kk:(kk-1)+NPPB));
%         MeanI=mean(Qrec(kk:(kk-1)+NPPB));
        if MeanI>=0.75
            ThisTag = ThisTag + 3;
        elseif (MeanI<0.75)&&(MeanI>=0.25)
            ThisTag = ThisTag + 4;
        elseif (MeanI<0.25)&&(MeanI>=-0.25)
            ThisTag = ThisTag + 1;
        elseif MeanI<-0.25
            ThisTag = ThisTag + 2;
        else
            ThisTag = ThisTag + 0;
        end
        switch ThisTag
            case 11
                DataRec = [DataRec 0 0];
                BitAnt  = [1 1];
            case 12
                DataRec = [DataRec 0 1];
                BitAnt  = [0 1];
            case 13
                DataRec = [DataRec 1 1];
                BitAnt  = [0 0];
            case 14
                DataRec = [DataRec 1 0];
                BitAnt  = [1 0];
            case 21
                DataRec = [DataRec 0 0];
                BitAnt  = [1 0];
            case 22
                DataRec = [DataRec 0 1];
                BitAnt  = [1 1];
            case 23
                DataRec = [DataRec 1 1];
                BitAnt  = [0 1];
            case 24
                DataRec = [DataRec 1 0];
                BitAnt  = [0 0];
            case 31
                DataRec = [DataRec 0 0];
                BitAnt  = [0 0];
            case 32
                DataRec = [DataRec 0 1];
                BitAnt  = [1 0];
            case 33
                DataRec = [DataRec 1 1];
                BitAnt  = [1 1];
            case 34
                DataRec = [DataRec 1 0];
                BitAnt  = [0 1];
            case 41
                DataRec = [DataRec 0 0];
                BitAnt  = [0 1];
            case 42
                DataRec = [DataRec 0 1];
                BitAnt  = [0 0];
            case 43
                DataRec = [DataRec 1 1];
                BitAnt  = [1 0];
            case 44
                DataRec = [DataRec 1 0];
                BitAnt  = [1 1];
            otherwise
        end
    end
    DataRecAux = rectpulse(DataRec,NPPB);
    figure;hold all;grid on;
    plot(tdqpsk,TxDataComp);
    plot(tdqpsk,DataRecAux);
    a=0;
else
    %%
    CoupConst = 3*pi/4;    %Couples Constant of the DI
    ExpCoef   = exp(-1j*1);%Exponential component from Ramaswami equation
    MutCoef   = 1;         %Coeficient to control polarity
    E0 = 0;
    %%
    % E11 = ExpCoef.*(cos(CoupConst).*Eout + MutCoef*1j*sin(CoupConst).*E0);        
    % E12 = ExpCoef.*(cos(CoupConst).*E0 + MutCoef*1j*sin(CoupConst).*Eout);
    E11 = Eout./2;
    E12 = Eout./2;

    TimDel = T;
    % aux1 = E12(t<=TimDel);
    % E12d = [E12(end-(length(aux1)-1):end) E12(1:end-(length(aux1)))];
    E12d = ifft(fft(E12).*exp(-1j*2*pi*TimDel.*f));

    PhaDel = 1*pi/4;
    E11d = E11.*exp(-1j*PhaDel);

    EI1 = E12d + E11d;
    EI2 = exp(-1j*pi).*EI1;
    % EI1 = ExpCoef.*(cos(CoupConst).*E11d + MutCoef*1j*sin(CoupConst).*E12d);        
    % EI2 = ExpCoef.*(cos(CoupConst).*E12d + MutCoef*1j*sin(CoupConst).*E11d);
    % EI1 = ExpCoef.*(cos(CoupConst).*E12d + MutCoef*1j*sin(CoupConst).*E11d);        
    % EI2 = ExpCoef.*(cos(CoupConst).*E11d + MutCoef*1j*sin(CoupConst).*E12d);
    %%
    % E11 = ExpCoef.*(cos(CoupConst).*Eout + MutCoef*1j*sin(CoupConst).*E0);        
    % E12 = ExpCoef.*(cos(CoupConst).*E0 + MutCoef*1j*sin(CoupConst).*Eout);
    E11 = Eout./2;
    E12 = Eout./2;

    TimDel = T;
    % aux1 = E12(t<=TimDel);
    % E12d = [E12(end-(length(aux1)-1):end) E12(1:end-(length(aux1)))];
    E12d = ifft(fft(E12).*exp(-1j*2*pi*TimDel.*f));

    PhaDel = -1*pi/4;
    E11d = E11.*exp(-1j*PhaDel);

    EQ1 = E12d + E11d;
    EQ2 = exp(-1j*pi).*EQ1;
    % EQ1 = ExpCoef.*(cos(CoupConst).*E11d + MutCoef*1j*sin(CoupConst).*E12d);        
    % EQ2 = ExpCoef.*(cos(CoupConst).*E12d + MutCoef*1j*sin(CoupConst).*E11d);
    % EQ1 = ExpCoef.*(cos(CoupConst).*E12d + MutCoef*1j*sin(CoupConst).*E11d);        
    % EQ2 = ExpCoef.*(cos(CoupConst).*E11d + MutCoef*1j*sin(CoupConst).*E12d);
    %%
    EoutI1Demd = EI1.*conj(EI1);%abs(EI1).^2;%Detection for one signal I at Uper arm
    EoutI2Demd = EI2.*conj(EI2);%abs(EI2).^2;%Detection for one signal I at Dow arm

    EoutQ1Demd = EQ1.*conj(EQ1);%abs(EQ1).^2;%Detection for one signal Q at Uper arm
    EoutQ2Demd = EQ2.*conj(EQ2);%abs(EQ2).^2;%Detection for one signal Q at Dow arm
    %%
    EoutA = EoutI2Demd;%Detection for one signal I at Uper arm
    EoutB = 0;%EoutI1Demd;%Detection for one signal I at Dow arm
    EoutC = EoutQ2Demd;%Detection for one signal Q at Uper arm
    EoutD = 0;%EoutQ1Demd;%Detection for one signal Q at Dow arm

    EoutIDemd  = EoutA - EoutB;
    EoutQDemd  = EoutC - EoutD;
    %%
    EoutI = EoutIDemd;
    EoutQ = EoutQDemd;
    % Normalazing the vectors
    EoutI = EoutI./max(abs(EoutI));
    EoutQ = EoutQ./max(abs(EoutQ));
    figure;
    hold all;
    grid on;
    plot(t,Isig./max(Isig))
    plot(t,abs(EoutI),t,real(EoutI),t,imag(EoutI))
    plot(t,Qsig./max(Qsig))
    plot(t,abs(EoutQ),t,real(EoutQ),t,imag(EoutQ))
end
a=1;