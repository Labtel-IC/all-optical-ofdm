%c
%c                                                       ..'Ž`'..'Ž`..'Ž`..                                                   
%c       File: DatumReceptionInputDAta File With the Main Script Parameters 
%c(Main fail to call each idividual parameters)
%c
%c     This main code is resposible to  load all necessary parameters of 
%c the MAIN script, just in this file changes should be done. Do not change
%c the MAIN file. 
%c
%c      
%c
%c                                           by P.Marciano LG
%c                                           28/10/2017
%c                                           pablorafael.mcx@gmail.com
%c 
%c     References:
%c@article{hillerkuss2010simple,
%c  title={Simple all-optical FFT scheme enabling Tbit/s real-time signal processing},
%c  author={Hillerkuss, D and Winter, M and Teschke, M and Marculescu, A and Li, J and Sigurdsson, G and Worms, K and Ezra, S Ben and Narkiss, N and Freude, W and others},
%c  journal={Optics express},
%c  volume={18},
%c  number={9},
%c  pages={9324--9340},
%c  year={2010},
%c  publisher={Optical Society of America}
%c}
%c@article{kim2002chirp,
%c  title={Chirp characteristics of dual-drive. Mach-Zehnder modulator with a finite DC extinction ratio},
%c  author={Kim, Hoon and Gnauck, Alan H},
%c  journal={IEEE Photonics Technology Letters},
%c  volume={14},
%c  number={3},
%c  pages={298--300},
%c  year={2002},
%c  publisher={IEEE}
%c}
%c
%c@phdthesis{togneri2005analise,
%c  title={An{\'a}lise de Sistemas de Multiplexa{\c{c}}{\~a}o por Subportadora-SCM},
%c  author={Togneri, Arnaldo Paterline},
%c  year={2005},
%c  school={UNIVERSIDADE FEDERAL DO ESP{\'I}RITO SANTO}
%c}
%c
%c@article{oliveiralarge,
%c  title={Large Signal Analysis of Mach-Zehnder Modulator Intensity Response in a Linear Dispersive Fiber},
%c  author={Oliveira, JMB and Salgado, HM and Rodrigues, MRD}
%c}
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                    G L O B A L  V A R I A B L E S                      %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c   
%c
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                   L I S T  O F  V A R I A B L E S                      %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c   
%c
%c
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                             S T R U C T S                              %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c   Here it will be added the functions that this code will call
%c   
%c
%c
%c
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Start of the Program                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    Control Variables of this program
close all;clear;clc;
OFF = 0;
ON = 1;
Ploting = ON;    %It is used to select an ploting.
PlotingThis = OFF; %It is used to select an especific ploting among others
Elucidation = OFF;%It is used to enable the OFFT example
UsePhotDiod = OFF;%To select or not the usage of the photodiod
BitFilBanWid  = 6.25e9; %Badnwidth for the bit conformation
BitFiltCenFre = 0;%7*12.5e9;    %Center frequency for the bit conformation
BitFilOrd     = 1;    %The order for the filter which will conform data
ThisCarr      = 4; %This variable will be used to select the carrier
%%
fc          = 12.5e9;%2^30;%center frequency of Electical signal ^10 kilo ^20 mega ^30 giga
fc2=fc;
f_max        = (2^2)*fc;                                                      %max frequency that matlab can reproduce
fsample      = (2^2)*f_max;                                                   %Sampling frequency
T            = 1/fc;                                                       %period of the filter
Ta           = 1/fsample;                                                       %period of the filter
Rb           = fc;                                                     %bit rate
% Rb           = fc;                                                         %bit rate
Tb           = 1/Rb;                                                       %time of one bit NRZ
NPPB         = 2^5;                                            %Number of Point Per Bit
NumberOf_T   = 2^10;                                                        %number of periods on the simulation 2^18 is the mï¿½ximum value for numberEyes = 2
FinalTime    = NumberOf_T*T;                                               %Final time of simulation
Nb           = floor(FinalTime/Tb);                                               %number of bits
% TotalSamples = FinalTime*fsample;                                          %Total Number of Samples
TotalSamples = NPPB*NumberOf_T*(T/Ta);                                                   %Total Number of Samples

t            = linspace(0,FinalTime,round(TotalSamples));                         %Time vector for simulation
f            = time2freq(t); 
tb=t;
t2=t;
f2=f;
fb=f;

A = exp(-1j*2*pi*fc2*t2);
B = exp(-1j*4*pi*fc2*t2);
C = exp(-1j*6*pi*fc2*t2);
D = exp(-1j*8*pi*fc2*t2);
E = exp(-1j*10*pi*fc2*t2);
F = exp(-1j*12*pi*fc2*t2);
G = exp(-1j*14*pi*fc2*t2);
H = exp(-1j*16*pi*fc2*t2);
I = exp(-1j*18*pi*fc2*t2);
J = exp(-1j*20*pi*fc2*t2);
K = exp(-1j*22*pi*fc2*t2);
L = exp(-1j*24*pi*fc2*t2);
M = exp(-1j*26*pi*fc2*t2);
N = exp(-1j*28*pi*fc2*t2);
O = exp(-1j*30*pi*fc2*t2);
P = exp(-1j*32*pi*fc2*t2);
Q = exp(-1j*34*pi*fc2*t2);
R = exp(-1j*36*pi*fc2*t2);
S = exp(-1j*38*pi*fc2*t2);
T = exp(-1j*40*pi*fc2*t2);
U = exp(-1j*42*pi*fc2*t2);
V = exp(-1j*44*pi*fc2*t2);
W = exp(-1j*46*pi*fc2*t2);
X = exp(-1j*48*pi*fc2*t2);
Y = exp(-1j*50*pi*fc2*t2);
Z = exp(-1j*52*pi*fc2*t2);
AA = exp(-1j*54*pi*fc2*t2);
AB = exp(-1j*56*pi*fc2*t2);
AC = exp(-1j*58*pi*fc2*t2);
AD = exp(-1j*60*pi*fc2*t2);
AE = exp(-1j*62*pi*fc2*t2);
AF = exp(-1j*64*pi*fc2*t2);
COMB = A+B+C+D+E+F+G+H+I+J+K+L+M+N+O+P+Q+R+S+T+U+V+W+X+Y+Z+AA+AB+AC+AD+AE+AF;
COMB = COMB + wgn(1,length(t2),-1*40);
figure;hold on
plot(f2,db(abs(fftshift(fft(COMB)./length(COMB)))))
EoutT = COMB;
EoutRec = COMB;
ON = 1;
OFF = 0;
%%   Control Variables for the All Optical FFT Example
% Time dellay used on the example
%Stage 1
Del0   = T/2;%Time delay at the frist stage
%Stage 2
Del1a  = T/4;%Time delay at the second stage
Del1b  = T/4;%Time delay at the second stage
%Stage 3
Del2aa = T/8;%Time delay at the final stage
Del2ab = T/8;%Time delay at the final stage
Del2ba = T/8;%Time delay at the final stage
Del2bb = T/8;%Time delay at the final stage

% Phase dellay used on the example given
%Stage 1
Pha0   = (0)*(pi/180); %Phase dalay at the first stage
%Stage 2
Pha1a  = (90)*(pi/180);%Phase dely at the second stage
Pha1b  = 0;            %Phase delay at the second stage
%Stage 3
Pha2aa = 3*pi/4;       %Phase delay at final stage
Pha2ab = pi/4;         %Phase delay at final stage
Pha2ba = pi/2;         %Phase delay at final stage
Pha2bb = 0;            %Phase delay at final stage

% Variables for the optical coupler
% it was fundamental to use a teoretical model for the couple. Just add or
% subtract results will not work. The model used at this script was taken
% from Optical Networks by Ramaswami
kl = (135)*pi/180;    %Coupling Constant. It is the argument of the function
MutCoef = 1;          %Multiplier, it is used to sete the polarization 
ExpCoef = exp(-1j*1); %Exponential that multiplies all frequency componente

%%    Constants to controu the Optical FFT

MaxNumStag = ceil(log2(8)); %Maximum Number of sub carriers. 32 is the maximum allowed, I did not adjust the the FFT to work with more than 16
if MaxNumStag>5
error(['MaxNumStag allowed value for MaxNumStag was achieved. Please'...
                                ' cange the number Maximum of carriers.']);
end
Count = ones(MaxNumStag,1);%Control Variable for adding correct variation
E0 = 0;                    %Second input, usualy zero
ActStag = 1;               %Starting point of Optical FFT

%%      Filters For Channel Selection
% Finaly, it was created filter responsible to better receive the signal
% and eliminate any non  linearities.O
SelecBand = 2*fc; %Bandwidth of the filter
CentFreq  = fc;   %Center Frequency
SelFilOrd = 5;    %rden for the selected Filter

[ SelecFilt3aaa ] = FiltroGaussiano(f,SelecBand,CentFreq,SelFilOrd);
% [ SelecFilt1 ] = Filtro_Retangular( SelecBand,CentFreq,fb);
SelecFilt3aaa = fftshift(SelecFilt3aaa);

[ SelecFilt3baa ] = FiltroGaussiano(f,SelecBand,2*CentFreq,SelFilOrd);
% [ SelecFilt2 ] = Filtro_Retangular( SelecBand,2*CentFreq,fb);
SelecFilt3baa = fftshift(SelecFilt3baa);

[ SelecFilt3aba ] = FiltroGaussiano(f,SelecBand,3*CentFreq,SelFilOrd);
% [ SelecFilt3 ] = Filtro_Retangular( SelecBand,3*CentFreq,fb);
SelecFilt3aba = fftshift(SelecFilt3aba);

[ SelecFilt3bba ] = FiltroGaussiano(f,SelecBand,4*CentFreq,SelFilOrd);
% [ SelecFilt4 ] = Filtro_Retangular( SelecBand,4*CentFreq,fb);
SelecFilt3bba = fftshift(SelecFilt3bba);

[ SelecFilt3aab ] = FiltroGaussiano(f,SelecBand,5*CentFreq,SelFilOrd);
% [ SelecFilt5 ] = Filtro_Retangular( SelecBand,5*CentFreq,fb);
SelecFilt3aab = fftshift(SelecFilt3aab);

[ SelecFilt3bab ] = FiltroGaussiano(f,SelecBand,6*CentFreq,SelFilOrd);
% [ SelecFilt6 ] = Filtro_Retangular( SelecBand,6*CentFreq,fb);
SelecFilt3bab = fftshift(SelecFilt3bab);

[ SelecFilt3abb ] = FiltroGaussiano(f,SelecBand,7*CentFreq,SelFilOrd);
% [ SelecFilt7 ] = Filtro_Retangular( SelecBand,7*CentFreq,fb);
SelecFilt3abb = fftshift(SelecFilt3abb);

[ SelecFilt3bbb ] = FiltroGaussiano(f,SelecBand,8*CentFreq,SelFilOrd);
% [ SelecFilt8 ] = Filtro_Retangular( SelecBand,8*CentFreq,fb);
SelecFilt3bbb = fftshift(SelecFilt3bbb);



%% Simple example of how the Optical FFT was implemented
% Them main code that implements an all optical FFT was previouly described
% at those following lines is described step by step this process.
%% First Stage
% At this point the scrip is implementing an FFT for 2 subcaries. The basic
% idea is to use an interferometre Delay to create an controled
% desconstrutive interferency between the input signal;
 %split
 % Initially the signal was splited in two different path. The signal that
 % passes within the uper arm will suffe no defformation.
Eout11 = ExpCoef.*(cos(kl).*EoutRec + MutCoef*1j*sin(kl).*0);
Eout12 = ExpCoef.*(cos(kl).*0 + MutCoef*1j*sin(kl).*EoutRec);
 %Dellay and change phase
 %Mean while the signal that passes through the Lower arm will be dellayed
 %and its phase will be changed.
aux1 = Eout12(t<=Del0);
Eout12d = [Eout12(end-(length(aux1)-1):end) Eout12(1:end-(length(aux1)))];
Eout12d = Eout12d.*exp(-1j*Pha0);
%Ploting results for qualitative analizes
if Ploting
    figure;
    plot(t,abs(Eout12),t,abs(Eout12d))
end
% Thus, the signal will be coupled again. This process was implemented with
% an teoretical device.
Eout1a = ExpCoef.*(cos(kl).*Eout11 + MutCoef*1j*sin(kl).*Eout12d);
Eout1b = ExpCoef.*(cos(kl).*Eout12d + MutCoef*1j*sin(kl).*Eout11);

%Ploting results for qualitative analizes
if Ploting
    figure;
    hold all;
    plot(f,db(abs((fftshift(fft(EoutT))))/length(EoutT)));
    plot(f,db(abs((fftshift(fft(Eout1a))))/length(Eout1a)));
    plot(f,db(abs((fftshift(fft(Eout1b))))/length(Eout1b)));
    axis([1e11 3e11 -370 10]);
end
%% Second Stage
% From the first stage two signals where generated. This stage will be
% exactly as the first. With the difference that now, there will be two
% delay interferometers. It is also important to mention that there will be
% variation on the input parameters of the DI such as different time and
% phase delay from the first stage
  %% Up arm
  % This is the first Delay Interferomentro within this stage.
 %split
 % Initially the signal was splited in two different path. The signal that
 % passes within the uper arm will suffe no defformation.
Eout21a = ExpCoef.*(cos(kl).*Eout1a + MutCoef*1j*sin(kl).*0);
Eout22a = ExpCoef.*(cos(kl).*0 + MutCoef*1j*sin(kl).*Eout1a);
 %Dellay and change phase
 %Mean while the signal that passes through the Lower arm will be dellayed
 %and its phase will be changed.

%Ploting results for qualitative analizes
if Ploting
    figure;
    hold all;
    plot(f,db(abs((fftshift(fft(EoutT))))/length(EoutT)));
    plot(f,db(abs((fftshift(fft(Eout21a))))/length(Eout21a)));
    plot(f,db(abs((fftshift(fft(Eout22a))))/length(Eout22a)));
    axis([1e11 3e11 -370 10]);
end

% Eout21a = (1/2).*Eout1a;
% Eout22a = (1/2).*Eout1a;
 %Dellay and change phase
aux2a = Eout22a(t<=Del1a);
Eout22da = [Eout22a(end-(length(aux2a)-1):end) Eout22a(1:end-(length(aux2a)))];
Eout22da = Eout22da.*exp(-1j*Pha1a);
%Ploting results for qualitative analizes
if Ploting
    figure;
    plot(t,abs(Eout22a),t,abs(Eout22da))
end
%Ploting results for qualitative analizes
if Ploting
    figure;
    hold all;
    plot(f,db(abs((fftshift(fft(EoutT))))/length(EoutT)));
    plot(f,db(abs((fftshift(fft(Eout21a))))/length(Eout21a)));
    plot(f,db(abs((fftshift(fft(Eout22da))))/length(Eout22da)));
    axis([1e11 3e11 -370 10]);
end
% Eout2a = Eout21a + Eout22da;
% Eout2aa = (1/2).*Eout2a;
% Eout2ab = (1/2).*Eout2a;
% kl =(90)*pi/180;
Eout2aa = ExpCoef.*(cos(kl).*Eout21a + MutCoef*1j*sin(kl).*Eout22da);
Eout2ab = ExpCoef.*(cos(kl).*Eout22da + MutCoef*1j*sin(kl).*Eout21a);

%Ploting results for qualitative analizes
if Ploting
    figure;
    hold all;
    plot(f,db(abs((fftshift(fft(EoutT))))/length(EoutT)));
    plot(f,db(abs((fftshift(fft(Eout2aa))))/length(Eout2aa)));
    plot(f,db(abs((fftshift(fft(Eout2ab))))/length(Eout2ab)));
    axis([1e11 3e11 -370 10]);
end
  %% Down arm
  %The second signal generate by the first stage will enter here
 %split
 %As previouly mentioned. The frist part is to split the input signal in
 %two parts that will travel through different paths.
Eout21b = ExpCoef.*(cos(kl).*Eout1b + MutCoef*1j*sin(kl).*0);
Eout22b = ExpCoef.*(cos(kl).*0 + MutCoef*1j*sin(kl).*Eout1b);
% Eout21b = (1/2).*Eout1b;
% Eout22b = (1/2).*Eout1b;
 %Dellay and change phase
aux2b = Eout22b(t<=Del1b);
Eout22db = [Eout22b(end-(length(aux2b)-1):end) Eout22b(1:end-(length(aux2b)))];
Eout22db = Eout22db.*exp(-1j*Pha1b);
%Ploting results for qualitative analizes
if Ploting
    figure;
    plot(t,abs(Eout22b),t,abs(Eout22db))
end
% Eout2b = Eout21b + Eout22db;
% Eout2ba = (1/2).*Eout2b;
% Eout2bb = (1/2).*Eout2b;
Eout2ba = ExpCoef.*(cos(kl).*Eout21b + MutCoef*1j*sin(kl).*Eout22db);
Eout2bb = ExpCoef.*(cos(kl).*Eout22db + MutCoef*1j*sin(kl).*Eout21b);

if Ploting
    figure;
    hold all;
    plot(f,db(abs((fftshift(fft(EoutT))))/length(EoutT)));
    plot(f,db(abs((fftshift(fft(Eout2ba))))/length(Eout2ba)));
    plot(f,db(abs((fftshift(fft(Eout2bb))))/length(Eout2bb)));
    axis([1e11 3e11 -370 10]);
end
%% Third Stage
% Finaly after the second stage, now we have 4 signals that was generated
% each one would filter one sub carrier.
  %% First Arm
 %split
Eout31aa = ExpCoef.*(cos(kl).*Eout2aa + MutCoef*1j*sin(kl).*0);
Eout32aa = ExpCoef.*(cos(kl).*0 + MutCoef*1j*sin(kl).*Eout2aa);
if Ploting
    figure;
    hold all;
    plot(f,db(abs((fftshift(fft(EoutT))))/length(EoutT)));
    plot(f,db(abs((fftshift(fft(Eout31aa))))/length(Eout31aa)));
    plot(f,db(abs((fftshift(fft(Eout32aa))))/length(Eout32aa)));
    axis([1e11 3e11 -370 10]);
end
% Eout31aa = (1/2).*Eout2aa;
% Eout32aa = (1/2).*Eout2aa;
 %Dellay and change phase
aux3aa = Eout32aa(t<=Del2aa);
Eout32daa = [Eout32aa(end-(length(aux3aa)-1):end) Eout32aa(1:end-(length(aux3aa)))];
Eout32daa = Eout32daa.*exp(-1j*Pha2aa);
if Ploting
    figure;
    plot(t,abs(Eout32aa),t,abs(Eout32daa))
end
% Eout3aa = Eout31aa + Eout32daa;
% Eout3aaa = (1/2).*Eout3aa;
% Eout3aab = (1/2).*Eout3aa;
Eout3aaa = ExpCoef.*(cos(kl).*Eout31aa + MutCoef*1j*sin(kl).*Eout32daa);
Eout3aab = ExpCoef.*(cos(kl).*Eout32daa + MutCoef*1j*sin(kl).*Eout31aa);

if Ploting
    figure;
    hold all;
    plot(f,db(abs((fftshift(fft(EoutT))))/length(EoutT)));
    plot(f,db(abs((fftshift(fft(Eout3aaa))))/length(Eout3aaa)));
    plot(f,db(abs((fftshift(fft(Eout3aab))))/length(Eout3aab)));
    axis([1e11 3e11 -370 10]);
end
  %% Second Arm
 %split
Eout31ab = ExpCoef.*(cos(kl).*Eout2ab + MutCoef*1j*sin(kl).*0);
Eout32ab = ExpCoef.*(cos(kl).*0 + MutCoef*1j*sin(kl).*Eout2ab);
% Eout31ab = (1/2).*Eout2ab;
% Eout32ab = (1/2).*Eout2ab;
 %Dellay and change phase
aux3ab = Eout32ab(t<=Del2ab);
Eout32dab = [Eout32ab(end-(length(aux3ab)-1):end) Eout32ab(1:end-(length(aux3ab)))];
Eout32dab = Eout32dab.*exp(-1j*Pha2ab);
if Ploting
    figure;
    plot(t,abs(Eout32ab),t,abs(Eout32dab))
end
% Eout3ab = Eout31ab + Eout32dab;
% Eout3aba = (1/2).*Eout3ab;
% Eout3abb = (1/2).*Eout3ab;
Eout3aba = ExpCoef.*(cos(kl).*Eout31ab + MutCoef*1j*sin(kl).*Eout32dab);
Eout3abb = ExpCoef.*(cos(kl).*Eout32dab + MutCoef*1j*sin(kl).*Eout31ab);

if Ploting
    figure;
    hold all;
    plot(f,db(abs((fftshift(fft(EoutT))))/length(EoutT)));
    plot(f,db(abs((fftshift(fft(Eout3aba))))/length(Eout3aba)));
    plot(f,db(abs((fftshift(fft(Eout3abb))))/length(Eout3abb)));
    axis([1e11 3e11 -370 10]);
end
  %% Third Arm
 %split
Eout31ba = ExpCoef.*(cos(kl).*Eout2ba + MutCoef*1j*sin(kl).*0);
Eout32ba = ExpCoef.*(cos(kl).*0 + MutCoef*1j*sin(kl).*Eout2ba);
% Eout31ba = (1/2).*Eout2ba;
% Eout32ba = (1/2).*Eout2ba;
 %Dellay and change phase
aux3ba = Eout32ba(t<=Del2ba);
Eout32dba = [Eout32ba(end-(length(aux3ba)-1):end) Eout32ba(1:end-(length(aux3ba)))];
Eout32dba = Eout32dba.*exp(-1j*Pha2ba);
if Ploting
    figure;
    plot(t,abs(Eout32ba),t,abs(Eout32dba))
end
% Eout3ba = Eout31ba + Eout32dba;
% Eout3baa = (1/2).*Eout3ba;
% Eout3bab = (1/2).*Eout3ba;
Eout3baa = ExpCoef.*(cos(kl).*Eout31ba + MutCoef*1j*sin(kl).*Eout32dba);
Eout3bab = ExpCoef.*(cos(kl).*Eout32dba + MutCoef*1j*sin(kl).*Eout31ba);

if Ploting
    figure;
    hold all;
    plot(f,db(abs((fftshift(fft(EoutT))))/length(EoutT)));
    plot(f,db(abs((fftshift(fft(Eout3baa))))/length(Eout3baa)));
    plot(f,db(abs((fftshift(fft(Eout3bab))))/length(Eout3bab)));
    axis([1e11 3e11 -370 10]);
end
  %% Forth Arm
 %split
Eout31bb = ExpCoef.*(cos(kl).*Eout2bb + MutCoef*1j*sin(kl).*0);
Eout32bb = ExpCoef.*(cos(kl).*0 + MutCoef*1j*sin(kl).*Eout2bb);
% Eout31bb = (1/2).*Eout2bb;
% Eout32bb = (1/2).*Eout2bb;
 %Dellay and change phase
aux3bb = Eout32bb(t<=Del2bb);
Eout32dbb = [Eout32bb(end-(length(aux3bb)-1):end) Eout32bb(1:end-(length(aux3bb)))];
Eout32dbb = Eout32dbb.*exp(-1j*Pha2bb);
if Ploting
    figure;
    plot(t,abs(Eout32bb),t,abs(Eout32dbb))
end
% Eout3bb = Eout31bb + Eout32dbb;
% Eout3bba = (1/2).*Eout3bb;
% Eout3bbb = (1/2).*Eout3bb;
Eout3bba = ExpCoef.*(cos(kl).*Eout31bb + MutCoef*1j*sin(kl).*Eout32dbb);
Eout3bbb = ExpCoef.*(cos(kl).*Eout32dbb + MutCoef*1j*sin(kl).*Eout31bb);

if Ploting
    figure;
    hold all;
    plot(f,db(abs((fftshift(fft(EoutT))))/length(EoutT)));
    plot(f,db(abs((fftshift(fft(Eout3bba))))/length(Eout3bba)));
    plot(f,db(abs((fftshift(fft(Eout3bbb))))/length(Eout3bbb)));
    axis([1e11 3e11 -370 10]);
end
% 
% aux = Eout3bbb(t<=GrupDel);
% Eout3bbb = [Eout3bbb(end-(length(aux3bb)-1):end) Eout3bbb(1:end-(length(aux3bb)))];
clear Eout11 Eout12 Eout12d Eout1 Eout21a Eout22a Eout22da Eout2a ...
Eout21b Eout22b Eout22db Eout2b Eout31aa Eout32aa Eout32daa Eout3aa ...
Eout31ab Eout32ab Eout32dab Eout3ab Eout31ba Eout32ba Eout32dba Eout3ba...
Eout31bb Eout32bb Eout32dbb Eout3bb;
%%
if Ploting
    figure;
    hold all;
    plot(f,db(abs((fftshift(fft(EoutT))))/length(EoutT)));
    plot(f,db(abs((fftshift(fft(Eout3aaa))))/length(Eout3aaa)));
    plot(f,db(abs((fftshift(fft(Eout3aab))))/length(Eout3aab)));
    plot(f,db(abs((fftshift(fft(Eout3aba))))/length(Eout3aba)));
    plot(f,db(abs((fftshift(fft(Eout3abb))))/length(Eout3abb)));
    plot(f,db(abs((fftshift(fft(Eout3baa))))/length(Eout3baa)));
    plot(f,db(abs((fftshift(fft(Eout3bab))))/length(Eout3bab)));
    plot(f,db(abs((fftshift(fft(Eout3bba))))/length(Eout3bba)));
    plot(f,db(abs((fftshift(fft(Eout3bbb))))/length(Eout3bbb)));
    axis([1e11 3e11 -370 10]);
end
Eout1 = ifft(fft(Eout3aaa).*SelecFilt3aaa);
Eout1Demd = abs(Eout1).^2;
% Eout1Demd = Eout1Demd-2;
% Eout1Demd = Eout1Demd-min(abs(Eout1Demd));
Eout1Demd = Eout1Demd./max(abs(Eout1Demd));
% Olho( Eout1Demd,T,NPPB );
if PlotingThis
    figure;
    hold all;
    plot(fb,db(abs((fftshift(fft(EoutT))))/length(EoutT)));
    plot(fb,db(abs((fftshift(fft(Eout3aaa))))/length(Eout3aaa)));
    plot(fb,db(abs((fftshift(fft(Eout1))))/length(Eout1)));
    plot(fb,db(abs(fftshift(SelecFilt3aaa))));
    axis([-4e10 6e10 -500 40]);
    a=1;
end
if Ploting
    figure;
    hold all;
    plot(tb,TxDataNrz);
    plot(tb,abs(Eout1Demd));
    a=1;
end

Eout4 = ifft(fft(Eout3bba).*SelecFilt3bba);
Eout4Demd = abs(Eout4).^2;
% Eout4Demd = Eout4Demd-min(abs(Eout4Demd));
Eout4Demd = Eout4Demd./max(abs(Eout4Demd));
% Olho( Eout4Demd,T,NPPB );
if PlotingThis
    figure;
    hold all;
    plot(fb,db(abs((fftshift(fft(EoutT))))/length(EoutT)));
    plot(fb,db(abs((fftshift(fft(Eout3bba))))/length(Eout3bba)));
    plot(fb,db(abs((fftshift(fft(Eout4))))/length(Eout4)));
    plot(fb,db(abs(fftshift(SelecFilt3bba))));
    axis([1e10 9e10 -500 40]);
    a=1;
end
if Ploting
    figure;
    hold all;
    plot(tb,TxDataNrz);
    plot(tb,abs(Eout4Demd));
    a=1;
end
a=1;