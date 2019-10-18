close all;clear;clc;


fc2           = 12.5e9;                                                     %center frequency of Electical signal ^10 kilo ^20 mega ^30 giga
Rb           = fc2;                                                         %bit rate
NPPB         = 2^7;                                                        %Number of Point Per Bit
fsample      = NPPB*Rb;                                                    %Sampling frequency
SymPer            = 1/Rb;                                                       %period of the filter
Ta           = 1/fsample;                                                  %period of the filter
% Tb           = 1/Rb;                                                     %time of one bit NRZ
NumberOf_T   = 2^13;                                                       %number of periods on the simulation 2^18 is the m�ximum value for numberEyes = 2
FinalTime    = NumberOf_T*SymPer;                                               %Final time of simulation
Nb           = floor(FinalTime/SymPer);                                         %number of bits
% TotalSamples = FinalTime*fsample;                                        %Total Number of Samples
TotalSamples = NumberOf_T*(SymPer/Ta);                                          %Total Number of Samples

t2            = linspace(0,FinalTime,round(TotalSamples));                         %Time vector for simulation
f2            = time2freq(t2); 



A = exp(1j*2*pi*fc2*t2);
B = exp(1j*4*pi*fc2*t2);
C = exp(1j*6*pi*fc2*t2);
D = exp(1j*8*pi*fc2*t2);
E = exp(1j*10*pi*fc2*t2);
F = exp(1j*12*pi*fc2*t2);
G = exp(1j*14*pi*fc2*t2);
H = exp(1j*16*pi*fc2*t2);
I = exp(1j*18*pi*fc2*t2);
J = exp(1j*20*pi*fc2*t2);
K = exp(1j*22*pi*fc2*t2);
L = exp(1j*24*pi*fc2*t2);
M = exp(1j*26*pi*fc2*t2);
N = exp(1j*28*pi*fc2*t2);
O = exp(1j*30*pi*fc2*t2);
P = exp(1j*32*pi*fc2*t2);
Q = exp(1j*34*pi*fc2*t2);
R = exp(1j*36*pi*fc2*t2);
S = exp(1j*38*pi*fc2*t2);
T = exp(1j*40*pi*fc2*t2);
U = exp(1j*42*pi*fc2*t2);
V = exp(1j*44*pi*fc2*t2);
W = exp(1j*46*pi*fc2*t2);
X = exp(1j*48*pi*fc2*t2);
Y = exp(1j*50*pi*fc2*t2);
Z = exp(1j*52*pi*fc2*t2);
AA = exp(1j*54*pi*fc2*t2);
AB = exp(1j*56*pi*fc2*t2);
AC = exp(1j*58*pi*fc2*t2);
AD = exp(1j*60*pi*fc2*t2);
AE = exp(1j*62*pi*fc2*t2);
AF = exp(1j*64*pi*fc2*t2);
COMB = A+B+C+D+E+F+G+H+I+J+K+L+M+N+O+P+Q+R+S+T+U+V+W+X+Y+Z+AA+AB+AC+AD+AE+AF;
COMB = COMB + wgn(1,length(t2),-1*40);
figure;hold on
plot(f2,db(abs(fftshift(fft(COMB)./length(COMB)))))
EoutT = COMB;
EoutRec = COMB;
ON = 1;
OFF = 0;
[Eout,~]=expfunc(t2,SymPer,5,EoutT);%,AngVec,ActStag,Count,E0)
figure;
hold all;
plot(f2,db(abs((fftshift(fft(EoutT))))/length(EoutT)));
for kk=1:size(Eout,1)
    plot(f2,db(abs((fftshift(fft(Eout(kk,:)))))/length(Eout(kk,:))));
end
axis([-0.1e11 (size(Eout,1))*1.1*fc2 -250 10]);
set(gcf,'units','normalized','outerposition',[0 0 1 1])
% DatumReception
a=1;