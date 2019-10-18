clear;close all;clc
%%        Parameters for the Main Program
% Basic specification for the simulation
fc   = 12.5e9;   % RF frequency [Hz]
% fc=f_RF;
n_sc   = 30;    % Number of required subcarriers
spare  = 10;    % Frequency spare [%]
N      = 2^19;  % Number of points in t and (of coarse) in f

%% Finding the apropriate frequency and sampling time
f_max2 = (1 + spare/100)*n_sc*(4*fc);                                      %Maximum frequency that this code 
                                                                           %will represent, according to Nyquist
df    = 2*f_max2/(N - 1);                                                  %Interval between samples frequencies
dt2   = ((N - 1)/N)/(2*f_max2);                                           %Time interval between samples
t_max = (N - 1)*dt2;                                                      %End time for the time vector
%
t2    = 0:dt2:(N - 1)*dt2;                                               %Time vector
f2    = time2freq(t2);                                                   %Frequency vector
%

% fc          = 12.5e9;%2^30;%center frequency of Electical signal ^10 kilo ^20 mega ^30 giga
f_max        = (2^2)*fc;                                                      %max frequency that matlab can reproduce
fsample      = (2^2)*f_max;                                                   %Sampling frequency
T            = 1/fc;                                                       %period of the filter
Ta           = 1/fsample;                                                       %period of the filter
Rb           = fc;                                                     %bit rate
% Rb           = fc;                                                         %bit rate
Tb           = 1/Rb;                                                       %time of one bit NRZ
NPPB         = 2^5;                                            %Number of Point Per Bit
NumberOf_T   = 2^10;                                                        %number of periods on the simulation 2^18 is the m�ximum value for numberEyes = 2
FinalTime    = NumberOf_T*T;                                               %Final time of simulation
Nb           = floor(FinalTime/Tb);                                               %number of bits
% TotalSamples = FinalTime*fsample;                                          %Total Number of Samples
TotalSamples = NPPB*NumberOf_T*(T/Ta);                                                   %Total Number of Samples

t            = linspace(0,FinalTime,round(TotalSamples));                         %Time vector for simulation
f            = time2freq(t); 

sinal1=sin(2*pi*fc*t);
sinal2=sin(2*pi*fc*t2);
figure
hold on
plot(t,sinal1)
plot(t2,sinal2)
a=1;