clear,close all,clc;
%%        Parameters for the Main Program
% Basic specification for the simulation
f_RF   = 1e9;   % RF frequency [Hz]
n_sc   = 30;    % Number of required subcarriers
spare  = 10;    % Frequency spare [%]
N      = 2^19;  % Number of points in t and (of coarse) in f

%% Finding the apropriate frequency and sampling time
f_max  = (1 + spare/100)*n_sc*(4*f_RF);                                    %Maximum frequency that this code 
                                                                           %will represent, according to Nyquist
df     = 2*f_max/(N - 1);                                                  %Interval between samples frequencies
dt     = ((N - 1)/N)/(2*f_max);                                            %Time interval between samples
t_max  = (N - 1)*dt;                                                       %End time for the time vector
%
t      = 0:dt:(N - 1)*dt;                                                  %Time vector
f      = time2freq(t);                                                     %Frequency vector
%%      Parameters of the Continous Wave Lenght
Eo_CW              = 1; %Amplitude of the CW signal
CW = ones(1,length(t));
CW=Eo_CW*CW;            %Continous Wave (laser)
Rad_f   = 2*pi*f_RF; %Radian frequency: Argument of the sinusoidal signal
%%
EleSig2 = sin(Rad_f*t);%sin(Rad_f*t + pi/2);           %The second arm will have a phase shift different of the first arm
EleSig1 = sin(Rad_f*t);                      %The second arm will have a phase shift different of the first arm
EleSig.U1t =EleSig1;                                   %The MZM model require that the voltages
EleSig.U2t =EleSig2;                                   %The MZM model require that the voltages
[Eout,~]=Mach_Zehnder_Modulator_DP(t,CW,EleSig,6666);