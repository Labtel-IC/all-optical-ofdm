
%%       General Parameter for the simulation
ON          =  1; %DO NOT CHANGE THIS PARAMETER
OFF         =  0; %DO NOT CHANGE THIS PARAMETER
Ploting     = ON;%Use On or Off to say if the figures should be plot or not
InfinitLoop = OFF;%Use on or off to say if the loop have output or not
Amplify     = ON;%use on or off to amplify or not the signal
Filtering   = ON;%Use on or off to filter or not the signal
AddCWToLoop = ON;%use on or off to add or not CW to the loop signal
AddNoise    = ON;%use on or off to add or not Noise to the loop signal
AxisOn      = OFF;%Use on or off to set the axis distances
AxisFig1    = [-2.5e10 2.5e12 -2 1];%
AxisFig2    = [-2.5e10 2.5e12 3 5];%
AxisFig3    = [-2.5e10 2.5e12 -0.1 1.2];%
CuplerGain  =  0.5; %The split ratio porpotion of the Coupler

%%        Parameters for the Main Program
% Basic specification for the simulation
fc   = 12.5e9;   % RF frequency [Hz]
% fc=f_RF;
n_sc   = 30;    % Number of required subcarriers
spare  = 10;    % Frequency spare [%]
N      = 2^19;  % Number of points in t and (of coarse) in f

%% Finding the apropriate frequency and sampling time
f_max = (1 + spare/100)*n_sc*(4*fc);                                      %Maximum frequency that this code 
                                                                           %will represent, according to Nyquist
df    = 2*f_max/(N - 1);                                                  %Interval between samples frequencies
dt2   = ((N - 1)/N)/(2*f_max);                                           %Time interval between samples
t_max = (N - 1)*dt2;                                                      %End time for the time vector
%
t2    = 0:dt2:(N - 1)*dt2;                                               %Time vector
f2    = time2freq(t2);                                                   %Frequency vector
%

% fc          = 12.5e9;%2^30;%center frequency of Electical signal ^10 kilo ^20 mega ^30 giga
f_max        = (2^2)*fc;                                                      %max frequency that matlab can reproduce
fsample      = (2^2)*f_max;                                                   %Sampling frequency
T            = 1/fc;                                                       %period of the filter
Tb           = T/2;
Ta           = 1/fsample;                                                       %period of the filter
% Rb           = fc/2;                                                     %bit rate
Rb           = fc;                                                         %bit rate
% Tb           = 1/Rb;                                                       %time of one bit NRZ
NPPB         = 2^5;                                            %Number of Point Per Bit
NumberOf_T   = 2^2;                                                        %number of periods on the simulation 2^18 is the m�ximum value for numberEyes = 2
FinalTime    = NumberOf_T*T;                                               %Final time of simulation
Nb           = FinalTime/Tb;                                               %number of bits
% TotalSamples = FinalTime*fsample;                                          %Total Number of Samples
TotalSamples = NPPB*NumberOf_T*(T/Ta);                                                   %Total Number of Samples

t            = linspace(0,FinalTime,round(TotalSamples));                         %Time vector for simulation
f            = time2freq(t); 
%%
dt = t(2)-t(1);
lt2=length(t2);lt1=length(t);
t2(end)*1e9-t(end)*1e9
lt2/lt1
dt2/dt


%%     Variable parameters for evaluation
Varm1_steps  = 1;%Number of different amplitudes on first arm of the MZM 
Varm2_steps  = 1;%Number of different amplitudes on second arm of the MZM 
Gain_steps   = 1;%Number of different gains used in the optical ring
phase_steps  = 1;%Number of different phase delays that will be used
phi_0_steps  = 1;%Number of different phase between arms of the DP-MZM
V1bias_steps = 1;%Number of different Bias voltage at the Up arm
V2bias_steps = 1;%Number of different Bias voltage at the Lower arm
%
%IMPORTANT: The RF amplitude Vpp must respect the device limits!
Varm1         = 1;%linspace(-5,5,Varm1_steps);                        %Variation of Amplitude for Arm 1
Varm2         = 1;%linspace(-5,5,Varm2_steps);                        %Variation of Amplitude for Arm 2
phase         = pi/2;% linspace(0,pi,phase_steps);                                %Variation of phase betwee Electrodes.;
phi_0_vet     = -0.5;%-2.05;%linspace(-7,7,phi_0_steps);                                %Variation of phase betwee arms.
V1bias        = 1.752;%-0.35;%linspace(-7,7,V1bias_steps);                               %Variation of bias voltage at the Up arm.
V2bias        = 1.052;%-1.2;%linspace(-7,7,V2bias_steps);                               %Variation of bias voltage at the Lower arm.
Gain_vet      = 1;%linspace(1,2,Gain_steps);                                  %Variation of Gain inside the optical ring


%%  MZM Parameters
MZ_Input_Data_DP;  

p = pwd;                    %Geting the current path of the main program
Local_Dp = [p '\input_files\'];%variable that shows where the input files for
Local = Local_Dp;

Make_MZ_Input_Files_DP;
% A = 2.*[-1 -1 -1 -1];
% B = 2.*[-1 -1 -1 -1];
% A = 2.*[1 1 1 1];
% B = 2.*[-1 -1 -1 -1];
% A = 2.*[-1 -1 -1 -1];
% B = 2.*[1 1 1 1];
% A = 2.*[1 1 1 1];
% B = 2.*[1 1 1 1];
% A = 2.*[-1 -1 1 1];
% B = 2.*[-1 -1 1 1];
A = (randi(2,1,Nb)-1);
B = (randi(2,1,Nb)-1);
% A = [0 1 1 0 1 0 1 1];
% B = [1 1 1 1 1 0 1 0];
A(A==0)=-1;
B(B==0)=-1;
A=1.*A;
B=1.*B;
C = cos(2*pi*fc*t);
SingA = rectpulse(A,TotalSamples/Nb);
SingB = rectpulse(B,TotalSamples/Nb);
% aux1 = SingB(t<=T/4);
% SingB = [SingB(end-(length(aux1)-1):end) SingB(1:end-(length(aux1)))];
aux1 = C(t<=T/4);
C2 = [C(end-(length(aux1)-1):end) C(1:end-(length(aux1)))];
% plot(t,C,t,C2);
% plot(t,C,t,SingA,t,SingB);
                                                     %the mzm will be saved
                                                     
                                                     
Ein = cos(2*pi*fc*t) + 1j*sin(2*pi*fc*t);
% U.U1t = SingA;
% U.U2t = SingB;
U.U1t = zeros(1,length(t));
U.U2t = zeros(1,length(t));
MZ_Input_File = 1;
Mach_Zehnder_Modulator_DP(t,Ein,U,MZ_Input_File,Local_Dp);