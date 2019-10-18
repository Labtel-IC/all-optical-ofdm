%c
%c                                                       ..'Ž`'..'Ž`..'Ž`..                                                   
%c       File: InputDataCharcAnaliticMmzDp 
%c
%c       This code is used to set up the input parameters that will be used
%c in the whole simualation.
%c
%c
%c
%c
%c
%c                                           by P.Marciano LG
%c                                           18/09/2017
%c                                           pablorafael.mcx@gmail.com
%c 
%c
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                    G L O B A L  V A R I A B L E S                      %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                   L I S T  O F  V A R I A B L E S                      %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c
%c
%c   
%c                     c
%c   f_RF        :This is the the Frequency of the signal times its order. 
%c                This is your alltogether frequency.                  [Hz]
%c   n_sc        :Number of subacaries that will be generated           [-]
%c   spare       :Engineering factor for better representation of the
%c                                                    maximum frequency [-]
%c   N           :Total number of points for the simulation             [-]
%c   f_max       :The last (highest) frequency to be represented      [GHz] 
%c   df          :Smallest frequency interval for the simulatio      [f/df]
%c   dt          :Smallest time interval for the simulation          [t/dt]
%c   t_max       :The last moment of the time vector                    [s] 
%c   t           :This is the time vector for your simulation, it 
%c                will run just as much periods of your signal you
%c                need and/or want to simulate.                        [s]
%c   f           :This is the frequency vector for your simulation, it 
%c                will run just as much periods of your signal you need 
%c                and/or want to simulate.                             [Hz]
%c   Po          :The power of the signal                               [w]
%c   Eo_CW       :This is the Amplitude of the signal                   [V]
%c   nsamp       :numpre of up sampling                                 [-]
%c   M           :Order of Gausian filter                               [-]
%c   m           :Modulation number                                     [-]
%c   k           :Number of bits per simble                             [-]
%c   Filter_Center_Freq :center frequency of the gausian filter       [GHz]
%c   FilterBand  :Pass band of the gausian filter                     [GHz]
%c   N_volta     :The number of times that the CW will recirculate      [-]
%c   EbN0        :Power of the bit to the noise flor        (Eb[dB]/No[dB])
%c   Varm1_steps :Number of different amplitudes on first arm of the 
%c                MZM                                                   [-] 
%c   Varm2_steps :Number of different amplitudes on second arm of the
%c                MZM                                                   [-] 
%c   Gain_steps  :Number of different gains used in the optical ring    [-]
%c   phase_steps :Variable responsible to store the number of possible
%c                values that the phase may assume through the whole 
%c                simulation.                                           [-]
%c   Vbias_steps :Variable responsible to store the number of possible
%c                values that the polatization voltage may assume through  
%c                the whole simulation.                                 [-]
%c   Varm1       :Variation of Amplitude for Arm 1                      [V]
%c   Varm2       :Variation of Amplitude for Arm 2                      [V]
%c   phase       :Variation of phase betwee arms.                     [phi]
%c   Vbias       :Stores the vector with all possible values for the
%c                polarization voltage.                                 [V]
%c   Gain_vet    :Variation of Gain inside the optical ring             [-]
%c   N_combination:Variable to calculate the percentage for completition of
%c                 the simulation                                       [-]
%c   Dif_Amp     :Variation of peaks power amplitude to be found   [Linear]
%c   ML          :Multiplier for the db measurement. 20 if the signal is a
%c                voltage 10 if the signal is the power in Watts       [dB]
%c   Rad_f       :Radian frequency: Argument of the sinusoidal
%c                signal                                        (2pi*f[Hz])
%c   Filtro      :Gausian filter to be used to limit the band of the OFC in
%c                the optical ring                                      [-]
%c   p           :Auxiliar variable to store the current directory      [-]
%c   Local       :It stores where the files for the MZM will be stored  [-]
%c
%c   
%c
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                             S T R U C T S                              %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c   Here it will be added the functions that this code will call
%c   
%c
%c Errase the folowing lines when this code is done
%c   MAIN_COMB_INPUT_DATA.
%c                    .Does not have input parameters
%c
%c
%c   The foloing functions are references they do not exist.
%c   
%c   Amplifier.
%c                   .id
%c                   .id_Node
%c                   .Gain
%c                   .Noise_Figure
%c                   .Product_Code
%c                   .type          : Must be 'electrical' or 'optical'
%c   Band.
%c   Optical_Channel.
%c                   .id
%c                   .id_Route
%c                   .lambda_used
%c                   .lambda.viable
%c   DCM.
%c   Demand.
%c                   .id
%c                   .id_Node_Source
%c                   .id_Node_Destination
%c                   .id_Optical_Channel
%c                   .type           : IMDD, Coherent, EOFDM, OOFDM 
%c                                     or FlexGrid.
%c                                     IMDD -  Intensity-Modulation 
%c                                             Direct Detection
%c                                     EOFDM - Electrically generated 
%c                                             optical OFDM
%c                                     OOFDM - Optically generated 
%c                                             optical OFDM
%c                   .value          : [Gbit/s]
%c                   
%c
%c   Fibre_Segment. 
%c   Interface.
%c   Link.
%c   Network.   
%c   Node.
%c   OADM_Terminal.
%c   Route.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Start of the Program                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%       General Parameter for the simulation
ON          =  1; %DO NOT CHANGE THIS PARAMETER
OFF         =  0; %DO NOT CHANGE THIS PARAMETER
InfinitLoop = ON;%Use on or off to say if the loop have output or not
Amplify     = ON;%use on or off to amplify or not the signal
Filtering   = ON;%Use on or off to filter or not the signal
AddCWToLoop = ON;%use on or off to add or not CW to the loop signal
AddNoise    = ON;%use on or off to add or not Noise to the loop signal
Ploting     = OFF;%Used say if it should plot or not the current results
PlotingThis = ON;%Used say if it should plot or not the current results
AxisOn      = OFF;%Use on or off to set the axis distances
AxisFig1    = [-2.5e10 2.5e12 -2 1];%
AxisFig2    = [-2.5e10 2.5e12 3 5];%
AxisFig3    = [-2.5e10 2.5e12 -0.1 1.2];%
CuplerGain  =  0.5; %The split ratio porpotion of the Coupler

%%        Parameters for the Main Program
% Basic specification for the simulation
fc   = 12.5e9;   % RF frequency [Hz]
% fc=f_RF;
n_sc   = 100;    % Number of required subcarriers
spare  = 28;    % Frequency spare [%]
N      = 2^19;  % Number of points in t and (of coarse) in f

%% Finding the apropriate frequency and sampling time
f_max = (1 + spare/100)*n_sc*(4*fc);                                      %Maximum frequency that this code 
                                                                           %will represent, according to Nyquist
df    = 2*f_max/(N - 1);                                                  %Interval between samples frequencies
dt2   = ((N - 1)/N)/(2*f_max);                                           %Time interval between samples
t_max = (N - 1)*dt2;                                                      %End time for the time vector
%
t1    = 0:dt2:(N - 1)*dt2;                                               %Time vector
f1    = time2freq(t1);                                                   %Frequency vector
%

% fc          = 12.5e9;%2^30;%center frequency of Electical signal ^10 kilo ^20 mega ^30 giga
f_max        = 32*fc;                                                      %max frequency that matlab can reproduce
fsample      = 32*f_max;                                                   %Sampling frequency
T            = 1/fc;                                                       %period of the filter
% Rb           = fc/2;                                                     %bit rate
Rb           = fc;                                                         %bit rate
Tb           = 1/Rb;                                                       %time of one bit NRZ
NumberOf_T   = 2^9;                                                        %number of periods on the simulation 2^18 is the mï¿½ximum value for numberEyes = 2
FinalTime    = NumberOf_T*T;                                               %Final time of simulation
Nb           = FinalTime/Tb;                                               %number of bits
TotalSamples = FinalTime*fsample;                                          %Total Number of Samples
NPPB         = TotalSamples/Nb;                                            %Number of Point Per Bit
t            = linspace(0,FinalTime,round(TotalSamples));                         %Time vector for simulation
f            = time2freq(t); 
%%
dt = t(2)-t(1);
lt1=length(t1);lt1=length(t);
t1(end)*1e9-t(end)*1e9
lt1/lt1
dt2/dt
%% Loop parametres
%
nsamp              = 1;       %numpre of up sampling
M                  = 5;       %Order of Gausian filter
m                  = 2;       %modulation number
k                  = log2(m); %number of bits per simble
Filter_Center_Freq = 0;       %center frequency of the gausian filter
FilterBand         = 4*f_max;%40e9;%pass band of the gausian filter
N_volta            = 100;      %The number of times that the CW will 
                              %recirculate
EbN0               = 40;      %power of the bit to the noise flor

%%      Parameters of the Continous Wave Lenght
Eo_CW              = 1; %Amplitude of the CW signal
CW = ones(1,length(t));
CW=Eo_CW*CW;            %Continous Wave (laser)

%%     Variable parameters for evaluation
Varm1_steps  = 10+1;%Number of different amplitudes on first arm of the MZM 
Varm2_steps  = 10+1;%Number of different amplitudes on second arm of the MZM 
Gain_steps   = 20+2;%Number of different gains used in the optical ring
phase_steps  = 20+2;%Number of different phase delays that will be used
phi_0_steps  = 20+1;%Number of different phase between arms of the DP-MZM
V1bias_steps = 20+1;%Number of different Bias voltage at the Up arm
V2bias_steps = 20+1;%Number of different Bias voltage at the Lower arm
%
%IMPORTANT: The RF amplitude Vpp must respect the device limits!
Varm1         = linspace(-5,5,Varm1_steps);                        %Variation of Amplitude for Arm 1
Varm2         = linspace(-5,5,Varm2_steps);                        %Variation of Amplitude for Arm 2
phase         = linspace(0,pi,phase_steps);                                %Variation of phase betwee Electrodes.;
phi_0_vet     = linspace(-7,7,phi_0_steps);                                %Variation of phase betwee arms.
V1bias        = linspace(-7,7,V1bias_steps);                               %Variation of bias voltage at the Up arm.
V2bias        = linspace(-7,7,V2bias_steps);                               %Variation of bias voltage at the Lower arm.
Gain_vet      = linspace(1,2,Gain_steps);                                  %Variation of Gain inside the optical ring
N_combination = phase_steps*V1bias_steps*V2bias_steps*Varm1_steps*...
                                        phi_0_steps*Gain_steps*Varm2_steps;%Variable to calculate the percentage for completition of the simulation
%
user = memory;                                                             %Variables to check if there will be enough memory to store all the possible combination
if N_combination>user.MaxPossibleArrayBytes*0.8          
    error(['The current simulation will overflow the maximum merrory' ... 
    ' available of the sistem. Please reduce the maximum number of' ...
                     ' possible combinations of this current simulation']);
end
%
Dif_Amp = 200;       %Variation of peaks power amplitude to be found
ML      = 20;        %Multiplier for the db measurement. 20 if the signal
                     %is a voltage 10 if the signal is the power in Watts
if (ML~=20) && (ML~=10)
    error(['ML - Multiplier must be 10 or 20. Other values will result'...
                                                 ' in wrong evaluations']);
end
Rad_f   = 2*pi*fc; %Radian frequency: Argument of the sinusoidal signal
%%                      Optical filter
[Filtro,FiltroDB] = FiltroGaussiano(f,FilterBand,Filter_Center_Freq,M);
Filtro = fftshift(Filtro);
%%  MZM Parameters
MZ_Input_Data_DP;                   %Loading the basic data to be saved
%% Files parameters
p = pwd;                                                                   %Geting the current path of the main program
Local = [p '\input_files\'];                                               %variable that shows where the input files for
                                                                           %the mzm will be saved
SFName = ['DpmzmRing' date];
LocSav = [p '\save_files\'];