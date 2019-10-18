%c
%c                                                       ..'Ž`'..'Ž`..'Ž`..                                                   
%c       File: GenerateOfcInputData
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
%c	 L		     :Comprimento do dispositivo                          [cm]
%c   U0          :Tensao de polarizacao                                [V]
%c   eletrodes	 :Numero de eletrodos (1,2). Definido pelo 
%c                numero de variaveis de entrada na funcao
%c                Mach_Zehnder_Modulator                                [-]
%c   U_pi1       :Tensao de chaveamento para 1 eletrodo                 [V]
%c   U_pi2       :Tensao de chaveamento para 2 eletrodos                [V]
%c   eta1        :Sensibilidade no caminho 1                        [1/V.m]
%c   eta2        :Sensibilidade no caminho 2                        [1/V.m]
%c   nopt	     :Indice de refracao optico                             [-]
%c   nel	     :Indice de refracao eletrico                           [-]
%c   alfa_ins	 :Perda por insercao                                   [dB]
%c   phi_0 		 :Constante de fase entre os dois caminhos              [-]
%c   alfa0		 :Perda condutiva                           [dB/cm.GHz^0.5]
%c   C           :Parametro de chirp                                    [-]
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
Ploting     = ON;%Use On or Off to say if the figures should be plot or not
InfinitLoop = OFF;%Use on or off to say if the loop have output or not
Amplify     = ON;%use on or off to amplify or not the signal
Filtering   = ON;%Use on or off to filter or not the signal
AddCWToLoop = ON;%use on or off to add or not CW to the loop signal
AddNoise    = OFF;%use on or off to add or not Noise to the loop signal
AxisOn      = OFF;%Use on or off to set the axis distances
AxisFig1    = [-2.5e10 2.5e11 -6 -2];%
AxisFig2    = [-2.5e10 2.5e11 0 2];%
AxisFig3    = [-2.5e10 2.5e11 -0.1 1.2];%
CuplerGain  =  0.5; %The split ratio porpotion of the Coupler
% VariablesToSave = [];

%%        Parameters for the Main Program
% Basic specification for the simulation
fc   = 12.5e9;   % RF frequency [Hz]
% fc=f_RF;
n_sc   = 30;    % Number of required subcarriers
spare  = 10;    % Frequency spare [%]
N      = 2^19;  % Number of points in t and (of coarse) in f

%% Finding the apropriate frequency and sampling time
f_max = (1 + spare/100)*n_sc*(4*fc);                                       %Maximum frequency that this code 
                                                                           %will represent, according to Nyquist
df    = 2*f_max/(N - 1);                                                   %Interval between samples frequencies
dt2   = ((N - 1)/N)/(2*f_max);                                             %Time interval between samples
t_max = (N - 1)*dt2;                                                       %End time for the time vector
%
t2    = 0:dt2:(N - 1)*dt2;                                                 %Time vector
f2    = time2freq(t2);                                                     %Frequency vector
%

fc           = 12.5e9;                                                     %center frequency of Electical signal ^10 kilo ^20 mega ^30 giga
Rb           = fc/(1e0);                                                         %bit rate
NumPor       = 2^7;                                                        %Number of carriers used
NumNyq       = 2^2;                                                        %Nyquest number to be at least two time the maximu frequency
fsample      = NumPor*NumNyq*fc;                                           %Sampling frequency
T            = 1/Rb;                                                       %period of the filter
Ta           = 1/fsample;                                                  %period of the filter
NPPB         = T/Ta;                                                       %Number of Point Per Bit
% Tb           = 1/Rb;                                                     %time of one bit NRZ
NtA          = 2^10;
NtF          =(2^19)/(NtA*NumNyq*NumPor);
NumberOf_T   = (2^1)*NtF*NtA;                                                    %number of periods on the simulation 2^18 is the mï¿½ximum value for numberEyes = 2
FinalTime    = NumberOf_T*T;                                               %Final time of simulation
Nb           = floor(FinalTime/T);                                         %number of bits
% TotalSamples = FinalTime*fsample;                                        %Total Number of Samples
TotalSamples = NumberOf_T*(T/Ta);                                          %Total Number of Samples

t            = linspace(0,FinalTime,round(TotalSamples));                         %Time vector for simulation
f            = time2freq(t); 
%%
% dt = t(2)-t(1);
% lt2=length(t2);lt1=length(t);
% t2(end)*1e9-t(end)*1e9
% lt2/lt1
% dt2/dt
%% Loop parametres
%
nsamp              = 1;       %numpre of up sampling
M                  = 5;       %Order of Gausian filter
m                  = 2;       %modulation number
k                  = log2(m); %number of bits per simble
if fc > 12.5e9
    Filter_Center_Freq = 70*fc;   %center frequency of the gausian filter
    FilterBand         = 240*fc;   %pass band of the gausian filter
else
    Filter_Center_Freq = NumPor*fc/2;   %center frequency of the gausian filter
    FilterBand         = 2.5*NumPor*fc;   %pass band of the gausian filter
%     Filter_Center_Freq = 60*fc;   %center frequency of the gausian filter
%     FilterBand         = 240*fc;   %pass band of the gausian filter
end
N_volta            = 400;     %The number of times that the CW will 
                              %recirculate
EbN0               = VarSNR-1;%power of the bit to the noise flor

%ASE power level
Pin        = 43.106;                                    %Input Power (dB)
Pin1       = 17.3195;                                   %Input Power (dB)
Pin2       = -8.6805;                                   %Input Power (dB)
SigRate    = 10^(EbN0/10);                              %SNR
F          = 10^(4.5/10);                               %Noise Factor (dB)
LightSpeed = 3e8;                                       %Ligth Speed (ms)
planck     = 6.626069e-34;                              %Plank constant
CarrLamb   = 1550e-9;                                   %Wave length (m)
% AsePower   = F*planck*(LightSpeed/CarrLamb)*36*fc*SigRate; %Ase Power(dB)
AsePower   = 10^((20-30)/10); %Ase Power(dB)
%%                      Optical filter
[Filtro,FiltroDB] = FiltroGaussiano(f,FilterBand,Filter_Center_Freq,M);
Filtro = fftshift(Filtro);

%%      Parameters of the Continous Wave Lenght
if fc > 12.5e9
    Eo_CW              = 10^((12.5-30)/10);    %Amplitude of the CW signal
else
    Eo_CW              = 10^((12-30)/10);    %Amplitude of the CW signal
end
CW                 = ones(1,length(t));
% CW                 = cos(2*pi*fc.*t);
% CW                 = exp(1j*2*pi*fc*t);
CW                 = Eo_CW*CW;          %Continous Wave (laser)
[~,PdBm] = MeasPower(CW);
a=1;
%%     Variable parameters for evaluation
Varm1_steps = 1;%Number of different amplitudes on first arm of the MZM 
Varm2_steps = 1;%Number of different amplitudes on second arm of the MZM 
Gain_steps  = 1;%Number of different gains used in the optical ring
phase_steps = 1;%Number of different phase delays that will be used
Vbias_steps = 1;%Number of different Bias voltage that will be used
%
%IMPORTANT: The RF amplitude Vpp must respect the device limits!
Varm1         = 1;%linspace(-3.8/2,3.8/2,Varm1_steps);                     %Variation of Amplitude for Arm 1
Varm2         = 1;%linspace(-3.8/2,3.8/2,Varm2_steps);                     %Variation of Amplitude for Arm 2
phase         = pi/2;%linspace(0,pi,phase_steps);                          %Variation of phase betwee arms. [pi/12 pi/4 pi/2 2*pi/3];%Phase Vector;
Vbias         = 3.80/2;%linspace(-7,7,Vbias_steps);                         %Variation of bias voltage on the simulation%[0.8 1.8 2.8 3.8];%Vbias Vector.
if InfinitLoop
    if PdBm>=15
        Gain_vet      = 30+10*log10(1.0223);
        SatGain       = 30+10*log10(7e5);
    elseif PdBm>=0
        Gain_vet      = 30+10*log10(2.049);%2.049;%2.177;%linspace(1,2,Gain_steps);                           %Variation of Gain inside the optical ring
        if EbN0==90
            SatGain   = 30+10*log10(4.6e5);
        elseif EbN0==70
            SatGain   = 30+10*log10(2.4e4);
        elseif EbN0==80
            SatGain   = 30+10*log10(7e4);
        else
            SatGain   = 30+10*log10(2.2e4);
        end
    elseif PdBm>=-7
        Gain_vet      = 30+10*log10(1.0224);%2.049;%2.177;%linspace(1,2,Gain_steps); 
        SatGain   = 30+10*log10(5e2);
    elseif PdBm>=-12
        Gain_vet      = 30+10*log10(1.0224);%2.049;%2.177;%linspace(1,2,Gain_steps); 
        SatGain   = 30+10*log10(1e2);
    elseif PdBm>=-60
        Gain_vet      = 30+10*log10(1.0224);%2.049;%2.177;%linspace(1,2,Gain_steps);                           %Variation of Gain inside the optical ring
        SatGain   = 30+10*log10(1e1);
    else
        Gain_vet      = 30+10*log10(1.0224);%2.049;%2.177;%linspace(1,2,Gain_steps);                           %Variation of Gain inside the optical ring
        SatGain   = 30+10*log10(1e0);
    end
else
    if AddNoise
        Gain_vet = 30+10*log10(2.0445);
        SatGain  = 30+10*log10(1e7);
    else
        if PdBm>=15
            Gain_vet = 30+10*log10(2.0445);
            SatGain  = 30+10*log10(8e4);
        else
            if fc > 12.5e9
                Gain_vet = 30+10*log10(2.0449);
                SatGain  = 30+10*log10(1e2);
            else
                Gain_vet = 30+10*log10(2.0449);
                SatGain  = 30+10*log10(1e2);
            end
        end
    end
end
% Vbias         = 3.78;%linspace(-7,7,Vbias_steps);                         %Variation of bias voltage on the simulation%[0.8 1.8 2.8 3.8];%Vbias Vector.
% Gain_vet      = 3.63;%linspace(1,2,Gain_steps);                           %Variation of Gain inside the optical ring
%
N_combination = phase_steps*Vbias_steps*Varm1_steps*Varm2_steps*Gain_steps;%Variable to calculate the percentage for completition of the simulation
%
Dif_Amp = 200;       %Variation of peaks power amplitude to be found
ML      = 20;        %Multiplier for the db measurement. 20 if the signal
                     %is a voltage 10 if the signal is the power in Watts
if (ML~=20) && (ML~=10)
    error(['ML - Multiplier must be 10 or 20. Other values will result'...
                                                 ' in wrong evaluations']);
end
Rad_f   = 2*pi*fc; %Radian frequency: Argument of the sinusoidal signal

user = memory;                                                             %Variables to check if there will be enough memory to store all the possible combination
if N_combination>user.MaxPossibleArrayBytes*0.8          
    error(['The current simulation will overflow the maximum merrory' ... 
    ' available of the sistem. Please reduce the maximum number of' ...
                     ' possible combinations of this current simulation']);
end
%%  MZM Parameters
MZ_Input_Data;                   %Loading the basic data to be saved
L        =  10.00;
% U0       =  2.390;
U_pi1    =   3.80;
U_pi2    =   3.80;
eta1     =  89.00;%1;%
eta2     = -12.70;%-1;%
nopt     =   2.17;
nel      =   2.60;
alfa_ins =   5.10;
phi_0    =   0.00;
alfa0    =   1.07;
% C        = (eta1+eta2)/(eta1-eta2);
% alfa0    = 10^(alfa0/20);
%% Files parameters
p = pwd;                    %Geting the current path of the main program
Local = [p '\input_files\'];%variable that shows where the input files for
                                                     %the mzm will be saved
SFName = ['OCS_fc_125_Car_' num2str(NumPor) '_Long_' date()];
% SFName = ['OCS_' num2str(EbN0)];
% SFName = ['OCS_' num2str(200+EbN0)];
% SFName = ['OCS_Pin_1'];
LocSav = [p '\save_files\'];