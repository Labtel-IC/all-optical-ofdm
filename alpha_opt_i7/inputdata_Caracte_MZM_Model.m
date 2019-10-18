%c
%c                                                       ..'Ž`'..'Ž`..'Ž`..                                                   
%c       File: inputdata_Caracte_MZM_Model 
%c
%c       This code is used to set up the input parameters that will be used
%c in the whole simualation.
%c
%c
%c
%c
%c
%c                                           by P.Marciano LG
%c                                           29/09/2017
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
%c                c
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
%c   phase_steps :Variable responsible to store the number of possible
%c                values that the phase may assume through the whole 
%c                simulation.                                           [-]
%c   Vbias_steps :Variable responsible to store the number of possible
%c                values that the polatization voltage may assume through  
%c                the whole simulation.                                 [-]
%c   Vbias       :Stores the vector with all possible values for the
%c                                                polarization voltage. [V]
%c   EleSig2     :The RF signal that will be inseted to the second arm of
%c                                                             the MZM [Hz]
%c   EleSig1     :The RF signal that will be inseted to the first arm of
%c                                                             the MZM [Hz]
%c   EleSig      :The struct RF signal that will be inseted to the MZM [Hz]
%c   AV1         :Configurable parameter to adjust the theoritical equation
%c                for the MZM output electrical field                   [-]
%c   p           :Auxiliar variable to store the current directory
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
%c
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Start of the Program                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%        Parameters for the Main Program
% Basic specification for the simulation
f_RF   = 1E9;   % RF frequency [Hz]
n_sc   = 10;    % Number of required subcarriers
spare  = 10;    % Frequency spare [%]
N      = 2^19;  % Number of points in t and (of coarse) in f

%% Finding the apropriate frequency and sampling time
f_max  = (1+spare/100)*n_sc*(4*f_RF); %Maximum frequency that this code 
                                      %will represent, according to Nyquist
df     = 2*f_max/(N - 1);             %Interval between samples frequencies
dt     = ((N - 1)/N)/(2*f_max);       %Time interval between samples
t_max  = (N - 1)*dt;                  %End time for the time vector
%
% t      = 0:dt:(N - 1)*dt;             %Time vector
% f      = time2freq(t);                %Frequency vector

fc2         = 1e9;%2^30;%center frequency of Electical signal ^10 kilo ^20 mega ^30 giga
f_max2      = 32*fc2;%max frequency that matlab can reproduce
fsample     = 16*f_max2;%Sampling frequency
T           = 1/fc2;%period of the filter
Rb          = fc2/2;%bit rate
Tb          = 1/Rb;%time of one bit NRZ
NumberOf_T  = 4*2^7;%number of periods on the simulation 2^18 is the mï¿½ximum value for numberEyes = 2
FinalTime   = NumberOf_T*T;%Final time of simulation
Nb          = FinalTime/Tb;%number of bits
TotalSamples= FinalTime*fsample;%Total Number of Samples
NPPB        = TotalSamples/Nb;%Number of Point Per Bit
t          = linspace(0,FinalTime,TotalSamples);%Time vector for simulation
f          = time2freq(t);  
%%      Parameters of the Continous Wave Lenght
Eo_CW              = 1; %Amplitude of the CW signal
CW = ones(1,length(t));
CW=Eo_CW*CW;            %Continous Wave (laser)
Po = mean(abs(CW).^2);  %Average Power of the CW
%%     Variable parameters for evaluation

phase_steps = 1;                    %Length of the vector that will store 
                                    %the fase variable of the MZM       [-]
Vbias_steps = 200+1;                %Length of the vector that will store 
                                    %the polarization voltage of the MZM[-]
Vbias = linspace(-7,7,Vbias_steps); %Vector to store the polarization 
                                    %voltage for the whole simulation [V]

EleSig2 = zeros(1,length(t));       %The signal applied at the first arm of
                                    %the MZM
EleSig1 = zeros(1,length(t));       %The second arm will have a phase shift 
                                    %different of the first arm
                                                                           
EleSig.U1t=EleSig1;                 %The MZM model require that the 
EleSig.U2t=EleSig2;                 %voltages on a variabel struct 

AV1 = 2.78;                            %Correction factor to corect the 
                                    %theoretical value with the 
                                    %experimental
                                    
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
%% Files parameters
p = pwd;                    %Geting the current path of the main program
Local = [p '\input_files\'];%variable that shows where the input files for
                                                     %the mzm will be saved
% CW = C;
% H = zeros(1,length(t));
% Powerdb = zeros(1,Vbias_steps);
% PowerW = zeros(1,Vbias_steps);
% PowerI = zeros(1,Vbias_steps);
% PowerIdb = zeros(1,Vbias_steps);