%%
ThisPath = pwd;                                                            %Path for the folter where the OFC is saved
UofcsLocal = '\..\save_files\';                                               %Files name for the location of the OFC
if exist('CurrentOCS','var')                                               %This was not been used, this code works when different OFCS need to be tested and it was called by another script
    OfcName = [ThisPath UofcsLocal 'OCS_' num2str(OcsToTest(CurrentOCS))];
else
    %     OfcName = [ThisPath UofcsLocal 'OCS_fc_125_Car_2_08-Jul-2018'];        %OFCS with 2 carriers NPPB 2^3 Nb 2^16 fc 12.5e9 EDFA Power 20 dB
    %     OfcName = [ThisPath UofcsLocal 'OCS_fc_125_Car_4_08-Jul-2018'];        %OFCS with 4 carriers NPPB 2^4 Nb 2^15 fc 12.5e9 EDFA Power 20 dB
    %     OfcName = [ThisPath UofcsLocal 'OCS_fc_125_Car_8_08-Jul-2018'];        %OFCS with 8 carriers NPPB 2^5 Nb 2^14 fc 12.5e9 EDFA Power 20 dB
    %     OfcName = [ThisPath UofcsLocal 'OCS_fc_125_Car_16_08-Jul-2018'];       %OFCS with 16 carriers NPPB 2^6 Nb 2^13 fc 12.5e9 EDFA Power 20 dB
    %     OfcName = [ThisPath UofcsLocal 'OCS_fc_125_Car_32_08-Jul-2018'];       %OFCS with 32 carriers NPPB 2^7 Nb 2^12 fc 12.5e9 EDFA Power 20 dB
    %     OfcName = [ThisPath UofcsLocal 'OCS_fc_125_Car_64_08-Jul-2018'];       %OFCS with 64 carriers NPPB 2^8 Nb 2^11 fc 12.5e9 EDFA Power 20 dB
    OfcName = [ThisPath UofcsLocal 'OCS_fc_125_Car_128_08-Jul-2018'];      %OFCS with 128 carriers NPPB 2^9 Nb 2^10 fc 12.5e9 EDFA Power 20 dB
    %     OfcName = [ThisPath UofcsLocal 'OCS_fc_125_Car_256_09-Jul-2018'];      %OFCS with 256 carriers NPPB 2^10 Nb 2^9 fc 12.5e9 EDFA Power 20 dB
    %     OfcName = [ThisPath UofcsLocal 'OCS_fc_125_16-May-2018'];              %OFCS with 128 carriers NPPB 2^9 Nb 2^10 fc 12.5e9 EDFA Power 20 db
    %     OfcName = [ThisPath UofcsLocal 'OCS_fc_125_Car_128_Long_12-Jul-2018']; %OFCS with 128 carriers NPPB 2^9 Nb 2^12 fc 12.5e9 EDFA Power 20 db
    %     OfcName = [ThisPath UofcsLocal 'OCS_fc_125_Car_128_Long_13-Jul-2018']; %OFCS with 128 carriers NPPB 2^9 Nb 2^11 fc 12.5e9 EDFA Power 20 db
end
load([OfcName '.mat']);                                                    %Loading a OFC from memory.
if ~exist('freqGHz','var')
    tps        = t/1E-12;
    freqTHz    = time2freq_lamb_2(tps);
    freqGHz    = freqTHz*1e-3;			                                  % Frequencia em GHz
    freqGHz    = -fftshift(freqGHz);
    freqGHz(1) = freqGHz(2);
end

UsingGpu      = 0;
CurTesPas     = 1;
if UsingGpu == 1
    fgpu = gpuArray(f);
end
%%                 Control Variable

%       General Parameter for the simulation
ON            =  1;                                                        %DO NOT CHANGE THIS PARAMETER
OFF           =  0;                                                        %DO NOT CHANGE THIS PARAMETER
Ploting       = OFF;                                                       %Use On or Off to say if the figures should be plot or not
PlotingThis   = ON;                                                        %Use On or Off to say if the figures should be plot or not
Selecting     = ON;                                                        %Use to filter or not the carriers to be transmited
ModAll        = OFF;                                                       %Set if all carriers will be modulate at once or not
AddCP         = ON;                                                        %Add the Cycle Prefix on the data stream

CurTesExtSiz  = length(1:CurTesPas);
SelecSetUP    = 1;                                                         %Select the simulated setup 1 for original 0 for IFFT/OFFT at receptor
SendingDowStr = 1;                                                         %Set 1 to send data downstream 0 to not send.
SigAten       = 0;                                                         %Set 1 to atenuate the OFCS signal. Set to zero to turn off attenuation
SnrRef        = 0;                                                         %Set 1 to select the own signal for the noise reference at receptor.
%Set 0 for the pilote carrier as noise reference at receptor.
FFTSplit      = 1;                                                         %Set 1 to use the OFFT to split the income signal on receptor. Set 0 to
%use filter to split the income signal at receptor.

Nsamp         = 2^9;%1;%                                                   %Number of samples per bit, it was necessary because the number of carriers here used
RecFilBanPas  = 1;                                                         %Filter after photodetector, 1 filter Before noise, 0 filter After noise.
SendingUp     = 1;                                                         %Select if the carrier upstream will be send. Set 1 to send, set 0 to not send.
SaveToScatter = 0;                                                         %Set 1 to save points to create scatter plot latter Set 0 to not save points (keep it zero to save data space)
TimeSys       = 1;                                                         %Select which time period will be accounted to the noise, set 1 to use T or set 0 to use t(2*NumAmosCP+NPPB).
%When CP is used Seting TimeSys to 1 gives higher SNR
PastTime      = 0;                                                         %Storage for the time flow on simulation
% Atenuation    = 10^(4.0);                                                %Set how much the OFCS will be attenuated
NumbOfTest    = 10;                                                      %The number of repetition that one test will be done
TestBundle    = 1;                                                         %Set if the scrip will be run once or many times
CurrentMedium = 0;                                                         %Sellection of channel 1 for Fiber other of B2B
PAM4Type      = [1 0];                                                     %Sellection the type of 4PAM that will be tested
GatingAtReception = 0;
UseOfdmElect  = 1;
ChanAwgn      = 0;
MinTesNumb    = 10;
MinSavTim     = 15;
%The main results were generated here using a simple algorithm that I
%develop, which seek the biggest opening on a given interval. The interval
%is found by using the eye diagram histogram. This variable will be used to
%verify whether the method is more efficient than the original. If so, it
%can be one explanation of the gain reduction needed at the receiver site.
%As the SNR value founded here are lower than what was expected.
CalcS         = 1;                                                         %Set 1 to enable the first eye openeing stimative
%There are many points where noise can be added. It was not defined yet the
%right place for it. But four points are being pointed as candidates. At
%the electrical signal that modulates the carrier #1, after the signal
%passes through the fiber #2, after the carrier was modulated #3 and just
%before the photodetector #4.
AddingNoiseE  = 0;                                                         %#1 Set 1 to add noise to the electrical signal
AddingNoiseF  = 0;                                                         %#2 Set 1 to add noise to the optical signal
AddingNoiseO  = 0;                                                         %#3 Set 1 to add noise to the optical signal
AddingNoiseP  = 0;                                                         %#4 Set 1 to add noise to the optical signal
ReceptorNoise = 0;
%Two other noises are not concidered here. One is just after the OFCS was
%loaded. The second is just after the photodetector.
SetCpSampZer  = 0;                                                         %Set 1 to use the samples on the Cp as zeros
PrintinEye    = 0;                                                         %Printing Eye Diagram of modulations

% At this process the decission level is taken from the eye diagram, which
% means one statistical tool trying to get a deterministic variable. Thus,
% if the implemented method was not able to get the right decission level
% it value will be used. The value stored will be the  one with less error.
% Thos values were achieved from simulations. They must be updated whenever
% this code be modified on the generation of the signals to be transmited.
DecLevDef3 = 0.6953;                                                       %Uper levels of PAM4 - between level 3 and level 4
DecLevDef2 = 0.4013;                                                       %Uper levels of PAM4 - between level 2 and level 3
DecLevDef1 = 0.1877;                                                       %Uper levels of PAM4 - between level 1 and level 2
LevDefDpqsk = 0.01;                                                        %Default level for the DPSK modulation



%This will be the variable that will select the type of modulation used by
%this script. So far, this the options goes from 1 to 4, where:
%       1   -->    OOK
%       2   -->    4PAM
%       3   -->    DQPSK
%       4   -->    DPSK
%       5   -->    OFDM Electrical
% ThisModula = 1;
UsedModula = 5;
UpStModula = 5;
if UsedModula ~= 5
    UseOfdmElect = 0;
end
%The following verification will determinate which type of test will be on
%running, a test bundle or a single test. In the section the user may
%change with variables will be loaded for testing.
%Test Bundle selected
ThisTest   = NumbOfTest;                                                   %The number of tests per type of modulation
ThisModula = UsedModula;                                                   %The modulations under test

%%          Variables for the Medium Fiber
%This part of the simulation it is important to chose before hand some
%general parameters, that will affect the whole simulation.

if ~exist('UsingAllCarriers','var')                                        %Those variables will just be seted if any other outsider script haven't call them yet
    IfftOrSum    = 1;                                                      %Selection if the modulated signal will be assembled by just a adder or by the OIFFT
    %Set 1 for the OFFT or set 0 for the adder
    if (UsedModula==5)&&(UseOfdmElect==1)
        NumCarr  = 2;                                                      %Total number of carriers (DownStream and UpStream) - This variable must be even
    else
        NumCarr  = 128;                                                      %Total number of carriers (DownStream and UpStream) - This variable must be even
    end
    lambdac      = 1550e-9;                                                %Central wave length of the signal throughout the fiber
    
    CarrPass     = 2;                                                      %The carrier up and down can be interleaved into odd and even or can be sequential.
    %Set 2 to make carier odd and even interleaved or set 1 to make sequential.
    RefCarr      = 1;                                                      %Setting the reference of the first carrier to be used.
    if CarrPass~=1&&CarrPass~=2
        error(['CarrPass can NOT be different from 1 or 2, please check'...
            'assigned value']);
    end
    if ~mod(CarrPass,2)                                                    %Setting parrameter when carriers were interleaved in odd and even order
        CarPosAux  = 1:NumCarr;                                            %Auxiliar variable to place each carrier possition accordingly
        CarrUsed   = ones(1,NumCarr);                                      %Setting flags on the carriers that will be used 1 for used 0 for unused
        CarrUsedDo = zeros(1,NumCarr);                                     %Starting downstream used carriers flags as all zero
        CarrUsedUp = zeros(1,NumCarr);                                     %Starting upstream used carriers flags as all zero
        CarrUsedDo(logical(mod(CarPosAux,2))) = 1;                         %Setting the downstream carrier as valid to be used
        CarrUsedUp(~mod(CarPosAux,2))         = 1;                         %Setting the upstream carrier as valid to be used
        %         CarrUsedDo(1) = 1;                                               %Setting the downstream carrier as valid to be used
        %         CarrUsedUp(2) = 1;                                               %Setting the upstream carrier as valid to be used
        InitCarrUp = 2;                                                    %Setting that the frist carrier of the upstream signal is at the 2 possition
        InitCarrDo = 1;                                                    %Setting that the frist carrier of the downstream signal is at the 1 possition
    else                                                                   %Setting parrameter when carriers were placed sequentialy
        CarrUsed   = ones(1,NumCarr);                                      %Setting flags on the carriers that will be used 1 for used 0 for unused
        %         CarrUsed   = zeros(1,NumCarr);                                   %Setting flags on the carriers that will be used 1 for used 0 for unused
        CarrUsedDo = CarrUsed;                                             %Starting downstream used carriers flags as all one
        CarrUsedDo(1) = 1;                                                 %Starting downstream used carriers flags as all one
        CarrUsedDo(1+end/2:end) = 0;                                       %Turnning off, on the downstream sequence, the upstream carriers
        CarrUsedUp = CarrUsed;                                             %Starting upstream used carriers flags as all ones
        CarrUsedUp(1+end/2) = 1;                                           %Starting upstream used carriers flags as all ones
        CarrUsedUp(1:end/2) = 0;                                           %Turnning off, on the upstream sequence, the downstream carriers
        InitCarrUp = (NumCarr/2) + 1;                                      %Setting that the frist carrier of the upstream signal is at the 1 possition
        InitCarrDo = 1;                                                    %Setting that the frist carrier of the downstream signal is at the 1 possition
    end
    % The optical FFT or the optical IFFT are very precise and needs that each
    % carrier be at the exactly and expected possition, because it is a
    % perriodic filter basicaly. Thus, if for some reason the first used
    % carrier is not the first carrier available they carriers used must be
    % mapped into a sequency that was expected by the MZIs clusters. If the
    % sequence incoming into the MZIs was not correct the FFT/IFFT will not
    % work.
    if ~mod(RefCarr,NumCarr)                                               %Luckly if the first carrier was a multiple of the number of used carriers
        ObsCarrPos = 1:NumCarr;                                            %nothing different needs to be done.The carriers order is straighforward
        ObsCarrUsed = [mod(RefCarr,NumCarr)+1:NumCarr 1:mod(RefCarr,...
            NumCarr)];                                          %As in this simulation all carrier will be used the observed carriers is equal
        %to the possition of the observed carriers
    else                                                                   %whereas, when the first carrier was not multiple of the number of carrier used
        ObsCarrPos  = [mod(RefCarr,NumCarr):NumCarr 1:mod(RefCarr,...
            NumCarr)-1];%The first carrier is what left from the divission Reference carrier by
        %the number of carriers used and this sequency goes untill the last carrier
        %selected. Then, the following carrier will be from 1 untill the divission result minus 1
        ObsCarrUsed = ObsCarrPos;                                          %As in this simulation all carrier will be used the observed carriers is equal
        %to the possition of the observed carriers
    end
%     SaveBerEbn0PreFix   = [pwd '\..\results\UsedModula_'];
%     SaveBerEbn0SuFix    = '_BerVsEbnONoOverSamp_1';
%     SaveBerEbn0FigSuFix = '_BerVsEbnONoOverSamp_fig_1.fig';
        SaveBerEbn0PreFix   = [pwd 'test0'];
        SaveBerEbn0SuFix    = '';
        SaveBerEbn0FigSuFix = '';
    
%     SavingAT     = [pwd '\..\results\' 'Ofdm4QamBerVsCarVsKm_' num2str(VetSnr) '_1'];
        SavingAT     = 'test_0';
end
%%    Simulation Parameters
SplitRatio  = 128;                                                         %Number of usser for the sistem
NRZPolarity = 1;                                                           %Polarity of NRZ coding, 0 for one-polar, 1 for bi-polar
Which4PAM   = PAM4Type(1);                                                 %Which 4PAM will be used, electrical "0" or optical "1"
MaxNumStagT = nextpow2(NumCarr);                                           %Next number, power of 2, by the given number of carriers
%%    Constants to controu the Optical FFT
% The optical FFT here implemented can take any number of carriers. The
% user must stay realistic because the number of actual devices used to
% implement this setup increases by the power of 2, thus higher orders of
% the FFT may not be feaseble to be implemented as the stability of each
% device may vary with temperature, aging and othter features.
MaxNumStagTx = nextpow2(NumCarr+1);                                        %Next number, power of 2, by the given number of carriers
MaxNumCarrTx = 2^MaxNumStagTx;                                             %With the number of stages the actual nubmer of carriers can be calculated%

MaxNumStag = nextpow2(NumCarr);                                            %Next number, power of 2, by the given number of carriers
MaxNumCarr = 2^MaxNumStag;                                                 %With the number of stages the actual nubmer of carriers can be calculated

%As result of the OFFT process the carriers get mingled at the output ports
%of this sistem. Thus, those vector were created as a map for the
%simulation. Where the possition of the vector is the output of the OFFT
%and the value within tha possition is the respective carrier filtered. The
%folloing ilustration better shows how each carrier get mingled in this
%process:
%               Ein : Is the input fild
%               MZI : Is the Mach-Zehnder Interferometre
%               Nº  : The number are each carrier from 1 to 32 within Ein
%
%                                                  /1
%                                             (MZI)
%                                            / 1   \17
%                                       (MZI)
%                                      /     \ 9   /9
%                                     / 1     (MZI)
%                                (MZI)             \25
%                                /    \ 5          /5
%                               /      \      (MZI)
%                              /        \    / 5   \21
%                             /         (MZI)
%                            /               \ 13  /13
%                           / 1               (MZI)
%                       (MZI)                      \29
%                       /   \ 3                    /3
%                      /     \                (MZI)
%                     /       \              / 3   \19
%                    /         \        (MZI)
%                   /           \      /     \ 11  /11
%                  /             \    / 3     (MZI)
%                 /              (MZI)             \27
%                /                    \ 7          /7
%               /                      \      (MZI)
%              /                        \    / 7   \23
%             /                          (MZI)
%            /                               \ 15  /15
%           /                                 (MZI)
%          / 1                                     \31
%Ein---(MZI)
%          \ 2                                     /2
%           \                                  (MZI)
%            \                               / 2   \18
%             \                         (MZI)
%              \                       /     \ 10  /10
%               \                     / 2     (MZI)
%                \               (MZI)             \26
%                 \              /    \            /6
%                  \            /      \ 6    (MZI)
%                   \          /        \    / 6   \22
%                    \        /         (MZI)
%                     \      /               \ 14  /14
%                      \    / 2               (MZI)
%                      (MZI)                       \30
%                           \ 4                    /4
%                            \                (MZI)
%                             \              / 4   \20
%                              \        (MZI)
%                               \      /     \ 12  /12
%                                \    / 4     (MZI)
%                                (MZI)             \28
%                                     \ 8          /8
%                                      \      (MZI)
%                                       \   / 8    \24
%                                       (MZI)
%                                           \ 16   /16
%                                             (MZI)
%                                                  \32
%Each MZI vertically aligned correspond to one Stage of the OFFT. The
%following verification creates the maping vector as previouly described.
%IMPORTANT! This mapping will only work if the first carrier is at 12.5 Ghz



%%      Parameters for the Cicle Prefix
% At this point, it is important to add the cycle prefix (CP) to mitigate
% the inter channel interference (ICI). For that there will be a reduction
% on the Rb. But to avoid complex modification on those scrips, the length
% of the period of one symbol will be altered accordingly with the size of
% the CP. For the whole simulation works it will be necessary some extra
% samples (StuffSamples) for matching the size of vector created. In real
% life, it would not be a problem because the signal would be continous not
% sampled.
if AddCP
    PerUtil         = 0.122;                                               %Percentage for the ultil period of symbol
    NumAmosCP       = ceil(PerUtil*NPPB);                                  %Half of the number of samples for the CP
    NewTotalSamples = (2*NumAmosCP+NPPB)*Nb;                               %New total samples of the simulation
    ExtraSamples    = NewTotalSamples - TotalSamples;                      %Extra samples created
    NumBitDesc      = ceil(ExtraSamples/(2*NumAmosCP+NPPB));               %Number of bits excided
    if mod(NumBitDesc,2)                                                   %Verify if the amount of bits to remove is odd
        NumBitDesc = NumBitDesc + 1;                                       %This simulation works just with even number of bit stream
    end                                                                    %Thus, adding one more unit to be removed is needed
    StuffSampels    = TotalSamples - (2*NumAmosCP+NPPB)*(Nb - NumBitDesc); %Number of samples for stuffing process
    Tcp             = t(2)*(2*NumAmosCP+NPPB);
else
    PerUtil         = 0;                                                   %Percentage for the ultil period of symbol
    NumAmosCP       = ceil(PerUtil*NPPB);                                  %Half of the number of samples for the CP
    NewTotalSamples = (2*NumAmosCP+NPPB)*Nb;                               %New total samples of the simulation
    ExtraSamples    = NewTotalSamples - TotalSamples;                      %Extra samples created
    NumBitDesc      = ceil(ExtraSamples/(2*NumAmosCP+NPPB));               %Number of bits excided
    if mod(NumBitDesc,2)                                                   %Verify if the amount of bits to remove is odd
        NumBitDesc = NumBitDesc + 1;                                       %This simulation works just with even number of bit stream
    end                                                                    %Thus, adding one more unit to be removed is needed
    StuffSampels    = TotalSamples - (2*NumAmosCP+NPPB)*(Nb - NumBitDesc); %Number of samples for stuffing process
    Tcp             = T;
end
%%     Generation the File for MZM Configuration
%%               MZM Parameters
%Every MZM hereafter used have a configuration file. Therefor, it is needed
%to creat the Input Data corectly to the aplication need. For the
%simplified MZM function there are already one file used to generate the
%COMB. Thus, a new file need to be created. This is why the Vbias_steps are
%seted to be 2 (the second configuration file) and it is also why the Vbias
%is a vector with two inputs. Because basiclay it will rewrite the first
%input file with the same information as before.

L          = 10;
U0         = 1.900;
U_pi1      = 3.8;
U_pi2      = 3.8;
eta1       = 89;% 89;
eta2       =-12.7;%-12.7;
%
nopt       =  2.17;
nel        =  2.60;
alfa_ins   =  5.1;
phi_0      =  0.0;
alfa0      =  0.55;

C     = (eta1+eta2)/(eta1-eta2);       	%  Parametro de chirp
alfa0 = 10^(alfa0/20); 			% [1/cm.GHz^0.5]

%At this point the file with the MZM will be created. The first file will
%have the polarization point at Vpi/2 because it was used to generate the
%UfOCS. The second file is used by the 4PAM and by the OOK modulation,
%because they are polarized at Vpi/2 as well. The third file will be used
%by the DPSK modulation as it is polarized at Vpi. The last one will vary
%in value for testing, therefore it can assume any value.
Vbias_steps = 5;                                                           %Seting the number of variation of the Input file
Vbias = [ U_pi1/2 U_pi1/2 U_pi1 0 -U_pi1/2];                               %Creating the Vbias that will be writen on the Inpu Data file
Local = [pwd '\input_files\'];                                             %Seting the location where this file will be stored.
if ~exist('MzFilesGenerated','var')
    S=1; %Variabel that give the name for the input file
    for Inc2=1:Vbias_steps%Secundary loop for the possible combinations
        U0               = Vbias(Inc2);%Variation of the Vbias
        [MZ_Input_Sdata] = Set_MZ_Input_Data_Simp(S,L,U0,U_pi1,U_pi2,eta1,...  %Functionresponsible to create the input files
            eta2,nopt,nel,alfa_ins,phi_0,alfa0,Local);
        S                = S+1;
    end
    %MzFilesGenerated = 1;
end
%%         Creating file for the DP-MZM
phi_0_steps  = 1;%Number of different phase between arms of the DP-MZM
V1bias_steps = 1;%Number of different Bias voltage at the Up arm
V2bias_steps = 1;%Number of different Bias voltage at the Lower arm
%
%IMPORTANT: The RF amplitude Vpp must respect the device limits!
phi_0_vet     = -2.05;%Variation of phase betwee arms.
V1bias        = 1.752;%Variation of bias voltage at the Up arm.
V2bias        = 1.052;%Variation of bias voltage at the Lower arm.

% Differently from the previouly Input Data file. This scrip does not have
% a previous configuration file for the MZM. Thus, it can create a new file
% to configure the MZM-IQ.
Local_Dp = [pwd '\input_files_dp\'];%Seting the location where this file will be stored.

% L          = 10;
% U10        = 1.550;
% U20        = 1.625;
% U_pi1      = 3.1;
% U_pi2      = 3.25;
% eta1       = 89;% 89;
% eta2       =-12.7;%-12.7;
% V1pi0      = 0.350;
% V2pi0      = 1.200;
% Vphi0      = 0.100;
% %
% nopt       =  2.17;
% nel        =  2.60;
% alfa_ins   =  5.1;
% phi_0      =  1.85;
% alfa0      =  0.55;
% 
% C     = (eta1+eta2)/(eta1-eta2);%Parametro de chirp
% alfa0 = 10^(alfa0/20);%[1/cm.GHz^0.5]
% if ~exist('MzFilesDpGenerated','var')
%     S=1;%Variabel that give the name for the input file
%     for Inc0=1:phi_0_steps
%         for Inc1=1:V1bias_steps%Secundary loop for the possible combinations
%             for Inc2=1:V2bias_steps%Secundary loop for the possible combinations
%                 phi_0 = phi_0_vet(Inc0);
%                 U10 = V1bias(Inc1);%Variation of the Vbias
%                 U20 = V2bias(Inc2);%Variation of the Vbias
%                 [MZ_Input_Sdata] = Set_MZ_Input_Data_DP(S,L,U10,U20,...%Functionresponsible to create the input files
%                     U_pi1,U_pi2,eta1,eta2,V1pi0,V2pi0,Vphi0,nopt,nel,...
%                     alfa_ins,phi_0,alfa0,Local_Dp);
%                 S = S+1;
%             end
%         end
%     end
%     %MzFilesDpGenerated = 1;
% end
% MZ_Input_File_Dp = 1;%Telling the script witch file will be used for modulation
%%           Control Variables
%This variable will be used to select the type of modulation used for the
%datum Enconding. Type 'DQPSK' for DQPSK encoding, '4PAM' for 4PAM
%encoding. Type 'NoN' if an encoding is not necessary.
switch UsedModula
    case 1
        Modulation = 'NoN';
    case 2
        Modulation = '4PAM';
    case 3
        Modulation = 'DQPSK';
    case 4
        Modulation = 'DPSK';
    case 5
        Modulation = 'OFDM';
    otherwise
        Modulation = 'NoN';
end

%%          Variables for the Medium Fiber

%This variable will be used to select the medium in which the information
%will travel Type 'B2B' for a Back-To-Back transmission or 'Fiber' for an
%optical Fiber
% Medium = 'B2B';%'Fiber';%
switch CurrentMedium
    case 1
        Medium       = 'Fiber';
    otherwise
        Medium       = 'B2B';
        FiberLength  = 0;                                                  %Length of the optical fiber [km]
end
%%          Variables for different types of Modulations
switch Modulation
    %The Electrica OFDM is based on some factor that needs to be calculated
    %before creating the OFDM signal.
    %The OFDM signal may be specified from different ways. This process
    %choose to start the OFDM modelation by seting a fixed number of FFT
    %points. This was riquered because the time vector of this simulation
    %was fixed. Thus, the symbols OFDM must fit withing the simulation
    %time.
    %####################OFDM PROJECT GOLD RULE################################
    %So far this simulation is limited by the vector length, which is given by:
    %                     Length = Nb*T/Ta
    %It is limitade by the numer of bits time the sample per bit (NPPB=T/Ta)
    %where T is the symbol period and Ta is the samping period. Thus the OFDM
    %frames are limited by:
    %
    %  !!GOLD RULE!! Nb*NPPB >= NumFram*NumFFT*NumPonPerElecCarr !!GOLD RULE!!
    %
    %Thus the Electrical OFDM project must attend this riquerment. Othets
    %constraints are the maximum simulation time which is Nb*T which is also
    %the minimal bandwidth allowed. The maximum bandwidth is 12.5 GHz
    %##########################################################################
    case 'OFDM'
        %% Parameters for the DownStream Recptor
        U0        = 1.9;
        BWD       = 1.0*fc;                                               %Badnwidth for the bit conformation
        CenFeq    = 0;                                                    %Center frequency for the bit conformation
        FiltOrd   = 0;                                                    %The order for the filter which will conform data
        %% Parameters for the DownStream Transmiter
        SelModFilt= 0;
        SelecGaus = 1;
        OfdMod    = 'qam';                                                %Chossing the type of modulation for the electrical subcarriers
        SelModTp  = 'BP';                                                  %Selecting with OFDM will be in Base Band "BP", Amplitude Modulation
        %with double side band "AM" or with single side band "AMSSB" when a
        %different option is choosen this script will implement phase modulation.
        UsingHermitian = 1;
        if ~UseOfdmElect
            UsingHermitian = 1;
        end
        ForcOvSam = 0;
        UsingParCar= 0;
        NuDmSuDi  = 5;
        DmtPas    = 2;
        DmtMinM   = 1;
        TapN      = 1;                                                     %Number of equalizer used
        ExtSam    = 1;                                                     %Set to the number of frames that will be used.
        NpEl      = 2^9;                                                  %Number of Electrical carriers
        NFFT      = NpEl/1;                                                %Size of the electical FFT
        OveT      = 2^9;                                                   %Variable to set the eletrical carrier period to be OveT time the symbol period
        M         = 4;                                                     %Modulation Level
        % Not all electrical carriers may be used. For those unused
        %carriers, zeros are placed in their locations. As a result, there
        %are free spaces to place the used carriers as the Zt variable
        %control how many zeros will be inserted on the FFT frame to offset
        %the used carriers from the beginning. The value of Zt is
        %approximately half of the number of the free carrier. Thus,
        %ControlZt variable controls the percentage that Zt will influence
        %the FFT offset. As a result, more flexibility was possible for
        %placing the used carriers on the available spectrum.
        ControlZt = 1;
        CpLe   = 1;                                                      %Set the percentage of electrical carriers that will be used
        if CpLe <= 0.5
            error(['CpLe can not be bellow 50%. More than half of ' ...
                'available carriers must be transmited']);
        end
        BW     = fc;                                                       %Signal available bandwidth
        %         NumFra = floor(BW*t(end));                                       %Number of frames
        %         if NFFT>Nb
        %             NFFT = Nb;
        %         end
        DatSiz = NFFT;                                                     %Set to the total number of electrical carriers to be used
        Te     = OveT*T;
        Tu     = Te;                                                       %Time that will actually be used to transmit information
        Tg     = 0.0*Tu;                                                     %Time for the guard band
        if ((NFFT*NpEl)==Nb)&&(Tg~=0)
            error('NFFT must be smaller than Nb to use some ciclic prefix');
        end
        Ts     = Tu + Tg;                                                  %Time of the OFDM symbol
        g      = Tu/Ts;                                                    %Percentage of band guard
        Dtf    = 1/Tu;                                                     %Chanal signaling ratio
        if UsingHermitian
            N      = NFFT/2 - 1;                                           %Number of carriers found
            Ns     = ceil((CpLe*(DatSiz))/2) - 1;                          %Number of carriers used
        else
            N      = NFFT;                                                 %Number of carriers found
            Ns     = ceil((CpLe*(DatSiz)));                                %Number of carriers used
            if mod(Ns,2)
                Ns = Ns - 1;
            end
        end
        DifNsN = ceil((N - Ns)/1);                                         %Number of carriers unused
        
        k      = log2(M);                                                  %Number of bits per symbol
        MZ_Input_File = 2;                                                 %Inputdata file to the MZM
        NuSaTs = Nb*NPPB;
        NuSaTu = ceil(g*NPPB);
        NuSaTg = NuSaTs - NuSaTu*Nb;
        OBw    = NFFT/Tu;                                                  %OFDM bandwidth
        Ofc    = 1.0*OBw;                                                  %Central frequency of out of base band transmission
        Ofs    = 2*Ofc;                                                    %Electrical sampling frequency
        OvSam = round(Te/t(2)/NFFT);                                       %Number of over samples for out of base band modulations
        if UsingHermitian
            NuPart = floor(Ns/NuDmSuDi);
            NuSuSe = floor(Ns/NuPart);
            DmtMve = ones(1,Ns);
            for ElCaM=NuSuSe:-1:2
                DmtMve(1+(NuSuSe-ElCaM)*NuPart:(NuSuSe-ElCaM+1)*NuPart)= 2^(ElCaM);%2^(ElCaM*DmtPas-1);
            end
        else
            NuPart = floor(Ns/NuDmSuDi);
            NuSuSe = floor(Ns/NuPart);
            DmtMve = ones(1,Ns/2);
            DmtMve(:) = 2^DmtMinM;
            for ElCaM=NuSuSe:-1:2
                DmtMve(1+(NuSuSe-ElCaM)*round(NuPart/2):(NuSuSe-ElCaM+1)*round(NuPart/2))= 2^(ElCaM);%2^(ElCaM*DmtPas-1);
            end
            DmtMve = [DmtMve fliplr(DmtMve)];
            %     figure;hold on;plot(1:Ns,DmtMve)
        end
        switch OfdMod
            case 'dpsk'
                DmtMve(:) = 2^DmtMinM;
            otherwise
                if NuDmSuDi==1
                    DmtMve(:) = 2^DmtMinM;
                else
                    DmtMve(1+(NuSuSe-ElCaM+1)*NuPart:end)= 2^DmtMinM;
                end
        end
        
        if UsingHermitian
            nRb    = (OBw/(2*(N+2)))*(Ns*mean(log2(DmtMve)));                                       %New Rb found
        else
            nRb    = (OBw/(N))*(Ns*mean(log2(DmtMve)));                                       %New Rb found
        end
        
        nRbAgg = nRb*NumCarr/2;                                            %New Aggregate Symbol ratio found
        %As not all electrical carries may be used, the user have the
        %choise to select where the OFDM signal will be centralized within
        %all carrier available. Zt set were the OFDM signal will start. If
        %Zt was set to zero, the OFDM signal will start at possition 1.
        %Thus, all unsued carriers will placed after it. If Zt was set to
        %10 the OFDM signal will start at electrical carrier 11. This is
        %important when using Hermitiam symetri. Zt will determinate how
        %further away the OFDM signal will be from the electrical carrier.
        if DifNsN > 51
            Zt  = 8;%floor(DifNsN/2);
        else
            if DifNsN < 1
                Zt = 0;
            else
                Zt  = floor(DifNsN/2);
            end
        end
        ZtC = floor(Zt*ControlZt);
        switch SelModTp
            case 'AM'
                NumFraPar = 1;
                FramSpac = 1;
                ForcOvSam = 0;
            case 'AMSSB'
                NumFraPar = 1;
                FramSpac = 1;
                ForcOvSam = 0;
            otherwise
                OBw    = NFFT/Tu;                                          %OFDM bandwidth
                Ofc    = 0*OBw;                                            %Central frequency of out of base band transmission
                Ofs    = 2*Ofc;                                            %Electrical sampling frequency
                nRb    = (OBw/(N+2))*(Ns*k);                               %New Rb found
                FramSpac = 2*OBw;
                if UsingParCar==1
                    NumFraPar = floor(Rb/((OBw+FramSpac)*1));
                    %                 NumFraPar = 3;
                    if (NumFraPar>1)&&(~mod(NumFraPar,2))
                        NumFraPar = NumFraPar - 1;
                    end
                    if NumFraPar>0
                        %                     OvSam  = round(Te*2*OBw*2^((NumFraPar-1)/2));                   %Number of over samples for out of base band modulations
                    else
                        %                     OvSam  = round(Te*2*OBw);
                        NumFraPar = 1;
                    end
                else
                    NumFraPar = 1;
                end
                OvSam = round(Te/t(2)/NFFT);
        end
        if ForcOvSam
            OvSam = ForcOvSam;
        end
        NumFra    = floor((Nb*NPPB)/((OvSam*NFFT)+(OvSam*NFFT*Tg/Tu)));
        NumFraTem = NumFra/NumFraPar;               %Number of frames
        %         if NumFra>1
        %             if mod(NumFra,2)
        %                 NumFra = NumFra - 1;
        %             end
        %         end
        NumFrU = floor((Nb*NPPB)/(OvSam*NFFT));                            %Number of frames
        NuAmOf = (OvSam*NumFra*NFFT);                                      %Total number of samples of the electrical OFDM signal
        NPOFEX = (NuAmOf/NumFra)*(Tg/Tu);                                  %Extra samples per OFDM frame. It is the number of CP
        NuAmTo = NuAmOf*(1+Tg/Tu);                                         %Total number of samples of the electrical OFDM signal added CP
        NPPOF = floor((Nb*NPPB)/NuAmTo);                                   %Over sampling for optical modulation transmission
        if (NuAmOf<=0)||(NPPOF<=0)
            error('OFDM paramters outof range. Please double check them.');
        end
        if NumFra/OBw > Nb*T
            error(['The simulation is exiding the maximum allowed time '...
                '.Please refrain the amount of serial frames used.']);
        end
        NPSTUf = Nb*NPPB-NPPOF*NuAmTo;                                     %Number of extra sample to make all vector to have the same length
        NumFrE = NumFrU - NumFra;                                          %Number of frames
        
        
        %         NPOFEX = floor(NumFrE*OvSam*NFFT/NumFra);
        %         NPPOF  = floor((2^19)/((NPOFEX*NumFra)+(OvSam*NumFra*NFFT)));
        %         NPSTUf = 2^19 - NPPOF*NumFra*OvSam*NFFT;
        %Variables for the Selection of carriers - Wave Divided.
        fin           = RefCarr*fc;                                        %The frequency of the first carrier of the optical COMB
        FBWD          = 1.0*fc;                                            %The filter BandWidth to split each individual carrier
        Order         = 5;                                                 %The order of the filter used to split each carrier
        a=a+1;
    case 'DPSK'
        
        U0            = 3.8;
        %% Parameters for the DownStream Recptor
        BWD2          = 1*fc;                                          %Badnwidth for the bit conformation
        CenFeq2       = 0;                                                 %Center frequency for the bit conformation
        FiltOrd2      = 1;                                                 %The order for the filter which will conform data
        SyncSymb      = [zeros(1,10) ones(1,4) zeros(1,10)];               %Generating the sync-word at the beginning of the data stream
        IniSyncPos    = 8*(2*NumAmosCP+NPPB);
        SyncPos       = 16*(2*NumAmosCP+NPPB);                             %Seting the possition where the sync-symble is
        SyncSymbEnd   = SyncSymb;                                          %Generating the sync-word at the end of the data stream
        SyncPeriod    = length(SyncSymb);                                  %Seting how many symble periods will be within the sync-word
        AuxSyncLength = SyncPeriod*(2*NumAmosCP+NPPB);                     %auxiliar variable to set where is the sync-symble
        tsync         = linspace(0,SyncPeriod*T,(2*NumAmosCP+NPPB)*...
            SyncPeriod);%create a time vetor for the sync-word to help on the ploting of this signal
        fsync         = time2freq(tsync);                                  %creating the the frequency vector of the sync-word for hereafter filtering
        BWD           = 2*fc;                                          %Badnwidth for the bit conformation
        CenFeq        = 0;                                                 %Center frequency for the bit conformation
        FiltOrd       = 1;                                                 %The order for the filter which will conform data
        BitFilBanWid  = 1*fc;                                              %Badnwidth for the bit conformation
        BitFiltCenFre = 0;                                                 %Center frequency for the bit conformation
        BitFilOrd     = 1;                                                 %The order for the filter which will conform data
        [BitFilt,~]   = FiltroGaussiano(fsync,BWD,CenFeq,FiltOrd);         %Creating filter for conformation of the input information
        BitFilt       = fftshift(BitFilt);
        PerDiv        = 8;
        
        SyncSymb      = rectpulse(SyncSymb,2*NumAmosCP+NPPB);
        SyncSymbEnd   = rectpulse(SyncSymbEnd,2*NumAmosCP+NPPB);
        SencondAdjust = OFF;                                               %Setting whether is needed to sync the information by the Syncronism Symble at the End
        LowSymPer     = 1;
        UpeSymPer     = 1;
        %% Parameters for the DownStream Transmiter
        %%      Variables for the DPSK Datum Modification
        %This variable allow the user to chose to send an random
        %information or just a high pulse within the midle of the datastram
        SendingData   = 1;                                                 %Selecting the nature of the information to be send
        NbDPSK        = Nb - NumBitDesc;                                   %The total number of bits to be transmited
        DatGai        = 1;                                                 %Seting the maximum range for the eletrical signal
        JusVal        = [zeros(1,10) ones(1,4) zeros(1,10)];               %Value that will be load into the justificaiton space
        JusValEnd     = JusVal;                                            %Value that will be load into the justificaiton space
        JusPos        = length(JusVal);                                    %Length of the justification symbol
        BWD           = 2*fc;                                              %Badnwidth for the bit conformation
        CenFeq        = 0;                                                 %Center frequency for the bit conformation
        FiltOrd       = 1;                                                 %The order for the filter which will conform data
        MZ_Input_File = 3;                                                 %Telling the script witch file will be used for modulation
        %Variables for the Selection of carriers - Wave Divided.
        fin           = RefCarr*fc;                                        %The frequency of the first carrier of the optical COMB
        FBWD          = 1.0*fc;                                            %The filter BandWidth to split each individual carrier
        Order         = 5;                                                 %The order of the filter used to split each carrier
    case 'DQPSK'
        %% Parameters for the DownStream Recptor
        BWD2          = 1.0*fc;                                        %Badnwidth for the bit conformation
        CenFeq2       = 0;                                                 %Center frequency for the bit conformation
        FiltOrd2      = 1;                                                 %The order for the filter which will conform data
        SyncSymb      = [zeros(1,10) ones(1,4) zeros(1,10)];               %Generating the sync-word at the beginning of the data stream
        IniSyncPos    = 8*(2*NumAmosCP+NPPB);
        SyncPos       = 16*(2*NumAmosCP+NPPB);                             %Seting the possition where the sync-symble is
        SyncSymbEnd   = SyncSymb;                                          %Generating the sync-word at the end of the data stream
        SyncPeriod    = length(SyncSymb);                                  %Seting how many symble periods will be within the sync-word
        AuxSyncLength = SyncPeriod*(2*NumAmosCP+NPPB);                     %auxiliar variable to set where is the sync-symble
        tsync         = linspace(0,SyncPeriod*T,(2*NumAmosCP+NPPB)*...
            SyncPeriod);%create a time vetor for the sync-word to help on the ploting of this signal
        fsync         = time2freq(tsync);                                  %creating the the frequency vector of the sync-word for hereafter filtering
        BWD           = 2*fc;                                          %Badnwidth for the bit conformation
        CenFeq        = 0;                                                 %Center frequency for the bit conformation
        FiltOrd       = 1;                                                 %The order for the filter which will conform data
        BitFilBanWid  = 1*fc;                                              %Badnwidth for the bit conformation
        BitFiltCenFre = 0;                                                 %Center frequency for the bit conformation
        BitFilOrd     = 1;                                                 %The order for the filter which will conform data
        [BitFilt,~]   = FiltroGaussiano(fsync,BWD,CenFeq,FiltOrd);         %Creating filter for conformation of the input information
        BitFilt       = fftshift(BitFilt);
        PerDiv        = 4;
        
        SyncSymb      = rectpulse(SyncSymb,2*NumAmosCP+NPPB);
        SyncSymbEnd   = rectpulse(SyncSymbEnd,2*NumAmosCP+NPPB);
        SencondAdjust = ON;                                                %Setting whether is needed to sync the information by the Syncronism Symble at the End
        LowSymPer     = 1;
        UpeSymPer     = 1;
        %% Parameters for the DownStream Transmiter
        %%      Variables for the DQPSK Datum Modification
        %This variable allow the user to chose to send an random
        %information or just a high pulse within the midle of the datastram
        SendingData = 1;                                                   %Selecting the nature of the information to be send
        NbDQPSK     = 2*(Nb - NumBitDesc);                                 %The total number of bits to be transmited
        JusVal      = [zeros(1,20) ones(1,8) zeros(1,20)];                 %Value that will be load into the justificaiton space
        JusValEnd   = JusVal;                                              %Value that will be load into the justificaiton space
        JusPos      = length(JusVal);                                      %Length of the justification symbol
        BWD         = 4*fc;                                                %Badnwidth for the bit conformation
        CenFeq      = 0;                                                   %Center frequency for the bit conformation
        FiltOrd     = 1;                                                   %The order for the filter which will conform data
        DatGai      = 1;                                                   %Gain to increase the modulatin signal amplitude
        V0          = 3.8;
        Vpi         = 3.8;
        %Variables for the Selection of carriers - Wave Divided.
        fin         = RefCarr*fc;                                      %The frequency of the first carrier of the optical COMB
        FBWD        = 1.0*fc;                                              %The filter BandWidth to split each individual carrier
        Order       = 5;                                                   %The order of the filter used to split each carrier
    case '4PAM'
        DecMod        = 1;  %Selec how the decission level will be measured
                            % 1 for measuring the zeros on the eye diagram (best)
                            % 0 for taking the meddle of the zeros vector
        U0 = 1.9;
        U0Up = -1.9;
        %% Parameters for the DownStream Recptor
        SyncSymb      = [zeros(1,10) ones(1,4) zeros(1,10)];               %Generating the sync-word at the beginning of the data stream
        IniSyncPos    = 8*(2*NumAmosCP+NPPB);
        SyncPos       = 16*(2*NumAmosCP+NPPB);                             %Seting the possition where the sync-symble is
        SyncSymbEnd   = SyncSymb;                                          %Generating the sync-word at the end of the data stream
        SyncPeriod    = length(SyncSymb);                                  %Seting how many symble periods will be within the sync-word
        AuxSyncLength = SyncPeriod*(2*NumAmosCP+NPPB);                     %auxiliar variable to set where is the sync-symble
        tsync         = linspace(0,SyncPeriod*T,(2*NumAmosCP+NPPB)*...
            SyncPeriod);%create a time vetor for the sync-word to help on the ploting of this signal
        fsync         = time2freq(tsync);                                  %creating the the frequency vector of the sync-word for hereafter filtering
        BWD           = 2*fc;                                          %Badnwidth for the bit conformation
        CenFeq        = 0;                                                 %Center frequency for the bit conformation
        FiltOrd       = 1;                                                 %The order for the filter which will conform data
        if PAM4Type(2)
            BWD2          = 1*fc;                                    %Badnwidth for the bit conformation
            CenFeq2       = 0;                                             %Center frequency for the bit conformation
            FiltOrd2      = 1;                                             %The order for the filter which will conform data
        else
            BWD2          = 1*fc;                                    %Badnwidth for the bit conformation
            CenFeq2       = 0;                                             %Center frequency for the bit conformation
            FiltOrd2      = 1;                                             %The order for the filter which will conform data
        end
        [BitFilt,~]   = FiltroGaussiano(fsync,BWD,CenFeq,FiltOrd);         %Creating filter for conformation of the input information
        BitFilt       = fftshift(BitFilt);
        SyncSymb      = rectpulse(SyncSymb,2*NumAmosCP+NPPB);              %Resampling the sync-word accordingly with the current system
        SyncSymbEnd   = rectpulse(SyncSymbEnd,2*NumAmosCP+NPPB);           %Resampling the sync-word accordingly with the current system
        IntervalStep  = 50;
        MinDist       = IntervalStep/6;
        %Because of reasons, the real world does not have components
        %capable instantly change the voltage level. Thus, to emulate this
        %behaviour the eletrical PAM signal pass through a gaussian filter
        %that will remove the frequencies of higher order, which will
        %result in small slope in the level variation
        SyncSymb      = ifft(fft(SyncSymb).*BitFilt);
        SyncSymbEnd   = ifft(fft(SyncSymbEnd).*BitFilt);
        SencondAdjust = OFF;                                               %Setting whether is needed to sync the information by the Syncronism Symble at the End
        LowSymPer     = 0.9;
        UpeSymPer     = 1.1;
        %% Parameters for the DownStream Transmiter
        %%      Variables for the 4PAM Datum Modification
        %This variable allow the user to chose to send an random
        %information or just a high pulse within the midle of the datastram
        SendingData   = 1;                                                 %Selecting the nature of the information to be send
        Nb4Pam        = 2*(Nb - NumBitDesc);                               %The total number of bits to be transmited
        %This variable controls if the electrical signal will vary from
        %0[V] to max[V] (signal unipolar) or from -max[V] to max[V]
        %(signal bipolar).
        Polirized     = 1;
        MaxAmp4PAM    = -3.8/4;                                            %Seting the maximum range for the eletrical signal
        JusVal        = [zeros(1,20) 1 0 1 0 1 0 1 0 zeros(1,20)];               %Value that will be load into the justificaiton space
        JusValEnd     = JusVal;                                            %Value that will be load into the justificaiton space
        JusPos        = length(JusVal);                                    %Length of the justification symbol
        U_pi2         = 3.8;                                               %Biar voltage of the selected device
        %There are three PAM4 techniques implemented in this simulation,
        %Electical PAM4, Optical PAM4 with IQ-MZM and Optical PAM4 with
        %DD-MZM. The last one is the modulation with the best result for
        %the current setup configuration so far. It is selected when
        %Which4PAM set 1 and ModSchem set 0.
        if Which4PAM
            ModSchem      = PAM4Type(2);
            BWD           = 1.0*fc;                                        %Badnwidth for the bit conformation
            CenFeq        = 0;                                             %Center frequency for the bit conformation
            FiltOrd       = 5;                                             %The order for the filter which will conform data
            Vbias         = 1.9;                                           %Biar voltage of the selected device
            if ModSchem
                MZ_Input_File = 2;                                         %Telling the script witch file will be used for modulation
                Vmax          = 1;
                Vmin          = -1;
            else
                MZ_Input_File = 2;                                         %Telling the script witch file will be used for modulation
                Vmax          = 1.25;
                Vmin          = 0;
            end
        else
            BWD           = 1.0*fc;                                        %Badnwidth for the bit conformation
            CenFeq        = 0;                                             %Center frequency for the bit conformation
            FiltOrd       = 1;                                             %The order for the filter which will conform data
            VPI           = 3.8;                                           %Characteristic voltage of the MZM - Vbias at VPI is the point of minimum
            Vbias         = 1.9;                                           %Biar voltage of the selected device
            MZ_Input_File = 2;                                             %Telling the script witch file will be used for modulation
        end
        %Variables for the Selection of carriers - Wave Divided.
        fin           = RefCarr*fc;                                        %The frequency of the first carrier of the optical COMB
        FBWD          = 1.0*fc;                                            %The filter BandWidth to split each individual carrier
        Order         = 5;                                                 %The order of the filter used to split each carrier
    otherwise
        %% Parameters for the DownStream Recptor
        SyncSymb      = [zeros(1,10) ones(1,4) zeros(1,10)];               %Generating the sync-word at the beginning of the data stream
        SyncPos       = 16*(2*NumAmosCP+NPPB);                             %Seting the possition where the sync-symble is
        IniSyncPos    = 8*(2*NumAmosCP+NPPB);
        SyncSymbEnd   = [zeros(1,10) ones(1,4) zeros(1,10)];               %Generating the sync-word at the end of the data stream
        SyncPeriod    = length(SyncSymb);                                  %Seting how many symble periods will be within the sync-word
        AuxSyncLength = SyncPeriod*(2*NumAmosCP+NPPB);                     %auxiliar variable to set where is the sync-symble
        tsync         = linspace(0,SyncPeriod*T,(2*NumAmosCP+NPPB)*...
            SyncPeriod);%create a time vetor for the sync-word to help on the ploting of this signal
        if NRZPolarity
            SyncSymb(SyncSymb==0) = -1;
            SyncSymbEnd(SyncSymbEnd==0) = -1;
        end
        fsync         = time2freq(tsync);                                  %creating the the frequency vector of the sync-word for hereafter filtering
        BWD           = 1*fc;                                        %Badnwidth for the bit conformation
        CenFeq        = 0;                                                 %Center frequency for the bit conformation
        FiltOrd       = 1;                                                 %The order for the filter which will conform data
        BitFilBanWid  = 1*fc;                                              %Badnwidth for the bit conformation
        BitFiltCenFre = 0;                                                 %Center frequency for the bit conformation
        BitFilOrd     = 1;                                                 %The order for the filter which will conform data
        [BitFilt,~]   = FiltroGaussiano(fsync,BitFilBanWid,BitFiltCenFre...
            ,BitFilOrd);%Creating filter for conformation of the input information
        BitFilt       = fftshift(BitFilt);
        SyncSymb      = rectpulse(SyncSymb,2*NumAmosCP+NPPB);              %Resampling the sync-word accordingly with the current system
        SyncSymbEnd   = rectpulse(SyncSymbEnd,2*NumAmosCP+NPPB);           %Resampling the sync-word accordingly with the current system
        
        %Because of reasons, the real world does not have components
        %capable instantly change the voltage level. Thus, to emulate this
        %behaviour the eletrical PAM signal pass through a gaussian filter
        %that will remove the frequencies of higher order, which will
        %result in small slope in the level variation
        SyncSymb      = ifft(fft(SyncSymb).*BitFilt);
        SyncSymbEnd   = ifft(fft(SyncSymbEnd).*BitFilt);
        SencondAdjust = OFF;                                               %Setting whether is needed to sync the information by the Syncronism Symble at the End
        LowSymPer     = 1;
        UpeSymPer     = 1;
        %% Parameters for the DownStream Transmiter
        %%      Variables for the OOK Datum Modification
        AdjusData     = 0;                                                 %Selecting the nature of the information to be send
        JusVal        = [zeros(1,10) ones(1,4) zeros(1,10)];               %Value that will be load into the justificaiton space
        JusValEnd     = [zeros(1,10) ones(1,4) zeros(1,10)];               %Value that will be load into the justificaiton space
        JusLen        = length(JusVal);                                    %Length for the information justification
        DatMin        = 0;                                                 %Minimum value that the Data can assume
        DatMax        = 1;                                                 %Maximum value that the DAta can assume
        NrzMin        = 0.96;%3.9;                                         %Minimum value for the NRZ bit
        NrzMax        = -0.96;%3.8;                                        %Maximum value for the NRZ bit
        IQAmp         = 1;                                                 %Amplitude of the signal IQ
        DatGai        = 1;                                                 %Gain to increase the modulatin signal amplitude
        MZ_Input_File = 5;                                                 %Telling the script witch file will be used for modulation
        %Variables for the Selection of carriers - Wave Divided.
        fin           = RefCarr*fc;                                        %The frequency of the first carrier of the optical COMB
        FBWD          = 1.0*fc;                                            %The filter BandWidth to split each individual carrier
        Order         = 5;                                                 %The order of the filter used to split each carrier
        BWD2          = 1*fc;
end
%%      Filter to Select the Carriers
[ SelecFilt ] = FiltroGaussiano(f,(NumCarr)*fc,((NumCarr+1)/2)*fc,11);
% [ SelecFilt ] = Filtro_Retangular( 34*fc,16*fc,f);
SelecFilt = fftshift(SelecFilt);
%%      Filter for the EVM measurement taking
EvmFilt = FiltroGaussiano(f,fc,0,1);
EvmFilt = fftshift(EvmFilt);

VetSnrIni     = 20;
VetSnrPass    = -10;
VerSnrEnd     = 20;
ThisPlotCont  = 1;
