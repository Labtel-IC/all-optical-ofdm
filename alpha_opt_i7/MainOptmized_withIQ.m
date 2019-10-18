%c
%c                                                       ..'Ž`'..'Ž`..'Ž`..
%c       File: Main File For the All Optical OFDM
%c(Main file to call each idividual parameters)
%c
%c     This main code is resposible to call and run the all the fuction to
%c related to this simulation, which is resposible to run a test for a
%c modulation type. For instance, simulate how a optical sinal modulated
%c with a PAM4 can be transmited in a PON system over km of SSMF.
%c Evaluating parameter such as eye opening, bit error ratio, etc. No
%c modification should be done here just to get the simulation running.
%c Changin in this file just for addint new features.
%c
%c
%c
%c                                           by P.Marciano LG
%c                                           28/10/2017
%c                                           Last UpDate
%c                                           09/05/2018
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This code was optimised at 15 of September of 2018. The principles here 
%followed for optimisation were retrieved from MatLab recommendations. 
%Those where: 

%1 - Except for function, place all script on the same file. The reason for 
%that is to make MatLab to “see” the whole script hence its optimization can 
%properly work. Otherwise, if part of the code is within other scrip MatLab 
%is not able to suggest a correction for that part of the code. Matlab sees 
%different scripts as different codes. 

%2 - Correct all warning given by MatLab. If one intends to create “mex” 
%function or using parallel computing, warning on code can be errors on 
%those approaches. 

%3 - Initialization of variables is needed to be done to save time by 
%avoiding resizing variables dynamically. MatLab allocates variables in the 
%continuous section of memory. If a variable change size over loop 
%interaction the code performance can be slow down as MatLab are looking 
%for free space to place this variable. 

%4 - Vectorization operation where it is possible. Matlab is optimized to 
%work with matrices thus, for instance, if it is possible to replace a “for” 
%loop by a matrix this change speeds up the coding process. 

%5 - Look at the function on the scrip, for instance, with the MatLab 
%profile function. Thus, optimizing the functions or script sections which 
%drag down the code speed. 

%6 - Transfer part (or all) processing to GPUs especially if the function 
%used were already on the list of MatLab function optimized to work on GPUs

%7 - Implement the slowest function or script section on “C” language.

%8 - Implement parallel processing where it is possible.
%%
close all;clc;clear;tic;a=1;
ThisPath = pwd;%Path for the folter where the OFC is saved
UofcsLocal = '\..\save_files\';%Files name for the location of the OFC
%Checking if this main file was not called by another script
%This was not been used, this code works when different OFCS need to be 
%tested and it was called by another script
if exist('CurrentOCS','var')
    OfcName = [ThisPath UofcsLocal 'OCS_' num2str(OcsToTest(CurrentOCS))];
else
    %OfcName = [ThisPath UofcsLocal 'OCS_fc_125_Car_2_08-Jul-2018'];        %OFCS with 2 carriers NPPB 2^3 Nb 2^16 fc 12.5e9 EDFA Power 20 dB
    %OfcName = [ThisPath UofcsLocal 'OCS_fc_125_Car_4_08-Jul-2018'];        %OFCS with 4 carriers NPPB 2^4 Nb 2^15 fc 12.5e9 EDFA Power 20 dB
    %OfcName = [ThisPath UofcsLocal 'OCS_fc_125_Car_8_08-Jul-2018'];        %OFCS with 8 carriers NPPB 2^5 Nb 2^14 fc 12.5e9 EDFA Power 20 dB
    %OfcName = [ThisPath UofcsLocal 'OCS_fc_125_Car_16_08-Jul-2018'];       %OFCS with 16 carriers NPPB 2^6 Nb 2^13 fc 12.5e9 EDFA Power 20 dB
    %OfcName = [ThisPath UofcsLocal 'OCS_fc_125_Car_32_08-Jul-2018'];       %OFCS with 32 carriers NPPB 2^7 Nb 2^12 fc 12.5e9 EDFA Power 20 dB
    %OfcName = [ThisPath UofcsLocal 'OCS_fc_125_Car_64_08-Jul-2018'];       %OFCS with 64 carriers NPPB 2^8 Nb 2^11 fc 12.5e9 EDFA Power 20 dB
    OfcName = [ThisPath UofcsLocal 'OCS_fc_125_Car_128_08-Jul-2018'];      %OFCS with 128 carriers NPPB 2^9 Nb 2^10 fc 12.5e9 EDFA Power 20 dB
    %OfcName = [ThisPath UofcsLocal 'OCS_fc_125_Car_256_09-Jul-2018'];      %OFCS with 256 carriers NPPB 2^10 Nb 2^9 fc 12.5e9 EDFA Power 20 dB
    %OfcName = [ThisPath UofcsLocal 'OCS_fc_125_16-May-2018'];              %OFCS with 128 carriers NPPB 2^9 Nb 2^10 fc 12.5e9 EDFA Power 20 db
    %OfcName = [ThisPath UofcsLocal 'OCS_fc_125_Car_128_Long_12-Jul-2018']; %OFCS with 128 carriers NPPB 2^9 Nb 2^12 fc 12.5e9 EDFA Power 20 db
    %OfcName = [ThisPath UofcsLocal 'OCS_fc_125_Car_128_Long_13-Jul-2018']; %OFCS with 128 carriers NPPB 2^9 Nb 2^11 fc 12.5e9 EDFA Power 20 db
end

%This code assumes that the OFC source was already generated. Thue the time
%and frequency vectors are loaded together with the OFC source.
load([OfcName '.mat']);%Loading a OFC from memory.

%Depending on the computer configurations some optimisation steps are not 
%possible. The UsingMex variable enables or disables the usage of C code 
%files. The CurTesPas control the vectorisation of the code. This variable 
%must be a divisor of the total number of simulation loop tests. The 
%CurTesSize marks the number of extra collum which means the interactions 
%calculated on a single interaction. The UsingGpu is just enabled when a 
%graphic video is present.

%Attention! The larger is the CurTesPas more RAM will be used. Thus, this 
%variable must be adjusted accordingly with the current PC configurations.

UsingMex      = 1;
CurTesPas     = 1;
CurTesSiz     = length(1:CurTesPas);
UseGpu        = 0;
if ~gpuDeviceCount()
    UsingGpu = 0;
else
    UsingGpu = UseGpu;
end
f = reshape(f,length(f),1);%Frequency Vector
if UsingGpu == 1
    fgpu = gpuArray(f);
else
    fgpu = 0;
end
f = repmat(f,1,CurTesSiz);
t = repmat(t.',1,CurTesSiz);%Time Vector
%%
%The MZM needs to have a frequency vector adjusted to work with the 
%wavelength used. Keeping in mind that the light here used is around 
%hundreds of terahertz. Thus, this frequency vector depends only on the 
%simulation original frequency vector. Thus, as the MZM is called several 
%times in this script, this frequency vector does not change. Thus, this 
%section was removed from the MZM function and calculated outside on the 
%script. When the MZM is called the new frequency vector is sent to the MZM 
%function as an input parameter.
if ~exist('freqGHz','var')
    freqGHz = zeros(size(t,1),size(t,2));
    for kk=1:CurTesSiz
        tps           = t(:,1)/1E-12;
        freqTHzA      = time2freq_lamb_2(tps);
        freqGHzA      = freqTHzA*1e-3;% Frequencia em GHz
        freqGHzA      = -fftshift(freqGHzA);
        freqGHzA(1)   = freqGHzA(2);
        freqGHz(:,kk) = freqGHzA;
    end
end
% fgpu = gpuArray(f);
%%                 Control Variable

%       General Parameter for the simulation
ON            =  1;                                                        %DO NOT CHANGE THIS PARAMETER
OFF           =  0;                                                        %DO NOT CHANGE THIS PARAMETER
Ploting       = ON;                                                       %Use On or Off to say if the figures should be plot or not
PlotingThis   = ON;                                                       %Use On or Off to say if the figures should be plot or not
PrintinEye    = ON;                                                       %Printing Eye Diagram of modulations
PlotingFD     = OFF;
%The main results were generated here using a simple algorithm that I
%develop, which seek the biggest opening on a given interval. The interval
%is found by using the eye diagram histogram. This variable will be used to
%verify whether the method is more efficient than the original. If so, it
%can be one explanation of the gain reduction needed at the receiver site.
%As the SNR value founded here are lower than what was expected.
CalcS         = OFF;                                                        %Set 1 to enable the first eye openeing stimative
%There are many points where noise can be added. It was not defined yet the
%right place for it. But four points are being pointed as candidates. At
%the electrical signal that modulates the carrier #1, after the signal
%passes through the fiber #2, after the carrier was modulated #3 and just
%before the photodetector #4.
AddingNoiseE  = OFF;                                                       %#1 Set 1 to add noise to the electrical signal
AddingNoiseF  = OFF;                                                       %#2 Set 1 to add noise to the optical signal
AddingNoiseO  = OFF;                                                       %#3 Set 1 to add noise to the optical signal
AddingNoiseP  = OFF;                                                       %#4 Set 1 to add noise to the optical signal
ReceptorNoise = OFF;
%Two other noises are not concidered here. One is just after the OFCS was
%loaded. The second is just after the photodetector.
SetCpSampZer  = OFF;                                                       %Set 1 to use the samples on the Cp as zeros
Selecting     = ON;                                                        %Use to filter or not the carriers to be transmited
ModAll        = OFF;                                                       %Set if all carriers will be modulate at once or not
AddCP         = ON;                                                        %Add the Cycle Prefix on the data stream
SelecSetUP    = ON;                                                        %Select the simulated setup 1 for original 0 for IFFT/OFFT at receptor
SendingDowStr = ON;                                                        %Set 1 to send data downstream 0 to not send.
SigAten       = OFF;                                                       %Set 1 to atenuate the OFCS signal. Set to zero to turn off attenuation
SnrRef        = OFF;                                                       %Set 1 to select the own signal for the noise reference at receptor.
%Set 0 for the pilote carrier as noise reference at receptor.
FFTSplit      = ON;                                                        %Set 1 to use the OFFT to split the income signal on receptor. Set 0 to
%use filter to split the income signal at receptor.
RecFilBanPas  = ON;                                                        %Filter after photodetector, 1 filter Before noise, 0 filter After noise.
SendingUp     = ON;                                                        %Select if the carrier upstream will be send. Set 1 to send, set 0 to not send.
SaveToScatter = OFF;                                                       %Set 1 to save points to create scatter plot latter Set 0 to not save points (keep it zero to save data space)
TimeSys       = ON;                                                        %Select which time period will be accounted to the noise, set 1 to use T or set 0 to use t(2*NumAmosCP+NPPB).
%When CP is used Seting TimeSys to 1 gives higher SNR
PastTime      = OFF;                                                       %Storage for the time flow on simulation
TestBundle    = ON;                                                        %Set if the scrip will be run once or many times
CurrentMedium = ON;                                                        %Sellection of channel 1 for Fiber other of B2B
UseOfdmElect  = OFF;
ChanAwgn      = OFF;

Nsamp         = 2^9;%1;%                                                   %Number of samples per bit, it was necessary because the number of carriers here used
% Atenuation    = 10^(4.0);                                                %Set how much the OFCS will be attenuated
NumbOfTest    = 2/CurTesPas;                                              %The number of repetition that one test will be done
PAM4Type      = [1 0];                                                     %Sellection the type of 4PAM that will be tested
MinTesNumb    = 10;
MinSavTim     = 15;

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
UpStModula = UsedModula;
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
        NumCarr  = 4;                                                    %Total number of carriers (DownStream and UpStream) - This variable must be even
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
        %CarrUsedDo(1) = 1;                                                 %Setting the downstream carrier as valid to be used
        %CarrUsedUp(2) = 1;                                                 %Setting the upstream carrier as valid to be used
        InitCarrUp = 2;                                                    %Setting that the frist carrier of the upstream signal is at the 2 possition
        InitCarrDo = 1;                                                    %Setting that the frist carrier of the downstream signal is at the 1 possition
    else                                                                   %Setting parrameter when carriers were placed sequentialy
        CarrUsed   = ones(1,NumCarr);                                      %Setting flags on the carriers that will be used 1 for used 0 for unused
        %CarrUsed   = zeros(1,NumCarr);                                     %Setting flags on the carriers that will be used 1 for used 0 for unused
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
            NumCarr)];                                                     %As in this simulation all carrier will be used the observed carriers is equal
        %to the possition of the observed carriers
    else                                                                   %whereas, when the first carrier was not multiple of the number of carrier used
        ObsCarrPos  = [mod(RefCarr,NumCarr):NumCarr 1:mod(RefCarr,...
            NumCarr)-1];%The first carrier is what left from the divission Reference carrier by
        %the number of carriers used and this sequency goes untill the last carrier
        %selected. Then, the following carrier will be from 1 untill the divission result minus 1
        ObsCarrUsed = ObsCarrPos;                                          %As in this simulation all carrier will be used the observed carriers is equal
        %to the possition of the observed carriers
    end
%         SaveBerEbn0PreFix   = [pwd '\..\results\UsedModula_'];
%         SaveBerEbn0SuFix    = '_BerVsEbnONoOverSamp_1';
%         SaveBerEbn0FigSuFix = '_BerVsEbnONoOverSamp_fig_1.fig';
    SaveBerEbn0PreFix   = [pwd 'test0'];
    SaveBerEbn0SuFix    = '';
    SaveBerEbn0FigSuFix = '';
    
%     SavingAT     = 'test_0';
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
        FreqCen   = 0;
        U0        = 1.9;                                                   %MZM Vbias
        BWD       = 1.0*fc;                                                %Badnwidth for the bit conformation
        CenFeq    = 0;                                                     %Center frequency for the bit conformation
        FiltOrd   = 0;                                                     %The order for the filter which will conform data
        %% Parameters for the DownStream Transmiter
        %It is important to mention that when SelModFilt is set to one it
        %causes attenuation on the OFDM signal higher frequencies. Thus
        %distorting the constellation diagram.
        FcIQ      = fc;
        UseIQ     = 1;
        SelModFilt= 0;                                                     %Set 1 to use filter to selec the OFDM signal from its allisings componentes
        SelecGaus = 1;                                                     %Select whether a gaussian or retangula filter will be used
        OfdMod    = 'qam';                                                 %Chossing the type of modulation for the electrical subcarriers
        SelModTp  = 'BP';                                                  %Selecting with OFDM will be in Base Band "BP", Amplitude Modulation
        %with double side band "AM" or with single side band "AMSSB" when a
        %different option is choosen this script will implement phase 
        %modulation.
        UsingHermitian = 0;                                                %Turn On or Off the usage of Hermitian Symetri
        if (~UseOfdmElect)&&(~UseIQ)                                                   %Whenever a non-electrical OFDM is used the Hermitian Symetri in needed on this script
            UsingHermitian = 1;
        end
        %It is not totally clear how the OFDM signal oversampling should 
        %be. The “ForcOvSam” variable force the OFDM signal oversampling to 
        %be exactly its value when it is different of zero. When this 
        %variable is set to zero the algorithm calculate the oversampling 
        %based on the signal bandwidth.
        ForcOvSam = 0;
        UsingParCar= 0;                                                    %This variable was create to select when OFDM agregation would be used or not
        %The main point of this script was to use an adaptative modulation 
        %or the discrete multitone (DMT) modulation. That is basically the 
        %OFDM signal but the electrical carriers have different modulation 
        %formats to adjust with the channel. For instance:
        %
        %    OFDM with NuDmSuDi = 1          OFDM with NuDmSuDi = 4
        %             M=4                    M=4
        %    ____________________           ______  M=3
        %    |                  |           |    |______  M=2
        %    |                  |           |          |______  M=1
        %    |                  |           |                |______
        %    |                  |           |                      |      
        NuDmSuDi  = 1;                                                     
        DmtPas    = 2;                                                     %Pass among OFDM division
        DmtMinM   = 2;                                                     %Minimum modulation level of the DMT signal
        TapN      = 1;                                                     %Number of equalizer used
        ExtSam    = 1;                                                     %Set to the number of frames that will be used.
        NpEl      = 2^9;                                                   %Number of Electrical carriers
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
        CpLe   = 1;                                                        %Set the percentage of electrical carriers that will be used
        if CpLe <= 0.5
            error(['CpLe can not be bellow 50%. More than half of ' ...
                'available carriers must be transmited']);
        end
        BW     = fc;                                                       %Signal available bandwidth
        %         NumFra = floor(BW*t(end));                               %Number of frames
        %         if NFFT>Nb
        %             NFFT = Nb;
        %         end
        DatSiz = NFFT;                                                     %Set to the total number of electrical carriers to be used
        Te     = OveT*T;                                                   %The period of electrical carrier
        Tu     = Te;                                                       %Time that will actually be used to transmit information
        Tg     = 0.0*Tu;                                                   %Time for the guard band
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
        BWD2          = 1*fc;%Badnwidth for the bit conformation
        CenFeq2       = 0;%Center frequency for the bit conformation
        FiltOrd2      = 1;%The order for the filter which will conform data
        SyncSymb      = [zeros(5,CurTesSiz);ones(2,CurTesSiz);zeros(5,CurTesSiz)];%Generating the sync-word at the beginning of the data stream
        IniSyncPos    = 4*(2*NumAmosCP+NPPB);
        SyncPos       = 8*(2*NumAmosCP+NPPB);%Seting the possition where the sync-symble is
        SyncSymbEnd   = SyncSymb;%Generating the sync-word at the end of the data stream
        SyncPeriod    = size(SyncSymb,1);%Seting how many symble periods will be within the sync-word
        AuxSyncLength = SyncPeriod*(2*NumAmosCP+NPPB);%auxiliar variable to set where is the sync-symble
        tsync         = linspace(0,SyncPeriod*T,(2*NumAmosCP+NPPB)*SyncPeriod);%create a time vetor for the sync-word to help on the ploting of this signal
        fsync         = time2freq(tsync);%creating the the frequency vector of the sync-word for hereafter filtering
        BWD           = 2*fc;%Badnwidth for the bit conformation
        CenFeq        = 0;%Center frequency for the bit conformation
        FiltOrd       = 1;%The order for the filter which will conform data
        BitFilBanWid  = 1*fc;%Badnwidth for the bit conformation
        BitFiltCenFre = 0;%Center frequency for the bit conformation
        BitFilOrd     = 1;%The order for the filter which will conform data
        [BitFilt,~]   = FiltroGaussiano(fsync,BWD,CenFeq,FiltOrd);%Creating filter for conformation of the input information
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
        SendingData   = 1;%Selecting the nature of the information to be send
        NbDPSK        = Nb - NumBitDesc;%The total number of bits to be transmited
        DatGai        = 1;%Seting the maximum range for the eletrical signal
        JusVal        = [zeros(5,CurTesSiz);ones(2,CurTesSiz);zeros(5,CurTesSiz)];%Value that will be load into the justificaiton space
        JusValEnd     = JusVal;%Value that will be load into the justificaiton space
        JusPos        = size(JusVal,1);%Length of the justification symbol
        BWD           = 2*fc;%Badnwidth for the bit conformation
        CenFeq        = 0;%Center frequency for the bit conformation
        FiltOrd       = 1;%The order for the filter which will conform data
        MZ_Input_File = 3;%Telling the script witch file will be used for modulation
        %Variables for the Selection of carriers - Wave Divided.
        fin           = RefCarr*fc;%The frequency of the first carrier of the optical COMB
        FBWD          = 1.0*fc;%The filter BandWidth to split each individual carrier
        Order         = 5;%The order of the filter used to split each carrier
    case 'DQPSK'
        %% Parameters for the DownStream Recptor
        BWD2          = 1.0*fc;%Badnwidth for the bit conformation
        CenFeq2       = 0;%Center frequency for the bit conformation
        FiltOrd2      = 1;%The order for the filter which will conform data
        SyncSymb      = [zeros(5,CurTesSiz);ones(2,CurTesSiz);zeros(5,CurTesSiz)];%Generating the sync-word at the beginning of the data stream
        IniSyncPos    = 4*(2*NumAmosCP+NPPB);
        SyncPos       = 8*(2*NumAmosCP+NPPB);%Seting the possition where the sync-symble is
        SyncSymbEnd   = SyncSymb;%Generating the sync-word at the end of the data stream
        SyncPeriod    = size(SyncSymb,1);%Seting how many symble periods will be within the sync-word
        AuxSyncLength = SyncPeriod*(2*NumAmosCP+NPPB);%auxiliar variable to set where is the sync-symble
        tsync         = linspace(0,SyncPeriod*T,(2*NumAmosCP+NPPB)*SyncPeriod);%create a time vetor for the sync-word to help on the ploting of this signal
        fsync         = time2freq(tsync);%creating the the frequency vector of the sync-word for hereafter filtering
        BWD           = 2*fc;%Badnwidth for the bit conformation
        CenFeq        = 0;%Center frequency for the bit conformation
        FiltOrd       = 1;%The order for the filter which will conform data
        BitFilBanWid  = 1*fc;%Badnwidth for the bit conformation
        BitFiltCenFre = 0;%Center frequency for the bit conformation
        BitFilOrd     = 1;%The order for the filter which will conform data
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
        SendingData = 1;%Selecting the nature of the information to be send
        NbDQPSK     = 2*(Nb - NumBitDesc);%The total number of bits to be transmited
        JusVal      = [zeros(10,CurTesSiz);ones(4,CurTesSiz);zeros(10,CurTesSiz)];%Value that will be load into the justificaiton space
        JusValEnd   = JusVal;%Value that will be load into the justificaiton space
        JusPos      = size(JusVal,1);%Length of the justification symbol
        BWD         = 4*fc;%Badnwidth for the bit conformation
        CenFeq      = 0;%Center frequency for the bit conformation
        FiltOrd     = 1;%The order for the filter which will conform data
        DatGai      = 1;%Gain to increase the modulatin signal amplitude
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
        SyncSymb      = [zeros(5,CurTesSiz);ones(2,CurTesSiz);zeros(5,CurTesSiz)];%Generating the sync-word at the beginning of the data stream
        SyncSymbSiz   = size(SyncSymb,1)*2;
        IniSyncPos    = 4*(2*NumAmosCP+NPPB);
        SyncPos       = 8*(2*NumAmosCP+NPPB);%Seting the possition where the sync-symble is
        SyncSymbEnd   = SyncSymb;%Generating the sync-word at the end of the data stream
        SyncPeriod    = size(SyncSymb,1);%Seting how many symble periods will be within the sync-word
        AuxSyncLength = SyncPeriod*(2*NumAmosCP+NPPB);%auxiliar variable to set where is the sync-symble
        tsync         = linspace(0,SyncPeriod*T,(2*NumAmosCP+NPPB)*SyncPeriod);%create a time vetor for the sync-word to help on the ploting of this signal
        fsync         = time2freq(tsync);%creating the the frequency vector of the sync-word for hereafter filtering
        BWD           = 2*fc;%Badnwidth for the bit conformation
        CenFeq        = 0;%Center frequency for the bit conformation
        FiltOrd       = 1;%The order for the filter which will conform data
        if PAM4Type(2)
            BWD2          = 1*fc;%Badnwidth for the bit conformation
            CenFeq2       = 0;%Center frequency for the bit conformation
            FiltOrd2      = 1;%The order for the filter which will conform data
        else
            BWD2          = 1*fc;%Badnwidth for the bit conformation
            CenFeq2       = 0;%Center frequency for the bit conformation
            FiltOrd2      = 1;%The order for the filter which will conform data
        end
        [BitFilt,~]   = FiltroGaussiano(fsync,BWD,CenFeq,FiltOrd);%Creating filter for conformation of the input information
        BitFilt       = fftshift(BitFilt);
        SyncSymb      = rectpulse(SyncSymb,2*NumAmosCP+NPPB);%Resampling the sync-word accordingly with the current system
        SyncSymbEnd   = rectpulse(SyncSymbEnd,2*NumAmosCP+NPPB);%Resampling the sync-word accordingly with the current system
        IntervalStep  = 50;
        MinDist       = IntervalStep/6;
        %Because of reasons, the real world does not have components
        %capable instantly change the voltage level. Thus, to emulate this
        %behaviour the eletrical PAM signal pass through a gaussian filter
        %that will remove the frequencies of higher order, which will
        %result in small slope in the level variation
        BitFilt       = repmat(BitFilt.',1,CurTesSiz);
        SyncSymb      = ifft(fft(SyncSymb).*BitFilt);
        SyncSymbEnd   = ifft(fft(SyncSymbEnd).*BitFilt);
        SencondAdjust = OFF;%Setting whether is needed to sync the information by the Syncronism Symble at the End
        LowSymPer     = 0.9;
        UpeSymPer     = 1.1;
        %% Parameters for the DownStream Transmiter
        %%      Variables for the 4PAM Datum Modification
        %This variable allow the user to chose to send an random
        %information or just a high pulse within the midle of the datastram
        SendingData   = 1;%Selecting the nature of the information to be send
        Nb4Pam        = 2*(Nb - NumBitDesc);%The total number of bits to be transmited
        %This variable controls if the electrical signal will vary from
        %0[V] to max[V] (signal unipolar) or from -max[V] to max[V]
        %(signal bipolar).
        Polirized     = 1;
        MaxAmp4PAM    = -3.8/4;%Seting the maximum range for the eletrical signal
        JusVal        = [zeros(10,CurTesSiz);repmat([1;0],2,CurTesSiz);zeros(10,CurTesSiz)];%Value that will be load into the justificaiton space
        JusValEnd     = JusVal;%Value that will be load into the justificaiton space
        JusPos        = size(JusVal,1);%Length of the justification symbol
        U_pi2         = 3.8;%Biar voltage of the selected device
        %There are three PAM4 techniques implemented in this simulation,
        %Electical PAM4, Optical PAM4 with IQ-MZM and Optical PAM4 with
        %DD-MZM. The last one is the modulation with the best result for
        %the current setup configuration so far. It is selected when
        %Which4PAM set 1 and ModSchem set 0.
        if Which4PAM
            ModSchem      = PAM4Type(2);
            BWD           = 1.0*fc;%Badnwidth for the bit conformation
            CenFeq        = 0;%Center frequency for the bit conformation
            FiltOrd       = 5;%The order for the filter which will conform data
            Vbias         = 1.9;%Biar voltage of the selected device
            if ModSchem
                MZ_Input_File = 2;%Telling the script witch file will be used for modulation
                Vmax          = 1;
                Vmin          = -1;
            else
                MZ_Input_File = 2;%Telling the script witch file will be used for modulation
                Vmax          = 1.25;
                Vmin          = 0;
            end
        else
            BWD           = 1.0*fc;%Badnwidth for the bit conformation
            CenFeq        = 0;%Center frequency for the bit conformation
            FiltOrd       = 1;%The order for the filter which will conform data
            VPI           = 3.8;%Characteristic voltage of the MZM - Vbias at VPI is the point of minimum
            Vbias         = 1.9;%Biar voltage of the selected device
            MZ_Input_File = 2;%Telling the script witch file will be used for modulation
        end
        %Variables for the Selection of carriers - Wave Divided.
        fin           = RefCarr*fc;%The frequency of the first carrier of the optical COMB
        FBWD          = 1.0*fc;%The filter BandWidth to split each individual carrier
        Order         = 5;%The order of the filter used to split each carrier
    otherwise
        %% Parameters for the DownStream Recptor
        SyncSymb      = [zeros(5,CurTesSiz);ones(2,CurTesSiz);zeros(5,CurTesSiz)];%Generating the sync-word at the beginning of the data stream
        SyncPos       = 8*(2*NumAmosCP+NPPB);%Seting the possition where the sync-symble is
        IniSyncPos    = 4*(2*NumAmosCP+NPPB);
        SyncSymbEnd   = [zeros(5,CurTesSiz);ones(2,CurTesSiz);zeros(5,CurTesSiz)];%Generating the sync-word at the end of the data stream
        SyncPeriod    = size(SyncSymb,1);%Seting how many symble periods will be within the sync-word
        AuxSyncLength = SyncPeriod*(2*NumAmosCP+NPPB);%auxiliar variable to set where is the sync-symble
        tsync         = linspace(0,SyncPeriod*T,(2*NumAmosCP+NPPB)*SyncPeriod);%create a time vetor for the sync-word to help on the ploting of this signal
        if NRZPolarity
            SyncSymb(SyncSymb==0) = -1;
            SyncSymbEnd(SyncSymbEnd==0) = -1;
        end
        fsync         = time2freq(tsync);%creating the the frequency vector of the sync-word for hereafter filtering
        BWD           = 1*fc;%Badnwidth for the bit conformation
        CenFeq        = 0;%Center frequency for the bit conformation
        FiltOrd       = 1;%The order for the filter which will conform data
        BitFilBanWid  = 1*fc;%Badnwidth for the bit conformation
        BitFiltCenFre = 0;%Center frequency for the bit conformation
        BitFilOrd     = 1;%The order for the filter which will conform data
        [BitFilt,~]   = FiltroGaussiano(fsync,BitFilBanWid,BitFiltCenFre,BitFilOrd);%Creating filter for conformation of the input information
        BitFilt       = fftshift(BitFilt);
        BitFilt       = repmat(BitFilt.',1,CurTesSiz);
        SyncSymb      = rectpulse(SyncSymb,2*NumAmosCP+NPPB);%Resampling the sync-word accordingly with the current system
        SyncSymbEnd   = rectpulse(SyncSymbEnd,2*NumAmosCP+NPPB);%Resampling the sync-word accordingly with the current system
        
        %Because of reasons, the real world does not have components
        %capable instantly change the voltage level. Thus, to emulate this
        %behaviour the eletrical PAM signal pass through a gaussian filter
        %that will remove the frequencies of higher order, which will
        %result in small slope in the level variation
        SyncSymb      = ifft(fft(SyncSymb).*BitFilt);
        SyncSymbEnd   = ifft(fft(SyncSymbEnd).*BitFilt);
        SencondAdjust = OFF;%Setting whether is needed to sync the information by the Syncronism Symble at the End
        LowSymPer     = 1;
        UpeSymPer     = 1;
        %% Parameters for the DownStream Transmiter
        %%      Variables for the OOK Datum Modification
        AdjusData     = 0;%Selecting the nature of the information to be send
        JusVal        = [zeros(5,CurTesSiz);ones(2,CurTesSiz);zeros(5,CurTesSiz)];%Value that will be load into the justificaiton space
        JusValEnd     = JusVal;%Value that will be load into the justificaiton space
        JusLen        = size(JusVal,1);%Length for the information justification
        DatMin        = 0;%Minimum value that the Data can assume
        DatMax        = 1;%Maximum value that the DAta can assume
        NrzMin        = 0.96;%3.9;%Minimum value for the NRZ bit
        NrzMax        = -0.96;%3.8;%Maximum value for the NRZ bit
        IQAmp         = 1;%Amplitude of the signal IQ
        DatGai        = 1;%Gain to increase the modulatin signal amplitude
        MZ_Input_File = 5;%Telling the script witch file will be used for modulation
        %Variables for the Selection of carriers - Wave Divided.
        fin           = RefCarr*fc;%The frequency of the first carrier of the optical COMB
        FBWD          = 1.0*fc;%The filter BandWidth to split each individual carrier
        Order         = 5;%The order of the filter used to split each carrier
        BWD2          = 1*fc;
end
%%      Filter to Select the Carriers
[ SelecFilt ] = FiltroGaussiano(f,(NumCarr)*fc,((NumCarr+1)/2)*fc,11);
% [ SelecFilt ] = Filtro_Retangular( 34*fc,16*fc,f);
SelecFilt = fftshift(SelecFilt);
%%      Filter for the EVM measurement taking
EvmFilt = FiltroGaussiano(f,fc,0,1);
EvmFilt = fftshift(EvmFilt);

VetSnrIni     = 20;                                                        %Start of the main Loop
VetSnrPass    = -10;                                                       %Main Loop pass
VerSnrEnd     = 10;                                                        %Stop of the main Loop
ThisPlotCont  = 1;                                                      

%This is the main loop, the VetSnr variable will store the EbNo that will 
%be used to test the system or the fiber length distance that will be 
%applied. It is not possible to control fiber length and EbN0 individually 
%because this script does not have a realistic receptor model. Thus, 
%evaluate the SNR effect over distance is not conclusive.

%clear MzFilesGenerated MzFilesDpGenerated
tic;                                                                       %Starting timer to check the simulation time
for VetSnr=VetSnrIni:VetSnrPass:VerSnrEnd
    clear BerOOK BerPAM4 BerDQPSK BerDPSK BerOOKS BerPAM4S BerDQPSKS ...
        BerDPSKS BerOFDM ChanelEqualizer FiberDelay
%     SavingAT     = [pwd '\..\results\' 'DqpskCp122FulSystKm_' num2str(VetSnr) '_2'];
    SavingAT     = [pwd 'test0'];
    %The initialisation of individual variables to avoid extra  
    %computational timing with a relocation of variables.
    FiberLength  = VetSnr;                                                 %Setting the length of the fiber in this simulation
    CarSNR       = VetSnr;                                                 %Set the snr to be applied at receptor.
    switch Modulation
        case 'OFDM'
            SNR = CarSNR + 10*log10(log2(max(DmtMve))) - 10*log10(OvSam);
            TxSymbAmos = zeros(NumbOfTest,NumCarr,NumFra*length(DmtMve));
        case 'DPSK' %Tested to fiberlength of 80 -> worked
            if TimeSys==1
                SNR = CarSNR + 10*log10(1) - 10*log10(1*T*BWD2);
            else
                SNR = CarSNR + 10*log10(1) - 10*log10(1*t(2*NumAmosCP+NPPB,1)*BWD2);
            end
            TxSymbAmos = zeros(NumbOfTest,NumCarr,Nb - NumBitDesc);
        case 'DQPSK'%Tested to fiberlength of 40 -> work in progress
            if TimeSys==1
                SNR = CarSNR + 10*log10(2) - 10*log10(1*T*BWD2);
            else
                SNR = CarSNR + 10*log10(2) - 10*log10(1*t(2*NumAmosCP+NPPB,1)*BWD2);
            end
            TxSymbAmos = zeros(NumbOfTest,NumCarr,Nb - NumBitDesc);
        case '4PAM'%Tested to fiberlength of 1 -> work in progress
            if TimeSys==1
                SNR = CarSNR + 10*log10(2) - 10*log10(1*T*BWD2);
            else
                SNR = CarSNR + 10*log10(2) - 10*log10(1*t(2*NumAmosCP+NPPB,1)*BWD2);
            end
            TxSymbAmos = zeros(NumbOfTest,NumCarr,Nb - NumBitDesc);
        otherwise %Tested to fiberlength of 80 -> worked
            if TimeSys==1
                SNR = CarSNR + 10*log10(1) - 10*log10(1*T*BWD2);
            else
                SNR = CarSNR + 10*log10(1) - 10*log10(1*t(2*NumAmosCP+NPPB,1)*BWD2);
            end
            TxSymbAmos = zeros(NumbOfTest,NumCarr,Nb - NumBitDesc);
    end
    %The pilot carrier must be close to the center of the OFDM frame for
    %probing the channel once the time delay varies from carrier to carrier.
    %The carries on the right side fo the pilot carrier will have a right shift
    %time delay whereas the carries on the left side of the pilot carrier will
    %have a left shit time delay, where the reference is the synchronization
    %symbol.
    % PilotCarrier = (NumCarr/2)*fc;
    %The main idea about the time delay compensation is that the receptor will
    %have previous knowledge about the length of the transmission link hence it
    %will posse the time delay of this link. Therefore, before any other
    %processing takes place first the income signal will have the time delay
    %adjusted by this information previously stored.
    if CurrentMedium==1
        if ~exist('FiberDelay','var')
            [FiberDelay,PowRef] = ProgationDelayT(RefCarr,NumCarr,t(:,1),lambdac,T,FiberLength,0,f(:,1),PlotingFD);
        end
    end
    
    % The electrical OFDM when modulated with QAM the signal will need a
    % equalizer, because the system implyes on phase rotation. Therefore, this
    % system needs correction before comparing the TX and RX signals. It is
    % also important to mention that the DPSK modulation format don't need
    % equalizer as the information is stored at the phase difference on
    % subsequent symbols phase rotaion will affect adjacent symbol in a very
    % similar way hence information will have very little sencibility to phase
    % rotation.
    % ChanelEqualizerP = FindChanelResponseP(OfdMod,M,NFFT,Ns,NumFra,MZ_Input_File,SelecFilt,Medium,DifNsN,fin,FBWD,Order,ExtSam);
    if exist('Ns','var')
        if UseOfdmElect ==1
            RxSigOfdmNoEq = zeros(NumCarr,Ns*NumFra*CurTesSiz);
            RxSigOfdm = zeros(NumCarr,Ns*NumFra*CurTesSiz);
        else
            RxSigOfdmNoEq = zeros(NumbOfTest,NumCarr,Ns*NumFra*CurTesSiz);
            RxSigOfdm = zeros(NumbOfTest,NumCarr,Ns*NumFra*CurTesSiz);
        end
        SigRecep4 = zeros(length(DmtMve),NumFra*CurTesSiz);
        
        TxSigOfdm = zeros(NumbOfTest,NumCarr,NumFra*CurTesSiz*length(DmtMve));
        BerToPlotOfdm = zeros(NumbOfTest,NumCarr,length(DmtMve)*NumFra*CurTesSiz);
    else
        RxSigOfdmNoEq = 0;
        SigRecep4 = 0;
        RxSigOfdm = 0;
        TxSigOfdm = 0;
        BerToPlotOfdm = 0;
    end
    EvmDB = zeros(NumbOfTest,1);
    EvmDBJ = zeros(NumbOfTest,1);
    EvmPer = zeros(NumbOfTest,1);
    EvmPerJ = zeros(NumbOfTest,1);
    EvmRms = zeros(NumbOfTest,1);
    EvmRmsJ = zeros(NumbOfTest,1);
    BerOFDM = zeros(NumbOfTest,1);
    CarrRecPowDo = zeros(NumbOfTest,NumCarr);
    CarrRecPowUp = zeros(NumbOfTest,NumCarr);
    VetElecPower = zeros(NumbOfTest,NumCarr,CurTesSiz);
    VetOptiPower = zeros(NumbOfTest,NumCarr,CurTesSiz);
    AberLevS = zeros(NumbOfTest,NumCarr);
    ValsLevS = zeros(NumbOfTest,NumCarr);
    EyeToPlot = zeros(NumbOfTest,Nb*NPPB);
    BerDPSK = zeros(NumbOfTest,NumCarr);
    BerDPSKS = zeros(NumbOfTest,NumCarr);
    AberLev = zeros(NumbOfTest,NumCarr);
    ValsLev = zeros(NumbOfTest,NumCarr);
    RxSymbAmos = zeros(NumbOfTest,NumCarr,Nb - NumBitDesc);
    VetElecPowerI = zeros(NumbOfTest,NumCarr,CurTesSiz);
    VetElecPowerQ = zeros(NumbOfTest,NumCarr,CurTesSiz);
    AberLevIS = zeros(NumbOfTest,NumCarr);
    ValsLevIS = zeros(NumbOfTest,NumCarr);
    AberLevQS = zeros(NumbOfTest,NumCarr);
    ValsLevQS = zeros(NumbOfTest,NumCarr);
    BerDQPSK = zeros(NumbOfTest,NumCarr);
    BerDQPSKS = zeros(NumbOfTest,NumCarr);
    AberLevI = zeros(NumbOfTest,NumCarr);
    ValsLevI = zeros(NumbOfTest,NumCarr);
    AberLevQ = zeros(NumbOfTest,NumCarr);
    ValsLevQ = zeros(NumbOfTest,NumCarr);
    ValsLev2 = zeros(3,NumCarr,NumbOfTest);
    Ber4PAM = zeros(NumbOfTest,NumCarr);
    Ber4PAMS = zeros(NumbOfTest,NumCarr);
    BerOOK = zeros(NumbOfTest,NumCarr);
    BerOOKS = zeros(NumbOfTest,NumCarr);
    VetElecPowerUp = zeros(NumbOfTest,NumCarr,CurTesSiz);
    VetOptiPowerUp = zeros(NumbOfTest,NumCarr,CurTesSiz);
    AberLevUp = zeros(NumbOfTest,NumCarr);
    ValsLevUp = zeros(NumbOfTest,NumCarr);
    SaveRxNotEq=[];
    SaveRxEq=[];
    CurrentModula=ThisModula;
    for CurrentTest=1:ThisTest
        %%               Datum Transmition
        % The second part of this scrip will modulate each given sub
        %carrier with its related information. It is important to mention
        %that the actual sub carriers were not selected yet. What is
        %presented so far is the raw OFC. The correct carriers need to be
        %selected or the user may take in account just those carriers that
        %matters for this simulation. At this moment the data will be
        %loaded, the carriers will be selected and the final result will be
        %transmited (Back-to-Back or within Optical Fiber). It is very
        %important to used the correct bit ratio, because it is deeply
        %related with the spacing between carriers. This scrips will assume
        %that the chosen OFC is suitble for the transmition data ratio.
        if (UseOfdmElect)&&(UsedModula==5)&&(NumCarr==2)
            %DatumElectricTransmission;TxSig = 0;
            TxData = zeros(NumFra*CurTesSiz,length(DmtMve));
            TxSymbAux = zeros(NumFra*CurTesSiz,length(DmtMve));
            TxDataMat = zeros(1,CurTesSiz*NumFra*length(DmtMve));
            %U1t = 0;
            UniDmtMve = unique(DmtMve);
            UniDmtMve = fliplr(UniDmtMve);
            DaSiAn = 0;
            DaSiPo = 0;
            for CarDmtM=1:length(UniDmtMve)
                M = UniDmtMve(CarDmtM);
                DaSiAu = sum(ismember(DmtMve,M));
                TxData          = randi([0 M-1],NumFra*CurTesSiz,DaSiAu);                %Generation random information
                TxDataMat(1,1+DaSiAn:DaSiAn + CurTesSiz*DaSiAu*NumFra) = TxData(:);%Saving data for future comparison
                DaSiAn = DaSiAn + CurTesSiz*DaSiAu*NumFra;
                switch OfdMod                                              %Sellect which modulation was used
                    case 'qam'
                        TxSymbAux(1:NumFra*CurTesSiz,1+DaSiPo:DaSiPo + DaSiAu)  = qammod(TxData,M);                    %Modulating information by QAM
                    otherwise
                        TxSymbAux(1:NumFra*CurTesSiz,1+DaSiPo:DaSiPo + DaSiAu)  = dpskmod(TxData,M);                   %Modulating information DPSK as default
                end
                DaSiPo = DaSiPo + DaSiAu;
            end
            TxSymb      = TxSymbAux.';                                %Converting from serial to paralel
            %EVMMatRef(1,:) = TxSymb(:);
            %Construction the Hermitian Simetry
            %The Hermitian simetry was composed of the signal its conjugate format and
            %zeros to fill empty spaces. The signal is formed as show below:
            %                    signal     midle  signal conjugated
            %                       |        |       |
            %                       V        V       V
            %            TxH = [0 a + jb 0 0 0 0 0 a - jb 0];
            %
            %I thought in many ways to form this signal, replacing the centred zeros by
            %copy of the signal bay just geting its information and adding
            %respectively. Also, by creating the signal with redundancy at the exact
            %length needed to make the hermitian singnal. But all tests end up with
            %similar results, no improvement was noticed.
            TxSymbConj      = conj(flipud(TxSymb));                    %Singnal conjugate format at reverse order
            if UsingHermitian
                TxSymbH = [zeros(1,NumFra*CurTesSiz);TxSymb;zeros(1,NumFra*CurTesSiz);TxSymbConj];%Hermitian Symetri
                TxSymbC = zeros(NFFT,NumFra*CurTesSiz);
                TxSymbC(1+Zt-ZtC:Zt-ZtC+size(TxSymbH,1)/2,:) = TxSymbH(1:size(TxSymbH,1)/2,:);
                TxSymbC(end-Zt+ZtC-(size(TxSymbH,1)/2)+1:end-Zt+ZtC,:) = TxSymbH(1+size(TxSymbH,1)/2:end,:);
            else
                TxSymbH = TxSymb;
                %     TxSymbC = [zeros((NFFT-size(TxSymbH,1))/2,NumFra*CurTesSiz);TxSymbH;zeros((NFFT-size(TxSymbH,1))/2,NumFra*CurTesSiz)];
                TxSymbC = [TxSymbH(1:size(TxSymbH,1)/2,:);zeros((NFFT-size(TxSymbH,1)),NumFra*CurTesSiz);TxSymbH(end+1-size(TxSymbH,1)/2:end,:)];
            end
            %Here was implemented the modulation through FFT process, it basica mix-up
            %all infomation following Fourier transform.
            TxSymbMod       = ifft(TxSymbC,NFFT);
            TxSymbMod  = rectpulse(TxSymbMod,OvSam);
            %Sellecting whether the OFDM signal will be transmited on
            %Base band or not
            tt = linspace(0,1*Te,size(TxSymbMod,1));
            ff = time2freq(tt);
            tta = repmat(tt.',1,NumFra*CurTesSiz);
            switch SelModTp
                case 'BP'                                              %Base Band transmission
                    %                     TxSymbModA  = rectpulse(TxSymbModA,OvSam);         %Over sampling
                    %                         TxSymbMod = TxSymbModA;
                    %                     TxSymbMod   = TxSymbMod;
                    if SelecGaus==1
                        [ ModFilt ] = FiltroGaussiano(ff,OBw/2,0,1);
                    else
                        [ ModFilt ] = Filtro_Retangular(OBw,0,ff);
                    end
                    ModFilt = fftshift(ModFilt);
                    if  (Ploting)
                        figure;hold all;
                        plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                    end
                    ModFilta = repmat(ModFilt.',1,NumFra);
                    if SelModFilt==1
                        TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
                    end
                    if  (Ploting)
                        %                             for jj=1:NumFra
                        %                                 ff2(:,jj) = ff;
                        %                             end
                        plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                        plot(ff,20*log10(abs(fftshift(ModFilt))));
                        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                    end
                case 'AM'                                              %Out of base band
                    %                     TxSymbModA  = rectpulse(TxSymbModA,OvSam);         %Over sampling
                    TxSymbMod   = TxSymbMod.*cos(2*pi*Ofc*tta);
                    if SelecGaus==1
                        [ ModFilt ] = FiltroGaussiano(ff,OBw,Ofc,1);
                    else
                        [ ModFilt ] = Filtro_Retangular( OBw,Ofc,ff);
                    end
                    ModFilt = fftshift(ModFilt);
                    if  (Ploting)
                        figure;hold all;
                        plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                    end
                    ModFilta = repmat(ModFilt.',1,NumFra);
                    if SelModFilt==1
                        TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
                    end
                    if  (Ploting)
                        plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                        plot(ff,20*log10(abs(fftshift(ModFilt))));
                        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                    end
                case 'AMSSB'                                           %Out of base band
                    %                     TxSymbModA  = rectpulse(TxSymbModA,OvSam);         %Over sampling
                    TxSymbMod   = TxSymbMod.*exp(1j*2*pi*Ofc*tta);
                    if SelecGaus==1
                        [ ModFilt ] = FiltroGaussiano(ff,OBw,Ofc,1);
                    else
                        [ ModFilt ] = Filtro_Retangular( OBw,Ofc,ff);
                    end
                    ModFilt = fftshift(ModFilt);
                    if  (Ploting)
                        figure;hold all;
                        plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                    end
                    ModFilta = repmat(ModFilt.',1,NumFra);
                    if SelModFilt==1
                        TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
                    end
                    if  (Ploting)
                        plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                        plot(ff,20*log10(abs(fftshift(ModFilt))));
                        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                        a=a+1;
                    end
                otherwise
                    %                     TxSymbModA  = rectpulse(TxSymbModA,OvSam);         %Over sampling
                    [TxSymbMod,tt] = modulate(TxSymbMod,Ofc,Ofs,'pm');
                    if SelecGaus==1
                        [ ModFilt ] = FiltroGaussiano(ff,3*Ofc,0,1);
                    else
                        [ ModFilt ] = Filtro_Retangular( 6*Ofc,0,ff);
                    end
                    ModFilt = fftshift(ModFilt);
                    if  (Ploting)
                        figure;hold all;
                        plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                    end
                    ModFilta = repmat(ModFilt.',1,NumFra);
                    TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
                    if  (Ploting)
                        plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                        plot(ff,20*log10(abs(fftshift(ModFilt))));
                        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                    end
            end
            %
            %                 figure;plot(TxSymbMod(:,1));
            
            %AuxTxSymbMod(1,:)  = TxSymbMod(:);
            TxSymbModA  = [TxSymbMod(1+end-NPOFEX:end,:);TxSymbMod];   %Adding Ciclyc prefix
            U1t  = rectpulse(TxSymbModA,NPPOF);                      %Over sampling
            
            TxSigA = reshape(U1t,OvSam*NumFra*NFFT,CurTesSiz);
            U1t = [TxSigA(1+end-NPSTUf:end,:);TxSigA];                         %Conforming vector sizes to the same length
            % CpPaFr = CpPaFr + 1;
            %Assigning the eletrical signal to one drive of the MZM - The aplitude of
            %the signal can be controlled by the variable DatGai, which can be
            %understood as an gain for the eletrical signal or an atenuation. The
            %second signal will be similar with the only difference a phase shift of pi
            NormFact = max(U1t)./0.95;
            NormFact = repmat(NormFact,size(U1t,1),1);
            U1t = (U1t./NormFact);
            
            % if  (PlotingThisThis)
            %     figure;hold all;grid on;
            %     plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
            %     axis([-25e9 37.5e9 -200 0]);
            %     set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
            % end
            
            
            if ChanAwgn==1
                [~,SigPower] = MeasPower(U1t);
                SigPower2 = SigPower-30;%20*log10(SigPowerI);
                %     SNR = CarSNR + 10*log10(log2(max(DmtMve))) - 10*log10(1*Ts*OBw);
                SNR = CarSNR + 10*log10(log2(max(DmtMve))) - 10*log10(OvSam);
                %     SNR
                U1t = awgn(U1t,SNR,SigPower2);
                %     EoutAux = TxSig;
            else
                %                 U1t = U1t;
            end
        else
            %% DatumTransmission;
            %%               Creating Data for transmission
            % Creating the information accordingly with the type of transmition
            
            %It is important to mention that an stream of data will be formed for each
            %individual carrier. For the last but not the least, the DownStream and the
            %UpStream carriers will be interleaved within the same transmition pass
            %band.For example, if the user choses carriers 1,3,5 and 7 as DownStream
            %carriers, the UpStream will be formed by carriers 2,4,6, and 8. This
            %manuver was suggested by Segatto to address the problem of ICI (Inter
            %Carrier Interference). The transmited signal still an OFDM once the
            %carriers from 1 to 8, for instance, still orthogonal to each other.
            %Therefore, the for loop to generate the information to be transmited will
            %have the passe of 2.
            
            if SendingDowStr==1
                
                %%  Selecting part of the UOFCS
                %The generation fo the Ultra Optical Flat Comb Source is not perfect, which
                %may drive the user to select a specific part of it for its actual use.
                %As a result, this condition verify if the user wants to select a part of
                %the UOFCS or use it as an whole.
                %There are two ways, implemented, to split the incoming OFCS. One by
                %filters e another by the optical FFT. But this step just need to be done
                %once in the simulation.
                
                %After a recent updated where the OIFFT was added to summup the modulated
                %signal, it is important to keep track of the right order of the carriers.
                %Because the OIFFT is very sensitive to the order of the income signals. As
                %it has a periodic response, if the correct signal was not placed at the
                %right port the response will destroy the OFDM signal.
                if ~exist('EoutTx','var')                                              %Verify is this step was previouly done
                    if SelecSetUP==1                                                      %Vefify is the split process will be by filter...
                        EoutTx=SelectEachCarrier(Eout,NumCarr,f(:,1).',fin,FBWD,Order,fc);    %This function is resposible to split each carrier
                        VetThisCarrTx = (RefCarr-1)+1:(RefCarr-1)+NumCarr;             %Keeping track of the right carrier sequence
                    else                                                               %... or by the OFFT
                        %As the OFFT has a periodic response the OFC needs to be constrained other
                        %whise carrier multiple carrier may interfir with other channels, This
                        %first selection was done with a broad pass-band filter. A higher OFFT
                        %order can also be used, although it may increase the computational time
                        %and my not be exactly feasible in the real world.
                        if Selecting==1
                            EoutA = ifft(fft(Eout).*(SelecFilt(:,1).'));
                            %                 PrintInfo(Ploting*1,f,EoutT);
                            %                 a=a+0;
                        else
                            EoutA = Eout;
                        end
                        if UsingGpu==1
                            EoutGpu = gpuArray(EoutA);
                            Tgpu = gpuArray(T);
                            MaxStagGpu = gpuArray(MaxNumStagTx);
                            [EoutAux1Gpu,~,VetThisCarrGpu]=OpticalFFTG(fgpu(:,1).',Tgpu,MaxStagGpu,EoutGpu);
                            EoutTx = gather(EoutAux1Gpu);
                            VetThisCarrTx = gather(VetThisCarrGpu);
                            clear EoutGpu EoutAux1Gpu VetThisCarrGpu Tgpu MaxStagGpu;
                        else
                            [EoutTx,~,VetThisCarrTx]=OpticalFFTN(f(:,1).',T,MaxNumStagTx,EoutA);  %This function is resposible to split each carrier
                        end
                        clear EoutA;
                    end
                    %         PrintInfo(Ploting*51,EoutTx,f);                              %Printing for qualitative analizes.
                    %         axis([min(f) max(f) -400 0]);
                    %         a=a+0;
                end
                switch Modulation
                    %The reason to implement the OFDM system was to work with ROF
                    %architecture. In that area the signal will not be converted to
                    %data before being transmited through optical fiber. The signal
                    %must be the modulator of a optical carrier. At LTE system, the
                    %signal received is OFDM electrical and it will also be the signal
                    %format for the 5G technology. Therefore, it was important to study
                    %how a OFDM signal would behaviour in the architecture here
                    %proposed.
                    case 'OFDM'
                        EoutModTem = zeros(NPPB*Nb,CurTesSiz,NumCarr);
                        TxDataMat = zeros(NumCarr,NumFra*length(DmtMve)*CurTesSiz);
                        for kk=InitCarrDo:CarrPass:NumCarr%Generating different data for each carrier
                            UniDmtMve = unique(DmtMve);
                            UniDmtMve = fliplr(UniDmtMve);
                            DaSiAn = 0;
                            DaSiPo = 0;
                            for CarDmtM=1:length(UniDmtMve)
                                M = UniDmtMve(CarDmtM);
                                DaSiAu = sum(ismember(DmtMve,M));
                                TxData = randi([0 M-1],NumFra*CurTesSiz,DaSiAu);%Generation random information
                                TxDataMat(kk,1+DaSiAn:DaSiAn + CurTesSiz*DaSiAu*NumFra) = TxData(:);%Saving data for future comparison
                                DaSiAn = DaSiAn + CurTesSiz*DaSiAu*NumFra;
                                switch OfdMod%Sellect which modulation was used
                                    case 'qam'
                                        TxSymbAux(1:NumFra*CurTesSiz,1+DaSiPo:DaSiPo + DaSiAu)  = qammod(TxData,M);%Modulating information by QAM
                                    otherwise
                                        TxSymbAux(1:NumFra*CurTesSiz,1+DaSiPo:DaSiPo + DaSiAu)  = dpskmod(TxData,M);%Modulating information DPSK as default
                                end
                                DaSiPo = DaSiPo + DaSiAu;
                            end
                            TxSymb = TxSymbAux.';
                            %Construction the Hermitian Simetry
                            TxSymbConj = conj(flipud(TxSymb));%Singnal conjugate format at reverse order
                            
                            %The Hermitian simetry was composed of the signal its conjugate format and
                            %zeros to fill empty spaces. The signal is formed as show below:
                            %                    signal     midle  signal conjugated
                            %                       |        |       |
                            %                       V        V       V
                            %            TxH = [0 a + jb 0 0 0 0 0 a - jb 0];
                            %
                            %I thought in many ways to form this signal, replacing the centred zeros by
                            %copy of the signal bay just geting its information and adding
                            %respectively. Also, by creating the signal with redundancy at the exact
                            %length needed to make the hermitian singnal. But all tests end up with
                            %similar results, no improvement was noticed.
                            if UsingHermitian
                                TxSymbH = [zeros(1,NumFra*CurTesSiz);TxSymb;zeros(1,NumFra*CurTesSiz);TxSymbConj];%Hermitian Symetri
                                TxSymbC = zeros(NFFT,NumFra*CurTesSiz);
                                TxSymbC(1+Zt-ZtC:Zt-ZtC+size(TxSymbH,1)/2,:) = TxSymbH(1:size(TxSymbH,1)/2,:);
                                TxSymbC(end-Zt+ZtC-(size(TxSymbH,1)/2)+1:end-Zt+ZtC,:) = TxSymbH(1+size(TxSymbH,1)/2:end,:);
                            else
                                TxSymbH = TxSymb;
                                %     TxSymbC = [zeros((NFFT-size(TxSymbH,1))/2,NumFra*CurTesSiz);TxSymbH;zeros((NFFT-size(TxSymbH,1))/2,NumFra*CurTesSiz)];
                                TxSymbC = [TxSymbH(1:size(TxSymbH,1)/2,:);zeros((NFFT-size(TxSymbH,1)),NumFra*CurTesSiz);TxSymbH(end+1-size(TxSymbH,1)/2:end,:)];
                            end
                            %Here was implemented the modulation through FFT process, it basica mix-up
                            %all infomation following Fourier transform.
                            TxSymbMod = ifft(TxSymbC,NFFT,1);
                            TxSymbMod = rectpulse(TxSymbMod,OvSam);
                            tt = linspace(0,1*Te,size(TxSymbMod,1));
                            ff = time2freq(tt);
                            tta = repmat(tt.',1,NumFra*CurTesSiz);
                            if UseIQ
                                RealTx = real(TxSymbMod);
                                ImagTx = imag(TxSymbMod);
                                TxSymbMod = RealTx.*cos(2*pi*FcIQ.*tta) + ImagTx.*sin(2*pi*FcIQ.*tta);
                            end
                            %Sellecting whether the OFDM signal will be transmited on
                            %Base band or not
                            switch SelModTp
                                case 'BP'%Base Band transmission
                                    if SelecGaus==1
                                        [ ModFilt ] = FiltroGaussiano(ff,OBw/2,0,1);%CarrOffSet*FramSpac,1);
                                    else
                                        [ ModFilt ] = Filtro_Retangular(OBw,0,ff);%CarrOffSet*FramSpac,ff);
                                    end
                                    ModFilt = fftshift(ModFilt);
                                    if (kk==1)&&(PlotingThis)
                                        figure;hold all;
                                        plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                                    end
                                    ModFilta = repmat(ModFilt.',1,NumFra*CurTesSiz);
                                    if SelModFilt==1
                                        TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
                                    end
                                    if (kk==1)&&(PlotingThis)
                                        plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                                        plot(ff,20*log10(abs(fftshift(ModFilt))));
                                        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                                    end
                                case 'AM'                                              %Out of base band
                                    TxSymbMod   = TxSymbMod.*cos(2*pi*Ofc*tta);
                                    if SelecGaus==1
                                        [ ModFilt ] = FiltroGaussiano(ff,OBw,Ofc,1);
                                    else
                                        [ ModFilt ] = Filtro_Retangular( OBw,Ofc,ff);
                                    end
                                    ModFilt = fftshift(ModFilt);
                                    if (kk==1)&&(PlotingThis)
                                        figure;hold all;
                                        plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                                    end
                                    ModFilta = repmat(ModFilt.',1,NumFra*CurTesSiz);
                                    if SelModFilt==1
                                        TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
                                    end
                                    if (kk==1)&&(PlotingThis)
                                        plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                                        plot(ff,20*log10(abs(fftshift(ModFilt))));
                                        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                                    end
                                case 'AMSSB'                                           %Out of base band
                                    TxSymbMod   = TxSymbMod.*exp(1j*2*pi*Ofc*tta);
                                    if SelecGaus==1
                                        [ ModFilt ] = FiltroGaussiano(ff,OBw,Ofc,1);
                                    else
                                        [ ModFilt ] = Filtro_Retangular( OBw,Ofc,ff);
                                    end
                                    ModFilt = fftshift(ModFilt);
                                    if (kk==1)&&(PlotingThis)
                                        figure;hold all;
                                        plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                                    end
                                    ModFilta = repmat(ModFilt.',1,NumFra*CurTesSiz);
                                    if SelModFilt==1
                                        TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
                                    end
                                    if (kk==1)&&(PlotingThis)
                                        plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                                        plot(ff,20*log10(abs(fftshift(ModFilt))));
                                        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                                        a=a+1;
                                    end
                                otherwise
                                    if SelecGaus==1
                                        [ ModFilt ] = FiltroGaussiano(ff,3*Ofc,0,1);
                                    else
                                        [ ModFilt ] = Filtro_Retangular( 6*Ofc,0,ff);
                                    end
                                    ModFilt = fftshift(ModFilt);
                                    if (kk==1)&&(Ploting)
                                        figure;hold all;
                                        plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                                    end
                                    ModFilta = repmat(ModFilt.',1,NumFra*CurTesSiz);
                                    TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
                                    if (kk==1)&&(Ploting)
                                        plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                                        plot(ff,20*log10(abs(fftshift(ModFilt))));
                                        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                                    end
                            end
                            TxSymbModA  = [TxSymbMod(1+end-NPOFEX:end,:);TxSymbMod];%Adding Ciclyc prefix
                            U1t  = rectpulse(TxSymbModA,NPPOF);%Over sampling
                            TxSigA = reshape(U1t,OvSam*NumFra*NFFT,CurTesSiz);
                            U1t = [TxSigA(1+end-NPSTUf:end,:);TxSigA];%Conforming vector sizes to the same length
                            
                            %Assigning the eletrical signal to one drive of the MZM - The aplitude of
                            %the signal can be controlled by the variable DatGai, which can be
                            %understood as an gain for the eletrical signal or an atenuation. The
                            %second signal will be similar with the only difference a phase shift of pi
                            NormFact = max(U1t)./0.95;
                            NormFact = repmat(NormFact,size(U1t,1),1);
                            U1t = (U1t./NormFact);%Normalizing the sinal to mach with the MZM espect to receive.
%                             U.U1t = U1t;
                            U2t = exp(-1j*pi).*U1t;
                            %As both signals will have the mostrly the same characteristics with the
                            %only difference the phase shift of 180 degress. The MZM-I will be working
                            %on the Push-Pull configuration. It is necessary to reduce the Chirp noise
                            %to zero.
                            EoutTxAux = repmat(EoutTx(VetThisCarrTx==(RefCarr-1+kk),:).',1,CurTesSiz);
                            if UsingGpu==1
                                Lgpu = gpuArray(L);
                                U0gpu = gpuArray(U0);
                                U_pi1gpu = gpuArray(U_pi1);
                                U_pi2gpu = gpuArray(U_pi2);
                                noptgpu = gpuArray(nopt);
                                nelgpu = gpuArray(nel);
                                Cgpu = gpuArray(C);
                                alfa0gpu = gpuArray(alfa0);
                                for jj=1:size(EoutTxAux,2)
                                    U1tgpu = gpuArray(U1t(:,jj));
                                    U2tgpu = gpuArray(U2t(:,jj));
                                    EoutTxGpu = gpuArray(EoutTxAux(:,jj));
                                    freqGHzGpu = gpuArray(freqGHz(:,jj));
                                    [EoutMod,~]=MZMGT(freqGHzGpu,EoutTxGpu,U1tgpu,U2tgpu,Lgpu,U0gpu,U_pi1gpu,U_pi2gpu,noptgpu,nelgpu,Cgpu,alfa0gpu);%Modulating individual carriers
                                    EoutModTem(:,jj,ObsCarrPos(kk)) = gather(EoutMod);%Adding the current Modulated output to the final OutPut
                                    clear U1tgpu U2tgpu EoutTxGpu freqGHzGpu EoutMod
                                end
                                clear Lgpu U0gpu U_pi1gpu U_pi2gpu noptgpu nelgpu Cgpu alfa0gpu;
                            else
                                U.U1t = U1t;
                                U.U2t = U2t;
                                [EoutMod,~]=MZM(freqGHz,EoutTxAux,U,L,U0,U_pi1,U_pi2,nopt,nel,C,alfa0);%Modulating individual carriers
                                EoutModTem(:,:,ObsCarrPos(kk)) = EoutMod;%Adding the current Modulated output to the final OutPut
                            end
                            if (kk==1)&&(PlotingThis)
                                figure;hold all;grid on;
                                plot(f(1:Nb*NPPB),20*log10(abs(fftshift(fft(EoutModTem(1:Nb*NPPB,1,ObsCarrPos(kk)))./length(EoutModTem(1:Nb*NPPB,1,ObsCarrPos(kk)))))));
                                axis([-25e9 37.5e9 -200 0]);
                                set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                            end
                            %This section adds the unmodulated upstream carrier on the optical OFDM
                            %signal
                            if (~mod(CarrPass,2))&&(SendingUp)
                                clear EoutModAux;
                                EoutModAux = EoutTx(VetThisCarrTx==(RefCarr-1+kk+1),:);
                                EoutMod = repmat(EoutModAux.',1,CurTesSiz);
                                EoutModTem(:,:,ObsCarrPos(kk+1)) = EoutMod;%Adding the current Unmodulated output to the final OutPut
                            end
                        end
                        %At this section individual and parralel carriers will be assemble to
                        %compose the optical OFDM signal. It can be done in two way, with a simple
                        %adder, which gives a small ICI clearance. Or by the optical IFFT, which
                        %gives a better ICI clearance.
                        if IfftOrSum==1                                                   %Use the OIFFT
                            if size(EoutModTem,3)>1                                    %Just if the signal has more than 1 dimension
                                if UsingGpu==1
                                    Tgpu = gpuArray(T);
                                    MaxStagGpu = gpuArray(MaxNumStagT);
                                    EoutMod = zeros(size(EoutModTem,1),size(EoutModTem,2));
                                    for jj=1:size(EoutModTem,2)
                                        EoutModTemGpu = gpuArray(EoutModTem(:,jj,:));
                                        [EoutAux1Gpu,~,~]=OpticalIFFTGT(fgpu,Tgpu,MaxStagGpu,EoutModTemGpu);
                                        EoutMod(:,jj) = gather(EoutAux1Gpu);
                                        clear EoutModTemGpu EoutAux1Gpu;
                                    end
                                    clear MaxStagGpu Tgpu;
                                else
                                    [EoutMod,~,~] = OpticalIFFTNT(f,T,MaxNumStagT,EoutModTem);
                                end
                            else
                                EoutMod = EoutModTem;                                  %otherwise just let the signal pass
                            end
                        else                                                           %Use a simple adder which don't display a better result
                            if size(EoutModTem,3)>1
                                EoutMod = sum(EoutModTem,3);
                            else
                                EoutMod = EoutModTem;
                            end
                        end
                    case 'DPSK'
                        EoutModTem = zeros(NPPB*Nb,CurTesSiz,NumCarr);
                        TxDataMat = zeros(NumCarr,NbDPSK*CurTesSiz);
                        %%            Generate the data DPSK
                        for kk=InitCarrDo:CarrPass:NumCarr%Generating different data for each carrier
                            if CarrUsedDo(kk)
                                if SendingData% Frist it is chosen to transmite a random information...
                                    TxData = (randi(2,NbDPSK,CurTesSiz)-1);%Creating the data stream to be transmited
                                    TxData(1:JusPos,:) = JusVal;%Adding the Justification simble at the begining of the stream to sincronize received data frame
                                    TxData(end-(JusPos-1):end,:) = JusValEnd;%Adding the Justification simble at the end of the stream to sincronize received data frame
                                else%... or just one high pulse for testing the channel
                                    TxData = zeros(NbDPSK,CurTesSiz);%Creating the base fo the data stream
                                    TxData(end/2,:) = 1;%Adding the code in the data stream to produce the highest level possible
                                end
                                TxDataMat(kk,:) = TxData(:);%Storring the transmited information for latter evaluation
                                %Once that we have the data to be transmited it needs to be converted to
                                %the eletrical signal that hereafter will modulate our carrier through the
                                %MZM-I
                                U1t = DpskEncodEqT(TxData);
                                U1t(U1t==1) =  1.9;
                                U1t(U1t==0) = -1.9;
                                U1t = rectpulse(U1t,NPPB);
                                % Adding CP to the data
                                if AddCP==1
                                    TxAuxB = reshape(U1t,NPPB,NbDPSK*CurTesSiz);
                                    if SetCpSampZer==1
                                        TxAuxA = [zeros(NumAmosCP,size(TxAuxB,2));TxAuxB;zeros(NumAmosCP,size(TxAuxB,2))];
                                    else
                                        %TxAuxA = [TxAuxB(1:NumAmosCP,:);TxAuxB;TxAuxB(end-(NumAmosCP-1):end,:)];
                                        TxAuxA = [flipud(TxAuxB(1:NumAmosCP,:));TxAuxB;flipud(TxAuxB(end-(NumAmosCP-1):end,:))];
                                    end
                                    TxAux = reshape(TxAuxA,(2*NumAmosCP+NPPB)*NbDPSK,CurTesSiz);
                                    U1t = [TxAux;TxAux(end-(StuffSampels-1):end,:)];
                                end
                                
                                %Because of reasons, the real world does not have components capable
                                %instantly change the voltage level. Thus, to emulate this behaviour the
                                %eletrical PAM signal pass through a gaussian filter that will remove the
                                %frequencies of higher order, which will result in small slope in the level
                                %variation
                                
                                %[BitFilt,~] = FiltroGaussiano(f,BWD,CenFeq,FiltOrd);   %Creating filter for conformation of the input information
                                %BitFilt     = fftshift(BitFilt);                       %Doing a shift on the vector for matching the transmited data
                                %TxSig      = ifft(fft(SigTx).*BitFilt);               %Conforming the information and Creating the modulation signal
                                %                                 U1t       = SigTx;                                   %Adding a gain to the eletrical signal
                                
                                %Assigning the eletrical signal to one drive of the MZM - The aplitude of
                                %the signal can be controlled by the variable DatGai, which can be
                                %understood as an gain for the eletrical signal or an atenuation. The
                                %second signal will be similar with the only difference a phase shift of pi
                                %U1t = U1t;
                                U2t = exp(-1j*pi).*U1t;
                                %As both signals will have the mostrly the same characteristics with the
                                %only difference the phase shift of 180 degress. The MZM-I will be working
                                %on the Push-Pull configuration. It is necessary to reduce the Chirp noise
                                %to zero.
                                EoutTxAux = repmat(EoutTx(VetThisCarrTx==(RefCarr-1+kk),:).',1,CurTesSiz);
                                if UsingGpu==1
                                    Lgpu = gpuArray(L);
                                    U0gpu = gpuArray(U0);
                                    U_pi1gpu = gpuArray(U_pi1);
                                    U_pi2gpu = gpuArray(U_pi2);
                                    noptgpu = gpuArray(nopt);
                                    nelgpu = gpuArray(nel);
                                    Cgpu = gpuArray(C);
                                    alfa0gpu = gpuArray(alfa0);
                                    for jj=1:size(EoutTxAux,2)
                                        U1tgpu = gpuArray(U1t(:,jj));
                                        U2tgpu = gpuArray(U2t(:,jj));
                                        EoutTxGpu = gpuArray(EoutTxAux(:,jj));
                                        freqGHzGpu = gpuArray(freqGHz(:,jj));
                                        [EoutMod,~]=MZMGT(freqGHzGpu,EoutTxGpu,U1tgpu,U2tgpu,Lgpu,U0gpu,U_pi1gpu,U_pi2gpu,noptgpu,nelgpu,Cgpu,alfa0gpu);%Modulating individual carriers
                                        EoutModTem(:,jj,ObsCarrPos(kk)) = gather(EoutMod);%Adding the current Modulated output to the final OutPut
                                        clear U1tgpu U2tgpu EoutTxGpu freqGHzGpu EoutMod
                                    end
                                    clear Lgpu U0gpu U_pi1gpu U_pi2gpu noptgpu nelgpu Cgpu alfa0gpu;
                                else
                                    U.U1t = U1t;
                                    U.U2t = U2t;
                                    [EoutMod,~]=MZM(freqGHz,EoutTxAux,U,L,U0,U_pi1,U_pi2,nopt,nel,C,alfa0);%Modulating individual carriers
                                    EoutModTem(:,:,ObsCarrPos(kk)) = EoutMod;%Adding the current Modulated output to the final OutPut
                                end
                                %[EoutMod,~]=MZM(freqGHz,EoutTxAux,U,L,U0,U_pi1,U_pi2,nopt,nel,C,alfa0);
                                %EoutModTem(:,:,ObsCarrPos(kk)) = EoutMod;%Adding the current Modulated output to the final OutPut
                                if (~mod(CarrPass,2))&&(SendingUp)
                                    clear EoutModAux;
                                    EoutMod = repmat(EoutTx(VetThisCarrTx==(RefCarr-1+kk+1),:).',1,CurTesSiz);
                                    EoutModTem(:,:,ObsCarrPos(kk+1)) = EoutMod;%Adding the current Modulated output to the final OutPut
                                end
                            else
                                EoutMod = EoutTx(VetThisCarrTx==(RefCarr-1+kk),:);
                                
                                EoutModTem(ObsCarrPos(kk),1:length(EoutMod)) = EoutMod;%Adding the current Modulated output to the final OutPut
                                if ~mod(CarrPass,2)
                                    clear EoutModAux;
                                    EoutMod = repmat(EoutTx(VetThisCarrTx==(RefCarr-1+kk+1),:).',1,CurTesSiz);
                                    EoutModTem(:,:,ObsCarrPos(kk+1)) = EoutMod;%Adding the current Unmodulated output to the final OutPut
                                end
                            end
                        end
                        if IfftOrSum==1
                            if size(EoutModTem,3)>1
                                if UsingGpu==1
                                    Tgpu = gpuArray(T);
                                    MaxStagGpu = gpuArray(MaxNumStagT);
                                    EoutMod = zeros(size(EoutModTem,1),size(EoutModTem,2));
                                    for jj=1:size(EoutModTem,2)
                                        EoutModTemGpu = gpuArray(EoutModTem(:,jj,:));
                                        [EoutAux1Gpu,~,~]=OpticalIFFTGT(fgpu,Tgpu,MaxStagGpu,EoutModTemGpu);
                                        EoutMod(:,jj) = gather(EoutAux1Gpu);
                                        clear EoutModTemGpu EoutAux1Gpu;
                                    end
                                    clear MaxStagGpu Tgpu;
                                else
                                    [EoutMod,~,~] = OpticalIFFTNT(f,T,MaxNumStagT,EoutModTem);
                                end
                            else
                                EoutMod = EoutModTem;
                            end
                        else
                            if size(EoutModTem,3)>1
                                EoutMod = sum(EoutModTem,3);
                            else
                                EoutMod = EoutModTem;
                            end
                        end
                    case 'DQPSK'
                        %%            Generate the data DQPSK
                        EoutModTem = zeros(NbDQPSK,CurTesSiz,NumCarr);
                        TxDataMat = zeros(2*NumCarr,(NbDQPSK/2)*CurTesSiz);
                        for kk=InitCarrDo:CarrPass:NumCarr%Generating different data for each carrier
                            if CarrUsedDo(kk)
                                if SendingData == 1% Frist it is chosen to transmite a random information...
                                    TxData = (randi(2,NbDQPSK,CurTesSiz)-1);%Creating the data stream to be transmited
                                    TxData(1:JusPos,:) = JusVal;%Adding the Justification simble at the begining of the stream to sincronize received data frame
                                    TxData(end-(JusPos-1):end,:) = JusValEnd;%Adding the Justification simble at the end of the stream to sincronize received data frame
                                    PreBits  = zeros(2,CurTesSiz);%Seting the start point for the DQPSK maping
                                    
                                else%... or just one high pulse for testing the channel
                                    TxData = zeros(NbDQPSK,CurTesSiz);%Creating the base fo the data stream
                                    TxData((end/2)-1:end/2,:) = 1;%Adding the code in the data stream to produce the highest level possible
                                    PreBits = zeros(2,CurTesSiz);
                                end
                                TxDataPos = linspace(1,NbDQPSK,NbDQPSK);%Auxiliar variable to split the TxData
                                [DataI,DataQ] = DqpskEncodEqT(TxData,PreBits);%Maping the income data to the DQPSK format
                                %Converting the I and Q components to the polirezed NRZ format
                                DataI(DataI==1) =  1.9;
                                DataI(DataI==0) = -1.9;
                                DataQ(DataQ==1) =  1.9;
                                DataQ(DataQ==0) = -1.9;
                                TxOdd = TxData(logical(mod(TxDataPos,2)),:);%Spliting the information of odd positions
                                TxEven = TxData(~(mod(TxDataPos,2)),:);%Spliting the information of even positions
                                TxDataMat(kk,:) = TxOdd(:);%Storring the transmited information for latter evaluation
                                TxDataMat(kk+NumCarr,:) = TxEven(:);%Storring the transmited information for latter evaluation
                                %Once that we have the data to be transmited it needs to be converted to
                                %the eletrical signal that hereafter will modulate our carrier through the
                                %MZM-I
                                SigI = rectpulse(DataI,NPPB);
                                SigQ = rectpulse(DataQ,NPPB);
                                
                                % Adding CP to the data
                                if AddCP
                                    TxAux1 = reshape(SigI,NPPB,CurTesSiz*(NbDQPSK/2));
                                    if SetCpSampZer==1
                                        TxAux2 = [zeros(NumAmosCP,size(TxAux1,2));TxAux1;zeros(NumAmosCP,size(TxAux1,2))];
                                    else
                                        %TxAux2 = [TxAux1(1:NumAmosCP,:);TxAux1;TxAux1(end-(NumAmosCP-1):end,:)];
                                        TxAux2 = [flipud(TxAux1(1:NumAmosCP,:));TxAux1;flipud(TxAux1(end-(NumAmosCP-1):end,:))];
                                    end
                                    TxAux1 = reshape(TxAux2,(2*NumAmosCP+NPPB)*(NbDQPSK/2),CurTesSiz);
                                    SigI = [TxAux1;TxAux1(end-(StuffSampels-1):end,:)];
                                    
                                    TxAux3 = reshape(SigQ,NPPB,CurTesSiz*(NbDQPSK/2));
                                    if SetCpSampZer==1
                                        TxAux4 = [zeros(NumAmosCP,size(TxAux3,2));TxAux3;zeros(NumAmosCP,size(TxAux3,2))];
                                    else
                                        %TxAux4 = [TxAux3(1:NumAmosCP,:);TxAux3;TxAux3(end-(NumAmosCP-1):end,:)];
                                        TxAux4 = [flipud(TxAux3(1:NumAmosCP,:));TxAux3;flipud(TxAux3(end-(NumAmosCP-1):end,:))];
                                    end
                                    TxAux5 = reshape(TxAux4,(2*NumAmosCP+NPPB)*(NbDQPSK/2),CurTesSiz);
                                    SigQ = [TxAux5;TxAux5(end-(StuffSampels-1):end,:)];
                                end
                                
                                %Because of reasons, the real world does not have components capable
                                %instantly change the voltage level. Thus, to emulate this behaviour the
                                %eletrical PAM signal pass through a gaussian filter that will remove the
                                %frequencies of higher order, which will result in small slope in the level
                                %variation
                                [BitFilt,~] = FiltroGaussiano(f,BWD,CenFeq,FiltOrd);   %Creating filter for conformation of the input information
                                BitFilt = fftshift(BitFilt);                           %Doing a shift on the vector for matching the transmited data
                                %TxSigI = SigI;                                         %Conforming the information and Creating the modulation signal
                                %TxSigQ = SigQ;                                         %Conforming the information and Creating the modulation signal
                                %TxSigI = ifft(fft(SigI).*BitFilt);                    %Conforming the information and Creating the modulation signal
                                %TxSigQ = ifft(fft(SigQ).*BitFilt);                    %Conforming the information and Creating the modulation signal
                                EoutTxAux = repmat(EoutTx(VetThisCarrTx==(RefCarr-1+kk),:).',1,CurTesSiz);
                                [EoutMod] = IqMod(EoutTxAux,SigI,SigQ,Vpi,V0);
                                EoutModTem(1:size(EoutMod,1),:,ObsCarrPos(kk)) = EoutMod;%Adding the current Modulated output to
                                %This section adds the unmodulated upstream carrier on the optical OFDM
                                %signal
                                if (~mod(CarrPass,2))&&(SendingUp)
                                    clear EoutModAux;
                                    EoutModAux = EoutTx(VetThisCarrTx==(RefCarr-1+kk+1),:);
                                    EoutMod = repmat(EoutModAux.',1,CurTesSiz);
                                    EoutModTem(:,:,ObsCarrPos(kk+1)) = EoutMod;%Adding the current Unmodulated output to the final OutPut
                                end
                                if (kk==1)&&(PlotingThis)
                                    figure;hold all;grid on;
                                    plot(f(1:Nb*NPPB,1),20*log10(abs(fftshift(fft(EoutMod(1:Nb*NPPB,1))./length(EoutMod(1:Nb*NPPB,1))))));
                                    axis([-25e9 37.5e9 -200 0]);
                                    set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                                end
                            else
                                EoutModAux = EoutTx(VetThisCarrTx==(RefCarr-1+kk),:);
                                EoutMod = repmat(EoutModAux.',1,CurTesSiz);
                                EoutModTem(ObsCarrPos(kk),1:length(EoutMod)) = ...
                                    EoutMod;%Adding the current Modulated output to the final OutPut
                                if ~mod(CarrPass,2)
                                    clear EoutModAux;
                                    EoutModAux = EoutTx(VetThisCarrTx==(RefCarr-1+kk+1),:);
                                    EoutMod = repmat(EoutModAux.',1,CurTesSiz);
                                    EoutModTem(:,:,ObsCarrPos(kk+1)) = EoutMod;%Adding the current Modulated output to the final OutPut
                                end
                            end
                        end
                        if IfftOrSum==1
                            if size(EoutModTem,3)>1
                                if UsingGpu==1
                                    Tgpu = gpuArray(T);
                                    MaxStagGpu = gpuArray(MaxNumStagT);
                                    EoutMod = zeros(size(EoutModTem,1),size(EoutModTem,2));
                                    for jj=1:size(EoutModTem,2)
                                        EoutModTemGpu = gpuArray(EoutModTem(:,jj,:));
                                        [EoutAux1Gpu,~,~]=OpticalIFFTGT(fgpu,Tgpu,MaxStagGpu,EoutModTemGpu);
                                        EoutMod(:,jj) = gather(EoutAux1Gpu);
                                        clear EoutModTemGpu EoutAux1Gpu;
                                    end
                                    clear MaxStagGpu Tgpu;
                                else
                                    [EoutMod,~,~] = OpticalIFFTNT(f,T,MaxNumStagT,EoutModTem);
                                end
                            else
                                EoutMod = EoutModTem;
                            end
                        else
                            if size(EoutModTem,3)>1
                                EoutMod = sum(EoutModTem,3);
                            else
                                EoutMod = EoutModTem;
                            end
                        end
                    case '4PAM'
                        %%        Generate the data 4PAM
                        EoutModTem = zeros(NPPB*Nb,CurTesSiz,NumCarr);
                        TxDataMat = zeros(NumCarr,Nb4Pam*CurTesSiz);
                        for kk=InitCarrDo:CarrPass:NumCarr%Generating different data for each carrier
                            if CarrUsedDo(kk)% Frist it is chosen to transmite a random information...
                                if SendingData==1
                                    TxData = (randi(2,Nb4Pam,CurTesSiz)-1);%Creating the data stream to be transmited
                                    TxData(1:JusPos,:) = JusVal;%Adding the Justification simble at the begining of the stream to sincronize received data frame
                                    TxData(end-(JusPos-1):end,:) = JusValEnd;%Adding the Justification simble at the end of the stream to sincronize received data frame
                                else%... or just one high pulse for testing the channel
                                    TxData = zeros(Nb4Pam,CurTesSiz);%Creating the base fo the data stream
                                    TxDataPos = linspace(1,Nb4Pam,Nb4Pam);%Creating a vector to auxilliate the addresing process
                                    TxData(~mod(TxDataPos,2),:)=1;%Adding the code in the data stream to produce the lowest level possible
                                    TxData(end/2 - 1,:) = 1;%Adding the code in the data stream to produce the highest level possible
                                end
                                TxDataMat(kk,:) = TxData(:);%Storring the transmited information for latter evaluation
                                %Once that we have the data to be transmited it needs to be
                                %converted to the eletrical signal that hereafter will modulate
                                %our carrier through the MZM-I
                                if Which4PAM==1%Selecting if the PAM will be done on Electrical or Optical domain
                                    %  Good for others kms
                                    [Phi1,Phi2] = Maping4PamIqT(TxData,Vmin,Vmax,ModSchem,FiberLength,SetCpSampZer);%Generating the eletrical signal for the optical PAM4
                                else
                                    [Phi1,Phi2] = Maping4PamT(TxData,VPI,Polirized,MaxAmp4PAM);%Generating the eletrical signal for the electrical PAM4
                                end
                                %The signal generated are not yet with the same number of samples as the
                                %OFCS loaded. These nexte lines do oversampling of the electrical signal.
                                TxSig1 = rectpulse(Phi1,NPPB);
                                TxSig2 = rectpulse(Phi2,NPPB);
                                %Thus, if it would be required to add cycle prefix the number of samples
                                %per symbol needs to change as well as some adjustments needs to be done
                                %for the new signal match in size with the size of the vector time. This
                                %problem just exist on simulation, at practice the main point is the
                                %syncronism of the signals.
                                % Adding CP to the data
                                if AddCP==1
                                    TxAux1 = reshape(TxSig1,NPPB,CurTesSiz*(Nb4Pam/2));
                                    if SetCpSampZer==1
                                        TxAux2 = [zeros(NumAmosCP,size(TxAux1,2));TxAux1;zeros(NumAmosCP,size(TxAux1,2))];
                                    else
                                        TxAux2 = [flipud(TxAux1(1:NumAmosCP,:));TxAux1;flipud(TxAux1(end-(NumAmosCP-1):end,:))];
                                    end
                                    TxAux3 = reshape(TxAux2,(2*NumAmosCP+NPPB)*(Nb4Pam/2),CurTesSiz);
                                    TxSig1 = [TxAux3;TxAux3(end-(StuffSampels-1):end,:)];
                                    TxAux4 = reshape(TxSig2,NPPB,CurTesSiz*(Nb4Pam/2));
                                    if SetCpSampZer==1
                                        TxAux5 = [zeros(NumAmosCP,size(TxAux4,2));TxAux4;zeros(NumAmosCP,size(TxAux4,2))];
                                    else
                                        TxAux5 = [flipud(TxAux4(1:NumAmosCP,:));TxAux4;flipud(TxAux4(end-(NumAmosCP-1):end,:))];
                                    end
                                    TxAux6 = reshape(TxAux5,(2*NumAmosCP+NPPB)*(Nb4Pam/2),CurTesSiz);
                                    TxSig2 = [TxAux6;TxAux6(end-(StuffSampels-1):end,:)];
                                end
                                %Because of reasons, the real world does not have components capable
                                %instantly change the voltage level. Thus, to emulate this behaviour the
                                %eletrical PAM signal pass through a gaussian filter that will remove the
                                %frequencies of higher order, which will result in small slope in the level
                                %variation
                                [BitFilt,~] = FiltroGaussiano(f,BWD,CenFeq,FiltOrd);   %Creating filter for conformation of the input information
                                BitFilt     = fftshift(BitFilt);                       %Doing a shift on the vector for matching the transmited data
                                %If the user choses to modulate all carriers at the same time there is no
                                %need of this script to generate data for each individual carrier.
                                %Therefore, this For loop can be halted
                                U1t = TxSig2;                                        %Assigning the electrical signal to one drive of the MZM
                                U2t = TxSig1;                                        %Assigning the electrical signal to another drive of the MZM
                                EoutTxAux = repmat(EoutTx(VetThisCarrTx==(RefCarr-1+kk),:).',1,CurTesSiz);
                                if UsingGpu==1
                                    Lgpu = gpuArray(L);
                                    U0gpu = gpuArray(U0);
                                    U_pi1gpu = gpuArray(U_pi1);
                                    U_pi2gpu = gpuArray(U_pi2);
                                    noptgpu = gpuArray(nopt);
                                    nelgpu = gpuArray(nel);
                                    Cgpu = gpuArray(C);
                                    alfa0gpu = gpuArray(alfa0);
                                    for jj=1:size(EoutTxAux,2)
                                        U1tgpu = gpuArray(U1t(:,jj));
                                        U2tgpu = gpuArray(U2t(:,jj));
                                        EoutTxGpu = gpuArray(EoutTxAux(:,jj));
                                        freqGHzGpu = gpuArray(freqGHz(:,jj));
                                        [EoutMod,~]=MZMGT(freqGHzGpu,EoutTxGpu,U1tgpu,U2tgpu,Lgpu,U0gpu,U_pi1gpu,U_pi2gpu,noptgpu,nelgpu,Cgpu,alfa0gpu);%Modulating individual carriers
                                        EoutModTem(:,jj,ObsCarrPos(kk)) = gather(EoutMod);%Adding the current Modulated output to the final OutPut
                                        clear U1tgpu U2tgpu EoutTxGpu freqGHzGpu EoutMod
                                    end
                                    clear Lgpu U0gpu U_pi1gpu U_pi2gpu noptgpu nelgpu Cgpu alfa0gpu;
                                else
                                    U.U1t = U1t;
                                    U.U2t = U2t;
                                    [EoutMod,~]=MZM(freqGHz,EoutTxAux,U,L,U0,U_pi1,U_pi2,nopt,nel,C,alfa0);%Modulating individual carriers
                                    EoutModTem(:,:,ObsCarrPos(kk)) = EoutMod;%Adding the current Modulated output to the final OutPut
                                end
                                %if Which4PAM==1
                                    %if ModSchem
                                        %[EoutMod] = MZM(freqGHz,EoutTxAux,U,L,U0,U_pi1,U_pi2,nopt,nel,C,alfa0);
                                    %else
                                        %[EoutMod] = MZM(freqGHz,EoutTxAux,U,L,U0,U_pi1,U_pi2,nopt,nel,C,alfa0);
                                    %end
                                %else
                                    %[EoutMod] = MZM(freqGHz,EoutTxAux,U,L,U0,U_pi1,U_pi2,nopt,nel,C,alfa0);
                                %end
                                %EoutModTem(:,:,ObsCarrPos(kk)) = EoutMod;%Adding the current Modulated output to
                                %This section adds the unmodulated upstream carrier on the optical OFDM
                                %signal
                                if (~mod(CarrPass,2))&&(SendingUp)
                                    clear EoutModAux;
                                    EoutModAux = EoutTx(VetThisCarrTx==(RefCarr-1+kk+1),:);
                                    EoutMod = repmat(EoutModAux.',1,CurTesSiz);
                                    EoutModTem(:,:,ObsCarrPos(kk+1)) = EoutMod;%Adding the current Unmodulated output to the final OutPut
                                end
                                if (kk==1)&&(PlotingThis)
                                    figure;hold all;grid on;
                                    plot(f(1:Nb*NPPB,1),20*log10(abs(fftshift(fft(EoutMod(1:Nb*NPPB,1))./length(EoutMod(1:Nb*NPPB,1))))));
                                    axis([-25e9 37.5e9 -200 0]);
                                    set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                                end
                            else
                                EoutModAux = EoutTx(VetThisCarrTx==(RefCarr-1+kk),:);
                                EoutMod = repmat(EoutModAux.',1,CurTesSiz);
                                EoutModTem(ObsCarrPos(kk),1:length(EoutMod)) = EoutMod;%Adding the current Modulated output to the final OutPut
                                if ~mod(CarrPass,2)
                                    clear EoutModAux;
                                    EoutModAux = EoutTx(VetThisCarrTx==(RefCarr-1+kk+1),:);
                                    EoutMod = repmat(EoutModAux.',1,CurTesSiz);
                                    EoutModTem(:,:,ObsCarrPos(kk+1)) = EoutMod;%Adding the current Modulated output to the final OutPut
                                end
                            end
                        end
                        if IfftOrSum==1
                            if size(EoutModTem,3)>1
                                if UsingGpu==1
                                    Tgpu = gpuArray(T);
                                    MaxStagGpu = gpuArray(MaxNumStagT);
                                    EoutMod = zeros(size(EoutModTem,1),size(EoutModTem,2));
                                    for jj=1:size(EoutModTem,2)
                                        EoutModTemGpu = gpuArray(EoutModTem(:,jj,:));
                                        [EoutAux1Gpu,~,~]=OpticalIFFTGT(fgpu,Tgpu,MaxStagGpu,EoutModTemGpu);
                                        EoutMod(:,jj) = gather(EoutAux1Gpu);
                                        clear EoutModTemGpu EoutAux1Gpu;
                                    end
                                    clear MaxStagGpu Tgpu;
                                else
                                    [EoutMod,~,~] = OpticalIFFTNT(f,T,MaxNumStagT,EoutModTem);
                                end
                            else
                                EoutMod = EoutModTem;
                            end
                            a=a+1;
                        else
                            if size(EoutModTem,3)>1
                                EoutMod = sum(EoutModTem,3);
                            else
                                EoutMod = EoutModTem;
                            end
                        end
                    otherwise
                        %%        Generate the data OOK
                        EoutModTem = zeros(NPPB*Nb,CurTesSiz,NumCarr);
                        TxDataMat = zeros(NumCarr,(Nb-NumBitDesc)*CurTesSiz);
                        for kk=InitCarrDo:CarrPass:NumCarr%Generating different data for each carrier
                            if CarrUsedDo(kk)% Frist it is chosen to transmite just one high pulse for
                                if AdjusData==1%testing the channel...
                                    TxData = zeros(Nb-NumBitDesc,CurTesSiz);%creating vector of zerros for correcting the delay caused by the transmission
                                    TxData(end/2,:) = 1;%Adding an strategical pulse for measuring the result.
                                    TxDataMat(kk,:) = TxData(:);%Storring the transmited information for latter evaluation
                                    if NRZPolarity==1
                                        TxData(TxData==0) = NrzMin;
                                        TxData(TxData==1) = NrzMax;
                                    end
                                else%... or just a random information
                                    TxData = (randi(2,Nb-NumBitDesc,CurTesSiz)-1);%Creating Random Information that will be loaded in each individual subcarrier
                                    TxData(1:JusLen,:) = JusVal;%Making the First 4 bits equal to zero
                                    TxData(end-(JusLen-1):end,:) = JusVal;%Making the Last 4 bits equal to zero
                                    TxDataMat(kk,:) = TxData(:);%Storring the transmited information for latter evaluation
                                    if NRZPolarity
                                        TxData(TxData==0) = NrzMin;
                                        TxData(TxData==1) = NrzMax;
                                    end
                                end
                                
                                %The signal generated are not yet with the same number of samples as the
                                %OFCS loaded. These nexte lines do oversampling
                                TxDataRes = rectpulse(TxData,NPPB);%Changing the length of the Data acordingly with the time vector
                                %Thus, if it would be required to add cycle prefix the number of samples
                                %per symbol needs to change as well as some adjustments needs to be done
                                %for the new signal match in size with the size of the vector time. This
                                %problem just exist on simulation, at practice the main point is the
                                %syncronism of the signals.
                                % Adding CP to the data
                                if AddCP==1
                                    TxAux = reshape(TxDataRes,NPPB,CurTesSiz*(Nb-NumBitDesc));
                                    if SetCpSampZer==1
                                        TxAux1 = [zeros(NumAmosCP,size(TxAux,2));TxAux;zeros(NumAmosCP,size(TxAux,2))];
                                    else
                                        TxAux1 = [flipud(TxAux(1:NumAmosCP,:));TxAux;flipud(TxAux(end-(NumAmosCP-1):end,:))];
                                    end
                                    TxAux2 = reshape(TxAux1,(2*NumAmosCP+NPPB)*(Nb-NumBitDesc),CurTesSiz);
                                    U1t = [TxAux2;TxAux2(end-(StuffSampels-1):end,:)];
                                end
                                
                                
                                %Because of reasons, the real world does not have components capable
                                %instantly change the voltage level. Thus, to emulate this behaviour the
                                %eletrical PAM signal pass through a gaussian filter that will remove the
                                %frequencies of higher order, which will result in small slope in the level
                                %variation
                                
                                [BitFilt,~] = FiltroGaussiano(f,BWD,CenFeq,FiltOrd);   %Creating filter for conformation of the input information
                                BitFilt = fftshift(BitFilt);                           %Doing a shift on the vector for matching the transmited data
                                %                   TxSig = ifft(fft(TxDataRes).*BitFilt);                 %Conforming the information and Creating the modulation signal
                                %U1t = TxDataRes;
                                %If the user choses to modulate all carriers at the same time there is no
                                %need of this script to generate data for each individual carrier.
                                %Therefore, this For loop can be halted
                                
                                
                                %Assigning the eletrical signal to one drive of the MZM - The aplitude of
                                %the signal can be controlled by the variable DatGai, which can be
                                %understood as an gain for the eletrical signal or an atenuation. The
                                %second signal will be similar with the only difference a phase shift of pi
                                %U.U1t = U1t;
                                U2t = exp(-1j*pi).*U1t;
                                %As both signals will have the mostrly the same characteristics with the
                                %only difference the phase shift of 180 degress. The MZM-I will be working
                                %on the Push-Pull configuration. It is necessary to reduce the Chirp noise
                                %to zero.
                                EoutTxAux = repmat(EoutTx(VetThisCarrTx==(RefCarr-1+kk),:).',1,CurTesSiz);
                                if UsingGpu==1
                                    Lgpu = gpuArray(L);
                                    U0gpu = gpuArray(U0);
                                    U_pi1gpu = gpuArray(U_pi1);
                                    U_pi2gpu = gpuArray(U_pi2);
                                    noptgpu = gpuArray(nopt);
                                    nelgpu = gpuArray(nel);
                                    Cgpu = gpuArray(C);
                                    alfa0gpu = gpuArray(alfa0);
                                    for jj=1:size(EoutTxAux,2)
                                        U1tgpu = gpuArray(U1t(:,jj));
                                        U2tgpu = gpuArray(U2t(:,jj));
                                        EoutTxGpu = gpuArray(EoutTxAux(:,jj));
                                        freqGHzGpu = gpuArray(freqGHz(:,jj));
                                        [EoutMod,~]=MZMGT(freqGHzGpu,EoutTxGpu,U1tgpu,U2tgpu,Lgpu,U0gpu,U_pi1gpu,U_pi2gpu,noptgpu,nelgpu,Cgpu,alfa0gpu);%Modulating individual carriers
                                        EoutModTem(:,jj,ObsCarrPos(kk)) = gather(EoutMod);%Adding the current Modulated output to the final OutPut
                                        clear U1tgpu U2tgpu EoutTxGpu freqGHzGpu EoutMod
                                    end
                                    clear Lgpu U0gpu U_pi1gpu U_pi2gpu noptgpu nelgpu Cgpu alfa0gpu;
                                else
                                    U.U1t = U1t;
                                    U.U2t = U2t;
                                    [EoutMod,~]=MZM(freqGHz,EoutTxAux,U,L,U0,U_pi1,U_pi2,nopt,nel,C,alfa0);%Modulating individual carriers
                                    EoutModTem(:,:,ObsCarrPos(kk)) = EoutMod;%Adding the current Modulated output to the final OutPut
                                end
                                %[EoutMod,~]=MZM(freqGHz,EoutTxAux,U,L,U0,U_pi1,U_pi2,nopt,nel,C,alfa0);
                                %EoutModTem(:,:,ObsCarrPos(kk)) = EoutMod;%Adding the current Modulated output to the final OutPut
                                if (~mod(CarrPass,2))&&(SendingUp)
                                    clear EoutModAux;
                                    EoutMod = repmat(EoutTx(VetThisCarrTx==(RefCarr-1+kk+1),:).',1,CurTesSiz);
                                    EoutModTem(:,:,ObsCarrPos(kk+1)) = EoutMod;%Adding the current Modulated output to the final OutPut
                                end
                                
                            else
                                EoutMod = repmat(EoutTx(VetThisCarrTx==(RefCarr-1+kk),:).',1,CurTesSiz);
                                EoutModTem(:,:,ObsCarrPos(kk)) = EoutMod;%Adding the current Modulated output to the final OutPut
                                if ~mod(CarrPass,2)
                                    clear EoutModAux;
                                    EoutMod = repmat(EoutTx(VetThisCarrTx==(RefCarr-1+kk+1),:).',1,CurTesSiz);
                                    EoutModTem(:,:,ObsCarrPos(kk+1)) = EoutMod;%Adding the current Modulated output to the final OutPut
                                end
                            end
                        end
                        if IfftOrSum==1
                            if size(EoutModTem,3)>1
                                if UsingGpu==1
                                    Tgpu = gpuArray(T);
                                    MaxStagGpu = gpuArray(MaxNumStagT);
                                    EoutMod = zeros(size(EoutModTem,1),size(EoutModTem,2));
                                    for jj=1:size(EoutModTem,2)
                                        EoutModTemGpu = gpuArray(EoutModTem(:,jj,:));
                                        [EoutAux1Gpu,~,~]=OpticalIFFTGT(fgpu,Tgpu,MaxStagGpu,EoutModTemGpu);
                                        EoutMod(:,jj) = gather(EoutAux1Gpu);
                                        clear EoutModTemGpu EoutAux1Gpu;
                                    end
                                    clear MaxStagGpu Tgpu;
                                else
                                    [EoutMod,~,~] = OpticalIFFTNT(f,T,MaxNumStagT,EoutModTem);
                                end
                            else
                                EoutMod = EoutModTem;
                            end
                        else
                            if size(EoutModTem,3)>1
                                EoutMod = sum(EoutModTem,3);
                            else
                                EoutMod = EoutModTem;
                            end
                        end
                end
            else
                %There are two ways, implemented, to split the incoming OFCS. One by
                %filters e another by the optical FFT. But this step just need to be done
                %once in the simulation.
                if ~exist('EoutTx','var')%Verify is this step was previouly done
                    if SelecSetUP==1%Vefify is the split process will be by filter...
                        EoutTx=SelectEachCarrier(Eout,NumCarr,f(:,1),fin,FBWD,Order,fc);%This function is resposible to split each carrier
                        VetThisCarrTx = (RefCarr-1)+1:(RefCarr-1)+NumCarr;%Keeping track of the right carrier sequence
                    else%... or by the OFFT
                        %As the OFFT has a periodic response the OFC needs to be constrained other
                        %whise carrier multiple carrier may interfir with other channels, This
                        %first selection was done with a broad pass-band filter. A higher OFFT
                        %order can also be used, although it may increase the computational time
                        %and my not be exactly feasible in the real world.
                        if Selecting==1
                            EoutA = ifft(fft(Eout).*SelecFilt(:,1));
                        else
                            EoutA = Eout;
                        end
                        if UsingGpu==1
                            EoutGpu = gpuArray(EoutA);
                            Tgpu = gpuArray(T);
                            MaxStagGpu = gpuArray(MaxNumStagTx);
                            [EoutAux1Gpu,~,VetThisCarrGpu]=OpticalFFTG(fgpu(:,1),Tgpu,MaxStagGpu,EoutGpu);
                            EoutTx = gather(EoutAux1Gpu);
                            VetThisCarrTx = gather(VetThisCarrGpu);
                            clear EoutGpu EoutAux1Gpu VetThisCarrGpu MaxStagGpu;
                        else
                            [EoutTx,~,VetThisCarrTx]=OpticalFFTN(f(:,1),T,MaxNumStagTx,EoutA);  %This function is resposible to split each carrier
                        end
                        clear EoutA;
                    end
                    %         PrintInfo(Ploting*51,EoutTx,f);                              %Printing for qualitative analizes.
                    %         axis([min(f) max(f) -400 0]);
                    %         a=a+0;
                end
            end
            %%   Transmission of the OFDM Symble through a channel
            % Having data stored and ready to be sent to end user. At the stage this
            % script is responsible to chose the medium where this signal will travel.
            % It may be withing an optical fiber or Back-toBack transmission.
            % [~,PdBm] = MeasPower(EoutMod)
            % [~,PdBm] = MeasPower(EoutMod)
            switch Medium
                case 'B2B'
                    EoutRec = EoutMod;
                case 'Fiber'
                    [EoutRec,~,~]=Fibra_Monomodo1(t,EoutMod,lambdac,T,FiberLength,0,f);
                    %         PrintInfo(Ploting*11,EoutRec.*conj(EoutRec)./max(abs(EoutRec...
                    %                                                  .*conj(EoutRec))),T,NPPB);
                otherwise
                    EoutRec = EoutMod;
            end
            
            EoutRec = EoutRec/SplitRatio;
            
        end
        if ~exist('ChanelEqualizer','var')
            if exist('OfdMod','var')
                
                
                ControlVars.FcIQ          = FcIQ;
                ControlVars.UseIQ         = UseIQ;
                ControlVars.SplitRatio    = SplitRatio;
                ControlVars.tta           = tta;
                ControlVars.ff            = ff;
                ControlVars.tt            = tt;
                ControlVars.ModFilta      = ModFilta;
                ControlVars.CurTesSiz     = CurTesSiz;
                ControlVars.fgpu          = fgpu;
                ControlVars.UsingGpu      = UsingGpu;
                ControlVars.ChanAwgn      = ChanAwgn;
                ControlVars.Tu            = Tu;
                ControlVars.ReceptorNoise = ReceptorNoise;
                ControlVars.ObsCarrUsed   = ObsCarrUsed;
                ControlVars.FFTSplit      = FFTSplit;
                ControlVars.IfftOrSum     = IfftOrSum;
                ControlVars.SendingUp     = SendingUp;
                ControlVars.ObsCarrPos    = ObsCarrPos;
                ControlVars.CarrPass      = CarrPass;
                ControlVars.RefCarr       = RefCarr;
                ControlVars.SelecSetUP    = SelecSetUP;
                ControlVars.NumCarr       = NumCarr;
                ControlVars.InitCarrUp    = InitCarrUp;
                ControlVars.InitCarrDo    = InitCarrDo;
                ControlVars.OfcName       = OfcName;
                ControlVars.alfa0         = alfa0;
                ControlVars.L             = L;
                ControlVars.U0            = U0;
                ControlVars.U_pi1         = U_pi1;
                ControlVars.U_pi2         = U_pi2;
                ControlVars.nopt          = nopt;
                ControlVars.nel           = nel;
                ControlVars.C             = C;
                ControlVars.SelModFilt    = SelModFilt;
                ControlVars.SNR           = SNR;
                ControlVars.UsingHermitian= UsingHermitian;
                ControlVars.VetSnr        = UsingHermitian;
                ControlVars.freqGHz       = freqGHz;
                ControlVars.NumFraPar     = NumFraPar;
                ControlVars.DmtMve        = DmtMve;
                ControlVars.ZtC           = ZtC;
                ControlVars.Te            = Te;
                ControlVars.OfdMod        = OfdMod;
                ControlVars.M             = M;
                ControlVars.NFFT          = NFFT;
                ControlVars.Ns            = Ns;
                ControlVars.NumFra        = NumFra;
                ControlVars.MZ_Input_File = MZ_Input_File;
                ControlVars.SelecFilt     = SelecFilt;
                ControlVars.DifNsN        = Medium;
                ControlVars.DifNsN        = DifNsN;
                ControlVars.fin           = fin;
                ControlVars.FBWD          = FBWD;
                ControlVars.Order         = Order;
                ControlVars.FibLen        = FiberLength;
                ControlVars.Ts            = Ts;
                ControlVars.OBw           = OBw;
                ControlVars.Zt            = Zt;
                ControlVars.OvSam         = OvSam;
                ControlVars.Ofc           = Ofc;
                ControlVars.Ofs           = Ofs;
                ControlVars.SelModTp      = SelModTp;
                ControlVars.NPOFEX        = NPOFEX;
                ControlVars.NPPOF         = NPPOF;
                ControlVars.NPSTUf        = NPSTUf;
                ControlVars.SelecGaus     = SelecGaus;
                switch OfdMod                                                              %Sellect which modulation was used
                    case 'qam'
                        ChanelEqualizer = 0;
                        for kk=1:TapN
                            if CurrentMedium==1
                                if (UseOfdmElect)&&(UsedModula==5)&&(NumCarr==2)
                                    ChanelEqualizer = (ChanelEqualizer + FindChanelElectricalResponse(ControlVars,FiberDelay));
                                else
                                    ChanelEqualizer = (ChanelEqualizer + FindChanelResponse(ControlVars,FiberDelay));
                                end
                            else
                                if (UseOfdmElect)&&(UsedModula==5)&&(NumCarr==2)
                                    ChanelEqualizer = (ChanelEqualizer + FindChanelElectricalResponse(ControlVars));
                                else
                                    ChanelEqualizer = (ChanelEqualizer + FindChanelResponse(ControlVars));
                                end
                            end
                        end
                        ChanelEqualizer = (ChanelEqualizer)./kk;
                    otherwise
                end
            end
        end
        %%              Datum Reception
        % The following step is to recover the information transmited
        %throughout an medium. It is the most important part of this work,
        %because it will deal with an All-Optical FFT implementation. It
        %aims to preciselly separate each sub carrier that compose the OFDM
        %Symble. It will apply methods previouly presented in the
        %literature. Looking to do all the FFT processing at the optical
        %domain.
        if (UseOfdmElect)&&(UsedModula==5)&&(NumCarr==2)
            %% DatumElectricReception%%
            
            Ix = U1t.*NormFact;
            %NormFact = max(Ix);
            %NormFact = repmat(NormFact,size(Ix,1),1);
            %Ix = Ix./NormFact;
            
            if PlotingThis
                figure;hold all;grid on;
                plot(f,20*log10(abs(fftshift(fft(Ix)./length(Ix)))));
                axis([-25e9 37.5e9 -200 0]);
                set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
            end
            Ix = Ix(1+NPSTUf:end);
            Ix = reshape(Ix,OvSam*NFFT,NumFra*CurTesSiz);
            Ix = intdump(Ix,NPPOF);                            %Downsampling the income signal
            SigRecepA = Ix(1+NPOFEX:end,:);
            
            switch SelModTp
                case 'BP'
                    %                     SigRecep  = intdump(SigRecepA,OvSam);
                    [BitFiltEle,~] = FiltroGaussiano(ff,OBw/2,0,1);        %Creating filter for selection of the received signal
                    %                     [ BitFiltEle ] = Filtro_Retangular(OBw,0,ff);
                    BitFiltEle = fftshift(BitFiltEle);
                    ModFilta = repmat(BitFiltEle.',1,NumFra);
                    SigRecepB   = SigRecepA;
                    if Ploting
                        figure;
                        hold all;
                        plot(ff,20*log10(abs(fftshift(fft(SigRecepB)./length(SigRecepB)))));
                        plot(ff,20*log10(abs(fftshift(BitFiltEle))));
                    end
                    if SelModFilt==1
                        SigRecepB = ifft(fft(SigRecepB).*ModFilta);
                    end
                    SigRecep  = intdump(SigRecepB,OvSam);
                    if Ploting
                        plot(ff,20*log10(abs(fftshift(fft(SigRecepB)./length(SigRecepB)))));
                        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                    end
                case 'AM'
                    [BitFiltEle,~] = FiltroGaussiano(ff,OBw/2,0,1);           %Creating filter for selection of the received signal
                    BitFiltEle = fftshift(BitFiltEle);
                    ModFilta = repmat(BitFiltEle.',1,NumFra);
                    SigRecepB   = SigRecepA.*cos(-2*pi*Ofc*tta);
                    if Ploting
                        figure;
                        hold all;
                        plot(ff,20*log10(abs(fftshift(fft(SigRecepB)./length(SigRecepB)))));
                        plot(ff,20*log10(abs(fftshift(BitFiltEle))));
                    end
                    if SelModFilt==1
                        SigRecepB = ifft(fft(SigRecepB).*ModFilta);
                    end
                    SigRecep  = intdump(SigRecepB,OvSam);
                    if Ploting
                        plot(ff,20*log10(abs(fftshift(fft(SigRecepB)./length(SigRecepB)))));
                        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                    end
                    
                case 'AMSSB'
                    [BitFiltEle,~] = FiltroGaussiano(ff,OBw/2,0,1);           %Creating filter for selection of the received signal
                    BitFiltEle = fftshift(BitFiltEle);
                    ModFilta = repmat(BitFiltEle.',1,NumFra);
                    SigRecepB   = SigRecepA.*exp(-1j*2*pi*Ofc*tta);
                    if Ploting
                        figure;
                        hold all;
                        plot(ff,20*log10(abs(fftshift(fft(SigRecepB)./length(SigRecepB)))));
                        plot(ff,20*log10(abs(fftshift(BitFiltEle))));
                    end
                    if SelModFilt==1
                        SigRecepB = ifft(fft(SigRecepB).*ModFilta);
                    end
                    SigRecep  = intdump(SigRecepB,OvSam);              %Over sampling
                    if Ploting
                        plot(ff,20*log10(abs(fftshift(fft(SigRecepB)./length(SigRecepB)))));
                        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                    end
                    %a=a+1;
                otherwise
                    SigRecepB = demod(SigRecepA,Ofc,Ofs,'pm');
                    SigRecep  = intdump(SigRecepB,OvSam);              %Over sampling
            end
            %             SigRecep = reshape(SigRecep,NumFra,NFFT);                      %Reshaping the signal when multiple frames were transmited
            SigRecep1 = fft(SigRecep,NFFT);                              %FFT demultiplexing process
            if UsingHermitian
                SigRecep2 = SigRecep1(2+Zt-ZtC:Ns+1+Zt-ZtC,:);                               %Taking the actual data
            else
                %     SigRecep2 = SigRecep1(1+(NFFT-size(TxSymbH,1))/2:size(TxSymbH,1)+(NFFT-size(TxSymbH,1))/2,:);
                SigRecep2 = [SigRecep1(1:size(TxSymbH,1)/2,:);SigRecep1(end+1-size(TxSymbH,1)/2:end,:)];
            end
            SigRecep3 = SigRecep2;
            %             SigRecep3 = SigRecep3(:).';                                    %Changing from parralel to serial
            %             SigRecep3 = SigRecep2.';
            if PlotingThis
                switch OfdMod                                                  %Modulation Tx signal for future comparison
                    case 'qam'
                        TxDataA = TxDataMat(1,:);
                        TxDataA = reshape(TxDataA,NumFra,length(TxDataA)/NumFra);
                        TxDataA = TxDataA.';
                        UniDmtMve = unique(DmtMve);
                        UniDmtMve = fliplr(UniDmtMve);
                        DaSiAn =  0;
                        for CarDmtM=1:length(UniDmtMve)
                            M = UniDmtMve(CarDmtM);
                            DaSiAu = sum(ismember(DmtMve,M));
                            TxSigToPlot(1+DaSiAn:DaSiAn + DaSiAu,:) = qammod(TxDataA(1+DaSiAn:DaSiAn + DaSiAu,:),M);%Saving data for future comparison
                            DaSiAn = DaSiAn + DaSiAu;
                            %TxSigToPlot(1+DaSiPo:DaSiPo + DaSiAu,:)  = qammod(TxDataA(1+DaSiAu*(CarDmtM-1):DaSiAu*CarDmtM,:),M);
                        end
                    otherwise
                        TxDataA = TxDataMat(1,:);
                        TxDataA = reshape(TxDataA,NumFra,length(TxDataA)/NumFra);
                        TxDataA = TxDataA.';
                        UniDmtMve = unique(DmtMve);
                        UniDmtMve = fliplr(UniDmtMve);
                        for CarDmtM=1:length(UniDmtMve)
                            M = UniDmtMve(CarDmtM);
                            DaSiAu = sum(ismember(DmtMve,M));
                            TxSigToPlot(1+DaSiAn:DaSiAn + DaSiAu,:) = dpskmod(TxDataA(1+DaSiAn:DaSiAn + DaSiAu,:),M);%Saving data for future comparison
                            DaSiAn = DaSiAn + DaSiAu;
                            %TxSigToPlot(1+DaSiAu*(CarDmtM-1):DaSiAu*CarDmtM,:)  = dpskmod(TxDataA(1+DaSiAu*(CarDmtM-1):DaSiAu*CarDmtM,:),M);
                        end
                end
                %         TxSigToPlot = dpskmod(TxDataMat(ThisCarr,2:end),M);
                if PlotingThis
                    figure;
                    txcolor = [0.2 0 1];
                    rxcolor = [1 0.4 0];
                    hold all;
                    plot(TxSigToPlot(:),'*','LineWidth',2,'color',txcolor);
                    plot(SigRecep3(:),'o','LineWidth',2,'color',rxcolor);
                    set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                end
            end
            %SaveRxNotEq = [SaveRxNotEq SigRecep3(:)];
            switch SelModTp
                case 'BP'
                    switch OfdMod
                        case 'qam'
                            RxSigOfdmNoEq(CurrentTest,:) = SigRecep3(:).';
                            clear SigRecep4;
                            UniDmtMve = unique(DmtMve);
                            UniDmtMve = fliplr(UniDmtMve);
                            DaSiAn = 0;
                            for CarDmtM=1:length(UniDmtMve)
                                M = UniDmtMve(CarDmtM);
                                DaSiAu = sum(ismember(DmtMve,M));
                                SigRecep4(1+DaSiAn:DaSiAn + DaSiAu,:)  = qamdemod(SigRecep3(1+DaSiAn:DaSiAn + DaSiAu,:),M);
                                DaSiAn = DaSiAn + DaSiAu;
                            end
                        otherwise
                            clear SigRecep4;
                            UniDmtMve = unique(DmtMve);
                            UniDmtMve = fliplr(UniDmtMve);
                            DaSiAn = 0;
                            for CarDmtM=1:length(UniDmtMve)
                                M = UniDmtMve(CarDmtM);
                                DaSiAu = sum(ismember(DmtMve,M));
                                SigRecep4(1+DaSiAn:DaSiAn + DaSiAu,:)  = dpskdemod(SigRecep3(1+DaSiAn:DaSiAn + DaSiAu,:),M);
                                DaSiAn = DaSiAn + DaSiAu;
                            end
                            %for CarDmtM=1:length(UniDmtMve)
                                %M = UniDmtMve(CarDmtM);
                                %DaSiAu = sum(ismember(DmtMve,M));
                                %SigRecep4(1+DaSiAu*(CarDmtM-1):DaSiAu*CarDmtM,:)  = dpskdemod(SigRecep3(1+DaSiAu*(CarDmtM-1):DaSiAu*CarDmtM,:),M);
                            %end
                    end
                case 'AMSSB'
                    switch OfdMod
                        case 'qam'
                            RxSigOfdmNoEq(CurrentTest,:) = SigRecep3(:).';
                            equa = ChanelEqualizer(1,:);
                            equ = reshape(equa,length(equa)/NumFra,NumFra);
                            if ~ChanAwgn
                                SigRecep3 = ((SigRecep3)./equ);%Fixing phase rotation and amplitude variation with the equalizer
                            end
                            clear SigRecep4;
                            UniDmtMve = unique(DmtMve);
                            UniDmtMve = fliplr(UniDmtMve);
                            DaSiAn = 0;
                            for CarDmtM=1:length(UniDmtMve)
                                M = UniDmtMve(CarDmtM);
                                DaSiAu = sum(ismember(DmtMve,M));
                                SigRecep4(1+DaSiAn:DaSiAn + DaSiAu,:)  = qamdemod(SigRecep3(1+DaSiAn:DaSiAn + DaSiAu,:),M);
                                DaSiAn = DaSiAn + DaSiAu;
                            end
                            %for CarDmtM=1:length(UniDmtMve)
                                %M = UniDmtMve(CarDmtM);
                                %DaSiAu = sum(ismember(DmtMve,M));
                                %SigRecep4(1+DaSiAu*(CarDmtM-1):DaSiAu*CarDmtM,:)  = qamdemod(SigRecep3(1+DaSiAu*(CarDmtM-1):DaSiAu*CarDmtM,:),M);
                            %end
                        otherwise
                            clear SigRecep4;
                            UniDmtMve = unique(DmtMve);
                            UniDmtMve = fliplr(UniDmtMve);
                            DaSiAn = 0;
                            for CarDmtM=1:length(UniDmtMve)
                                M = UniDmtMve(CarDmtM);
                                DaSiAu = sum(ismember(DmtMve,M));
                                SigRecep4(1+DaSiAn:DaSiAn + DaSiAu,:)  = dpskdemod(SigRecep3(1+DaSiAn:DaSiAn + DaSiAu,:),M);
                                DaSiAn = DaSiAn + DaSiAu;
                            end
                            %for CarDmtM=1:length(UniDmtMve)
                                %M = UniDmtMve(CarDmtM);
                                %DaSiAu = sum(ismember(DmtMve,M));
                                %SigRecep4(1+DaSiAu*(CarDmtM-1):DaSiAu*CarDmtM,:)  = dpskdemod(SigRecep3(1+DaSiAu*(CarDmtM-1):DaSiAu*CarDmtM,:),M);
                            %end
                    end
                otherwise
                    switch OfdMod
                        case 'qam'
                            RxSigOfdmNoEq(CurrentTest,:) = SigRecep3(:).';
                            equa = ChanelEqualizer(1,:);
                            equ = reshape(equa,length(equa)/NumFra,NumFra);
                            SigRecep3 = ((SigRecep3)./equ);%Fixing phase rotation and amplitude variation with the equalizer
                            clear SigRecep4;
                            UniDmtMve = unique(DmtMve);
                            UniDmtMve = fliplr(UniDmtMve);
                            DaSiAn = 0;
                            for CarDmtM=1:length(UniDmtMve)
                                M = UniDmtMve(CarDmtM);
                                DaSiAu = sum(ismember(DmtMve,M));
                                SigRecep4(1+DaSiAn:DaSiAn + DaSiAu,:)  = qamdemod(SigRecep3(1+DaSiAn:DaSiAn + DaSiAu,:),M);
                                DaSiAn = DaSiAn + DaSiAu;
                            end
                            %for CarDmtM=1:length(UniDmtMve)
                                %M = UniDmtMve(CarDmtM);
                                %DaSiAu = sum(ismember(DmtMve,M));
                                %SigRecep4(1+DaSiAu*(CarDmtM-1):DaSiAu*CarDmtM,:)  = qamdemod(SigRecep3(1+DaSiAu*(CarDmtM-1):DaSiAu*CarDmtM,:),M);
                            %end
                        otherwise
                            clear SigRecep4;
                            UniDmtMve = unique(DmtMve);
                            UniDmtMve = fliplr(UniDmtMve);
                            DaSiAn = 0;
                            for CarDmtM=1:length(UniDmtMve)
                                M = UniDmtMve(CarDmtM);
                                DaSiAu = sum(ismember(DmtMve,M));
                                SigRecep4(1+DaSiAn:DaSiAn + DaSiAu,:)  = dpskdemod(SigRecep3(1+DaSiAn:DaSiAn + DaSiAu,:),M);
                                DaSiAn = DaSiAn + DaSiAu;
                            end
                            %for CarDmtM=1:length(UniDmtMve)
                                %M = UniDmtMve(CarDmtM);
                                %DaSiAu = sum(ismember(DmtMve,M));
                                %SigRecep4(1+DaSiAu*(CarDmtM-1):DaSiAu*CarDmtM,:)  = dpskdemod(SigRecep3(1+DaSiAu*(CarDmtM-1):DaSiAu*CarDmtM,:),M);
                            %end
                    end
            end
            %EVMMatRec(1,:) = SigRecep3(:);%Taking just the middle samples as references
            %[EvmDB(CurrentTest,1), EvmPer(CurrentTest,1), EvmRms(CurrentTest,1) ] = ...
            %    EvmCalc( EvmMatRef(1,:),SigRecep3(:).' );
            %[EvmDBJ(CurrentTest,1),EvmPerJ(CurrentTest,1),EvmRmsJ(CurrentTest,1)] = ...
            %    evm1(M,OfdMod,EvmMatRef(1,:),SigRecep3(:).');
            %         SigRecep4 = dpskdemod(SigRecep3,M);
            SigRecep4 = SigRecep4.';
            SigRecep4 = SigRecep4(:).';
            RxSigOfdm(CurrentTest,:) = SigRecep4;
            %SaveRxEq = [SaveRxEq SigRecep4];
            if PlotingThis
                figure(111);
                txcolor = [0.2 0 1];
                rxcolor = [1 0.4 0];
                hold all;
                plot(TxSigToPlot(:),'*','LineWidth',2,'color',txcolor);
                plot(SigRecep3(:),'o','LineWidth',2,'color',rxcolor);
                set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
            end
            if PlotingThis
                figure;
                hold all;
                plot(TxDataMat(1,:));
                plot(SigRecep4);
                set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
            end
            SymbErr = ~(TxDataMat(1,:)==SigRecep4);%Measuring system bit error ration
            DmtMvep = repmat(DmtMve.',1,NumFra*CurTesSiz);
            DmtKvep = log2(DmtMvep(:).');
            %             BerOFDM(CurrentTest,ThisCarr) = sum(SymbErr)/length(SymbErr);
            BerOFDM(CurrentTest,1) = sum(SymbErr.*DmtKvep)/sum(DmtKvep);
            % BerOFDM
            a=a+0;
            close all;
        else
            if SendingDowStr==1
                %% DatumReception;
                if FFTSplit==1
                    if UsingGpu==1
                        Tgpu = gpuArray(T);
                        MaxStagGpu = gpuArray(MaxNumStag);
                        EoutAux1 = zeros(size(EoutRec,1),size(EoutRec,2),length(ObsCarrUsed));
                        for jj=1:size(EoutRec,2)
                            EoutGpu = gpuArray(EoutRec(:,jj));
                            [EoutAux1Gpu,~,VetThisCarrGpu]=OpticalFFTGT(fgpu,Tgpu,MaxStagGpu,EoutGpu);
                            EoutAux1(:,jj,:) = gather(EoutAux1Gpu);
                            VetThisCarr = gather(VetThisCarrGpu);
                            clear EoutGpu EoutAux1Gpu VetThisCarrGpu;
                        end
                        clear Tgpu MaxStagGpu;
                    else
                        [EoutAux1,~,VetThisCarr]=OpticalFFTNT(f,T,MaxNumStag,EoutRec);
                    end
                else
                    VetThisCarr = ObsCarrPos;
                    if SendingUp==1%&&(NumCarr>2)
                        EoutAux1    = SelectEachCarrier(EoutRec,NumCarr,f,fin,1.0*fc,11,fc);
                    else
                        EoutAux1(1,:)    = EoutRec;
                        EoutAux1(2,:)    = 0;
                    end
                end
                for kk=1:size(EoutAux1,3)
                    [~,CarrRecPowDo(CurrentTest,kk)] = MeasPower(EoutAux1(:,:,VetThisCarr==ObsCarrUsed(kk)),t);
                end
                if Ploting
                    figure;hold all;
                    for kk=1:size(EoutAux1,3)
                        plot(f(:,1),20*log10(abs(fftshift(fft(EoutAux1(:,1,VetThisCarr==kk)./length(EoutAux1(:,1,VetThisCarr==kk)))))));
                    end
                    a=1;
                end
                %%  Ploting some results for qualitative analizes
                % figure; plot(f,db(abs(fftshift(fft(EoutRec)./length(EoutRec))))); set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                % PrintInfo(Ploting*12,f,db(abs((fftshift(fft(EoutRec))))/length(EoutRec)),EoutAux1,RefCarr*fc);
                %%      Recovering the Data
                %This is basicaly the final step, so far, in this process. Here, the
                %transmited signal will be received and processed to recover the actual
                %data within it.
                switch Modulation
                    case 'OFDM'
                        for ThisCarr=InitCarrDo:CarrPass:NumCarr
                            Ix = EoutAux1(:,:,VetThisCarr==ObsCarrUsed(ThisCarr));
                            switch Medium
                                case 'Fiber'
                                    Ix = ifft(fft(Ix).*exp(1j*2*pi*f*(...
                                        FiberDelay(ThisCarr)*Ta)));
                                otherwise
                            end
                            %             switch OfdMod                                                %Modulation Tx signal for future comparison
                            %                 case 'qam'
                            %                     sn=((1/2)*(((1/2)^MaxNumStag)-1))/((1/2)-1);
                            %                     DiffPos = round((sn/(2/2^IfftOrSum))*T/Ta);
                            %                     EoutAux = [EoutAux(DiffPos+1:end) EoutAux(1:DiffPos)];
                            %                 otherwise
                            %             end
                            %The current incoming signal is them converted from the optical
                            %domain to the eletrical domain with the help of an photo
                            %detector.
                            Ix = Ix.*conj(Ix);
                            %             Ix = Ix.*NormFact;
                            [BitFilt,~] = FiltroGaussiano(f,BWD/2,CenFeq,FiltOrd);           %Creating filter for selection of the received signal
                            BitFilt = fftshift(BitFilt);                                   %Shifting the filter for matching the received signal
                            Ix = ifft(fft(Ix).*BitFilt);                                   %Filtering the signal for better detection
                            Ix = Ix - min(min(Ix));                                             %Removing any off-set that may exist
                            Ix = Ix.*NormFact;                                         
                            NorAux = max(Ix);
                            NorAux = repmat(NorAux,size(Ix,1),1);
                            Ix = (Ix./NorAux);%Normalizing the eletrical signal (amplifying what is needed)
                            
                            if ReceptorNoise==1                                               %Verify whether the noise is to be added or not
                                [~,SigPower] = MeasPower(Ix);
                                SigPower2 = SigPower-30;%20*log10(SigPowerI);
                                Ix = awgn(Ix,SNR,SigPower2);
                            end
                            if Ploting
                                figure;hold all;grid on;plot(f,20*log10(abs(fftshift(fft(Ix)./length(Ix)))));axis([-25e9 37.5e9 -200 0]);set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                            end
                            %Ix = reshape(Ix,Nb*NPPB,CurTesSiz);
                            Ix = Ix(1+NPSTUf:end,:);
                            Ix = reshape(Ix,OvSam*NFFT,NumFra*CurTesSiz);
                            %Ix = reshape(Ix,length(Ix)/NumFra,NumFra);
                            Ix = intdump(Ix,NPPOF);                                        %Downsampling the income signal
                            SigRecepA = Ix(1+NPOFEX:end,:);
                            if UseIQ
                                SigReal = SigRecepA.*cos(2*pi*FcIQ.*tta);
                                SigImag = SigRecepA.*sin(2*pi*FcIQ.*tta);
                                SigRecepA = SigReal + 1j.*SigImag;
                            end
                            %CpPaFr = 1;
                            %for CarrOffSet=-1*(NumFraPar-1)/2:1:(NumFraPar-1)/2
                            switch SelModTp
                                case 'BP'
                                    %                     SigRecep  = intdump(SigRecepA,OvSam);
                                    if SelecGaus==1
                                        [ BitFiltEle ] = FiltroGaussiano(ff,OBw/2,0,1);%CarrOffSet*FramSpac,1);
                                    else
                                        [ BitFiltEle ] = Filtro_Retangular(OBw,0,ff);%CarrOffSet*FramSpac,ff);
                                    end
                                    BitFiltEle = fftshift(BitFiltEle);
                                    ModFilta = repmat(BitFiltEle.',1,NumFra*CurTesSiz);
                                    SigRecepB   = SigRecepA;%.*exp(-1j*2*pi*CarrOffSet*FramSpac*tta);
                                    if Ploting
                                        figure;
                                        hold all;
                                        plot(ff,20*log10(abs(fftshift(fft(SigRecepB)./length(SigRecepB)))));
                                        plot(ff,20*log10(abs(fftshift(BitFiltEle))));
                                    end
                                    if SelModFilt==1
                                        SigRecepB = ifft(fft(SigRecepB).*ModFilta);
                                    end
                                    SigRecep  = intdump(SigRecepB,OvSam);
                                    if Ploting
                                        plot(ff,20*log10(abs(fftshift(fft(SigRecepB)./length(SigRecepB)))));
                                        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                                    end
                                case 'AM'
                                    [BitFiltEle,~] = FiltroGaussiano(ff,OBw/2,0,1);           %Creating filter for selection of the received signal
                                    BitFiltEle = fftshift(BitFiltEle);
                                    ModFilta = repmat(BitFiltEle.',1,NumFra*CurTesSiz);
                                    SigRecepB   = SigRecepA.*cos(-2*pi*Ofc*tta);
                                    if Ploting
                                        figure;
                                        hold all;
                                        plot(ff,20*log10(abs(fftshift(fft(SigRecepB)./length(SigRecepB)))));
                                        plot(ff,20*log10(abs(fftshift(BitFiltEle))));
                                    end
                                    if SelModFilt==1
                                        SigRecepB = ifft(fft(SigRecepB).*ModFilta);
                                    end
                                    SigRecep  = intdump(SigRecepB,OvSam);
                                    if Ploting
                                        plot(ff,20*log10(abs(fftshift(fft(SigRecepB)./length(SigRecepB)))));
                                        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                                    end
                                    %                         x_ofdm = x_ofdm_ssb;
                                case 'AMSSB'
                                    [BitFiltEle,~] = FiltroGaussiano(ff,OBw/2,0,1);           %Creating filter for selection of the received signal
                                    BitFiltEle = fftshift(BitFiltEle);
                                    ModFilta = repmat(BitFiltEle.',1,NumFra*CurTesSiz);
                                    SigRecepB   = SigRecepA.*exp(-1j*2*pi*Ofc*tta);
                                    if Ploting
                                        figure;
                                        hold all;
                                        plot(ff,20*log10(abs(fftshift(fft(SigRecepB)./length(SigRecepB)))));
                                        plot(ff,20*log10(abs(fftshift(BitFiltEle))));
                                    end
                                    if SelModFilt==1
                                        SigRecepB = ifft(fft(SigRecepB).*ModFilta);
                                    end
                                    SigRecep  = intdump(SigRecepB,OvSam);              %Over sampling
                                    if Ploting
                                        plot(ff,20*log10(abs(fftshift(fft(SigRecepB)./length(SigRecepB)))));
                                        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                                    end
                                otherwise
                                    SigRecepB = demod(SigRecepA,Ofc,Ofs,'pm');
                                    SigRecep  = intdump(SigRecepB,OvSam);              %Over sampling
                            end
                            %             SigRecep = reshape(SigRecep,NumFra,NFFT);                      %Reshaping the signal when multiple frames were transmited
                            SigRecep1 = fft(SigRecep,NFFT);                              %FFT demultiplexing process
                            %             SigRecep2 = SigRecep1(2+Zt-ZtC:Ns+1+Zt-ZtC,:);                               %Taking the actual data
                            if UsingHermitian
                                SigRecep2 = SigRecep1(2+Zt-ZtC:Ns+1+Zt-ZtC,:);                               %Taking the actual data
                            else
                                %     SigRecep2 = SigRecep1(1+(NFFT-size(TxSymbH,1))/2:size(TxSymbH,1)+(NFFT-size(TxSymbH,1))/2,:);
                                SigRecep2 = [SigRecep1(1:size(TxSymbH,1)/2,:);SigRecep1(end+1-size(TxSymbH,1)/2:end,:)];
                            end
                            SigRecep3 = SigRecep2;
                            %             SigRecep3 = SigRecep3(:).';                                    %Changing from parralel to serial
                            %             SigRecep3 = SigRecep2.';
                            if PlotingThis
                                switch OfdMod                                                  %Modulation Tx signal for future comparison
                                    case 'qam'
                                        TxDataA = TxDataMat(ThisCarr,:);
                                        TxDataA = reshape(TxDataA,NumFra*CurTesSiz,length(TxDataA)/(NumFra*CurTesSiz));
                                        TxDataA = TxDataA.';
                                        UniDmtMve = unique(DmtMve);
                                        UniDmtMve = fliplr(UniDmtMve);
                                        DaSiAn =  0;
                                        for CarDmtM=1:length(UniDmtMve)
                                            M = UniDmtMve(CarDmtM);
                                            DaSiAu = sum(ismember(DmtMve,M));
                                            TxSigToPlot(1+DaSiAn:DaSiAn + DaSiAu,:) = qammod(TxDataA(1+DaSiAn:DaSiAn + DaSiAu,:),M);%Saving data for future comparison
                                            DaSiAn = DaSiAn + DaSiAu;
                                            %TxSigToPlot(1+DaSiAu*(CarDmtM-1):DaSiAu*CarDmtM,:)  = qammod(TxDataA(1+DaSiAu*(CarDmtM-1):DaSiAu*CarDmtM,:),M);
                                        end
                                    otherwise
                                        TxDataA = TxDataMat(ThisCarr,:);
                                        TxDataA = reshape(TxDataA,NumFra*CurTesSiz,length(TxDataA)/(NumFra*CurTesSiz));
                                        TxDataA = TxDataA.';
                                        UniDmtMve = unique(DmtMve);
                                        UniDmtMve = fliplr(UniDmtMve);
                                        DaSiAn =  0;
                                        for CarDmtM=1:length(UniDmtMve)
                                            M = UniDmtMve(CarDmtM);
                                            DaSiAu = sum(ismember(DmtMve,M));
                                            TxSigToPlot(1+DaSiAn:DaSiAn + DaSiAu,:) = dpskmod(TxDataA(1+DaSiAn:DaSiAn + DaSiAu,:),M);%Saving data for future comparison
                                            DaSiAn = DaSiAn + DaSiAu;
                                            %TxSigToPlot(1+DaSiAu*(CarDmtM-1):DaSiAu*CarDmtM,:)  = dpskmod(TxDataA(1+DaSiAu*(CarDmtM-1):DaSiAu*CarDmtM,:),M);
                                        end
                                end
                                %         TxSigToPlot = dpskmod(TxDataMat(ThisCarr,2:end),M);
                                if PlotingThis
                                    figure;
                                    txcolor = [0.2 0 1];
                                    rxcolor = [1 0.4 0];
                                    hold all;
                                    plot(SigRecep3(:),'o','LineWidth',2,'color',rxcolor);
                                    plot(TxSigToPlot(:),'*','LineWidth',2,'color',txcolor);
                                    set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                                end
                            end
                            
                            switch OfdMod
                                case 'qam'
                                    RxSigOfdmNoEq(CurrentTest,ThisCarr,:) = SigRecep3(:).';
                                    equa = ChanelEqualizer(ObsCarrUsed(ThisCarr),:);
                                    equ = reshape(equa,length(equa)/(NumFra*CurTesSiz),NumFra*CurTesSiz);
                                    SigRecep3 = ((SigRecep3)./equ);%Fixing phase rotation and amplitude variation with the equalizer
                                    clear SigRecep4;
                                    UniDmtMve = unique(DmtMve);
                                    UniDmtMve = fliplr(UniDmtMve);
                                    DaSiAn = 0;
                                    for CarDmtM=1:length(UniDmtMve)
                                        M = UniDmtMve(CarDmtM);
                                        DaSiAu = sum(ismember(DmtMve,M));
                                        SigRecep4(1+DaSiAn:DaSiAn + DaSiAu,:)  = qamdemod(SigRecep3(1+DaSiAn:DaSiAn + DaSiAu,:),M);
                                        DaSiAn = DaSiAn + DaSiAu;
                                    end
                                    %for CarDmtM=1:length(UniDmtMve)
                                        %M = UniDmtMve(CarDmtM);
                                        %DaSiAu = sum(ismember(DmtMve,M));
                                        %SigRecep4(1+DaSiAu*(CarDmtM-1):DaSiAu*CarDmtM,:)  = qamdemod(SigRecep3(1+DaSiAu*(CarDmtM-1):DaSiAu*CarDmtM,:),M);
                                    %end
                                    %                     SigRecep4  = qamdemod(SigRecep3,M);                    %Modulating information
                                otherwise
                                    clear SigRecep4;
                                    UniDmtMve = unique(DmtMve);
                                    UniDmtMve = fliplr(UniDmtMve);
                                    DaSiAn = 0;
                                    for CarDmtM=1:length(UniDmtMve)
                                        M = UniDmtMve(CarDmtM);
                                        DaSiAu = sum(ismember(DmtMve,M));
                                        SigRecep4(1+DaSiAn:DaSiAn + DaSiAu,:)  = dpskdemod(SigRecep3(1+DaSiAn:DaSiAn + DaSiAu,:),M);
                                        DaSiAn = DaSiAn + DaSiAu;
                                    end
                                    %for CarDmtM=1:length(UniDmtMve)
                                        %M = UniDmtMve(CarDmtM);
                                        %DaSiAu = sum(ismember(DmtMve,M));
                                        %SigRecep4(1+DaSiAu*(CarDmtM-1):DaSiAu*CarDmtM,:)  = dpskdemod(SigRecep3(1+DaSiAu*(CarDmtM-1):DaSiAu*CarDmtM,:),M);
                                    %end                             %Serializing the signal
                            end
                            SigRecep4 = SigRecep4.';
                            SigRecep4 = SigRecep4(:).';
                            RxSigOfdm(CurrentTest,ThisCarr,:) = SigRecep4;
                            if PlotingThis
                                figure;
                                txcolor = [0.2 0 1];
                                rxcolor = [1 0.4 0];
                                hold all;
                                plot(SigRecep3(:),'o','LineWidth',2,'color',rxcolor);
                                plot(TxSigToPlot(:),'*','LineWidth',2,'color',txcolor);
                                set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                            end
                            if PlotingThis
                                figure;
                                hold all;
                                plot(TxDataMat(ThisCarr,:));
                                plot(SigRecep4);
                                SampPos = 1:length(SigRecep4);
                                plot(SampPos(~(TxDataMat(ThisCarr,:)==SigRecep4)),...
                                    SigRecep4(~(TxDataMat(ThisCarr,:)==SigRecep4)),'ko','MarkerFaceColor',[0 0 0]);
                                set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                            end
                            SymbErr = ~(TxDataMat(ThisCarr,:)==SigRecep4);%Measuring system bit error ration
                            %CpPaFr = CpPaFr + 1;
                            %end
                            DmtMvep = repmat(DmtMve,NumFra*CurTesSiz,1);
                            DmtKvep = log2(DmtMvep(:).');
                            BerToPlotOfdm(CurrentTest,ThisCarr,:) = SymbErr.*DmtKvep;
                            %             BerOFDM(CurrentTest,ThisCarr) = sum(SymbErr)/length(SymbErr);
                            BerOFDM(CurrentTest,ThisCarr) = sum(SymbErr.*DmtKvep)/sum(DmtKvep);
                            close all;
                        end
                    case 'DPSK'
                        %%                   Receiver DPSK
                        for ThisCarr=InitCarrDo:CarrPass:NumCarr                                    %For each carrier the same process of reception will be used.
                            if CarrUsedDo(ThisCarr)
                                %%  Reception Process: Handling the income optical field
                                %At the first step the income field will be selected from the
                                %optical FFT Output with the help of the maping vector
                                %(previously described).
                                Ix = EoutAux1(:,:,VetThisCarr==ObsCarrUsed(ThisCarr));
                                %%            Fiber Time Delay Compensation
                                switch Medium
                                    case 'Fiber'
                                        Ix = ifft(fft(Ix).*exp(1j*2*pi*f*(...
                                            FiberDelay(ThisCarr)*Ta)));
                                    otherwise
                                end
                                %             EoutAux = EoutRec;
                                %%  Plot for Qualitative analizes
                                %PrintInfo(Ploting*13,t,abs(EoutAux));
                                %             PrintInfo(Ploting*14,f,20*log10(abs(fftshift(fft(EoutAux)./...
                                %                                                        length(EoutAux)))));
                                
                                %%                   synchronizing
                                %For the reception process work properly, it is needed to
                                %sincronized the recieved signal or the sampling process will
                                %not work.
                                
                                %At this first moment a small part of the income signal needs
                                %be analized for the sincronism process. As we are looking for
                                %the actual synchronis symbol, which is the bigest phase shift
                                %of the income signal. The next delay interferometer convert
                                %this phase shift to a amplitude variation.
                                
                                %The phase delay is not needed because the optical field will
                                %be analized as a whole as it is composed only by the real part
                                PhaDel          = 0;
                                %Remember that at DPSK modulation the information is stored at
                                %the phase difference of the current symbol with the previous
                                %one hence this time delay is needed to analyse each symbol by
                                %analizes of its interaction of the current symbel with the
                                %previous one.
                                TimDel          = T;
                                if UsingGpu==1
                                    TimeDelGpu = gpuArray(TimDel);
                                    PhaDelGpu  = gpuArray(PhaDel);
                                    ESync1 = zeros(size(Ix,1),size(Ix,2));
                                    ESync2 = zeros(size(Ix,1),size(Ix,2));
                                    for jj=1:size(Ix,2)
                                        IxGpu = gpuArray(Ix(:,jj));
                                        [ESync1Gpu,ESync2Gpu] = DelayInterfG(fgpu,TimeDelGpu,PhaDelGpu,IxGpu);
                                        ESync1(:,jj) = gather(ESync1Gpu);
                                        ESync2(:,jj) = gather(ESync2Gpu);
                                        clear IxGpu ESync1Gpu ESync2Gpu;
                                    end
                                    clear TimeDelGpu PhaDelGpu;
                                else
                                    [ESync1,ESync2] = DelayInterfExp(f,TimDel,PhaDel,Ix);        %Coverting phase shift to amplitude variation
                                end
                                %The second moment is to transfer this information from the
                                %optical domain to the eletrical domain for an eletronic
                                %processing. It is done with the help of an photo diode. The
                                %configuration here used is an balanced receiver as the output
                                %of the Delay Interferometer has two signals resulting from
                                %the constructive and destructive signal interaction.
                                ESync1 = ESync1.*conj(ESync1);
                                %             ESync1 = ifft(fft(ESync1)./fft(PulseResp(VetThisCarr==ThisCarr,:)));
                                ESync2 = ESync2.*conj(ESync2);
                                %             ESync2 = ifft(fft(ESync2)./fft(PulseResp(VetThisCarr==ThisCarr,:)));
                                
                                %%           Creating the Reception Filter
                                
                                
                                [BitFilt,~] = FiltroGaussiano(f,BWD2,CenFeq2,FiltOrd2);%Creating filter for selection of the received signal
                                BitFilt = fftshift(BitFilt);%Shifting the filter for matching the received signal
                                
                                %%
                                Esync  = ESync2-ESync1;
                                
                                if RecFilBanPas==1
                                    Esync  = ifft(fft(Esync).*BitFilt);%Filter is used to remove higher order components
                                    Emean  = mean(Esync);
                                    Emean  = repmat(Emean,size(Esync,1),1);
                                    Esync  = Esync-Emean;
                                    Enorm  = max(Esync);
                                    Enorm  = repmat(Enorm,size(Esync,1),1);
                                    Esync  = Esync./Enorm;
                                end
                                
                                %%         Adding Noise to Received Signal
                                %Just after the photo diod is where the noise must be added as
                                %it is the main constrain of the system. Once the signal is
                                %converted from optical to electrical and it digital values
                                %were recevered, there is no point to add noise, because the
                                %information was already received, it just needs decoding.
                                %The noise here added will be gaussian basically it represents
                                %the reception sensibility. In another words, how much the
                                %signal must be about the noise for an acceptable recovering.
                                
                                if ReceptorNoise==1%Verify whether the noise is to be added or not
                                    [~,SigPower] = MeasPower(Esync);
                                    SigPower2 = SigPower-30;%20*log10(SigPowerI);
                                    SNR2 = CarSNR + 10*log10(1) - 10*log10(Nsamp) + 10*log10(10^0.36);
                                    Esync = awgn(Esync,SNR,SigPower2);
                                end
                                if ~RecFilBanPas
                                    Esync  = ifft(fft(Esync).*BitFilt);%Filter is used to remove higher order components
                                    Emean  = mean(Esync);
                                    Emean  = repmat(Emean,size(Esync,1),1);
                                    Esync  = Esync-Emean;
                                    Enorm  = max(Esync);
                                    Enorm  = repmat(Enorm,size(Esync,1),1);
                                    Esync  = Esync./Enorm;
                                end
                                
                                %EoutAuxF = 20*log10(abs(fftshift(fft(Ix)./length(Ix))));
                                %for jj=1:size(Ix,2)
                                    %VetElecPower(CurrentTest,ThisCarr,jj)= MeasPower(Ix(:,jj));
                                    %[VetOptiPower(CurrentTest,ThisCarr,jj),~]= findpeaks(EoutAuxF(:,jj),'SortStr','descend','NPeaks',1);
                                %end
								if UsingMex == 1
									[BerDPSK(CurrentTest,ThisCarr),ValsLev(CurrentTest,ThisCarr),AberLev(CurrentTest,ThisCarr),Esync,InterAB,LocMaxAB,MaxValAB,SeqFinAB,TxData,Data] = RecDowDepsk_mex(MaxNumStag,IfftOrSum,T,Ta,Esync,AddCP,StuffSampels,NumAmosCP,NPPB,NbDPSK,CurTesSiz,SyncPeriod,ThisCarr,TxDataMat,LevDefDpqsk);
								else
									[BerDPSK(CurrentTest,ThisCarr),ValsLev(CurrentTest,ThisCarr),AberLev(CurrentTest,ThisCarr),Esync,InterAB,LocMaxAB,MaxValAB,SeqFinAB,TxData,Data] = RecDowDepsk(MaxNumStag,IfftOrSum,T,Ta,Esync,AddCP,StuffSampels,NumAmosCP,NPPB,NbDPSK,CurTesSiz,SyncPeriod,ThisCarr,TxDataMat,LevDefDpqsk);
								end
                                if CalcS==1
									if UsingMex == 1
										[~,~,EyeOpenI,EyeOpenHighI,EyeOpenLowI] = Olho_mex(Esync,T,NPPB,0,1);
									else
										[~,~,EyeOpenI,EyeOpenHighI,EyeOpenLowI] = Olho(Esync,T,NPPB,0,1);
									end
                                    AberLevS  = EyeOpenI;
                                    ValsLevS  = EyeOpenLowI + EyeOpenI/2;          %Comparison between the Transmited and received and counting the differences
                                    ThisDataSize = NPPB/2:NPPB:size(Esync,2);
                                    DataS = zeros(1,length(ThisDataSize));                                                  %Initialization of the vector that will store the income data
                                    for kk=1:length(ThisDataSize)%length(Esync(ThisDataSize))                                    %The comparison process will be made for each symbol period
                                        %MeanOfData = mean(EoutI((kk-1)+SymLocI(1)));
                                        MeanOfData = mean(Esync((kk-1)*NPPB+NPPB/2));
                                        if MeanOfData > EyeOpenLowI + EyeOpenI/2%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                                            DataS(kk) = 1;                                 %is 1
                                        else                                                       %If it is the lowest level the incoming data
                                            DataS(kk) = 0;                                 %is 0
                                        end
                                    end
                                    %       Calculating the Bit Error Ratio (BER)
                                    %The final process here is to count the number of wrongdoings
                                    %of this whole process upon the transmited data for
                                    %quantitative analizes

                                    TxDataA = TxDataMat(ThisCarr,:);
                                    TxDataB = reshape(TxDataA,NbDPSK,CurTesSiz);
                                    TxDataC = TxDataB(1+SyncPeriod:end-SyncPeriod,:);
                                    TxData  = reshape(TxDataC,1,size(TxDataC,1)*size(TxDataC,2));
%                                     if CalcS
                                    %DataS = DataS(1+SyncPeriod:end-SyncPeriod);
                                    BitErrS = sum(xor(TxData,DataS));%Comparison between the Transmited and received and counting the differences
                                    BerDPSKS(CurrentTest,ThisCarr) = BitErrS/((NbDPSK-(2*SyncPeriod))*CurTesSiz);%Calculating the ration of wrong bits and the total number of bits transmited
%                                     end
                                end
                                %% Ploting the result for qualitative analizes
                                if PrintinEye==1
                                    PrintInfo(PlotingThis*17,Esync,T,NPPB);
                                    hold on;
                                    plot(t((NPPB)/2),InterAB(round(SeqFinAB(LocMaxAB)-MaxValAB)),'k^');
                                    plot(t((NPPB)/2),InterAB(SeqFinAB(LocMaxAB)),'kv');
                                    plot(t((NPPB)/2),InterAB(round(SeqFinAB(LocMaxAB)-MaxValAB/2)),'k*');
                                    if CalcS==1
                                        plot(t((NPPB)/2),EyeOpenLowI,'mx');
                                        plot(t((NPPB)/2),EyeOpenHighI,'mo');
                                        plot(t((NPPB)/2),EyeOpenLowI + EyeOpenI/2,'md');
                                    end
                                    set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                                    %figure;
                                    %hold all;
                                    %plotpos = zeros(1,length(IxAux1));
                                    %plot(IxAux1,plotpos,'o','color',[1 0.4 0]);
                                    %plot(EvmMatRef(ObsCarrUsed(ThisCarr),:),plotpos,'b*');
                                    %set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                                    %a=a+0;
                                end
                                if PlotingThis
                                    figure;
                                    hold all;
                                    plot(TxData);
                                    plot(Data);
                                    set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                                    display(BerDPSK);
                                end
                                close all;
                                %                         BerDPSK(CurrentTest,ThisCarr+1) = 1;
                            end
                        end
                    case 'DQPSK'
                        
                        %%                   Receiver DQPSK
                        for ThisCarr=InitCarrDo:CarrPass:NumCarr                                           %For each carrier the same process of reception will be used.
                            if CarrUsedDo(ThisCarr)
                                %%  Reception Process: Handling the income optical field
                                %At the first step the income field will be selected from the
                                %optical FFT Output with the help of the maping vector
                                %(previously described).
                                Ix = EoutAux1(:,:,VetThisCarr==ObsCarrUsed(ThisCarr));
                                %%            Fiber Time Delay Compensation
                                switch Medium
                                    case 'Fiber'
                                        Ix = ifft(fft(Ix).*exp(1j*2*pi*f*(FiberDelay(ThisCarr)*Ta)));
                                    otherwise
                                end
                                %             EoutAux = EoutRec;
                                %%  Plot for Qualitative analizes
                                %PrintInfo(Ploting*20,t,abs(EoutAux));
                                %             PrintInfo(Ploting*21,f,20*log10(abs(fftshift(fft(EoutAux)./...
                                %                                                        length(EoutAux)))));
                                
                                %%                   synchronizing
                                %%               Recovering the Information
                                %Once the Income signal was synchronized it is possible to
                                %split its components in phase and in quadrature. It is very
                                %important to understand here that those components will carry
                                %the actual data transmited. From it we can't recover the  I
                                %and Q eletrical signals that were used to modulate our
                                %carrier. Those signals get mingled together to composed the
                                %actual data. Therefore, when we pass the income signal
                                %throughout the receptor we will recover the TxData generated
                                %and transmited. It is important that it get very clear here
                                %because I toke two months to understand it, no one could
                                %explain it well to me and those articles about it make this
                                %information unclear. Thus, to sheed a light on this concept,
                                %it is important to know that the Data is encoded in the I and
                                %Q eletrical components that will get mixed together to
                                %modulate the optical carrier. When this optical signal passes
                                %throughout the MZ-Interferometer the acutal data (TX) will be
                                %recover, not the I and Q component. But the data at the odd
                                %possition (encoded in I) will be at the real component (in
                                %phase), whereas the data at the even possition (encoded in Q)
                                %will be at the imaginary component(in quadrature). The
                                %MZ-Interferometer will be resposible to take apart the real
                                %and imaginary components of the income optical field.
                                
                                %                                 E_rec3 = Ix./max(abs(Ix));                           %Normalizing income signal
                                %For the interferometric process  take in account just the real
                                %component it is needed a phase delay of 45° degrees;
                                PhaDel = 1*pi/4;
                                %Remember that at DQPSK modulation the information is stored at
                                %the difference of phase between the current and previous
                                %symbol hence this time delay of one symbol period is needed.
                                TimDel = T;
                                if UsingGpu==1
                                    TimeDelGpu = gpuArray(TimDel);
                                    PhaDelGpu  = gpuArray(PhaDel);
                                    EoutA = zeros(size(Ix,1),size(Ix,2));
                                    EoutB = zeros(size(Ix,1),size(Ix,2));
                                    for jj=1:size(Ix,2)
                                        IxGpu = gpuArray(Ix(:,jj));
                                        [ESync1Gpu,ESync2Gpu] = DelayInterfG(fgpu,TimeDelGpu,PhaDelGpu,IxGpu);
                                        EoutA(:,jj) = gather(ESync1Gpu);
                                        EoutB(:,jj) = gather(ESync2Gpu);
                                        clear IxGpu ESync1Gpu ESync2Gpu;
                                    end
                                    clear TimeDelGpu PhaDelGpu;
                                else
                                    [EoutA,EoutB] = DelayInterfExp(f,TimDel,PhaDel,Ix);        %Coverting phase shift to amplitude variation
                                end
                                %[EoutA,EoutB] = DelayInterfExp(t,TimDel,PhaDel,Ix);          %Coverting phase shift to amplitude variation
                                %For the interferometric process  take in account just the real
                                %component it is needed a phase delay of -45° degrees;
                                PhaDel = -1*pi/4;
                                %Remember that at DQPSK modulation the information is stored at
                                %the difference of phase between the current and previous
                                %symbol hence this time delay of one symbol period is needed.
                                TimDel = T;
                                if UsingGpu==1
                                    TimeDelGpu = gpuArray(TimDel);
                                    PhaDelGpu  = gpuArray(PhaDel);
                                    EoutC = zeros(size(Ix,1),size(Ix,2));
                                    EoutD = zeros(size(Ix,1),size(Ix,2));
                                    for jj=1:size(Ix,2)
                                        IxGpu = gpuArray(Ix(:,jj));
                                        [ESync1Gpu,ESync2Gpu] = DelayInterfG(fgpu,TimeDelGpu,PhaDelGpu,IxGpu);
                                        EoutC(:,jj) = gather(ESync1Gpu);
                                        EoutD(:,jj) = gather(ESync2Gpu);
                                        clear IxGpu ESync1Gpu ESync2Gpu;
                                    end
                                    clear TimeDelGpu PhaDelGpu;
                                else
                                    [EoutC,EoutD] = DelayInterfExp(f,TimDel,PhaDel,Ix);        %Coverting phase shift to amplitude variation
                                end
                                %[EoutC,EoutD] = DelayInterfExp(t,TimDel,PhaDel,Ix);          %Coverting phase shift to amplitude variation
                                
                                %The second moment is to transfer this information from the
                                %optical domain to the eletrical domain for an eletronic
                                %processing. It is done with the help of an photo diode.
                                EoutA = EoutA.*conj(EoutA);
                                EoutB = EoutB.*conj(EoutB);
                                EoutC = EoutC.*conj(EoutC);
                                EoutD = EoutD.*conj(EoutD);
                                
                                %The process with the photo diode is self-coherent, which means
                                %the rusult will be a component of the signal centered in f=0
                                %and another component centered at f=2*fc (frenquecy central).
                                %Therefore, to remove the higher order component a low pass
                                %filter will be used.
                                %%           Creating the Reception Filter
                                [BitFilt,~] = FiltroGaussiano(f,BWD2,CenFeq2,FiltOrd2);        %Creating filter for selection of the received signal
                                BitFilt = fftshift(BitFilt);                                   %Shifting the filter for matching the received signal
                                %%
                                
                                if RecFilBanPas
                                    EoutA = ifft(fft(EoutA).*BitFilt);
                                    EoutB = ifft(fft(EoutB).*BitFilt);
                                    EoutC = ifft(fft(EoutC).*BitFilt);
                                    EoutD = ifft(fft(EoutD).*BitFilt);
                                end
                                
                                % The configuration here used is an balanced receiver as the
                                %output of the Delay Interferometer has two signals resulting
                                %from the constructive and destructive signal interaction.
                                %%         Adding Noise to Received Signal
                                %Just after the photo diod is where the noise must be added as
                                %it is the main constrain of the system. Once the signal is
                                %converted from optical to electrical and it digital values
                                %were recevered, there is no point to add noise, because the
                                %information was already received, it just needs decoding.
                                %The noise here added will be gaussian basically it represents
                                %the reception sensibility. In another words, how much the
                                %signal must be about the noise for an acceptable recovering.
                                
                                if ReceptorNoise==1                                               %Verify whether the noise is to be added or not
                                    [~,SigPowerA] = MeasPower(EoutA);                          %it own received signal or
                                    [~,SigPowerB] = MeasPower(EoutB);
                                    [~,SigPowerC] = MeasPower(EoutC);                          %it own received signal or
                                    [~,SigPowerD] = MeasPower(EoutD);
                                    SNR2 = CarSNR + 10*log10(2) - 10*log10(Nsamp) + 10*log10(10^0.3);
                                    SigPowerA2 = SigPowerA-30;%20*log10(SigPowerI);
                                    SigPowerB2 = SigPowerB-30;%20*log10(SigPowerQ);
                                    SigPowerC2 = SigPowerC-30;%20*log10(SigPowerI);
                                    SigPowerD2 = SigPowerD-30;%20*log10(SigPowerQ);
                                    EoutA = awgn(EoutA,SNR,SigPowerA2);
                                    EoutB = awgn(EoutB,SNR,SigPowerB2);
                                    EoutC = awgn(EoutC,SNR,SigPowerC2);
                                    EoutD = awgn(EoutD,SNR,SigPowerD2);
                                end
                                
                                
                                if ~RecFilBanPas
                                    EoutA = ifft(fft(EoutA).*BitFilt);
                                    EoutB = ifft(fft(EoutB).*BitFilt);
                                    EoutC = ifft(fft(EoutC).*BitFilt);
                                    EoutD = ifft(fft(EoutD).*BitFilt);
                                end
								if UsingMex==1
									[AberLevI(CurrentTest,ThisCarr),ValsLevI(CurrentTest,ThisCarr),AberLevQ(CurrentTest,ThisCarr),ValsLevQ(CurrentTest,ThisCarr),BerDQPSK(CurrentTest,ThisCarr),EoutI,EoutQ,TxDataOdd,TxDataEven,DataOdd,DataEven,LocMaxAB,MaxValAB,SeqFinAB,QMaxValAB,QLocMaxAB,QSeqFinAB,QInterAB,InterAB,IxAuxAB,QIxAuxAB]=RecDowDqpsk_mex(EoutA,EoutB,EoutC,EoutD,T,Ta,IfftOrSum,MaxNumStag,StuffSampels,NbDQPSK,CurTesSiz,NumAmosCP,NPPB,SyncPeriod,TxDataMat,ThisCarr,LevDefDpqsk,AddCP,NumCarr);
                                else
									[AberLevI(CurrentTest,ThisCarr),ValsLevI(CurrentTest,ThisCarr),AberLevQ(CurrentTest,ThisCarr),ValsLevQ(CurrentTest,ThisCarr),BerDQPSK(CurrentTest,ThisCarr),EoutI,EoutQ,TxDataOdd,TxDataEven,DataOdd,DataEven,LocMaxAB,MaxValAB,SeqFinAB,QMaxValAB,QLocMaxAB,QSeqFinAB,QInterAB,InterAB,IxAuxAB,QIxAuxAB]=RecDowDqpsk(EoutA,EoutB,EoutC,EoutD,T,Ta,IfftOrSum,MaxNumStag,StuffSampels,NbDQPSK,CurTesSiz,NumAmosCP,NPPB,SyncPeriod,TxDataMat,ThisCarr,LevDefDpqsk,AddCP,NumCarr);
								end
								%%
                                %             berpos = 1:2:size(BerDQPSK,2);
                                %             BerDQPSK(size(BerDQPSK,1),berpos)
                                if CalcS==1
                                    if UsingMex==1
                                        [~,~,EyeOpenI,EyeOpenHighI,EyeOpenLowI] = Olho_mex(EoutI,T,NPPB,0,1);
                                        [~,~,EyeOpenQ,EyeOpenHighQ,EyeOpenLowQ] = Olho_mex(EoutQ,T,NPPB,0,1);
                                    else
                                        [~,~,EyeOpenI,EyeOpenHighI,EyeOpenLowI] = Olho(EoutI,T,NPPB,0,1);
                                        [~,~,EyeOpenQ,EyeOpenHighQ,EyeOpenLowQ] = Olho(EoutQ,T,NPPB,0,1);
                                    end
									AberLevIS(CurrentTest,ThisCarr)  = EyeOpenI;
                                    ValsLevIS(CurrentTest,ThisCarr)  = EyeOpenLowI + EyeOpenI/2;
                                    AberLevQS(CurrentTest,ThisCarr)  = EyeOpenQ;
                                    ValsLevQS(CurrentTest,ThisCarr)  = EyeOpenLowQ + EyeOpenQ/2;
                                    
                                    ThisDataSize = NPPB/2:NPPB:length(EoutI);
                                    DataOddS = zeros(1,length(ThisDataSize));                                                  %Initialization of the vector that will store the income data
                                    for kk=1:length(ThisDataSize)%length(EoutI(ThisDataSize))                                    %The comparison process will be made for each symbol period
                                        %MeanOfData = mean(EoutI((kk-1)+SymLocI(1)));
                                        MeanOfData = mean(EoutI((kk-1)*NPPB+NPPB/2));
                                        if MeanOfData > EyeOpenLowI + EyeOpenI/2%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                                            DataOddS(kk) = 1;                                 %is 1
                                        else                                                       %If it is the lowest level the incoming data
                                            DataOddS(kk) = 0;                                 %is 0
                                        end
                                        
                                    end
                                    ThisDataSize = NPPB/2:NPPB:length(EoutQ);
                                    DataEvenS = zeros(1,length(ThisDataSize));                                                 %Initialization of the vector that will store the income data
                                    for kk=1:length(ThisDataSize)%length(EoutQ(ThisDataSize))                                    %The comparison process will be made for each symbol period
                                        %                 MeanOfData = mean(EoutQ((kk-1)+SymLocI(1)));
                                        MeanOfData = mean(EoutQ((kk-1)*NPPB+NPPB/2));
                                        if MeanOfData > EyeOpenLowQ + EyeOpenQ/2%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                                            DataEvenS(kk) = 1;                                 %is 1
                                        else                                                       %If it is the lowest level the incoming data
                                            DataEvenS(kk) = 0;                                 %is 0
                                        end
                                        
                                    end
                                    
                                    BitErrOddS    = sum(xor(TxDataOdd,DataOddS));                      %Comparison between the Transmited and received and counting the differences
                                    BitErrEvenS   = sum(xor(TxDataEven,DataEvenS));                    %Comparison between the Transmited and received and counting the differences
                                    
                                    BerDQPSKS(CurrentTest,ThisCarr) = (BitErrOddS+BitErrEvenS)/...
                                        ((NbDQPSK)-(2*SyncPeriod));%Calculating the ration of wrong bits and the total number of bits transmited
                                end
                                if PrintinEye==1
                                    PrintInfo(PlotingThis*24,EoutI,T,NPPB);
                                    hold on;
                                    plot(t((NPPB)/2),InterAB(round(SeqFinAB(LocMaxAB)-MaxValAB)),'k^');
                                    plot(t((NPPB)/2),InterAB(SeqFinAB(LocMaxAB)),'kv');
                                    plot(t((NPPB)/2),InterAB(round(SeqFinAB(LocMaxAB)-MaxValAB/2)),'kd');
                                    %PrintInfo(Ploting*25,Ui,T,NPPB);
                                    PrintInfo(Ploting*26,EoutQ,T,NPPB);
                                    hold on;
                                    plot(t((NPPB)/2),QInterAB(round(QSeqFinAB(QLocMaxAB)-QMaxValAB)),'k^');
                                    plot(t((NPPB)/2),QInterAB(QSeqFinAB(QLocMaxAB)),'kv');
                                    plot(t((NPPB)/2),QInterAB(round(QSeqFinAB(QLocMaxAB)-QMaxValAB/2)),'kd');
                                    if CalcS
                                        plot(t((NPPB)/2),EyeOpenLowI,'mx');
                                        plot(t((NPPB)/2),EyeOpenHighI,'mo');
                                        plot(t((NPPB)/2),EyeOpenLowI + EyeOpenI/2,'md');
                                        plot(t((NPPB)/2),EyeOpenLowQ,'cx');
                                        plot(t((NPPB)/2),EyeOpenHighQ,'co');
                                        plot(t((NPPB)/2),EyeOpenLowQ + EyeOpenQ/2,'cd');
                                    end
                                    set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                                    figure;
                                    hold all;
                                    plotpos = zeros(1,length(TxDataOdd));
                                    plot(IxAuxAB,QIxAuxAB,'o','color',[1 0.4 0]);
                                    TxDataOddP = TxDataOdd;
                                    TxDataOddP(TxDataOddP==0) = -1;
                                    TxDataEvenP = TxDataEven;
                                    TxDataEvenP(TxDataEvenP==0) = -1;
                                    plot(TxDataOddP,TxDataEvenP,'b*');
                                    set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                                    %a=a+0;
                                end
                                %             a=a+6;
                                if PlotingThis
                                    figure;
                                    hold all;
                                    plot([TxDataOdd TxDataEven]);
                                    plot([DataOdd DataEven]);
                                    set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                                    display(BerDQPSK);
                                end
                                close all;
                                %                 BerDQPSK(CurrentTest,ThisCarr+1) = 1;
                            end
                        end
                    case '4PAM'
                        
                        %%              Receiver 4PAM
                        for ThisCarr=InitCarrDo:CarrPass:NumCarr                                           %For each carrier the same process of reception will be used.
                            if CarrUsedDo(ThisCarr)
                                %%  Reception Process: Handling the income optical field
                                %At the first step the income field will be selected from the
                                %optical FFT Output with the help of the maping vector
                                %(previously described).
                                Ix = EoutAux1(:,:,VetThisCarr==ObsCarrUsed(ThisCarr));
                                %%            Fiber Time Delay Compensation
                                switch Medium
                                    case 'Fiber'
                                        Ix = ifft(fft(Ix).*exp(1j*2*pi*f*(FiberDelay(ThisCarr)*Ta)));
                                    otherwise
                                end
                                %The current incoming signal is them converted from the optical
                                %domain to the eletrical domain with the help of an photo
                                %detector.
                                Ix = Ix.*conj(Ix);
                                %The process with the photo diode is self-coherent, which means
                                %the rusult will be a component of the signal centered in f=0
                                %and another component centered at f=2*fc (frenquecy central).
                                %Therefore, to remove the higher order component a low pass
                                %filter will be used.
                                %             switch Medium
                                %                 case 'Fiber'
                                %%           Creating the Reception Filter
                                %             taux = t(1:length(Ix));
                                %             faux = time2freq(taux);
                                [BitFilt,~] = FiltroGaussiano(f,BWD2,CenFeq2,FiltOrd2);        %Creating filter for selection of the received signal
                                BitFilt = fftshift(BitFilt);                                   %Shifting the filter for matching the received signal
                                if RecFilBanPas==1
                                    Ix = ifft(fft(Ix).*BitFilt);
                                    Ix(1:4*(2*NumAmosCP+NPPB),:) = repmat(Ix(5*(2*NumAmosCP+NPPB),:),size(Ix(1:4*(2*NumAmosCP+NPPB),:),1),1);
                                    Ix(end+1-4*(2*NumAmosCP+NPPB):end,:) = repmat(Ix(5*(2*NumAmosCP+NPPB),:),size(Ix(1:4*(2*NumAmosCP+NPPB),:),1),1);
                                    IxMin = min(Ix);
                                    IxMin = repmat(IxMin,size(Ix,1),1);
                                    Ix = Ix - IxMin;                                         %Removing the DC component from them Eletrical signal received
                                    IxMax = max(Ix);
                                    IxMax = repmat(IxMax,size(Ix,1),1);
                                    Ix = Ix./IxMax;                                     %Normalizing the eletrical signal (amplifying what is needed)
                                end
                                %                 figure; plot(f,db(abs(fftshift(fft(Ix)./length(Ix))))); set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                                %                 otherwise
                                %             end
                                
                                
                                %%         Adding Noise to Received Signal
                                %Just after the photo diod is where the noise must be added as
                                %it is the main constrain of the system. Once the signal is
                                %converted from optical to electrical and it digital values
                                %were recevered, there is no point to add noise, because the
                                %information was already received, it just needs decoding.
                                %The noise here added will be gaussian basically it represents
                                %the reception sensibility. In another words, how much the
                                %signal must be about the noise for an acceptable recovering.
                                if ReceptorNoise==1                                               %Verify whether the noise is to be added or not
                                    [~,SigPower] = MeasPower(Ix);                              %it own received signal or
                                    SigPower2 = SigPower-30;%10*log10(SigPower);
                                    Ix = awgn(Ix,SNR,'measured');
                                end
                                
                                if ~RecFilBanPas
                                    Ix = ifft(fft(Ix).*BitFilt);
                                    Ix(1:4*(2*NumAmosCP+NPPB),:) = repmat(Ix(5*(2*NumAmosCP+NPPB),:),size(Ix(1:4*(2*NumAmosCP+NPPB),:),1),1);
                                    Ix(end+1-4*(2*NumAmosCP+NPPB):end,:) = repmat(Ix(5*(2*NumAmosCP+NPPB),:),size(Ix(1:4*(2*NumAmosCP+NPPB),:),1),1);
                                    IxMin = min(Ix);
                                    IxMin = repmat(IxMin,size(Ix,1),1);
                                    Ix = Ix - IxMin;                                         %Removing the DC component from them Eletrical signal received
                                    IxMax = max(Ix);
                                    IxMax = repmat(IxMax,size(Ix,1),1);
                                    Ix = Ix./IxMax;                                    %Normalizing the eletrical signal (amplifying what is needed)
                                end
                                
                                %                 [BitFilt,~] = FiltroGaussiano(f,200*BWD2,CenFeq2,FiltOrd2);        %Creating filter for selection of the received signal
                                %                 BitFilt = fftshift(BitFilt);                                   %Shifting the filter for matching the received signal
                                %                 Ix = ifft(fft(Ix).*BitFilt);
                                
                                %EoutAuxF = 20*log10(abs(fftshift(fft(Ix)./length(Ix))));
                                %for jj=1:size(Ix,2)
                                    %VetElecPower(CurrentTest,ThisCarr,jj)= MeasPower(Ix(:,jj));
                                    %[VetOptiPower(CurrentTest,ThisCarr,jj),~]= findpeaks(EoutAuxF(:,jj),'SortStr','descend','NPeaks',1);
                                %end
                                %For the reception process work properly, it is needed to
                                %sincronized the recieved signal or the sampling process will
                                %not work.
                                %AuxSync = (Ix(IniSyncPos:SyncPos));                            %Selecting the sync-word within the received signal
                                %AuxSync1 = AuxSync;
                                %This synchronization process is based on the mean expected
                                %value. That means, the information mean value should be within
                                %one period of symbol. Thus, the mean value of the received
                                %signal is acquired and compare of the known sync-word to
                                %verify if this mean value is at the right possition
                                %AuxSync1(AuxSync1<0)               = 0;                        %To keep the mean value above zero anything under is neglected
                                %AuxSync1(AuxSync1>=mean(AuxSync1)) = 1;                        %Adding a flag to the first sample of the received mean value
                                %AuxSync1(AuxSync1<mean(AuxSync1))  = -1;                       %All the others samples at set to the lowest level
                                %AuxSync2                           = SyncSymb(IniSyncPos:...
                                %SyncPos);%Selecting the sync-word within the known signal
                                %AuxSync2(AuxSync2>=mean(AuxSync2)) = 1;                        %Adding a flag to the first sample of the known mean value
                                %AuxSync2(AuxSync2<mean(AuxSync2))  = -1;                       %All the others samples at set to the lowest level
                                
                                %PosToSyn  = find(ismember(AuxSync1,1));                        %Finding where is the location of the first sample to synchronize
                                %PosSyn = find(ismember(AuxSync2,1));                           %Finding where is the location of the first sample to synchronize
                                
                                %AuxSyncCorr = PosToSyn(round(end/2))-PosSyn(round(end/2));
								if UsingMex==1
									[LevDec1,LevDec2,LevDec3,Ix,Ber4PAM(CurrentTest,ThisCarr),DecLevDef3,DecLevDef2,DecLevDef1,TxData,IxRecDef,IxRec,AberLev(1,CurrentTest,ThisCarr),AberLev(2,CurrentTest,ThisCarr),AberLev(3,CurrentTest,ThisCarr),ValsLev(1,CurrentTest,ThisCarr),ValsLev(2,CurrentTest,ThisCarr),ValsLev(3,CurrentTest,ThisCarr),ValsLev2(1,CurrentTest,ThisCarr),ValsLev2(2,CurrentTest,ThisCarr),ValsLev2(3,CurrentTest,ThisCarr),InterAB,InterCD,InterEF,SeqFinAB,SeqFinCD,SeqFinEF,LocMaxAB,LocMaxCD,LocMaxEF,MaxValAB,MaxValCD,MaxValEF,Levels] = RecDowPam4_mex(Ix,T,Ta,MaxNumStag,StuffSampels,NumAmosCP,NPPB,CurTesSiz,Nb4Pam,IntervalStep,MinDist,DecLevDef1,DecLevDef2,DecLevDef3,TxDataMat,ThisCarr,IfftOrSum,AddCP,SyncPeriod,DecMod);
                                else
									[LevDec1,LevDec2,LevDec3,Ix,Ber4PAM(CurrentTest,ThisCarr),DecLevDef3,DecLevDef2,DecLevDef1,TxData,IxRecDef,IxRec,AberLev(1,CurrentTest,ThisCarr),AberLev(2,CurrentTest,ThisCarr),AberLev(3,CurrentTest,ThisCarr),ValsLev(1,CurrentTest,ThisCarr),ValsLev(2,CurrentTest,ThisCarr),ValsLev(3,CurrentTest,ThisCarr),ValsLev2(1,CurrentTest,ThisCarr),ValsLev2(2,CurrentTest,ThisCarr),ValsLev2(3,CurrentTest,ThisCarr),InterAB,InterCD,InterEF,SeqFinAB,SeqFinCD,SeqFinEF,LocMaxAB,LocMaxCD,LocMaxEF,MaxValAB,MaxValCD,MaxValEF,Levels] = RecDowPam4(Ix,T,Ta,MaxNumStag,StuffSampels,NumAmosCP,NPPB,CurTesSiz,Nb4Pam,IntervalStep,MinDist,DecLevDef1,DecLevDef2,DecLevDef3,TxDataMat,ThisCarr,IfftOrSum,AddCP,SyncPeriod,DecMod);
								end
								%if BitErrS<BitErr
                                    %Ber4PAM(CurrentTest,ThisCarr) = BitErrS/((Nb4Pam-(4*SyncPeriod))*CurTesSiz);%Calculating the ration of wrong bits and the total number of bits transmited
                                %else
                                    %Ber4PAM(CurrentTest,ThisCarr) = BitErr/((Nb4Pam-(4*SyncPeriod))*CurTesSiz);%Calculating the ration of wrong bits and the total number of bits transmited
                                %end
                                %%           Ploting for Qualitative Analizes
                                %PrintInfo(Ploting*41,t(length(TxDataMat(ThisCarr,:))),Nb4Pam/2,TxDataMat(ThisCarr,:),IxRec);
                                %%       Calculating the Bit Error Ratio (BER)
                                %The final process here is to count the number of wrongdoings
                                %of this whole process upon the transmited data for
                                %quantitative analizes
                                
                                if CalcS==1
                                    limiar3 = (4/6)*abs(max(Ix));%(1/2)*abs(max(Ix)-min(Ix))/3;
                                    limiar2 = (2.5/6)*abs(max(Ix));%((1/2)*abs(max(Ix)-min(Ix))/3) + abs(max(Ix)-min(Ix))/3;
                                    limiar1 = (1/6)*abs(max(Ix));%((1/2)*abs(max(Ix)-min(Ix))/3) + 2*(abs(max(Ix)-min(Ix))/3);
                                    limiarPos1 = Levels(1);%abs(min(Ix));%Interval>limiar1;
                                    limiarPos2 = Levels(2);%0.33*(abs(max(Ix))-abs(min(Ix)));%(Interval<=limiar1)&(Interval>limiar2);
                                    limiarPos3 = Levels(3);%0.63*(abs(max(Ix))-abs(min(Ix)));%(Interval<=limiar2)&(Interval>limiar3);
                                    limiarPos4 = Levels(4);%abs(max(Ix));%Interval<=limiar3;
                                    
                                    EyeSymMa1 = reshape(Ix,NPPB,((Nb4Pam/2)-SyncSymbSiz)*CurTesSiz);
                                    Hi1 = zeros(1,size(EyeSymMa1,1));
                                    Hi2 = zeros(1,size(EyeSymMa1,1));
                                    Hi3 = zeros(1,size(EyeSymMa1,1));
                                    Lo1 = zeros(1,size(EyeSymMa1,1));
                                    Lo2 = zeros(1,size(EyeSymMa1,1));
                                    Lo3 = zeros(1,size(EyeSymMa1,1));
                                    LevHi1 = zeros(1,size(EyeSymMa1,1));
                                    LevHi2 = zeros(1,size(EyeSymMa1,1));
                                    LevHi3 = zeros(1,size(EyeSymMa1,1));
                                    LevLo1 = zeros(1,size(EyeSymMa1,1));
                                    LevLo2 = zeros(1,size(EyeSymMa1,1));
                                    LevLo3 = zeros(1,size(EyeSymMa1,1));
                                    EyeAb1 = zeros(1,size(EyeSymMa1,1));
                                    EyeAb2 = zeros(1,size(EyeSymMa1,1));
                                    EyeAb3 = zeros(1,size(EyeSymMa1,1));
                                    for kk = 1:size(EyeSymMa1,1)
                                        EyeHi3    = find((EyeSymMa1(kk,:)<limiarPos4)&(EyeSymMa1(kk,:)>limiar3));
                                        EyeHi2    = find((EyeSymMa1(kk,:)<limiarPos3)&(EyeSymMa1(kk,:)>limiar2));
                                        EyeHi1    = find((EyeSymMa1(kk,:)<limiarPos2)&(EyeSymMa1(kk,:)>limiar1));
                                        
                                        EyeLo3    = find((EyeSymMa1(kk,:)<limiar3)&(EyeSymMa1(kk,:)>limiarPos3));
                                        EyeLo2    = find((EyeSymMa1(kk,:)<limiar2)&(EyeSymMa1(kk,:)>limiarPos2));
                                        EyeLo1    = find((EyeSymMa1(kk,:)<limiar1)&(EyeSymMa1(kk,:)>limiarPos1));
                                        
                                        Hi1(kk)   = mean((EyeSymMa1(kk,EyeHi1)));
                                        Hi2(kk)   = mean((EyeSymMa1(kk,EyeHi2)));
                                        Hi3(kk)   = mean((EyeSymMa1(kk,EyeHi3)));
                                        
                                        Lo1(kk)   = mean((EyeSymMa1(kk,EyeLo1)));
                                        Lo2(kk)   = mean((EyeSymMa1(kk,EyeLo2)));
                                        Lo3(kk)   = mean((EyeSymMa1(kk,EyeLo3)));
                                        
                                        LevHi1(kk)= mean(Hi1) - std(EyeSymMa1(kk,EyeHi1));
                                        LevHi2(kk)= mean(Hi2) - std(EyeSymMa1(kk,EyeHi2));
                                        LevHi3(kk)= mean(Hi3) - std(EyeSymMa1(kk,EyeHi3));
                                        
                                        LevLo1(kk)= mean(Lo1) + std(EyeSymMa1(kk,EyeLo1));
                                        LevLo2(kk)= mean(Lo2) + std(EyeSymMa1(kk,EyeLo2));
                                        LevLo3(kk)= mean(Lo3) + std(EyeSymMa1(kk,EyeLo3));
                                        
                                        EyeAb1(kk)= LevHi1(kk) - LevLo1(kk);
                                        EyeAb2(kk)= LevHi2(kk) - LevLo2(kk);
                                        EyeAb3(kk)= LevHi3(kk) - LevLo3(kk);
                                    end
                                    EyeAbertura1 = mean(EyeAb1);
                                    EyeAbertura2 = mean(EyeAb2);
                                    EyeAbertura3 = mean(EyeAb3);
                                    ThisDeciLev1 = mean(LevLo1) + EyeAbertura1/2;
                                    ThisDeciLev2 = mean(LevLo2) + EyeAbertura2/2;
                                    ThisDeciLev3 = mean(LevLo3) + EyeAbertura3/2;
                                    AberLevS(1,ThisCarr,CurrentTest)  = EyeAbertura3;
                                    AberLevS(2,ThisCarr,CurrentTest)  = EyeAbertura2;
                                    AberLevS(3,ThisCarr,CurrentTest)  = EyeAbertura1;
                                    ValsLevS(1,ThisCarr,CurrentTest)  = ThisDeciLev3;
                                    ValsLevS(2,ThisCarr,CurrentTest)  = ThisDeciLev2;
                                    ValsLevS(3,ThisCarr,CurrentTest)  = ThisDeciLev1;

                                    %plot(t(NPPB/2),ThisDeciLev3,'b*');plot(t(NPPB/2),mean(LevHi3),'bv');plot(t(NPPB/2),mean(LevLo3),'b^');
                                    %plot(t(NPPB/2),ThisDeciLev2,'k*');plot(t(NPPB/2),mean(LevHi2),'kv');plot(t(NPPB/2),mean(LevLo2),'k^');
                                    %plot(t(NPPB/2),ThisDeciLev1,'g*');plot(t(NPPB/2),mean(LevHi1),'gv');plot(t(NPPB/2),mean(LevLo1),'g^');
                                    
                                    ThisDataSize = NPPB/2:NPPB:length(Ix);
                                    IxRecDeS = zeros(1,2*size(ThisDataSize,2));%Initialization of the vector that will store the income data
                                    ContBit2  = 1;
                                    for kk=1:length(ThisDataSize)%NPPB:length(Ix)                                       %The comparison process will be made for each symbol period
                                        %                 midaux = round(mean(SymLoc(1:round(end/2))));
                                        midaux = NPPB/2;%SymLoc(1);
                                        aux1 = Ix((kk-1)*NPPB+midaux+1);     %An small portion of the income signal is take for evaluation
                                        MeanRec = mean(aux1);                                      %Measuring the avarage value of the samples taken
                                        if MeanRec <= ThisDeciLev1                                   %If it is the lowest level the incoming data
                                            IxRecDeS(ContBit2) = 0;
                                            ContBit2  = 1 + ContBit2;
                                            IxRecDeS(ContBit2) = 0;
                                            ContBit2  = 1 + ContBit2;
                                            %IxRecDeS = [IxRecDeS 0 0];                             %is 01 (1)
                                        elseif (MeanRec <= ThisDeciLev2)&&(MeanRec > ThisDeciLev1)     %If it is the second level the incoming data
                                            IxRecDeS(ContBit2) = 0;
                                            ContBit2  = 1 + ContBit2;
                                            IxRecDeS(ContBit2) = 1;
                                            ContBit2  = 1 + ContBit2;
                                            %IxRecDeS = [IxRecDeS 0 1];                             %is 00 (0)
                                        elseif (MeanRec <= ThisDeciLev3)&&(MeanRec > ThisDeciLev2)     %If it is the tird level the incoming data
                                            IxRecDeS(ContBit2) = 1;
                                            ContBit2  = 1 + ContBit2;
                                            IxRecDeS(ContBit2) = 1;
                                            ContBit2  = 1 + ContBit2;
                                            %IxRecDeS = [IxRecDeS 1 1];                             %is 10 (2)
                                        elseif MeanRec > ThisDeciLev3                                %If it is the uper level the incoming data
                                            IxRecDeS(ContBit2) = 1;
                                            ContBit2  = 1 + ContBit2;
                                            IxRecDeS(ContBit2) = 0;
                                            ContBit2  = 1 + ContBit2;
                                            %IxRecDeS = [IxRecDeS 1 0];                             %is 11 (3)
                                        else                                                       %If for some misteriose reason neither of previous verification were sucedded
                                            IxRecDeS(ContBit2) = 0;
                                            ContBit2  = 1 + ContBit2;
                                            IxRecDeS(ContBit2) = 0;
                                            ContBit2  = 1 + ContBit2;
                                            %IxRecDeS = [IxRecDeS 0 0];                             %by default the current data is set to be 00 (0)
                                        end
                                    end
                                    
                                    BitErrS = sum(xor(TxData,IxRecDeS));%Comparison between the Transmited and received and counting the differences
                                    Ber4PAMS(CurrentTest,ThisCarr) = BitErrS/((Nb4Pam-(4*SyncPeriod))*CurTesSiz);%Calculating the ration of wrong bits and the total number of bits transmited
                               
                                end
                                
                                %%           Ploting for Qualitative Analizes
                                %if ThisCarr==126
                                %EyeToPlot(CurrentTest,1:length(Esync)) = Esync;
                                %save(['savingforplotingeye' num2str(VetSnr)],'EyeToPlot','VetSnr','IfftOrSum','UsedModula','T','NPPB');
                                %end
                                if PrintinEye==1
                                    PrintInfo(PlotingThis*36,Ix,T,NPPB);
                                    hold all;
                                    plot(t(NPPB/2),LevDec1,'b*');plot(t(NPPB/2),InterEF(SeqFinEF(LocMaxEF)-1),'bv');plot(t(NPPB/2),InterEF(round(SeqFinEF(LocMaxEF)-MaxValEF+1)),'b^');
                                    plot(t(NPPB/2),LevDec2,'gd');plot(t(NPPB/2),InterCD(SeqFinCD(LocMaxCD)-1),'gv');plot(t(NPPB/2),InterCD(round(SeqFinCD(LocMaxCD)-MaxValCD+1)),'g^');
                                    plot(t(NPPB/2),LevDec3,'kd');plot(t(NPPB/2),InterAB(SeqFinAB(LocMaxAB)-1),'kv');plot(t(NPPB/2),InterAB(round(SeqFinAB(LocMaxAB)-MaxValAB+1)),'k^');
                                  
                                    %figure;
                                    %hold all;
                                    %plotpos = zeros(1,length(IxAux));
                                    %plot(IxAux,plotpos,'o','color',[1 0.4 0]);
                                    %plot(EvmMatRef(ObsCarrUsed(ThisCarr),:),plotpos,'b*');
                                    set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                                    %a=a+1;
                                    %                     drawnow;
                                    %                     pause(0.1);
                                end
                                if PlotingThis
                                    figure;
                                    hold all;
                                    plot(TxData);
                                    plot(IxRec);
                                    set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                                    display(Ber4PAM);
                                end
                                close all;
                                %                             Ber4PAM(CurrentTest,ThisCarr+1) = 1;
                                
                            end
                        end
                    otherwise
                        %%              Receiver OOK
                        %         AberLev(CurrentTest,:) = zeros(1,NumCarr);
                        for ThisCarr=InitCarrDo:CarrPass:NumCarr                                           %For each carrier the same process of reception will be used.
                            if CarrUsedDo(ThisCarr)
                                %%  Reception Process: Handling the income optical field
                                %At the first step the income field will be selected from the
                                %optical FFT Output with the help of the maping vector
                                %(previously described).
                                Ix = EoutAux1(:,:,VetThisCarr==ObsCarrUsed(ThisCarr));
                                %%            Fiber Time Delay Compensation
                                switch Medium
                                    case 'Fiber'
                                        Ix = ifft(fft(Ix).*exp(1j*2*pi*f*(...
                                            FiberDelay((ThisCarr))*Ta)));
                                    otherwise
                                end
                                
                                %The current incoming signal is them converted from the optical
                                %domain to the eletrical domain with the help of an photo
                                %detector.
                                Ix =Ix.*conj(Ix);
                                %%           Creating the Reception Filter
                                [BitFilt,~] = FiltroGaussiano(f,BWD,CenFeq,FiltOrd);%Creating filter for selection of the received signal
                                BitFilt = fftshift(BitFilt);%Shifting the filter for matching the received signal
                                %%
                                if RecFilBanPas==1
                                    Ix = ifft(fft(Ix).*BitFilt);
                                    IxMin = min(Ix);
                                    IxMin = repmat(IxMin,size(Ix,1),1);
                                    Ix = Ix - IxMin;%Removing the DC component from them Eletrical signal received
                                    IxMax = max(Ix);
                                    IxMax = repmat(IxMax,size(Ix,1),1);
                                    Ix = Ix./IxMax;
                                    %Ix = Ix - min(Ix);%Removing the DC component from them Eletrical signal received
                                    %Ix = Ix./max(abs(Ix));%Normalizing the eletrical signal (amplifying what is needed)
                                end
                                
                                if ReceptorNoise==1
                                    [~,SigPower] = MeasPower(Ix);
                                    SigPower2 = SigPower-30;%10*log10(SigPower);
                                    Ix = awgn(Ix,SNR,SigPower2);
                                end
                                if ~RecFilBanPas
                                    Ix = ifft(fft(Ix).*BitFilt);
                                    IxMin = min(Ix);
                                    IxMin = repmat(IxMin,size(Ix,1),1);
                                    Ix = Ix - IxMin;%Removing the DC component from them Eletrical signal received
                                    IxMax = max(Ix);
                                    IxMax = repmat(IxMax,size(Ix,1),1);
                                    Ix = Ix./IxMax;
                                end
                                
                                %EoutAuxF = 20*log10(abs(fftshift(fft(Ix)./length(Ix))));
                                %for jj=1:size(Ix,2)
                                    %VetElecPower(CurrentTest,ThisCarr,jj)= MeasPower(Ix(:,jj));
                                    %[VetOptiPower(CurrentTest,ThisCarr,jj),~]= findpeaks(EoutAuxF(:,jj),'SortStr','descend','NPeaks',1);
                                %end
                                %% Ploting the result for qualitative analizes
                                %             PrintInfo(Ploting*42,f,20.*log10(abs(fftshift(fft(Ix)./...
                                %                                                             length(Ix)))));
                                %             PrintInfo(Ploting*43,t,TxSigMat(ThisCarr,:)./max(TxSigMat(ThisCarr,:)),Ix);
                                %             PrintInfo(Ploting*44,Ix,T,NPPB);
                                %             PrintInfo(Ploting*45,t(IniSyncPos:SyncPos),Ix(IniSyncPos:...
                                %                                     SyncPos),SyncSymb(IniSyncPos:SyncPos));
                                %%        synchronizing
                                %For the reception process work properly, it is needed to
                                %sincronized the recieved signal or the sampling process will
                                %not work.
                                %AuxSync = (Ix(IniSyncPos:SyncPos));                            %Selecting the sync-word within the received signal
                                %AuxSync1 = AuxSync;
                                %This synchronization process is based on the mean expected
                                %value. That means, the information mean value should be within
                                %one period of symbol. Thus, the mean value of the received
                                %signal is acquired and compare of the known sync-word to
                                %verify if this mean value is at the right possition
                                %AuxSync1(AuxSync1<0) = 0;                                      %To keep the mean value above zero anything under is neglected
                                %AuxSync1(AuxSync1>=mean(AuxSync1)) = 1;                        %Adding a flag to the first sample of the received mean value
                                %AuxSync1(AuxSync1<mean(AuxSync1)) = -1;                        %All the others samples at set to the lowest level
                                %AuxSync2 = SyncSymb(IniSyncPos:SyncPos);                       %Selecting the sync-word within the known signal
                                %AuxSync2(AuxSync2>=mean(AuxSync2)) = 1;                        %Adding a flag to the first sample of the known mean value
                                %AuxSync2(AuxSync2<mean(AuxSync2)) = -1;                        %All the others samples at set to the lowest level
                                %PosToSyn  = find(ismember(AuxSync1,1));
                                %PosSyn    = find(ismember(AuxSync2,1));
                                
                                %AuxSyncCorr = PosToSyn(round(end/2))-PosSyn(round(end/2));
                                sn =((1/2)*(((1/2)^MaxNumStag)-1))/((1/2)-1);
                                AuxSyncCorr = round((sn/(2/2^IfftOrSum))*T/Ta);
                                Ix = [Ix(AuxSyncCorr+1:end,:);Ix(1:AuxSyncCorr,:)];
                                %The difference between the PossitionTosynchronize and
                                %Possitionsynchronized will be used to correct the time
                                %shifting on the transmition and reception process.
                                
                                %if AuxSyncCorr>=0%If the difference is positive, left-shift...
                                %                 Ix = ifft(fft(Ix).*exp(1j*2*pi*f*(AuxSyncCorr*Ta)));     %Shift based on time change
                                %Ix = [Ix(AuxSyncCorr+1:end) Ix(1:AuxSyncCorr)];            %Shift based on sampling sliding
                                %else%... but if the difference is negative, right-shift
                                %                 Ix = ifft(fft(Ix).*exp(-1j*2*pi*f*(AuxSyncCorr*Ta)));    %Shift based on time change
                                %Ix = [Ix(end+AuxSyncCorr+1:end) Ix(1:end + AuxSyncCorr)];  %Shift based on sampling sliding
                                %end
                                
                                
                                %%        synchronizing
                                %if SencondAdjust
                                %For some reason that we could not understand sometimes the
                                %time (sampling) sliding of the signal is not equal
                                %throught the data stream. Thus, the second part of the
                                %synchronism process will be turn ON or OFF according to the
                                %user's will.
                                %AuxSyncEnd = (Ix(end-(SyncPos+1*abs(AuxSyncCorr))+1:end-...
                                %(2*NPPB+1*abs(AuxSyncCorr))));%Selecting the sync-word within the received signal
                                %SyncSymbEndAux = SyncSymbEnd(end-SyncPos+1:end-2*NPPB);    %Selecting the sync-word within the known signal
                                %AuxSyncEnd1 = AuxSyncEnd;
                                %AuxSyncEnd1(AuxSyncEnd1<0) = 0;                            %To keep the mean value above zero anything under is neglected
                                %AuxSyncEnd1(AuxSyncEnd1>=mean(AuxSyncEnd1)) = 1;           %Adding a flag to the first sample of the received mean value
                                %AuxSyncEnd1(AuxSyncEnd1<mean(AuxSyncEnd1)) = -1;           %All the others samples at set to the lowest level
                                %AuxSyncEnd2 = SyncSymbEndAux;
                                %AuxSyncEnd2(AuxSyncEnd2>=mean(AuxSyncEnd2)) = 1;           %Adding a flag to the first sample of the known mean value
                                %AuxSyncEnd2(AuxSyncEnd2<mean(AuxSyncEnd2)) = -1;           %All the others samples at set to the lowest level
                                %PosToSynEnd  = find(ismember(AuxSyncEnd1,1));
                                %PosSynEnd = find(ismember(AuxSyncEnd2,1));
                                %AuxSyncCorrEnd = PosToSynEnd(round(end/2))-PosSynEnd(...
                                %round(end/2));
                                %The difference between the PossitionTosynchronize and
                                %Possitionsynchronized will be used to correct the time
                                %shifting on the transmition and reception process.
                                %if AuxSyncCorrEnd>=0%If possitive difference, left-shift...
                                %Ix = ifft(fft(Ix).*exp(1j*2*pi*f*(AuxSyncCorrEnd*...
                                %Ta)));%Shift based on time change
                                %Ix = [Ix(AuxSyncCorrEnd+1:end) Ix(1:AuxSyncCorrEnd)];
                                %else%... but if the difference is negative, right-shift
                                %Ix = ifft(fft(Ix).*exp(-1j*2*pi*f*(AuxSyncCorrEnd*...
                                %Ta)));%Shift based on time change
                                %Ix = [Ix(end-AuxSyncCorrEnd+1:end) Ix(1:end - ...
                                %AuxSyncCorrEnd)];%Shift based on sampling sliding
                                %end
                                %end
                                
                                %% Ploting the result for qualitative analizes
                                
                                %                             PrintInfo(Ploting*45,t(IniSyncPos:SyncPos),Ix(IniSyncPos:SyncPos),SyncSymb(IniSyncPos:SyncPos));
                                %                             PrintInfo(Ploting*46,t(1:length(Ix)),Ix);
                                %%                  Removing CP
                                if AddCP
                                    IxAux = Ix(1:end - StuffSampels,:);
                                    IxAux = reshape(IxAux,(2*NumAmosCP+NPPB),(Nb-NumBitDesc)*CurTesSiz);
                                    IxAux = IxAux(1+NumAmosCP:end-NumAmosCP,:);
                                    IxAux = reshape(IxAux,NPPB*(Nb-NumBitDesc),CurTesSiz);
                                    Ix    = IxAux;
                                end
                                %% Taking the sampling the EVM meassurement
                                %                 evm = comm.EVM('MaximumEVMOutputPort',true,...
                                %                 'XPercentileEVMOutputPort',true, 'XPercentileValue',90,...
                                %                 'SymbolCountOutputPort',true);
                                %PosAuxEout = NPPB/2:NPPB:length(Ix);                       %Varriable respossible to take just the samples at the middle of the symbol
                                %IxAux      = Ix(PosAuxEout);                               %Normalizing the reference
                                %a=a+0;
                                %RxSymbAmos(CurrentTest,ObsCarrPos==ThisCarr,:) = IxAux;%RxSymbAmosUp = [];
                                %EvmMatRec(ObsCarrPos==ThisCarr,:) = IxAux;                       %Taking just the middle samples as references
                                %[EvmDB(CurrentTest,ThisCarr), EvmPer(CurrentTest,ThisCarr), EvmRms(CurrentTest,ThisCarr) ] = EvmCalc( EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux );
                                %[EvmDBJ(CurrentTest,ThisCarr),EvmPerJ(CurrentTest,ThisCarr),EvmRmsJ(CurrentTest,ThisCarr)] = evm1(2,'pam',EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux);
                                %##########################################################################
                                %if (CalcS==1)||(PrintinEye==1)
                                    %IxS = Ix.';
                                %end
                                Ix = Ix(1+SyncPeriod*NPPB:end-SyncPeriod*NPPB,:);
                                Ix = reshape(Ix,1,size(Ix,1)*size(Ix,2));
                                if UsingMex==1
                                    [~,~,EyeOpen,EyeOpenHigh,EyeOpenLow] = Olho_mex((Ix),T,NPPB,0,1);
                                    [AberLev(CurrentTest,ThisCarr),ValsLev(CurrentTest,ThisCarr),Inter,SeqFin,LocMax,MaxVal,TxData,EoutCorrD,EoutCorr2,EoutCorr,BerOOK(CurrentTest,ThisCarr),LocMax2,SeqFin2,MaxVal2,Inter2] = RedDowOok_mex(Ix,TxDataMat,NPPB,ThisCarr,Nb,NumBitDesc,SyncPeriod,CurTesSiz,EyeOpenLow,EyeOpen);
                                else
                                    [~,~,EyeOpen,EyeOpenHigh,EyeOpenLow] = Olho((Ix),T,NPPB,0,1);
                                    [AberLev(CurrentTest,ThisCarr),ValsLev(CurrentTest,ThisCarr),Inter,SeqFin,LocMax,MaxVal,TxData,EoutCorrD,EoutCorr2,EoutCorr,BerOOK(CurrentTest,ThisCarr),LocMax2,SeqFin2,MaxVal2,Inter2] = RedDowOok(Ix,TxDataMat,NPPB,ThisCarr,Nb,NumBitDesc,SyncPeriod,CurTesSiz,EyeOpenLow,EyeOpen);
                                end
                                if CalcS==1
                                    AberLevS(CurrentTest,ThisCarr)= EyeOpen;
                                    ValsLevS(CurrentTest,ThisCarr)= EyeOpenLow + EyeOpen/2;
                                    
                                    ThisDataSize = NPPB/2:NPPB:length(Ix);
                                    EoutCorrS = zeros(1,length(ThisDataSize));                                                 %Initialization of the vector that will store the income data
                                    for kk=1:length(ThisDataSize)%length(Ix(ThisDataSize))                                       %The comparison process will be made for each symbol period
                                        %An small portion of the income signal is take for
                                        %evaluation by measuring the avarage value of the samples
                                        %taken
                                        %                 CalcMean = mean((Ix((kk-1)+SymLoc(1))));
                                        CalcMeanS = mean((Ix((kk-1)*NPPB+NPPB/2)));
                                        %Verifying the interval for each symbol received.
                                        if CalcMeanS >= EyeOpenLow + EyeOpen/2%EyeOpenLow+EyeOpen/2                        %If it is the uper level the incoming data
                                            EoutCorrS(kk) = 1;                               %is 1
                                        else                                                       %If it is the lowest level the incoming data
                                            EoutCorrS(kk) = 0;                               %is 0
                                        end
                                    end
                                    BitErrS = sum(xor(TxData,EoutCorrS));%Comparison between the Transmited and received and counting the differences
                                    BerOOKS(CurrentTest,ThisCarr) = BitErrS/(((Nb-NumBitDesc)-2*SyncPeriod)*CurTesSiz);%Calculating the ration of wrong bits and the total number of bits transmited
                                end
                                
                                %% Ploting the result for qualitative analizes
                                
                                %if ThisCarr==126
                                %EyeToPlot(CurrentTest,1:length(Ix)) = Ix;
                                %save(['savingforplotingeye' num2str(VetSnr)],'EyeToPlot','VetSnr','IfftOrSum','UsedModula','T','NPPB');
                                %end
                                if PrintinEye==1
                                    PrintInfo(PlotingThis*47,(Ix),T,NPPB);
                                    hold on;
                                    plot(t((NPPB)/2),Inter(round(SeqFin(LocMax)-MaxVal)),'k^');
                                    plot(t((NPPB)/2),Inter(SeqFin(LocMax)),'kv');
                                    plot(t((NPPB)/2),Inter(round(SeqFin(LocMax)-MaxVal/2)),'kd');
                                    plot(t((NPPB)/2),Inter2(round(SeqFin2(LocMax2)-MaxVal2)),'b^');
                                    plot(t((NPPB)/2),Inter2(SeqFin2(LocMax2)),'bv');
                                    plot(t((NPPB)/2),Inter2(round(SeqFin2(LocMax2)-MaxVal2/2)),'b*');
                                    if CalcS
                                        plot(t((NPPB)/2),EyeOpenLow,'mx');
                                        plot(t((NPPB)/2),EyeOpenHigh,'mo');
                                        plot(t((NPPB)/2),EyeOpenLow + EyeOpen/2,'md');
                                    end
                                    set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                                    %figure;
                                    %hold all;
                                    %plotpos = zeros(1,length(IxAux));
                                    %plot(IxAux,plotpos,'o','color',[1 0.4 0]);
                                    %plot(EvmMatRef(ObsCarrUsed(ThisCarr),:),plotpos,'b*');
                                    %set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                                    %a=a+0;
                                end
                                if PlotingThis
                                    figure;
                                    hold all;
                                    plot(TxData);
                                    plot(EoutCorr);
                                    set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                                    display(BerOOK);
                                end
                                close all;
                                %                             BerOOK(CurrentTest,ThisCarr+1) = 1;
                                
                            end
                            %         end
                        end
                end
                
            else
                %DatumReceptionInputData;
                if UsingGpu==1
                    Tgpu = gpuArray(T);
                    MaxStagGpu = gpuArray(MaxNumStag);
                    EoutAux1 = zeros(size(EoutRec,1),size(EoutRec,2),length(ObsCarrUsed));
                    for jj=1:size(EoutRec,2)
                        EoutGpu = gpuArray(EoutRec(:,jj));
                        [EoutAux1Gpu,~,VetThisCarrGpu]=OpticalFFTG(fgpu,Tgpu,MaxStagGpu,EoutGpu);
                        EoutAux1(:,jj,:) = gather(EoutAux1Gpu);
                        VetThisCarr = gather(VetThisCarrGpu);
                        clear EoutGpu EoutAux1Gpu VetThisCarrGpu;
                    end
                    clear Tgpu MaxStagGpu;
                else
                    [EoutAux1,~,VetThisCarr]=OpticalFFTN(f,T,MaxNumStag,EoutRec);
                end
            end
        end
        if (UseOfdmElect)&&(UsedModula==5)&&(NumCarr==2)
        else
            if SendingUp
                %%              Datum Upstream
                %This section is resposible to generate and to send the
                %information from the user to the OLT.
                %% DatumUpstream;
                clear TxData TxDataMat EoutAux1 VetThisCarr
                if UsingGpu==1
                    Tgpu = gpuArray(T);
                    MaxStagGpu = gpuArray(MaxNumStag); 
                    EoutAux1 = zeros(size(EoutRec,1),size(EoutRec,2),length(ObsCarrUsed));
                    for jj=1:size(EoutRec,2)
                        EoutGpu = gpuArray(EoutRec(:,jj));
                        [EoutAux1Gpu,~,VetThisCarrGpu]=OpticalFFTGT(fgpu,Tgpu,MaxStagGpu,EoutGpu);
                        EoutAux1(:,jj,:) = gather(EoutAux1Gpu);
                        VetThisCarr = gather(VetThisCarrGpu);
                        clear EoutGpu EoutAux1Gpu VetThisCarrGpu;
                    end
                    clear Tgpu MaxStagGpu;
                else
                    [EoutAux1,~,VetThisCarr]=OpticalFFTNT(f,T,MaxNumStag,EoutRec);
                end
                clear EoutRec
                if Ploting
                    figure;hold all;
                    for kk=1:size(EoutAux1,3)
                        plot(f(:,1),20*log10(abs(fftshift(fft(EoutAux1(:,1,VetThisCarr==kk)./length(EoutAux1(:,1,VetThisCarr==kk)))))));
                    end
                    a=1;
                end
                % PrintInfo(Ploting*12,f,db(abs((fftshift(fft(EoutRec))))/length(EoutRec)),EoutAux1,RefCarr*fc);
                %%               Creating Data for transmission
                % Creating the information accordingly with the type of transmition
                
                %It is important to mention that an stream of data will be formed for each
                %individual carrier. For the last but not the least, the DownStream and the
                %UpStream carriers will be interleaved within the same transmition pass
                %band.For example, if the user choses carriers 1,3,5 and 7 as DownStream
                %carriers, the UpStream will be formed by carriers 2,4,6, and 8. This
                %manuver was suggested by Segatto to address the problem of ICI (Inter
                %Carrier Interference). The transmited signal still an OFDM once the
                %carriers from 1 to 8, for instance, still orthogonal to each other.
                %Therefore, the for loop to generate the information to be transmited will
                %have the passe of 2.
                
                switch Modulation
                    case 'OFDM'
                        EoutModTem = zeros(NPPB*Nb,CurTesSiz,NumCarr);
                        TxDataMat = zeros(NumCarr,NumFra*length(DmtMve)*CurTesSiz);
                        for kk=InitCarrUp:CarrPass:NumCarr%Generating different data for each carrier
                            UniDmtMve = unique(DmtMve);
                            UniDmtMve = fliplr(UniDmtMve);
                            DaSiAn = 0;
                            DaSiPo = 0;
                            for CarDmtM=1:length(UniDmtMve)
                                M = UniDmtMve(CarDmtM);
                                DaSiAu = sum(ismember(DmtMve,M));
                                TxData = randi([0 M-1],NumFra*CurTesSiz,DaSiAu);%Generation random information
                                TxDataMat(kk,1+DaSiAn:DaSiAn + CurTesSiz*DaSiAu*NumFra) = TxData(:);%Saving data for future comparison
                                DaSiAn = DaSiAn + CurTesSiz*DaSiAu*NumFra;
                                switch OfdMod%Sellect which modulation was used
                                    case 'qam'
                                        TxSymbAux(1:NumFra*CurTesSiz,1+DaSiPo:DaSiPo + DaSiAu)  = qammod(TxData,M);%Modulating information by QAM
                                    otherwise
                                        TxSymbAux(1:NumFra*CurTesSiz,1+DaSiPo:DaSiPo + DaSiAu)  = dpskmod(TxData,M);%Modulating information DPSK as default
                                end
                                DaSiPo = DaSiPo + DaSiAu;
                            end
                            TxSymb = TxSymbAux.';
                            TxSymbConj      = conj(flipud(TxSymb));                    %Singnal conjugate format at reverse order
                            %The Hermitian simetry was composed of the signal its conjugate format and
                            %zeros to fill empty spaces. The signal is formed as show below:
                            %                    signal     midle  signal conjugated
                            %                       |        |       |
                            %                       V        V       V
                            %            TxH = [0 a + jb 0 0 0 0 0 a - jb 0];
                            %
                            %I thought in many ways to form this signal, replacing the centred zeros by
                            %copy of the signal bay just geting its information and adding
                            %respectively. Also, by creating the signal with redundancy at the exact
                            %length needed to make the hermitian singnal. But all tests end up with
                            %similar results, no improvement was noticed.
                            if UsingHermitian
                                TxSymbH = [zeros(1,NumFra*CurTesSiz);TxSymb;zeros(1,NumFra*CurTesSiz);TxSymbConj];%Hermitian Symetri
                                TxSymbC = zeros(NFFT,NumFra*CurTesSiz);
                                TxSymbC(1+Zt-ZtC:Zt-ZtC+size(TxSymbH,1)/2,:) = TxSymbH(1:size(TxSymbH,1)/2,:);
                                TxSymbC(end-Zt+ZtC-(size(TxSymbH,1)/2)+1:end-Zt+ZtC,:) = TxSymbH(1+size(TxSymbH,1)/2:end,:);
                            else
                                TxSymbH = TxSymb;
                                %     TxSymbC = [zeros((NFFT-size(TxSymbH,1))/2,NumFra*CurTesSiz);TxSymbH;zeros((NFFT-size(TxSymbH,1))/2,NumFra*CurTesSiz)];
                                TxSymbC = [TxSymbH(1:size(TxSymbH,1)/2,:);zeros((NFFT-size(TxSymbH,1)),NumFra*CurTesSiz);TxSymbH(end+1-size(TxSymbH,1)/2:end,:)];
                            end
                            %                 plot(TxSymbC);
                            %Here was implemented the modulation through FFT process, it basica mix-up
                            %all infomation following Fourier transform.
                            TxSymbMod = ifft(TxSymbC,NFFT,1);
                            TxSymbMod = rectpulse(TxSymbMod,OvSam);
                            tt = linspace(0,1*Te,size(TxSymbMod,1));
                            ff = time2freq(tt);
                            tta = repmat(tt.',1,NumFra*CurTesSiz);
                            switch SelModTp
                                case 'BP'%Base Band transmission
                                    if SelecGaus==1
                                        [ ModFilt ] = FiltroGaussiano(ff,OBw/2,0,1);%CarrOffSet*FramSpac,1);
                                    else
                                        [ ModFilt ] = Filtro_Retangular(OBw,0,ff);%CarrOffSet*FramSpac,ff);
                                    end
                                    ModFilt = fftshift(ModFilt);
                                    if (kk==1)&&(PlotingThis)
                                        figure;hold all;
                                        plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                                    end
                                    ModFilta = repmat(ModFilt.',1,NumFra*CurTesSiz);
                                    if SelModFilt==1
                                        TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
                                    end
                                    if (kk==1)&&(PlotingThis)
                                        plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                                        plot(ff,20*log10(abs(fftshift(ModFilt))));
                                        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                                    end
                                case 'AM'                                              %Out of base band
                                    TxSymbMod   = TxSymbMod.*cos(2*pi*Ofc*tta);
                                    if SelecGaus==1
                                        [ ModFilt ] = FiltroGaussiano(ff,OBw,Ofc,1);
                                    else
                                        [ ModFilt ] = Filtro_Retangular( OBw,Ofc,ff);
                                    end
                                    ModFilt = fftshift(ModFilt);
                                    if (kk==1)&&(PlotingThis)
                                        figure;hold all;
                                        plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                                    end
                                    ModFilta = repmat(ModFilt.',1,NumFra*CurTesSiz);
                                    if SelModFilt==1
                                        TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
                                    end
                                    if (kk==1)&&(PlotingThis)
                                        plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                                        plot(ff,20*log10(abs(fftshift(ModFilt))));
                                        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                                    end
                                case 'AMSSB'                                           %Out of base band
                                    TxSymbMod   = TxSymbMod.*exp(1j*2*pi*Ofc*tta);
                                    if SelecGaus==1
                                        [ ModFilt ] = FiltroGaussiano(ff,OBw,Ofc,1);
                                    else
                                        [ ModFilt ] = Filtro_Retangular( OBw,Ofc,ff);
                                    end
                                    ModFilt = fftshift(ModFilt);
                                    if (kk==1)&&(PlotingThis)
                                        figure;hold all;
                                        plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                                    end
                                    ModFilta = repmat(ModFilt.',1,NumFra*CurTesSiz);
                                    if SelModFilt==1
                                        TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
                                    end
                                    if (kk==1)&&(PlotingThis)
                                        plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                                        plot(ff,20*log10(abs(fftshift(ModFilt))));
                                        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                                        a=a+1;
                                    end
                                otherwise
                                    if SelecGaus==1
                                        [ ModFilt ] = FiltroGaussiano(ff,3*Ofc,0,1);
                                    else
                                        [ ModFilt ] = Filtro_Retangular( 6*Ofc,0,ff);
                                    end
                                    ModFilt = fftshift(ModFilt);
                                    if (kk==1)&&(Ploting)
                                        figure;hold all;
                                        plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                                    end
                                    ModFilta = repmat(ModFilt.',1,NumFra*CurTesSiz);
                                    TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
                                    if (kk==1)&&(Ploting)
                                        plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                                        plot(ff,20*log10(abs(fftshift(ModFilt))));
                                        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                                    end
                            end
                            TxSymbModA  = [TxSymbMod(1+end-NPOFEX:end,:);TxSymbMod];%Adding Ciclyc prefix
                            U1t  = rectpulse(TxSymbModA,NPPOF);%Over sampling
                            TxSigA = reshape(U1t,OvSam*NumFra*NFFT,CurTesSiz);
                            U1t = [TxSigA(1+end-NPSTUf:end,:);TxSigA];%Conforming vector sizes to the same length
                            
                            %Assigning the eletrical signal to one drive of the MZM - The aplitude of
                            %the signal can be controlled by the variable DatGai, which can be
                            %understood as an gain for the eletrical signal or an atenuation. The
                            %second signal will be similar with the only difference a phase shift of pi
                            NormFact = max(U1t)./0.95;
                            NormFact = repmat(NormFact,size(U1t,1),1);
                            U1t = (U1t./NormFact);%Normalizing the sinal to mach with the MZM espect to receive.
                            %U.U1t = U1t;
                            U2t = exp(-1j*pi).*U1t;
                            %As both signals will have the mostrly the same characteristics with the
                            %only difference the phase shift of 180 degress. The MZM-I will be working
                            %on the Push-Pull configuration. It is necessary to reduce the Chirp noise
                            %to zero.
                            EoutTxAux = EoutAux1(:,:,VetThisCarr==(RefCarr-1+kk));
                            if UsingGpu==1
                                Lgpu = gpuArray(L);
                                U0gpu = gpuArray(U0);
                                U_pi1gpu = gpuArray(U_pi1);
                                U_pi2gpu = gpuArray(U_pi2);
                                noptgpu = gpuArray(nopt);
                                nelgpu = gpuArray(nel);
                                Cgpu = gpuArray(C);
                                alfa0gpu = gpuArray(alfa0);
                                for jj=1:size(EoutTxAux,2)
                                    U1tgpu = gpuArray(U1t(:,jj));
                                    U2tgpu = gpuArray(U2t(:,jj));
                                    EoutTxGpu = gpuArray(EoutTxAux(:,jj));
                                    freqGHzGpu = gpuArray(freqGHz(:,jj));
                                    [EoutMod,~]=MZMGT(freqGHzGpu,EoutTxGpu,U1tgpu,U2tgpu,Lgpu,U0gpu,U_pi1gpu,U_pi2gpu,noptgpu,nelgpu,Cgpu,alfa0gpu);%Modulating individual carriers
                                    EoutModTem(:,jj,ObsCarrPos(kk)) = gather(EoutMod);%Adding the current Modulated output to the final OutPut
                                    clear U1tgpu U2tgpu EoutTxGpu freqGHzGpu EoutMod
                                end
                                clear Lgpu U0gpu U_pi1gpu U_pi2gpu noptgpu nelgpu Cgpu alfa0gpu;
                            else
                                U.U1t = U1t;
                                U.U2t = U2t;
                                [EoutMod,~]=MZM(freqGHz,EoutTxAux,U,L,U0,U_pi1,U_pi2,nopt,nel,C,alfa0);%Modulating individual carriers
                                EoutModTem(:,:,ObsCarrPos(kk)) = EoutMod;%Adding the current Modulated output to the final OutPut
                            end
                            %[EoutMod,~]=MZM(freqGHz,EoutTxAux,U,L,U0,U_pi1,U_pi2,nopt,nel,C,alfa0);%Modulating individual carriers
                            %EoutModTem(:,:,ObsCarrPos(kk)) = EoutMod;%Adding the current Modulated output to the final OutPut
                            if (kk==2)&&(PlotingThis)
                                figure;hold all;grid on;
                                plot(f(:,1),20*log10(abs(fftshift(fft(EoutMod(:,1))./length(EoutMod(:,1))))));
                                axis([20e9 110e9 -200 0]);
                                set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                            end
                        end
                        if IfftOrSum==1                                                   %Use the OIFFT
                            if size(EoutModTem,3)>1                                    %Just if the signal has more than 1 dimension
                                if UsingGpu==1
                                    Tgpu = gpuArray(T);
                                    MaxStagGpu = gpuArray(MaxNumStagT);
                                    EoutMod = zeros(size(EoutModTem,1),size(EoutModTem,2));
                                    for jj=1:size(EoutModTem,2)
                                        EoutModTemGpu = gpuArray(EoutModTem(:,jj,:));
                                        [EoutAux1Gpu,~,~]=OpticalIFFTGT(fgpu,Tgpu,MaxStagGpu,EoutModTemGpu);
                                        EoutMod(:,jj) = gather(EoutAux1Gpu);
                                        clear EoutModTemGpu EoutAux1Gpu;
                                    end
                                    clear MaxStagGpu Tgpu;
                                else
                                    [EoutMod,~,~] = OpticalIFFTNT(f,T,MaxNumStagT,EoutModTem);
                                end
                            else
                                EoutMod = EoutModTem;                                  %otherwise just let the signal pass
                            end
                        else                                                           %Use a simple adder which don't display a better result
                            if size(EoutModTem,3)>1
                                EoutMod = sum(EoutModTem,3);
                            else
                                EoutMod = EoutModTem;
                            end
                        end
                    case 'DPSK'
                        %%            Generate the data DPSK
                        EoutModTem = zeros(NPPB*Nb,CurTesSiz,NumCarr);
                        TxDataMat = zeros(NumCarr,NbDPSK*CurTesSiz);
                        % 						EoutModTem = 0;
                        for kk=InitCarrUp:CarrPass:NumCarr%Generating different data for each carrier
                            if CarrUsedUp(kk)
                                % Frist it is chosen to transmite a random information...
                                TxData = (randi(2,NbDPSK,CurTesSiz)-1);%Creating the data stream to be transmited
                                TxData(1:JusPos,:) = JusVal;%Adding the Justification simble at the begining of the stream to sincronize received data frame
                                TxData(end-(JusPos-1):end,:) = JusValEnd;%Adding the Justification simble at the end of the stream to sincronize received data frame
                                TxDataMat(kk,:) = TxData(:);%Storring the transmited information for latter evaluation
                                %Once that we have the data to be transmited it needs to be converted to
                                %the eletrical signal that hereafter will modulate our carrier through the
                                %MZM-I
                                U1t = DpskEncodEqT(TxData);
                                U1t(U1t==1) =  1.9;
                                U1t(U1t==0) = -1.9;
                                U1t = rectpulse(U1t,NPPB);
                                % Adding CP to the data
                                if AddCP==1
                                    TxAuxB = reshape(U1t,NPPB,NbDPSK*CurTesSiz);
                                    if SetCpSampZer==1
                                        TxAuxA = [zeros(NumAmosCP,size(TxAuxB,2));TxAuxB;zeros(NumAmosCP,size(TxAuxB,2))];
                                    else
                                        %TxAuxA = [TxAuxB(1:NumAmosCP,:);TxAuxB;TxAuxB(end-(NumAmosCP-1):end,:)];
                                        TxAuxA = [flipud(TxAuxB(1:NumAmosCP,:));TxAuxB;flipud(TxAuxB(end-(NumAmosCP-1):end,:))];
                                    end
                                    TxAux = reshape(TxAuxA,(2*NumAmosCP+NPPB)*NbDPSK,CurTesSiz);
                                    U1t = [TxAux;TxAux(end-(StuffSampels-1):end,:)];
                                end
                                
                                %Because of reasons, the real world does not have components
                                %capable instantly change the voltage level. Thus, to emulate
                                %this behaviour the eletrical PAM signal pass through a
                                %gaussian filter that will remove the frequencies of higher
                                %order, which will result in small slope in the level variation
                                
                                [BitFilt,~] = FiltroGaussiano(f,BWD,CenFeq,FiltOrd);           %Creating filter for conformation of the input information
                                BitFilt     = fftshift(BitFilt);                               %Doing a shift on the vector for matching the transmited data
                                %                 TxSig       = ifft(fft(SigTx).*BitFilt);                       %Conforming the information and Creating the modulation signal
                                %                                 U1t       = U1t;                                   %Adding a gain to the eletrical signal
                                %             PrintInfo(Ploting*2,TxSig,T,NPPB);
                                %                 TxSigMat(kk,:) = TxSig;                                        %Storing the eletrical signal for further evaluation
                                %                 a=a+4;
                                
                                EoutTxAux = EoutAux1(:,:,VetThisCarr==ObsCarrUsed(kk));
                                %Assigning the eletrical signal to one drive of the MZM -
                                %The aplitude of the signal can be controlled by the
                                %variable DatGai, which can be understood as an gain for
                                %the eletrical signal or an atenuation. The second signal
                                %will be similar with the only difference a phase shift of
                                %pi.
                                %U.U1t  = U1t;
                                U2t  = exp(-1j*pi).*U1t;
                                %As both signals will have the mostrly the same
                                %characteristics with the only difference the phase shift
                                %of 180 degress. The MZM-I will be working on the Push-Pull
                                %configuration. It is necessary to reduce the Chirp noise
                                %to zero.
                                if UsingGpu==1
                                    Lgpu = gpuArray(L);
                                    U0gpu = gpuArray(U0);
                                    U_pi1gpu = gpuArray(U_pi1);
                                    U_pi2gpu = gpuArray(U_pi2);
                                    noptgpu = gpuArray(nopt);
                                    nelgpu = gpuArray(nel);
                                    Cgpu = gpuArray(C);
                                    alfa0gpu = gpuArray(alfa0);
                                    for jj=1:size(EoutTxAux,2)
                                        U1tgpu = gpuArray(U1t(:,jj));
                                        U2tgpu = gpuArray(U2t(:,jj));
                                        EoutTxGpu = gpuArray(EoutTxAux(:,jj));
                                        freqGHzGpu = gpuArray(freqGHz(:,jj));
                                        [EoutMod,~]=MZMGT(freqGHzGpu,EoutTxGpu,U1tgpu,U2tgpu,Lgpu,U0gpu,U_pi1gpu,U_pi2gpu,noptgpu,nelgpu,Cgpu,alfa0gpu);%Modulating individual carriers
                                        EoutModTem(:,jj,ObsCarrPos(kk)) = gather(EoutMod);%Adding the current Modulated output to the final OutPut
                                        clear U1tgpu U2tgpu EoutTxGpu freqGHzGpu EoutMod
                                    end
                                    clear Lgpu U0gpu U_pi1gpu U_pi2gpu noptgpu nelgpu Cgpu alfa0gpu;
                                else
                                    U.U1t = U1t;
                                    U.U2t = U2t;
                                    [EoutMod,~]=MZM(freqGHz,EoutTxAux,U,L,U0,U_pi1,U_pi2,nopt,nel,C,alfa0);%Modulating individual carriers
                                    EoutModTem(:,:,ObsCarrPos(kk)) = EoutMod;%Adding the current Modulated output to the final OutPut
                                end
                                %[EoutMod,~]=MZM(freqGHz,EoutMod,U,L,U0,U_pi1,U_pi2,nopt,nel,C,alfa0);    %Modulating individual carriers
                                %EoutModTem(:,:,ObsCarrPos(kk)) = EoutMod;%Adding the current Modulated output to the final OutPut
                            else
                                EoutMod = EoutAux1(:,:,VetThisCarr==ObsCarrPos(kk));
                                EoutModTem(:,:,ObsCarrPos(kk)) = EoutMod;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                            end
                        end
                        if IfftOrSum==1
                            if size(EoutModTem,3)>1
                                if UsingGpu==1
                                    Tgpu = gpuArray(T);
                                    MaxStagGpu = gpuArray(MaxNumStagT);
                                    EoutMod = zeros(size(EoutModTem,1),size(EoutModTem,2));
                                    for jj=1:size(EoutModTem,2)
                                        EoutModTemGpu = gpuArray(EoutModTem(:,jj,:));
                                        [EoutAux1Gpu,~,~]=OpticalIFFTGT(fgpu,Tgpu,MaxStagGpu,EoutModTemGpu);
                                        EoutMod(:,jj) = gather(EoutAux1Gpu);
                                        clear EoutModTemGpu EoutAux1Gpu;
                                    end
                                    clear MaxStagGpu Tgpu;
                                else
                                    [EoutMod,~,~] = OpticalIFFTNT(f,T,MaxNumStagT,EoutModTem);
                                end
                            else
                                EoutMod = EoutModTem;
                            end
                        else
                            EoutMod = sum(EoutModTem,3);
                        end
                    case 'DQPSK'
                        %%            Generate the data DQPSK
                        EoutModTem = zeros(NPPB*Nb,CurTesSiz,NumCarr);
                        TxDataMat = zeros(2*NumCarr,(NbDQPSK/2)*CurTesSiz);
                        for kk=InitCarrUp:CarrPass:NumCarr%Generating different data for each carrier
                            if CarrUsedUp(kk)
                                % Frist it is chosen to transmite a random information...
                                TxData = (randi(2,NbDQPSK,CurTesSiz)-1);%Creating the data stream to be transmited
                                TxData(1:JusPos,:) = JusVal;%Adding the Justification simble at the begining of the stream to sincronize received data frame
                                TxData(end-(JusPos-1):end,:) = JusValEnd;%Adding the Justification simble at the end of the stream to sincronize received data frame
                                PreBits  = zeros(2,CurTesSiz);%Seting the start point for the DQPSK maping
                                TxDataPos = linspace(1,NbDQPSK,NbDQPSK);%Auxiliar variable to split the TxData
                                [DataI,DataQ] = DqpskEncodEqT(TxData,PreBits);%Maping the income data to the DQPSK format
                                %Converting the I and Q components to the polirezed NRZ format
                                DataI(DataI==1) =  1.9;
                                DataI(DataI==0) = -1.9;
                                DataQ(DataQ==1) =  1.9;
                                DataQ(DataQ==0) = -1.9;
                                TxOdd = TxData(logical(mod(TxDataPos,2)),:);%Spliting the information of odd positions
                                TxEven = TxData(~(mod(TxDataPos,2)),:);%Spliting the information of even positions
                                TxDataMat(kk,:) = TxOdd(:);%Storring the transmited information for latter evaluation
                                TxDataMat(kk+NumCarr,:) = TxEven(:);%Storring the transmited information for latter evaluation
                                %Once that we have the data to be transmited it needs to be converted to
                                %the eletrical signal that hereafter will modulate our carrier through the
                                %MZM-I
                                SigI = rectpulse(DataI,NPPB);
                                SigQ = rectpulse(DataQ,NPPB);
                                
                                % Adding CP to the data
                                if AddCP
                                    TxAux1 = reshape(SigI,NPPB,CurTesSiz*(NbDQPSK/2));
                                    if SetCpSampZer==1
                                        TxAux2 = [zeros(NumAmosCP,size(TxAux1,2));TxAux1;zeros(NumAmosCP,size(TxAux1,2))];
                                    else
                                        %TxAux2 = [TxAux1(1:NumAmosCP,:);TxAux1;TxAux1(end-(NumAmosCP-1):end,:)];
                                        TxAux2 = [flipud(TxAux1(1:NumAmosCP,:));TxAux1;flipud(TxAux1(end-(NumAmosCP-1):end,:))];
                                    end
                                    TxAux1 = reshape(TxAux2,(2*NumAmosCP+NPPB)*(NbDQPSK/2),CurTesSiz);
                                    SigI = [TxAux1;TxAux1(end-(StuffSampels-1):end,:)];
                                    
                                    TxAux3 = reshape(SigQ,NPPB,CurTesSiz*(NbDQPSK/2));
                                    if SetCpSampZer==1
                                        TxAux4 = [zeros(NumAmosCP,size(TxAux3,2));TxAux3;zeros(NumAmosCP,size(TxAux3,2))];
                                    else
                                        %TxAux4 = [TxAux3(1:NumAmosCP,:);TxAux3;TxAux3(end-(NumAmosCP-1):end,:)];
                                        TxAux4 = [flipud(TxAux3(1:NumAmosCP,:));TxAux3;flipud(TxAux3(end-(NumAmosCP-1):end,:))];
                                    end
                                    TxAux5 = reshape(TxAux4,(2*NumAmosCP+NPPB)*(NbDQPSK/2),CurTesSiz);
                                    SigQ = [TxAux5;TxAux5(end-(StuffSampels-1):end,:)];
                                end
                                
                                %Because of reasons, the real world does not have components
                                %capable instantly change the voltage level. Thus, to emulate
                                %this behaviour the eletrical PAM signal pass through a
                                %gaussian filter that will remove the frequencies of higher
                                %order, which will result in small slope in the level variation
                                
                                [BitFilt,~] = FiltroGaussiano(f,BWD,CenFeq,FiltOrd);           %Creating filter for conformation of the input information
                                BitFilt = fftshift(BitFilt);                                   %Doing a shift on the vector for matching the transmited data
                                %                 TxSigI = ifft(fft(SigI).*BitFilt);                             %Conforming the information and Creating the modulation signal
                                %                 TxSigQ = ifft(fft(SigQ).*BitFilt);                             %Conforming the information and Creating the modulation signal
                                %TxSigI = SigI;                             %Conforming the information and Creating the modulation signal
                                %TxSigQ = SigQ;                             %Conforming the information and Creating the modulation signal
                                %             PrintInfo(Ploting*3,TxSigI,T,NPPB,TxSigQ);
                                %                 TxSigMat(kk,:)         = TxSigI;                               %Storing the eletrical signal for further evaluation
                                %                 TxSigMat(kk+NumCarr,:) = TxSigQ;                               %Storing the eletrical signal for further evaluation
                                EoutMod = EoutAux1(:,:,VetThisCarr==ObsCarrUsed(kk));
                                [EoutMod] = IqMod(EoutMod,SigI,SigQ,Vpi,V0);
                                %                 a=a+4;
                                EoutModTem(:,:,ObsCarrPos(kk)) = EoutMod;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                            else
                                EoutMod = EoutAux1(:,:,VetThisCarr==ObsCarrPos(kk));
                                EoutModTem(:,:,ObsCarrPos(kk)) = EoutMod;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                            end
                        end
                        if IfftOrSum==1
                            if size(EoutModTem,3)>1
                                if UsingGpu==1
                                    Tgpu = gpuArray(T);
                                    MaxStagGpu = gpuArray(MaxNumStagT);
                                    EoutMod = zeros(size(EoutModTem,1),size(EoutModTem,2));
                                    for jj=1:size(EoutModTem,2)
                                        EoutModTemGpu = gpuArray(EoutModTem(:,jj,:));
                                        [EoutAux1Gpu,~,~]=OpticalIFFTGT(fgpu,Tgpu,MaxStagGpu,EoutModTemGpu);
                                        EoutMod(:,jj) = gather(EoutAux1Gpu);
                                        clear EoutModTemGpu EoutAux1Gpu;
                                    end
                                    clear MaxStagGpu Tgpu;
                                else
                                    [EoutMod,~,~] = OpticalIFFTNT(f,T,MaxNumStagT,EoutModTem);
                                end
                            else
                                EoutMod = EoutModTem;
                            end
                        else
                            EoutMod = sum(EoutModTem,3);
                        end
                    case '4PAM'
                        %%        Generate the data 4PAM
                        EoutModTem = zeros(NPPB*Nb,CurTesSiz,NumCarr);
                        TxDataMat = zeros(NumCarr,Nb4Pam*CurTesSiz);
                        for kk=InitCarrUp:CarrPass:NumCarr%Generating different data for each carrier
                            % Frist it is chosen to transmite a random information...
                            if CarrUsedUp(kk)
                                TxData = (randi(2,Nb4Pam,CurTesSiz)-1);%Creating the data stream to be transmited
                                TxData(1:JusPos,:) = JusVal;%Adding the Justification simble at the begining of the stream to sincronize received data frame
                                TxData(end-(JusPos-1):end,:) = JusValEnd;%Adding the Justification simble at the end of the stream to sincronize received data frame
                                TxDataMat(kk,:) = TxData(:);                                      %Storring the transmited information for latter evaluation
                                %Once that we have the data to be transmited it needs to be
                                %converted to the eletrical signal that hereafter will modulate
                                %our carrier through the MZM-I
                                
                                %Selecting if the PAM will be done on Electrical or Optical
                                %domain
                                if Which4PAM==1%Selecting if the PAM will be done on Electrical or Optical domain
                                    %  Good for others kms
                                    [Phi1,Phi2] = Maping4PamIqT(TxData,Vmin,Vmax,ModSchem,FiberLength,SetCpSampZer);%Generating the eletrical signal for the optical PAM4
                                else
                                    [Phi1,Phi2] = Maping4PamT(TxData,VPI,Polirized,MaxAmp4PAM);%Generating the eletrical signal for the electrical PAM4
                                end
                                %The signal generated are not yet with the same number of
                                %samples as the OFCS loaded. These nexte lines do oversampling
                                %of the electrical signal.
                                TxSig1 = rectpulse(Phi1,NPPB);
                                TxSig2 = rectpulse(Phi2,NPPB);
                                
                                %Thus, if it would be required to add cycle prefix the number
                                %of samples per symbol needs to change as well as some
                                %adjustments needs to be done for the new signal match in size
                                %with the size of the vector time. This problem just exist on
                                %simulation, at practice the main point is the syncronism of
                                %the signals.
                                % Adding CP to the data
                                if AddCP==1
                                    TxAux1 = reshape(TxSig1,NPPB,CurTesSiz*(Nb4Pam/2));
                                    if SetCpSampZer==1
                                        TxAux2 = [zeros(NumAmosCP,size(TxAux1,2));TxAux1;zeros(NumAmosCP,size(TxAux1,2))];
                                    else
                                        TxAux2 = [flipud(TxAux1(1:NumAmosCP,:));TxAux1;flipud(TxAux1(end-(NumAmosCP-1):end,:))];
                                    end
                                    TxAux3 = reshape(TxAux2,(2*NumAmosCP+NPPB)*(Nb4Pam/2),CurTesSiz);
                                    TxSig1 = [TxAux3;TxAux3(end-(StuffSampels-1):end,:)];
                                    TxAux4 = reshape(TxSig2,NPPB,CurTesSiz*(Nb4Pam/2));
                                    if SetCpSampZer==1
                                        TxAux5 = [zeros(NumAmosCP,size(TxAux4,2));TxAux4;zeros(NumAmosCP,size(TxAux4,2))];
                                    else
                                        TxAux5 = [flipud(TxAux4(1:NumAmosCP,:));TxAux4;flipud(TxAux4(end-(NumAmosCP-1):end,:))];
                                    end
                                    TxAux6 = reshape(TxAux5,(2*NumAmosCP+NPPB)*(Nb4Pam/2),CurTesSiz);
                                    TxSig2 = [TxAux6;TxAux6(end-(StuffSampels-1):end,:)];
                                end
                                %Because of reasons, the real world does not have components
                                %capable instantly change the voltage level. Thus, to emulate
                                %this behaviour the eletrical PAM signal pass through a
                                %gaussian filter that will remove the frequencies of higher
                                %order, which will result in small slope in the level variation
                                
                                [BitFilt,~]            = FiltroGaussiano(f,BWD,CenFeq,FiltOrd);%Creating filter for conformation of the input information
                                BitFilt                = fftshift(BitFilt);                    %Doing a shift on the vector for matching the transmited data
                                %                 TxSig1                 = ifft(fft(TxSig1).*BitFilt);           %Conforming the information and Creating the modulation signal
                                %                 TxSig2                 = ifft(fft(TxSig2).*BitFilt);           %Conforming the information and Creating the modulation signal
                                %             TxSigMat(kk,:)         = TxSig1;                               %Storing the eletrical signal for further evaluation
                                %             TxSigMat(kk+NumCarr,:) = TxSig2;                               %Storing the eletrical signal for further evaluation
                                
                                EoutTxAux = EoutAux1(:,:,VetThisCarr==ObsCarrUsed(kk));
                                U1t  = TxSig1;                                               %Assigning the electrical signal to one drive of the MZM
                                U2t  = TxSig2;                                               %Assigning the electrical signal to another drive of the MZM
                                
                                %to zero.
                                if UsingGpu==1
                                    Lgpu = gpuArray(L);
                                    U0gpu = gpuArray(U0Up);
                                    U_pi1gpu = gpuArray(U_pi1);
                                    U_pi2gpu = gpuArray(U_pi2);
                                    noptgpu = gpuArray(nopt);
                                    nelgpu = gpuArray(nel);
                                    Cgpu = gpuArray(C);
                                    alfa0gpu = gpuArray(alfa0);
                                    for jj=1:size(EoutTxAux,2)
                                        U1tgpu = gpuArray(U1t(:,jj));
                                        U2tgpu = gpuArray(U2t(:,jj));
                                        EoutTxGpu = gpuArray(EoutTxAux(:,jj));
                                        freqGHzGpu = gpuArray(freqGHz(:,jj));
                                        [EoutMod,~]=MZMGT(freqGHzGpu,EoutTxGpu,U1tgpu,U2tgpu,Lgpu,U0gpu,U_pi1gpu,U_pi2gpu,noptgpu,nelgpu,Cgpu,alfa0gpu);%Modulating individual carriers
                                        EoutModTem(:,jj,ObsCarrPos(kk)) = gather(EoutMod);%Adding the current Modulated output to the final OutPut
                                        clear U1tgpu U2tgpu EoutTxGpu freqGHzGpu EoutMod
                                    end
                                    clear Lgpu U0gpu U_pi1gpu U_pi2gpu noptgpu nelgpu Cgpu alfa0gpu;
                                else
                                    U.U1t = U1t;
                                    U.U2t = U2t;
                                    [EoutMod,~]=MZM(freqGHz,EoutTxAux,U,L,U0Up,U_pi1,U_pi2,nopt,nel,C,alfa0);%Modulating individual carriers
                                    EoutModTem(:,:,ObsCarrPos(kk)) = EoutMod;%Adding the current Modulated output to the final OutPut
                                end
                                %if Which4PAM
                                    %if ModSchem
                                        %                         [EoutModAux] = IqMod4Pam (EoutT,U.U1t,U.U2t,U_pi2,Vbias);
                                        %[EoutMod] = MZM(freqGHz,EoutMod,U,L,U0Up,U_pi1,U_pi2,nopt,nel,C,alfa0);
                                    %else
                                        %[EoutMod] = MZM(freqGHz,EoutMod,U,L,U0Up,U_pi1,U_pi2,nopt,nel,C,alfa0);
                                    %end
                                %else
                                    %[EoutMod] = MZM(freqGHz,EoutMod,U,L,U0Up,U_pi1,U_pi2,nopt,nel,C,alfa0);
                                %end
                                %                 EoutMod = EoutMod + EoutModAux;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                                %EoutModTem(:,:,ObsCarrPos(kk)) = EoutMod;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                            else
                                EoutMod = EoutAux1(:,:,VetThisCarr==ObsCarrPos(kk));
                                EoutModTem(:,:,ObsCarrPos(kk)) = EoutMod;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                            end
                        end
                        if IfftOrSum==1
                            if size(EoutModTem,3)>1
                                if UsingGpu==1
                                    Tgpu = gpuArray(T);
                                    MaxStagGpu = gpuArray(MaxNumStagT);
                                    EoutMod = zeros(size(EoutModTem,1),size(EoutModTem,2));
                                    for jj=1:size(EoutModTem,2)
                                        EoutModTemGpu = gpuArray(EoutModTem(:,jj,:));
                                        [EoutAux1Gpu,~,~]=OpticalIFFTGT(fgpu,Tgpu,MaxStagGpu,EoutModTemGpu);
                                        EoutMod(:,jj) = gather(EoutAux1Gpu);
                                        clear EoutModTemGpu EoutAux1Gpu;
                                    end
                                    clear MaxStagGpu Tgpu;
                                else
                                    [EoutMod,~,~] = OpticalIFFTNT(f,T,MaxNumStagT,EoutModTem);
                                end
                            else
                                EoutMod = EoutModTem;
                            end
                        else
                            EoutMod = sum(EoutModTem,3);
                        end
                    otherwise
                        %%        Generate the data OOK
                        EoutModTem = zeros(NPPB*Nb,CurTesSiz,NumCarr);
                        TxDataMat = zeros(NumCarr,(Nb-NumBitDesc)*CurTesSiz);
                        for kk=InitCarrUp:CarrPass:NumCarr%Generating different data for each carrier
                            % Frist it is chosen to transmite just one high pulse for
                            %testing the channel...
                            if CarrUsedUp(kk)
                                TxData = (randi(2,Nb-NumBitDesc,CurTesSiz)-1);%Creating Random Information that will be loaded in each individual subcarrier
                                TxData(1:JusLen,:) = JusVal;%Making the First 4 bits equal to zero
                                TxData(end-(JusLen-1):end,:) = JusVal;%Making the Last 4 bits equal to zero
                                TxDataMat(kk,:) = TxData(:);%Storring the transmited information for latter evaluation
                                if NRZPolarity
                                    TxData(TxData==0) = NrzMin;
                                    TxData(TxData==1) = NrzMax;
                                end
                                %The signal generated are not yet with the same number of
                                %samples as the OFCS loaded. These nexte lines do oversampling
                                TxDataRes = rectpulse(TxData,NPPB);                            %Changing the length of the Data acordingly with the time vector
                                
                                %Thus, if it would be required to add cycle prefix the number
                                %of samples per symbol needs to change as well as some
                                %adjustments needs to be done for the new signal match in size
                                %with the size of the vector time. This problem just exist on
                                %simulation, at practice the main point is the syncronism of
                                %the signals.
                                % Adding CP to the data
                                if AddCP==1
                                    TxAux = reshape(TxDataRes,NPPB,CurTesSiz*(Nb-NumBitDesc));
                                    if SetCpSampZer==1
                                        TxAux1 = [zeros(NumAmosCP,size(TxAux,2));TxAux;zeros(NumAmosCP,size(TxAux,2))];
                                    else
                                        TxAux1 = [flipud(TxAux(1:NumAmosCP,:));TxAux;flipud(TxAux(end-(NumAmosCP-1):end,:))];
                                    end
                                    TxAux2 = reshape(TxAux1,(2*NumAmosCP+NPPB)*(Nb-NumBitDesc),CurTesSiz);
                                    U1t = [TxAux2;TxAux2(end-(StuffSampels-1):end,:)];
                                end
                                %Because of reasons, the real world does not have components
                                %capable instantly change the voltage level. Thus, to emulate
                                %this behaviour the eletrical PAM signal pass through a
                                %gaussian filter that will remove the frequencies of higher
                                %order, which will result in small slope in the level variation
                                
                                [BitFilt,~] = FiltroGaussiano(f,BWD,CenFeq,FiltOrd);           %Creating filter for conformation of the input information
                                BitFilt = fftshift(BitFilt);                                   %Doing a shift on the vector for matching the transmited data
                                %                 TxSig = ifft(fft(TxDataRes).*BitFilt);                         %Conforming the information and Creating the modulation signal
                                %U1t = TxDataRes;                         %Conforming the information and Creating the modulation signal
                                %             TxSigMat(kk,:) = TxSig;                                        %Storring the transmited information for latter evaluation
                                
                                EoutTxAux = EoutAux1(:,:,VetThisCarr==ObsCarrUsed(kk));
                                %Assigning the eletrical signal to one drive of the MZM -
                                %The aplitude of the signal can be controlled by the
                                %variable DatGai, which can be understood as an gain for
                                %the eletrical signal or an atenuation. The second signal
                                %will be similar with the only difference a phase shift of
                                %pi.
                                U1t  = DatGai.*U1t;
                                U2t  = DatGai.*exp(-1j*pi).*U1t;
                                %As both signals will have the mostrly the same
                                %characteristics with the only difference the phase shift
                                %of 180 degress. The MZM-I will be working on the Push-Pull
                                %configuration. It is necessary to reduce the Chirp noise
                                %to zero.
                                
                                if UsingGpu==1
                                    Lgpu = gpuArray(L);
                                    U0gpu = gpuArray(U0);
                                    U_pi1gpu = gpuArray(U_pi1);
                                    U_pi2gpu = gpuArray(U_pi2);
                                    noptgpu = gpuArray(nopt);
                                    nelgpu = gpuArray(nel);
                                    Cgpu = gpuArray(C);
                                    alfa0gpu = gpuArray(alfa0);
                                    for jj=1:size(EoutTxAux,2)
                                        U1tgpu = gpuArray(U1t(:,jj));
                                        U2tgpu = gpuArray(U2t(:,jj));
                                        EoutTxGpu = gpuArray(EoutTxAux(:,jj));
                                        freqGHzGpu = gpuArray(freqGHz(:,jj));
                                        [EoutMod,~]=MZMGT(freqGHzGpu,EoutTxGpu,U1tgpu,U2tgpu,Lgpu,U0gpu,U_pi1gpu,U_pi2gpu,noptgpu,nelgpu,Cgpu,alfa0gpu);%Modulating individual carriers
                                        EoutModTem(:,jj,ObsCarrPos(kk)) = gather(EoutMod);%Adding the current Modulated output to the final OutPut
                                        clear U1tgpu U2tgpu EoutTxGpu freqGHzGpu EoutMod
                                    end
                                    clear Lgpu U0gpu U_pi1gpu U_pi2gpu noptgpu nelgpu Cgpu alfa0gpu;
                                else
                                    U.U1t = U1t;
                                    U.U2t = U2t;
                                    [EoutMod,~]=MZM(freqGHz,EoutTxAux,U,L,U0,U_pi1,U_pi2,nopt,nel,C,alfa0);%Modulating individual carriers
                                    EoutModTem(:,:,ObsCarrPos(kk)) = EoutMod;%Adding the current Modulated output to the final OutPut
                                end
                                %[EoutMod,~]=MZM(freqGHz,EoutMod,U,L,U0,U_pi1,U_pi2,nopt,nel,C,alfa0);    %Modulating individual carriers
                                %                 EoutMod = EoutMod + EoutModAux;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                                %EoutModTem(:,:,ObsCarrPos(kk)) = EoutMod;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                                %% Taking the sampling the EVM meassurement
                                %PosAuxEout = (2*NumAmosCP+NPPB)/2:(2*NumAmosCP+NPPB)...
                                %:(length(EoutMod)-StuffSampels);%Varriable respossible to take just the samples at the middle of the symbol
                                %Ix         = EoutMod.*conj(EoutMod);             %Recovering the signal that will be transmited
                                %Ix         = ifft(fft(Ix).*EvmFilt);                   %Removing higher signal generated by the receiving process
                                %Ix         = Ix - min(Ix);
                                %Ix         = Ix./max(abs(Ix));                         %Normalizing the reference
                                %IxAux      = Ix(PosAuxEout);                           %Normalizing the reference
                                %a=a+0;
                                %TxSymbAmos(CurrentTest,ObsCarrPos==kk,:) = IxAux;%RxSymbAmosUp = [];
                                %EvmMatRef(ObsCarrPos==kk,:) = IxAux;                   %Taking just the middle samples as references
                            else
                                EoutMod = EoutAux1(:,:,VetThisCarr==ObsCarrPos(kk));
                                EoutModTem(:,:,ObsCarrPos(kk)) = EoutMod;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                            end
                        end
                        if IfftOrSum==1
                            if size(EoutModTem,3)>1
                                if UsingGpu==1
                                    Tgpu = gpuArray(T);
                                    MaxStagGpu = gpuArray(MaxNumStagT);
                                    EoutMod = zeros(size(EoutModTem,1),size(EoutModTem,2));
                                    for jj=1:size(EoutModTem,2)
                                        EoutModTemGpu = gpuArray(EoutModTem(:,jj,:));
                                        [EoutAux1Gpu,~,~]=OpticalIFFTGT(fgpu,Tgpu,MaxStagGpu,EoutModTemGpu);
                                        EoutMod(:,jj) = gather(EoutAux1Gpu);
                                        clear EoutModTemGpu EoutAux1Gpu;
                                    end
                                    clear MaxStagGpu Tgpu;
                                else
                                    [EoutMod,~,~] = OpticalIFFTNT(f,T,MaxNumStagT,EoutModTem);
                                end
                            else
                                EoutMod = EoutModTem;
                            end
                        else
                            EoutMod = sum(EoutModTem,3);
                        end
                end
                
                %%   Transmission of the OFDM Symble through a channel
                % Having data stored and ready to be sent to end user. At the stage this
                % script is responsible to chose the medium where this signal will travel.
                % It may be withing an optical fiber or Back-toBack transmission.
                switch Medium
                    case 'B2B'
                        EoutRec = EoutMod;
                    case 'Fiber'
                        [EoutRec,~,~]=Fibra_Monomodo1(t,EoutMod,lambdac,T,FiberLength,0,f);
                        %         PrintInfo(Ploting*11,EoutRec.*conj(EoutRec)./max(abs(EoutRec.*...
                        %                                                    conj(EoutRec))),T,NPPB);
                    otherwise
                        EoutRec = EoutMod;
                end
                
                %%                Datum Upstream Reception
                %This section is responsible to receive the information from
                %users
                %% DatumUpstrRec;
                
                %%            Recovering the signal
                % At this point the FFT will be inplemented. The received signal need to be
                % given as an parameter. This following function recursively implement the
                % FFT Operation. One may need to investigate how would work a receptor with
                % filters insted of a OFFT. The second implementation is an array of
                % filters.
                if FFTSplit==1
                    if UsingGpu==1
                        Tgpu = gpuArray(T);
                        MaxStagGpu = gpuArray(MaxNumStag);
                        EoutAux1 = zeros(size(EoutRec,1),size(EoutRec,2),length(ObsCarrUsed));
                        for jj=1:size(EoutRec,2)
                            EoutGpu = gpuArray(EoutRec(:,jj));
                            [EoutAux1Gpu,~,VetThisCarrGpu]=OpticalFFTGT(fgpu,Tgpu,MaxStagGpu,EoutGpu);
                            EoutAux1(:,jj,:) = gather(EoutAux1Gpu);
                            VetThisCarr = gather(VetThisCarrGpu);
                            clear EoutGpu EoutAux1Gpu VetThisCarrGpu;
                        end
                        clear MaxStagGpu Tgpu;
                    else
                        [EoutAux1,~,VetThisCarr]=OpticalFFTNT(f,T,MaxNumStag,EoutRec);
                    end
                else
                    VetThisCarr = ObsCarrPos;
                    EoutAux1    = SelectEachCarrier(EoutRec,NumCarr,f,fin,1.0*fc,5,fc);
                end
                for kk=1:size(EoutAux1,3)
                    [~,CarrRecPowUp(CurrentTest,kk)] = MeasPower(EoutAux1(:,:,VetThisCarr==ObsCarrUsed(kk)),t);
                end
                %%  Ploting some results for qualitative analizes
                % PrintInfo(Ploting*12,f,db(abs((fftshift(fft(EoutRec))))/length(EoutRec)),EoutAux1,RefCarr*fc);
                %%      Recovering the Data
                %This is basicaly the final step, so far, in this process. Here, the
                %transmited signal will be received and processed to recover the actual
                %data within it.
                switch Modulation
                    case 'OFDM'
                        for ThisCarr=InitCarrUp:CarrPass:NumCarr
                            Ix = EoutAux1(:,:,VetThisCarr==ObsCarrUsed(ThisCarr));
                            switch Medium
                                case 'Fiber'
                                    Ix = ifft(fft(Ix).*exp(1j*2*pi*f*(...
                                        FiberDelay(ThisCarr)*Ta)));
                                otherwise
                            end
                            %             switch OfdMod                                                  %Modulation Tx signal for future comparison
                            %                 case 'qam'
                            %                     sn=((1/2)*(((1/2)^MaxNumStag)-1))/((1/2)-1);
                            %                     DiffPos = round((sn/(2/2^IfftOrSum))*T/Ta);
                            %                     EoutAux = [EoutAux(DiffPos+1:end) EoutAux(1:DiffPos)];
                            %                 otherwise
                            %             end
                            %The current incoming signal is them converted from the optical
                            %domain to the eletrical domain with the help of an photo
                            %detector.
                            Ix =Ix.*conj(Ix);
                            [BitFilt,~] = FiltroGaussiano(f,BWD,CenFeq,FiltOrd);%Creating filter for selection of the received signal
                            BitFilt = fftshift(BitFilt);%Shifting the filter for matching the received signal
                            Ix = ifft(fft(Ix).*BitFilt);%Filtering the signal for better detection
                            Ix = Ix - min(min(Ix));%Removing any off-set that may exist
                            Ix = Ix.*NormFact;
                            NorAux = max(Ix);
                            NorAux = repmat(NorAux,size(Ix,1),1);
                            Ix = (Ix./NorAux);%Normalizing the eletrical signal (amplifying what is needed)
                            
                            if ReceptorNoise==1%Verify whether the noise is to be added or not
                                [~,SigPower] = MeasPower(Ix);
                                SigPower2 = SigPower-30;%20*log10(SigPowerI);
                                Ix = awgn(Ix,SNR,SigPower2);
                            end
                            if PlotingThis
                                figure;hold all;grid on;plot(f,20*log10(abs(fftshift(fft(Ix)./length(Ix)))));axis([-25e9 37.5e9 -200 0]);set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                            end
                            Ix = Ix(1+NPSTUf:end,:);
                            Ix = reshape(Ix,OvSam*NFFT,NumFra*CurTesSiz);
                            Ix = intdump(Ix,NPPOF);%Downsampling the income signal
                            SigRecepA = Ix(1+NPOFEX:end,:);
                            
                            switch SelModTp
                                case 'BP'
                                    if SelecGaus==1
                                        [ BitFiltEle ] = FiltroGaussiano(ff,OBw/2,0,1);%CarrOffSet*FramSpac,1);
                                    else
                                        [ BitFiltEle ] = Filtro_Retangular(OBw,0,ff);%CarrOffSet*FramSpac,ff);
                                    end
                                    BitFiltEle = fftshift(BitFiltEle);
                                    ModFilta = repmat(BitFiltEle.',1,NumFra*CurTesSiz);
                                    SigRecepB   = SigRecepA;%.*exp(-1j*2*pi*CarrOffSet*FramSpac*tta);
                                    if Ploting
                                        figure;
                                        hold all;
                                        plot(ff,20*log10(abs(fftshift(fft(SigRecepB)./length(SigRecepB)))));
                                        plot(ff,20*log10(abs(fftshift(BitFiltEle))));
                                    end
                                    if SelModFilt==1
                                        SigRecepB = ifft(fft(SigRecepB).*ModFilta);
                                    end
                                    SigRecep  = intdump(SigRecepB,OvSam);
                                    if Ploting
                                        plot(ff,20*log10(abs(fftshift(fft(SigRecepB)./length(SigRecepB)))));
                                        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                                    end
                                case 'AM'
                                    [BitFiltEle,~] = FiltroGaussiano(ff,OBw/2,0,1);           %Creating filter for selection of the received signal
                                    BitFiltEle = fftshift(BitFiltEle);
                                    ModFilta = repmat(BitFiltEle.',1,NumFra*CurTesSiz);
                                    SigRecepB   = SigRecepA.*cos(-2*pi*Ofc*tta);
                                    if Ploting
                                        figure;
                                        hold all;
                                        plot(ff,20*log10(abs(fftshift(fft(SigRecepB)./length(SigRecepB)))));
                                        plot(ff,20*log10(abs(fftshift(BitFiltEle))));
                                    end
                                    if SelModFilt==1
                                        SigRecepB = ifft(fft(SigRecepB).*ModFilta);
                                    end
                                    SigRecep  = intdump(SigRecepB,OvSam);
                                    if Ploting
                                        plot(ff,20*log10(abs(fftshift(fft(SigRecepB)./length(SigRecepB)))));
                                        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                                    end
                                    %                         x_ofdm = x_ofdm_ssb;
                                case 'AMSSB'
                                    [BitFiltEle,~] = FiltroGaussiano(ff,OBw/2,0,1);           %Creating filter for selection of the received signal
                                    BitFiltEle = fftshift(BitFiltEle);
                                    ModFilta = repmat(BitFiltEle.',1,NumFra*CurTesSiz);
                                    SigRecepB   = SigRecepA.*exp(-1j*2*pi*Ofc*tta);
                                    if Ploting
                                        figure;
                                        hold all;
                                        plot(ff,20*log10(abs(fftshift(fft(SigRecepB)./length(SigRecepB)))));
                                        plot(ff,20*log10(abs(fftshift(BitFiltEle))));
                                    end
                                    if SelModFilt==1
                                        SigRecepB = ifft(fft(SigRecepB).*ModFilta);
                                    end
                                    SigRecep  = intdump(SigRecepB,OvSam);              %Over sampling
                                    if Ploting
                                        plot(ff,20*log10(abs(fftshift(fft(SigRecepB)./length(SigRecepB)))));
                                        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                                    end
                                otherwise
                                    SigRecepB = demod(SigRecepA,Ofc,Ofs,'pm');
                                    SigRecep  = intdump(SigRecepB,OvSam);              %Over sampling
                            end
                            %             SigRecep = reshape(SigRecep,NumFra,NFFT);                      %Reshaping the signal when multiple frames were transmited
                            SigRecep1 = fft(SigRecep,NFFT);                              %FFT demultiplexing process
                            %             SigRecep2 = SigRecep1(2+Zt-ZtC:Ns+1+Zt-ZtC,:);                               %Taking the actual data
                            if UsingHermitian
                                SigRecep2 = SigRecep1(2+Zt-ZtC:Ns+1+Zt-ZtC,:);                               %Taking the actual data
                            else
                                %     SigRecep2 = SigRecep1(1+(NFFT-size(TxSymbH,1))/2:size(TxSymbH,1)+(NFFT-size(TxSymbH,1))/2,:);
                                SigRecep2 = [SigRecep1(1:size(TxSymbH,1)/2,:);SigRecep1(end+1-size(TxSymbH,1)/2:end,:)];
                            end
                            SigRecep3 = SigRecep2;
                            %             SigRecep3 = SigRecep3(:).';                                    %Changing from parralel to serial
                            %             SigRecep3 = SigRecep2.';
                            if PlotingThis
                                switch OfdMod                                                  %Modulation Tx signal for future comparison
                                    case 'qam'
                                        TxDataA = TxDataMat(ThisCarr,:);
                                        TxDataA = reshape(TxDataA,NumFra*CurTesSiz,length(TxDataA)/(NumFra*CurTesSiz));
                                        TxDataA = TxDataA.';
                                        UniDmtMve = unique(DmtMve);
                                        UniDmtMve = fliplr(UniDmtMve);
                                        DaSiAn =  0;
                                        for CarDmtM=1:length(UniDmtMve)
                                            M = UniDmtMve(CarDmtM);
                                            DaSiAu = sum(ismember(DmtMve,M));
                                            TxSigToPlot(1+DaSiAn:DaSiAn + DaSiAu,:) = qammod(TxDataA(1+DaSiAn:DaSiAn + DaSiAu,:),M);%Saving data for future comparison
                                            DaSiAn = DaSiAn + DaSiAu;
                                            %TxSigToPlot(1+DaSiAu*(CarDmtM-1):DaSiAu*CarDmtM,:)  = qammod(TxDataA(1+DaSiAu*(CarDmtM-1):DaSiAu*CarDmtM,:),M);
                                        end
                                    otherwise
                                        TxDataA = TxDataMat(ThisCarr,:);
                                        TxDataA = reshape(TxDataA,NumFra*CurTesSiz,length(TxDataA)/(NumFra*CurTesSiz));
                                        TxDataA = TxDataA.';
                                        UniDmtMve = unique(DmtMve);
                                        UniDmtMve = fliplr(UniDmtMve);
                                        DaSiAn =  0;
                                        for CarDmtM=1:length(UniDmtMve)
                                            M = UniDmtMve(CarDmtM);
                                            DaSiAu = sum(ismember(DmtMve,M));
                                            TxSigToPlot(1+DaSiAn:DaSiAn + DaSiAu,:) = dpskmod(TxDataA(1+DaSiAn:DaSiAn + DaSiAu,:),M);%Saving data for future comparison
                                            DaSiAn = DaSiAn + DaSiAu;
                                            %TxSigToPlot(1+DaSiAu*(CarDmtM-1):DaSiAu*CarDmtM,:)  = dpskmod(TxDataA(1+DaSiAu*(CarDmtM-1):DaSiAu*CarDmtM,:),M);
                                        end
                                end
                                %         TxSigToPlot = dpskmod(TxDataMat(ThisCarr,2:end),M);
                                if PlotingThis
                                    figure;
                                    txcolor = [0.2 0 1];
                                    rxcolor = [1 0.4 0];
                                    hold all;
                                    plot(SigRecep3(:),'o','LineWidth',2,'color',rxcolor);
                                    plot(TxSigToPlot(:),'*','LineWidth',2,'color',txcolor);
                                    set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                                end
                            end
                            
                            switch OfdMod
                                case 'qam'
                                    RxSigOfdmNoEq(CurrentTest,ThisCarr,:) = SigRecep3(:).';
                                    equa = ChanelEqualizer(ObsCarrUsed(ThisCarr),:);
                                    equ = reshape(equa,length(equa)/(NumFra*CurTesSiz),NumFra*CurTesSiz);
                                    SigRecep3 = ((SigRecep3)./equ);%Fixing phase rotation and amplitude variation with the equalizer
                                    clear SigRecep4;
                                    UniDmtMve = unique(DmtMve);
                                    UniDmtMve = fliplr(UniDmtMve);
                                    DaSiAn = 0;
                                    for CarDmtM=1:length(UniDmtMve)
                                        M = UniDmtMve(CarDmtM);
                                        DaSiAu = sum(ismember(DmtMve,M));
                                        SigRecep4(1+DaSiAn:DaSiAn + DaSiAu,:)  = qamdemod(SigRecep3(1+DaSiAn:DaSiAn + DaSiAu,:),M);
                                        DaSiAn = DaSiAn + DaSiAu;
                                    end
                                    %for CarDmtM=1:length(UniDmtMve)
                                        %M = UniDmtMve(CarDmtM);
                                        %DaSiAu = sum(ismember(DmtMve,M));
                                        %SigRecep4(1+DaSiAu*(CarDmtM-1):DaSiAu*CarDmtM,:)  = qamdemod(SigRecep3(1+DaSiAu*(CarDmtM-1):DaSiAu*CarDmtM,:),M);
                                    %end
                                    %                     SigRecep4  = qamdemod(SigRecep3,M);                    %Modulating information
                                otherwise
                                    %                     SigRecep3 = reshape(SigRecep3,NumFra,NFFT);            %Reshaping the signal if multiple carriers were transmited
                                    SigRecep4 = dpskdemod(SigRecep3,M);                    %Modulating information
                                    %                     SigRecep4 = SigRecep4(:);                              %Serializing the signal
                            end
                            SigRecep4 = SigRecep4.';
                            SigRecep4 = SigRecep4(:).';
                            RxSigOfdm(CurrentTest,ThisCarr,:) = SigRecep4;
                            if PlotingThis
                                figure;
                                txcolor = [0.2 0 1];
                                rxcolor = [1 0.4 0];
                                hold all;
                                plot(SigRecep3(:),'o','LineWidth',2,'color',rxcolor);
                                plot(TxSigToPlot(:),'*','LineWidth',2,'color',txcolor);
                                set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                            end
                            if PlotingThis
                                figure;
                                hold all;
                                plot(TxDataMat(ThisCarr,:));
                                plot(SigRecep4);
                                SampPos = 1:length(SigRecep4);
                                plot(SampPos(~(TxDataMat(ThisCarr,:)==SigRecep4)),...
                                    SigRecep4(~(TxDataMat(ThisCarr,:)==SigRecep4)),'ko','MarkerFaceColor',[0 0 0]);
                                set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                            end
                            SymbErr = ~(TxDataMat(ThisCarr,:)==SigRecep4);%Measuring system bit error ration
                            DmtMvep = repmat(DmtMve,NumFra*CurTesSiz,1);
                            DmtKvep = log2(DmtMvep(:).');
                            
                            BerToPlotOfdm(CurrentTest,ThisCarr,:) = SymbErr.*DmtKvep;
                            %             BerOFDM(CurrentTest,ThisCarr) = sum(SymbErr)/length(SymbErr);
                            BerOFDM(CurrentTest,ThisCarr) = sum(SymbErr.*DmtKvep)/sum(DmtKvep);
                            a=a+0;
                            close all;
                        end
                    case 'DPSK'
                        %%                   Receiver DPSK
                        for ThisCarr=InitCarrUp:CarrPass:NumCarr                                    %For each carrier the same process of reception will be used.
                            if CarrUsedUp(ThisCarr)
                                %%  Reception Process: Handling the income optical field
                                %At the first step the income field will be selected from the
                                %optical FFT Output with the help of the maping vector
                                %(previously described).
                                % 								EoutAux = EoutAux1(VetThisCarr==ObsCarrUsed(ThisCarr),:);
                                Ix = EoutAux1(:,:,VetThisCarr==ObsCarrUsed(ThisCarr));
                                %%            Fiber Time Delay Compensation
                                switch Medium
                                    case 'Fiber'
                                        Ix = ifft(fft(Ix).*exp(1j*2*pi*f*(FiberDelay(ThisCarr)*Ta)));
                                    otherwise
                                end
                                %             EoutAux = EoutRec;
                                %%  Plot for Qualitative analizes
                                %PrintInfo(Ploting*13,t,abs(EoutAux));
                                %             PrintInfo(Ploting*14,f,20*log10(abs(fftshift(fft(EoutAux)./...
                                %                                                        length(EoutAux)))));
                                
                                %%                   synchronizing
                                %For the reception process work properly, it is needed to
                                %sincronized the recieved signal or the sampling process will
                                %not work.
                                
                                %At this first moment a small part of the income signal needs
                                %be analized for the sincronism process. As we are looking for
                                %the actual synchronis symbol, which is the bigest phase shift
                                %of the income signal. The next delay interferometer convert
                                %this phase shift to a amplitude variation.
                                
                                %The phase delay is not needed because the optical field will
                                %be analized as a whole as it is composed only by the real part
                                PhaDel          = 0;
                                %Remember that at DPSK modulation the information is stored at
                                %the phase difference of the current symbol with the previous
                                %one hence this time delay is needed to analyse each symbol by
                                %analizes of its interaction of the current symbel with the
                                %previous one.
                                TimDel          = T;
                                if UsingGpu==1
                                    TimeDelGpu = gpuArray(TimDel);
                                    PhaDelGpu  = gpuArray(PhaDel);
                                    ESync1 = zeros(size(Ix,1),size(Ix,2));
                                    ESync2 = zeros(size(Ix,1),size(Ix,2));
                                    for jj=1:size(Ix,2)
                                        IxGpu = gpuArray(Ix(:,jj));
                                        [ESync1Gpu,ESync2Gpu] = DelayInterfG(fgpu,TimeDelGpu,PhaDelGpu,IxGpu);
                                        ESync1(:,jj) = gather(ESync1Gpu);
                                        ESync2(:,jj) = gather(ESync2Gpu);
                                        clear IxGpu ESync1Gpu ESync2Gpu;
                                    end
                                    clear TimeDelGpu PhaDelGpu;
                                else
                                    [ESync1,ESync2] = DelayInterfExp(f,TimDel,PhaDel,Ix);        %Coverting phase shift to amplitude variation
                                end
                                %[ESync1,ESync2] = DelayInterfExp(t,TimDel,PhaDel,Ix);        %Coverting phase shift to amplitude variation
                                
                                %The second moment is to transfer this information from the
                                %optical domain to the eletrical domain for an eletronic
                                %processing. It is done with the help of an photo diode. The
                                %configuration here used is an balanced receiver as the output
                                %of the Delay Interferometer has two signals resulting from
                                %the constructive and destructive signal interaction.
                                ESync1 = ESync1.*conj(ESync1);
                                
                                ESync2 = ESync2.*conj(ESync2);
                                
                                %%           Creating the Reception Filter
                                
                                
                                [BitFilt,~] = FiltroGaussiano(f,BWD2,CenFeq2,FiltOrd2);        %Creating filter for selection of the received signal
                                BitFilt = fftshift(BitFilt);                                   %Shifting the filter for matching the received signal
                                
                                %%
                                Esync  = ESync2-ESync1;
                                
                                if RecFilBanPas==1
                                    Esync  = ifft(fft(Esync).*BitFilt);%Filter is used to remove higher order components
                                    Emean  = mean(Esync);
                                    Emean  = repmat(Emean,size(Esync,1),1);
                                    Esync  = Esync-Emean;
                                    Enorm  = max(Esync);
                                    Enorm  = repmat(Enorm,size(Esync,1),1);
                                    Esync  = Esync./Enorm;
                                end
                                %%         Adding Noise to Received Signal
                                %Just after the photo diod is where the noise must be added as
                                %it is the main constrain of the system. Once the signal is
                                %converted from optical to electrical and it digital values
                                %were recevered, there is no point to add noise, because the
                                %information was already received, it just needs decoding.
                                %The noise here added will be gaussian basically it represents
                                %the reception sensibility. In another words, how much the
                                %signal must be about the noise for an acceptable recovering.
                                if ReceptorNoise==1                                               %Verify whether the noise is to be added or not
                                    [~,SigPower] = MeasPower(Esync);
                                    SigPower2 = SigPower-30;%20*log10(SigPowerI);
                                    SNR2 = CarSNR + 10*log10(1) - 10*log10(Nsamp) + 10*log10(10^0.36);
                                    if TimeSys==1
                                        SNR = CarSNR + 10*log10(1) - 10*log10(1*T*BWD2);
                                    else
                                        SNR = CarSNR + 10*log10(1) - 10*log10(1*t(2*NumAmosCP+NPPB)*BWD2);
                                    end
                                    Esync = awgn(Esync,SNR,SigPower2);
                                    %                     Esync = ifft(fft(Esync).*BitFiltN);                             %Filter is used to remove higher order components
                                end
                                if ~RecFilBanPas
                                    Esync  = ifft(fft(Esync).*BitFilt);%Filter is used to remove higher order components
                                    Emean  = mean(Esync);
                                    Emean  = repmat(Emean,size(Esync,1),1);
                                    Esync  = Esync-Emean;
                                    Enorm  = max(Esync);
                                    Enorm  = repmat(Enorm,size(Esync,1),1);
                                    Esync  = Esync./Enorm;
                                end
                                
                                
%                                 EoutAuxF = 20*log10(abs(fftshift(fft(Ix)./length(Ix))));
%                                 for jj=1:size(Ix,2)
%                                     VetElecPower(CurrentTest,ThisCarr,jj)= MeasPower(Ix(:,jj));
%                                     [VetOptiPower(CurrentTest,ThisCarr,jj),~]= findpeaks(EoutAuxF(:,jj),'SortStr','descend','NPeaks',1);
%                                 end
								if UsingMex==1
									[BerDPSK(CurrentTest,ThisCarr),ValsLev(CurrentTest,ThisCarr),AberLev(CurrentTest,ThisCarr),Esync,InterAB,LocMaxAB,MaxValAB,SeqFinAB,TxData,Data] = RecDowDepsk_mex(MaxNumStag,IfftOrSum,T,Ta,Esync,AddCP,StuffSampels,NumAmosCP,NPPB,NbDPSK,CurTesSiz,SyncPeriod,ThisCarr,TxDataMat,LevDefDpqsk);
								else
									[BerDPSK(CurrentTest,ThisCarr),ValsLev(CurrentTest,ThisCarr),AberLev(CurrentTest,ThisCarr),Esync,InterAB,LocMaxAB,MaxValAB,SeqFinAB,TxData,Data] = RecDowDepsk(MaxNumStag,IfftOrSum,T,Ta,Esync,AddCP,StuffSampels,NumAmosCP,NPPB,NbDPSK,CurTesSiz,SyncPeriod,ThisCarr,TxDataMat,LevDefDpqsk);
								end
                                if CalcS==1
									if UsingMex==1
										[~,~,EyeOpenI,EyeOpenHighI,EyeOpenLowI] = Olho_mex(Esync,T,NPPB,0,1);
                                    else
										[~,~,EyeOpenI,EyeOpenHighI,EyeOpenLowI] = Olho(Esync,T,NPPB,0,1);
									end
									AberLevS  = EyeOpenI;
                                    ValsLevS  = EyeOpenLowI + EyeOpenI/2;          %Comparison between the Transmited and received and counting the differences
                                    ThisDataSize = NPPB/2:NPPB:size(Esync,2);
                                    DataS = zeros(1,length(ThisDataSize));                                                  %Initialization of the vector that will store the income data
                                    for kk=1:length(ThisDataSize)%length(Esync(ThisDataSize))                                    %The comparison process will be made for each symbol period
                                        %MeanOfData = mean(EoutI((kk-1)+SymLocI(1)));
                                        MeanOfData = mean(Esync((kk-1)*NPPB+NPPB/2));
                                        if MeanOfData > EyeOpenLowI + EyeOpenI/2%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                                            DataS(kk) = 1;                                 %is 1
                                        else                                                       %If it is the lowest level the incoming data
                                            DataS(kk) = 0;                                 %is 0
                                        end
                                    end
                                    %       Calculating the Bit Error Ratio (BER)
                                    %The final process here is to count the number of wrongdoings
                                    %of this whole process upon the transmited data for
                                    %quantitative analizes

                                    TxDataA = TxDataMat(ThisCarr,:);
                                    TxDataB = reshape(TxDataA,NbDPSK,CurTesSiz);
                                    TxDataC = TxDataB(1+SyncPeriod:end-SyncPeriod,:);
                                    TxData  = reshape(TxDataC,1,size(TxDataC,1)*size(TxDataC,2));
%                                     if CalcS
                                    %DataS = DataS(1+SyncPeriod:end-SyncPeriod);
                                    BitErrS = sum(xor(TxData,DataS));%Comparison between the Transmited and received and counting the differences
                                    BerDPSKS(CurrentTest,ThisCarr) = BitErrS/((NbDPSK-(2*SyncPeriod))*CurTesSiz);%Calculating the ration of wrong bits and the total number of bits transmited
%                                     end
                                end
                                %% Ploting the result for qualitative analizes
                                if PrintinEye==1
                                    PrintInfo(PlotingThis*17,Esync,T,NPPB);
                                    hold on;
                                    plot(t((NPPB)/2),InterAB(round(SeqFinAB(LocMaxAB)-MaxValAB)),'k^');
                                    plot(t((NPPB)/2),InterAB(SeqFinAB(LocMaxAB)),'kv');
                                    plot(t((NPPB)/2),InterAB(round(SeqFinAB(LocMaxAB)-MaxValAB/2)),'k*');
                                    if CalcS==1
                                        plot(t((NPPB)/2),EyeOpenLowI,'mx');
                                        plot(t((NPPB)/2),EyeOpenHighI,'mo');
                                        plot(t((NPPB)/2),EyeOpenLowI + EyeOpenI/2,'md');
                                    end
                                    set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                                    %figure;
                                    %hold all;
                                    %plotpos = zeros(1,length(IxAux1));
                                    %plot(IxAux1,plotpos,'o','color',[1 0.4 0]);
                                    %plot(EvmMatRef(ObsCarrUsed(ThisCarr),:),plotpos,'b*');
                                    %set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                                    %a=a+0;
                                end
                                if PlotingThis
                                    figure;
                                    hold all;
                                    plot(TxData);
                                    plot(Data);
                                    set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                                    display(BerDPSK);
                                end
                                close all;
                            end
                        end
                    case 'DQPSK'
                        
                        %%                   Receiver DQPSK
                        for ThisCarr=InitCarrUp:CarrPass:NumCarr                                           %For each carrier the same process of reception will be used.
                            if CarrUsedUp(ThisCarr)
                                %%  Reception Process: Handling the income optical field
                                %At the first step the income field will be selected from the
                                %optical FFT Output with the help of the maping vector
                                %(previously described).
                                Ix = EoutAux1(:,:,VetThisCarr==ObsCarrUsed(ThisCarr));
                                if AddingNoiseP==1
                                    [~,sigpower] = MeasPower(Ix);
                                    Ix      = awgn(Ix,osnrp,sigpower-30);
                                    [~,sigpower] = MeasPower(Ix);
                                end
                                %%            Fiber Time Delay Compensation
                                switch Medium
                                    case 'Fiber'
                                        Ix = ifft(fft(Ix).*exp(1j*2*pi*f*(...
                                            FiberDelay(ThisCarr)*Ta)));
                                    otherwise
                                end
                                %%               Recovering the Information
                                %Once the Income signal was synchronized it is possible to
                                %split its components in phase and in quadrature. It is very
                                %important to understand here that those components will carry
                                %the actual data transmited. From it we can't recover the  I
                                %and Q eletrical signals that were used to modulate our
                                %carrier. Those signals get mingled together to composed the
                                %actual data. Therefore, when we pass the income signal
                                %throughout the receptor we will recover the TxData generated
                                %and transmited. It is important that it get very clear here
                                %because I toke two months to understand it, no one could
                                %explain it well to me and those articles about it make this
                                %information unclear. Thus, to sheed a light on this concept,
                                %it is important to know that the Data is encoded in the I and
                                %Q eletrical components that will get mixed together to
                                %modulate the optical carrier. When this optical signal passes
                                %throughout the MZ-Interferometer the acutal data (TX) will be
                                %recover, not the I and Q component. But the data at the odd
                                %possition (encoded in I) will be at the real component (in
                                %phase), whereas the data at the even possition (encoded in Q)
                                %will be at the imaginary component(in quadrature). The
                                %MZ-Interferometer will be resposible to take apart the real
                                %and imaginary components of the income optical field.
                                
                                %                                 E_rec3 = Ix./max(abs(Ix));                           %Normalizing income signal
                                %For the interferometric process  take in account just the real
                                %component it is needed a phase delay of 45° degrees;
                                PhaDel = 1*pi/4;
                                %Remember that at DQPSK modulation the information is stored at
                                %the difference of phase between the current and previous
                                %symbol hence this time delay of one symbol period is needed.
                                TimDel = T;
                                if UsingGpu==1
                                    TimeDelGpu = gpuArray(TimDel);
                                    PhaDelGpu  = gpuArray(PhaDel);
                                    EoutA = zeros(size(Ix,1),size(Ix,2));
                                    EoutB = zeros(size(Ix,1),size(Ix,2));
                                    for jj=1:size(Ix,2)
                                        IxGpu = gpuArray(Ix(:,jj));
                                        [ESync1Gpu,ESync2Gpu] = DelayInterfG(fgpu,TimeDelGpu,PhaDelGpu,IxGpu);
                                        EoutA(:,jj) = gather(ESync1Gpu);
                                        EoutB(:,jj) = gather(ESync2Gpu);
                                        clear IxGpu ESync1Gpu ESync2Gpu;
                                    end
                                    clear TimeDelGpu PhaDelGpu;
                                else
                                    [EoutA,EoutB] = DelayInterfExp(f,TimDel,PhaDel,Ix);        %Coverting phase shift to amplitude variation
                                end
                                %[EoutA,EoutB] = DelayInterfExp(t,TimDel,PhaDel,Ix);          %Coverting phase shift to amplitude variation
                                %For the interferometric process  take in account just the real
                                %component it is needed a phase delay of -45° degrees;
                                PhaDel = -1*pi/4;
                                %Remember that at DQPSK modulation the information is stored at
                                %the difference of phase between the current and previous
                                %symbol hence this time delay of one symbol period is needed.
                                TimDel = T;
                                if UsingGpu==1
                                    TimeDelGpu = gpuArray(TimDel);
                                    PhaDelGpu  = gpuArray(PhaDel);
                                    EoutC = zeros(size(Ix,1),size(Ix,2));
                                    EoutD = zeros(size(Ix,1),size(Ix,2));
                                    for jj=1:size(Ix,2)
                                        IxGpu = gpuArray(Ix(:,jj));
                                        [ESync1Gpu,ESync2Gpu] = DelayInterfG(fgpu,TimeDelGpu,PhaDelGpu,IxGpu);
                                        EoutC(:,jj) = gather(ESync1Gpu);
                                        EoutD(:,jj) = gather(ESync2Gpu);
                                        clear IxGpu ESync1Gpu ESync2Gpu;
                                    end
                                    clear TimeDelGpu PhaDelGpu;
                                else
                                    [EoutC,EoutD] = DelayInterfExp(f,TimDel,PhaDel,Ix);        %Coverting phase shift to amplitude variation
                                end
                                %[EoutC,EoutD] = DelayInterfExp(t,TimDel,PhaDel,Ix);          %Coverting phase shift to amplitude variation
                                
                                %The second moment is to transfer this information from the
                                %optical domain to the eletrical domain for an eletronic
                                %processing. It is done with the help of an photo diode.
                                EoutA = EoutA.*conj(EoutA);
                                EoutB = EoutB.*conj(EoutB);
                                EoutC = EoutC.*conj(EoutC);
                                EoutD = EoutD.*conj(EoutD);
                                
                                %The process with the photo diode is self-coherent, which means
                                %the rusult will be a component of the signal centered in f=0
                                %and another component centered at f=2*fc (frenquecy central).
                                %Therefore, to remove the higher order component a low pass
                                %filter will be used.
                                %%           Creating the Reception Filter
                                [BitFilt,~] = FiltroGaussiano(f,BWD2,CenFeq2,FiltOrd2);        %Creating filter for selection of the received signal
                                BitFilt = fftshift(BitFilt);                                   %Shifting the filter for matching the received signal
                                %%
                                
                                if RecFilBanPas
                                    EoutA = ifft(fft(EoutA).*BitFilt);
                                    EoutB = ifft(fft(EoutB).*BitFilt);
                                    EoutC = ifft(fft(EoutC).*BitFilt);
                                    EoutD = ifft(fft(EoutD).*BitFilt);
                                end
                                
                                % The configuration here used is an balanced receiver as the
                                %output of the Delay Interferometer has two signals resulting
                                %from the constructive and destructive signal interaction.
                                %%         Adding Noise to Received Signal
                                %Just after the photo diod is where the noise must be added as
                                %it is the main constrain of the system. Once the signal is
                                %converted from optical to electrical and it digital values
                                %were recevered, there is no point to add noise, because the
                                %information was already received, it just needs decoding.
                                %The noise here added will be gaussian basically it represents
                                %the reception sensibility. In another words, how much the
                                %signal must be about the noise for an acceptable recovering.
                                
                                if ReceptorNoise==1                                               %Verify whether the noise is to be added or not
                                    [~,SigPowerA] = MeasPower(EoutA);                          %it own received signal or
                                    [~,SigPowerB] = MeasPower(EoutB);
                                    [~,SigPowerC] = MeasPower(EoutC);                          %it own received signal or
                                    [~,SigPowerD] = MeasPower(EoutD);
                                    SNR2 = CarSNR + 10*log10(2) - 10*log10(Nsamp) + 10*log10(10^0.3);
                                    
                                    if TimeSys==1
                                        SNR = CarSNR + 10*log10(2) - 10*log10(T*BWD2);
                                    else
                                        SNR = CarSNR + 10*log10(2) - 10*log10(1*t(2*NumAmosCP+NPPB)*BWD2);
                                    end
                                    SigPowerA2 = SigPowerA-30;%20*log10(SigPowerI);
                                    SigPowerB2 = SigPowerB-30;%20*log10(SigPowerQ);
                                    SigPowerC2 = SigPowerC-30;%20*log10(SigPowerI);
                                    SigPowerD2 = SigPowerD-30;%20*log10(SigPowerQ);
                                    EoutA = awgn(EoutA,SNR,SigPowerA2);
                                    EoutB = awgn(EoutB,SNR,SigPowerB2);
                                    EoutC = awgn(EoutC,SNR,SigPowerC2);
                                    EoutD = awgn(EoutD,SNR,SigPowerD2);
                                end
                                
                                
                                if ~RecFilBanPas
                                    EoutA = ifft(fft(EoutA).*BitFilt);
                                    EoutB = ifft(fft(EoutB).*BitFilt);
                                    EoutC = ifft(fft(EoutC).*BitFilt);
                                    EoutD = ifft(fft(EoutD).*BitFilt);
                                end
                                if UsingMex==1
									[AberLevI(CurrentTest,ThisCarr),ValsLevI(CurrentTest,ThisCarr),AberLevQ(CurrentTest,ThisCarr),ValsLevQ(CurrentTest,ThisCarr),BerDQPSK(CurrentTest,ThisCarr),EoutI,EoutQ,TxDataOdd,TxDataEven,DataOdd,DataEven,LocMaxAB,MaxValAB,SeqFinAB,QMaxValAB,QLocMaxAB,QSeqFinAB,QInterAB,InterAB,IxAuxAB,QIxAuxAB]=RecDowDqpsk_mex(EoutA,EoutB,EoutC,EoutD,T,Ta,IfftOrSum,MaxNumStag,StuffSampels,NbDQPSK,CurTesSiz,NumAmosCP,NPPB,SyncPeriod,TxDataMat,ThisCarr,LevDefDpqsk,AddCP,NumCarr);
                                else
									[AberLevI(CurrentTest,ThisCarr),ValsLevI(CurrentTest,ThisCarr),AberLevQ(CurrentTest,ThisCarr),ValsLevQ(CurrentTest,ThisCarr),BerDQPSK(CurrentTest,ThisCarr),EoutI,EoutQ,TxDataOdd,TxDataEven,DataOdd,DataEven,LocMaxAB,MaxValAB,SeqFinAB,QMaxValAB,QLocMaxAB,QSeqFinAB,QInterAB,InterAB,IxAuxAB,QIxAuxAB]=RecDowDqpsk(EoutA,EoutB,EoutC,EoutD,T,Ta,IfftOrSum,MaxNumStag,StuffSampels,NbDQPSK,CurTesSiz,NumAmosCP,NPPB,SyncPeriod,TxDataMat,ThisCarr,LevDefDpqsk,AddCP,NumCarr);
                                end
								%%
                                %             berpos = 1:2:size(BerDQPSK,2);
                                %             BerDQPSK(size(BerDQPSK,1),berpos)
                                if CalcS==1
                                    if UsingMex==1
                                        [~,~,EyeOpenI,EyeOpenHighI,EyeOpenLowI] = Olho_mex(EoutI,T,NPPB,0,1);
                                        [~,~,EyeOpenQ,EyeOpenHighQ,EyeOpenLowQ] = Olho_mex(EoutQ,T,NPPB,0,1);
                                    else
                                        [~,~,EyeOpenI,EyeOpenHighI,EyeOpenLowI] = Olho(EoutI,T,NPPB,0,1);
                                        [~,~,EyeOpenQ,EyeOpenHighQ,EyeOpenLowQ] = Olho(EoutQ,T,NPPB,0,1);
                                    end
									AberLevIS(CurrentTest,ThisCarr)  = EyeOpenI;
                                    ValsLevIS(CurrentTest,ThisCarr)  = EyeOpenLowI + EyeOpenI/2;
                                    AberLevQS(CurrentTest,ThisCarr)  = EyeOpenQ;
                                    ValsLevQS(CurrentTest,ThisCarr)  = EyeOpenLowQ + EyeOpenQ/2;
                                    
                                    ThisDataSize = NPPB/2:NPPB:length(EoutI);
                                    DataOddS = zeros(1,length(ThisDataSize));                                                  %Initialization of the vector that will store the income data
                                    for kk=1:length(ThisDataSize)%length(EoutI(ThisDataSize))                                    %The comparison process will be made for each symbol period
                                        %MeanOfData = mean(EoutI((kk-1)+SymLocI(1)));
                                        MeanOfData = mean(EoutI((kk-1)*NPPB+NPPB/2));
                                        if MeanOfData > EyeOpenLowI + EyeOpenI/2%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                                            DataOddS(kk) = 1;                                 %is 1
                                        else                                                       %If it is the lowest level the incoming data
                                            DataOddS(kk) = 0;                                 %is 0
                                        end
                                        
                                    end
                                    ThisDataSize = NPPB/2:NPPB:length(EoutQ);
                                    DataEvenS = zeros(1,length(ThisDataSize));                                                 %Initialization of the vector that will store the income data
                                    for kk=1:length(ThisDataSize)%length(EoutQ(ThisDataSize))                                    %The comparison process will be made for each symbol period
                                        %                 MeanOfData = mean(EoutQ((kk-1)+SymLocI(1)));
                                        MeanOfData = mean(EoutQ((kk-1)*NPPB+NPPB/2));
                                        if MeanOfData > EyeOpenLowQ + EyeOpenQ/2%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                                            DataEvenS(kk) = 1;                                 %is 1
                                        else                                                       %If it is the lowest level the incoming data
                                            DataEvenS(kk) = 0;                                 %is 0
                                        end
                                        
                                    end
                                    
                                    BitErrOddS    = sum(xor(TxDataOdd,DataOddS));                      %Comparison between the Transmited and received and counting the differences
                                    BitErrEvenS   = sum(xor(TxDataEven,DataEvenS));                    %Comparison between the Transmited and received and counting the differences
                                    
                                    BerDQPSKS(CurrentTest,ThisCarr) = (BitErrOddS+BitErrEvenS)/...
                                        ((NbDQPSK)-(2*SyncPeriod));%Calculating the ration of wrong bits and the total number of bits transmited
                                end
%                                 if (ThisCarr==2)&&(Ploting==1)
%                                     IxAux1  = [EoutI(1:end - StuffSampels,:);EoutQ(1:end - StuffSampels,:)];
%                                     Olho(IxAux1(:),Tcp,(2*NumAmosCP+NPPB),1,4);
%                                 end
                                if PrintinEye==1
                                    PrintInfo(PlotingThis*24,EoutI,T,NPPB);
                                    hold on;
                                    plot(t((NPPB)/2),InterAB(round(SeqFinAB(LocMaxAB)-MaxValAB)),'k^');
                                    plot(t((NPPB)/2),InterAB(SeqFinAB(LocMaxAB)),'kv');
                                    plot(t((NPPB)/2),InterAB(round(SeqFinAB(LocMaxAB)-MaxValAB/2)),'kd');
                                    %PrintInfo(Ploting*25,Ui,T,NPPB);
                                    PrintInfo(Ploting*26,EoutQ,T,NPPB);
                                    hold on;
                                    plot(t((NPPB)/2),QInterAB(round(QSeqFinAB(QLocMaxAB)-QMaxValAB)),'k^');
                                    plot(t((NPPB)/2),QInterAB(QSeqFinAB(QLocMaxAB)),'kv');
                                    plot(t((NPPB)/2),QInterAB(round(QSeqFinAB(QLocMaxAB)-QMaxValAB/2)),'kd');
                                    if CalcS
                                        plot(t((NPPB)/2),EyeOpenLowI,'mx');
                                        plot(t((NPPB)/2),EyeOpenHighI,'mo');
                                        plot(t((NPPB)/2),EyeOpenLowI + EyeOpenI/2,'md');
                                        plot(t((NPPB)/2),EyeOpenLowQ,'cx');
                                        plot(t((NPPB)/2),EyeOpenHighQ,'co');
                                        plot(t((NPPB)/2),EyeOpenLowQ + EyeOpenQ/2,'cd');
                                    end
                                    set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                                    figure;
                                    hold all;
                                    plotpos = zeros(1,length(TxDataOdd));
                                    plot(IxAuxAB,QIxAuxAB,'o','color',[1 0.4 0]);
                                    TxDataOddP = TxDataOdd;
                                    TxDataOddP(TxDataOddP==0) = -1;
                                    TxDataEvenP = TxDataEven;
                                    TxDataEvenP(TxDataEvenP==0) = -1;
                                    plot(TxDataOddP,TxDataEvenP,'b*');
                                    set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                                    %a=a+0;
                                end
                                %             a=a+6;
                                if PlotingThis
                                    figure;
                                    hold all;
                                    plot([TxDataOdd TxDataEven]);
                                    plot([DataOdd DataEven]);
                                    set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                                    display(BerDQPSK);
                                end
                                close all;
                            end
                        end
                    case '4PAM'
                        
                        %%              Receiver 4PAM
                        for ThisCarr=InitCarrUp:CarrPass:NumCarr                                           %For each carrier the same process of reception will be used.
                            if CarrUsedUp(ThisCarr)
                                %%  Reception Process: Handling the income optical field
                                %At the first step the income field will be selected from the
                                %optical FFT Output with the help of the maping vector
                                %(previously described).
                                Ix = EoutAux1(:,:,VetThisCarr==ObsCarrUsed(ThisCarr));
                                %%            Fiber Time Delay Compensation
                                switch Medium
                                    case 'Fiber'
                                        Ix = ifft(fft(Ix).*exp(1j*2*pi*f*(FiberDelay((ThisCarr))*Ta)));
                                    otherwise
                                end
                                %The current incoming signal is them converted from the optical
                                %domain to the eletrical domain with the help of an photo
                                %detector.
                                Ix = Ix.*conj(Ix);
                                %The process with the photo diode is self-coherent, which means
                                %the rusult will be a component of the signal centered in f=0
                                %and another component centered at f=2*fc (frenquecy central).
                                %Therefore, to remove the higher order component a low pass
                                %filter will be used.
                                %             switch Medium
                                %                 case 'Fiber'
                                %%           Creating the Reception Filter
                                %             taux = t(1:length(Ix));
                                %             faux = time2freq(taux);
                                [BitFilt,~] = FiltroGaussiano(f,BWD2,CenFeq2,FiltOrd2);        %Creating filter for selection of the received signal
                                BitFilt = fftshift(BitFilt);                                   %Shifting the filter for matching the received signal
                                if RecFilBanPas==1
                                    Ix = ifft(fft(Ix).*BitFilt);
                                    Ix(1:4*(2*NumAmosCP+NPPB),:) = repmat(Ix(5*(2*NumAmosCP+NPPB),:),size(Ix(1:4*(2*NumAmosCP+NPPB),:),1),1);
                                    Ix(end+1-4*(2*NumAmosCP+NPPB):end,:) = repmat(Ix(5*(2*NumAmosCP+NPPB),:),size(Ix(1:4*(2*NumAmosCP+NPPB),:),1),1);
                                    IxMin = min(Ix);
                                    IxMin = repmat(IxMin,size(Ix,1),1);
                                    Ix = Ix - IxMin;                                         %Removing the DC component from them Eletrical signal received
                                    IxMax = max(Ix);
                                    IxMax = repmat(IxMax,size(Ix,1),1);
                                    Ix = Ix./IxMax;                                     %Normalizing the eletrical signal (amplifying what is needed)
                                end
                                %                 otherwise
                                %             end
                                
                                
                                %%         Adding Noise to Received Signal
                                %Just after the photo diod is where the noise must be added as
                                %it is the main constrain of the system. Once the signal is
                                %converted from optical to electrical and it digital values
                                %were recevered, there is no point to add noise, because the
                                %information was already received, it just needs decoding.
                                %The noise here added will be gaussian basically it represents
                                %the reception sensibility. In another words, how much the
                                %signal must be about the noise for an acceptable recovering.
                                if ReceptorNoise==1                                               %Verify whether the noise is to be added or not
                                    [~,SigPower] = MeasPower(Ix);                              %it own received signal or
                                    SigPower2 = SigPower-30;%10*log10(SigPower);
                                    Ix = awgn(Ix,SNR,SigPower2);
                                end
                                
                                if ~RecFilBanPas
                                    Ix = ifft(fft(Ix).*BitFilt);
                                    Ix(1:4*(2*NumAmosCP+NPPB),:) = repmat(Ix(5*(2*NumAmosCP+NPPB),:),size(Ix(1:4*(2*NumAmosCP+NPPB),:),1),1);
                                    Ix(end+1-4*(2*NumAmosCP+NPPB):end,:) = repmat(Ix(5*(2*NumAmosCP+NPPB),:),size(Ix(1:4*(2*NumAmosCP+NPPB),:),1),1);
                                    IxMin = min(Ix);
                                    IxMin = repmat(IxMin,size(Ix,1),1);
                                    Ix = Ix - IxMin;                                         %Removing the DC component from them Eletrical signal received
                                    IxMax = max(Ix);
                                    IxMax = repmat(IxMax,size(Ix,1),1);
                                    Ix = Ix./IxMax;                                    %Normalizing the eletrical signal (amplifying what is needed)
                                end
                                
                                %                 [BitFilt,~] = FiltroGaussiano(f,200*BWD2,CenFeq2,FiltOrd2);        %Creating filter for selection of the received signal
                                %                 BitFilt = fftshift(BitFilt);                                   %Shifting the filter for matching the received signal
                                %                 Ix = ifft(fft(Ix).*BitFilt);
                                
                                %EoutAuxF = 20*log10(abs(fftshift(fft(Ix)./length(Ix))));
                                %for jj=1:size(Ix,2)
                                    %VetElecPower(CurrentTest,ThisCarr,jj)= MeasPower(Ix(:,jj));
                                    %[VetOptiPower(CurrentTest,ThisCarr,jj),~]= findpeaks(EoutAuxF(:,jj),'SortStr','descend','NPeaks',1);
                                %end
                                %For the reception process work properly, it is needed to
                                %sincronized the recieved signal or the sampling process will
                                %not work.
                                %AuxSync = (Ix(IniSyncPos:SyncPos));                            %Selecting the sync-word within the received signal
                                %AuxSync1 = AuxSync;
                                %This synchronization process is based on the mean expected
                                %value. That means, the information mean value should be within
                                %one period of symbol. Thus, the mean value of the received
                                %signal is acquired and compare of the known sync-word to
                                %verify if this mean value is at the right possition
                                %AuxSync1(AuxSync1<0)               = 0;                        %To keep the mean value above zero anything under is neglected
                                %AuxSync1(AuxSync1>=mean(AuxSync1)) = 1;                        %Adding a flag to the first sample of the received mean value
                                %AuxSync1(AuxSync1<mean(AuxSync1))  = -1;                       %All the others samples at set to the lowest level
                                %AuxSync2                           = SyncSymb(IniSyncPos:...
                                %SyncPos);%Selecting the sync-word within the known signal
                                %AuxSync2(AuxSync2>=mean(AuxSync2)) = 1;                        %Adding a flag to the first sample of the known mean value
                                %AuxSync2(AuxSync2<mean(AuxSync2))  = -1;                       %All the others samples at set to the lowest level
                                
                                %PosToSyn  = find(ismember(AuxSync1,1));                        %Finding where is the location of the first sample to synchronize
                                %PosSyn = find(ismember(AuxSync2,1));                           %Finding where is the location of the first sample to synchronize
                                
                                %AuxSyncCorr = PosToSyn(round(end/2))-PosSyn(round(end/2));
                                if UsingMex==1
									[LevDec1,LevDec2,LevDec3,Ix,Ber4PAM(CurrentTest,ThisCarr),DecLevDef3,DecLevDef2,DecLevDef1,TxData,IxRecDef,IxRec,AberLev(1,CurrentTest,ThisCarr),AberLev(2,CurrentTest,ThisCarr),AberLev(3,CurrentTest,ThisCarr),ValsLev(1,CurrentTest,ThisCarr),ValsLev(2,CurrentTest,ThisCarr),ValsLev(3,CurrentTest,ThisCarr),ValsLev2(1,CurrentTest,ThisCarr),ValsLev2(2,CurrentTest,ThisCarr),ValsLev2(3,CurrentTest,ThisCarr),InterAB,InterCD,InterEF,SeqFinAB,SeqFinCD,SeqFinEF,LocMaxAB,LocMaxCD,LocMaxEF,MaxValAB,MaxValCD,MaxValEF,Levels] = RecDowPam4_mex(Ix,T,Ta,MaxNumStag,StuffSampels,NumAmosCP,NPPB,CurTesSiz,Nb4Pam,IntervalStep,MinDist,DecLevDef1,DecLevDef2,DecLevDef3,TxDataMat,ThisCarr,IfftOrSum,AddCP,SyncPeriod,DecMod);
                                else
									[LevDec1,LevDec2,LevDec3,Ix,Ber4PAM(CurrentTest,ThisCarr),DecLevDef3,DecLevDef2,DecLevDef1,TxData,IxRecDef,IxRec,AberLev(1,CurrentTest,ThisCarr),AberLev(2,CurrentTest,ThisCarr),AberLev(3,CurrentTest,ThisCarr),ValsLev(1,CurrentTest,ThisCarr),ValsLev(2,CurrentTest,ThisCarr),ValsLev(3,CurrentTest,ThisCarr),ValsLev2(1,CurrentTest,ThisCarr),ValsLev2(2,CurrentTest,ThisCarr),ValsLev2(3,CurrentTest,ThisCarr),InterAB,InterCD,InterEF,SeqFinAB,SeqFinCD,SeqFinEF,LocMaxAB,LocMaxCD,LocMaxEF,MaxValAB,MaxValCD,MaxValEF,Levels] = RecDowPam4(Ix,T,Ta,MaxNumStag,StuffSampels,NumAmosCP,NPPB,CurTesSiz,Nb4Pam,IntervalStep,MinDist,DecLevDef1,DecLevDef2,DecLevDef3,TxDataMat,ThisCarr,IfftOrSum,AddCP,SyncPeriod,DecMod);
                                end
								%if BitErrS<BitErr
                                    %Ber4PAM(CurrentTest,ThisCarr) = BitErrS/((Nb4Pam-(4*SyncPeriod))*CurTesSiz);%Calculating the ration of wrong bits and the total number of bits transmited
                                %else
                                    %Ber4PAM(CurrentTest,ThisCarr) = BitErr/((Nb4Pam-(4*SyncPeriod))*CurTesSiz);%Calculating the ration of wrong bits and the total number of bits transmited
                                %end
                                %%           Ploting for Qualitative Analizes
                                %PrintInfo(Ploting*41,t(length(TxDataMat(ThisCarr,:))),Nb4Pam/2,TxDataMat(ThisCarr,:),IxRec);
                                %%       Calculating the Bit Error Ratio (BER)
                                %The final process here is to count the number of wrongdoings
                                %of this whole process upon the transmited data for
                                %quantitative analizes
                                
                                if CalcS==1
                                    limiar3 = (4/6)*abs(max(Ix));%(1/2)*abs(max(Ix)-min(Ix))/3;
                                    limiar2 = (2.5/6)*abs(max(Ix));%((1/2)*abs(max(Ix)-min(Ix))/3) + abs(max(Ix)-min(Ix))/3;
                                    limiar1 = (1/6)*abs(max(Ix));%((1/2)*abs(max(Ix)-min(Ix))/3) + 2*(abs(max(Ix)-min(Ix))/3);
                                    limiarPos1 = Levels(1);%abs(min(Ix));%Interval>limiar1;
                                    limiarPos2 = Levels(2);%0.33*(abs(max(Ix))-abs(min(Ix)));%(Interval<=limiar1)&(Interval>limiar2);
                                    limiarPos3 = Levels(3);%0.63*(abs(max(Ix))-abs(min(Ix)));%(Interval<=limiar2)&(Interval>limiar3);
                                    limiarPos4 = Levels(4);%abs(max(Ix));%Interval<=limiar3;
                                    
                                    EyeSymMa1 = reshape(Ix,NPPB,((Nb4Pam/2)-SyncSymbSiz)*CurTesSiz);
                                    Hi1 = zeros(1,size(EyeSymMa1,1));
                                    Hi2 = zeros(1,size(EyeSymMa1,1));
                                    Hi3 = zeros(1,size(EyeSymMa1,1));
                                    Lo1 = zeros(1,size(EyeSymMa1,1));
                                    Lo2 = zeros(1,size(EyeSymMa1,1));
                                    Lo3 = zeros(1,size(EyeSymMa1,1));
                                    LevHi1 = zeros(1,size(EyeSymMa1,1));
                                    LevHi2 = zeros(1,size(EyeSymMa1,1));
                                    LevHi3 = zeros(1,size(EyeSymMa1,1));
                                    LevLo1 = zeros(1,size(EyeSymMa1,1));
                                    LevLo2 = zeros(1,size(EyeSymMa1,1));
                                    LevLo3 = zeros(1,size(EyeSymMa1,1));
                                    EyeAb1 = zeros(1,size(EyeSymMa1,1));
                                    EyeAb2 = zeros(1,size(EyeSymMa1,1));
                                    EyeAb3 = zeros(1,size(EyeSymMa1,1));
                                    for kk = 1:size(EyeSymMa1,1)
                                        EyeHi3    = find((EyeSymMa1(kk,:)<limiarPos4)&(EyeSymMa1(kk,:)>limiar3));
                                        EyeHi2    = find((EyeSymMa1(kk,:)<limiarPos3)&(EyeSymMa1(kk,:)>limiar2));
                                        EyeHi1    = find((EyeSymMa1(kk,:)<limiarPos2)&(EyeSymMa1(kk,:)>limiar1));
                                        
                                        EyeLo3    = find((EyeSymMa1(kk,:)<limiar3)&(EyeSymMa1(kk,:)>limiarPos3));
                                        EyeLo2    = find((EyeSymMa1(kk,:)<limiar2)&(EyeSymMa1(kk,:)>limiarPos2));
                                        EyeLo1    = find((EyeSymMa1(kk,:)<limiar1)&(EyeSymMa1(kk,:)>limiarPos1));
                                        
                                        Hi1(kk)   = mean((EyeSymMa1(kk,EyeHi1)));
                                        Hi2(kk)   = mean((EyeSymMa1(kk,EyeHi2)));
                                        Hi3(kk)   = mean((EyeSymMa1(kk,EyeHi3)));
                                        
                                        Lo1(kk)   = mean((EyeSymMa1(kk,EyeLo1)));
                                        Lo2(kk)   = mean((EyeSymMa1(kk,EyeLo2)));
                                        Lo3(kk)   = mean((EyeSymMa1(kk,EyeLo3)));
                                        
                                        LevHi1(kk)= mean(Hi1) - std(EyeSymMa1(kk,EyeHi1));
                                        LevHi2(kk)= mean(Hi2) - std(EyeSymMa1(kk,EyeHi2));
                                        LevHi3(kk)= mean(Hi3) - std(EyeSymMa1(kk,EyeHi3));
                                        
                                        LevLo1(kk)= mean(Lo1) + std(EyeSymMa1(kk,EyeLo1));
                                        LevLo2(kk)= mean(Lo2) + std(EyeSymMa1(kk,EyeLo2));
                                        LevLo3(kk)= mean(Lo3) + std(EyeSymMa1(kk,EyeLo3));
                                        
                                        EyeAb1(kk)= LevHi1(kk) - LevLo1(kk);
                                        EyeAb2(kk)= LevHi2(kk) - LevLo2(kk);
                                        EyeAb3(kk)= LevHi3(kk) - LevLo3(kk);
                                    end
                                    EyeAbertura1 = mean(EyeAb1);
                                    EyeAbertura2 = mean(EyeAb2);
                                    EyeAbertura3 = mean(EyeAb3);
                                    ThisDeciLev1 = mean(LevLo1) + EyeAbertura1/2;
                                    ThisDeciLev2 = mean(LevLo2) + EyeAbertura2/2;
                                    ThisDeciLev3 = mean(LevLo3) + EyeAbertura3/2;
                                    AberLevS(1,ThisCarr,CurrentTest)  = EyeAbertura3;
                                    AberLevS(2,ThisCarr,CurrentTest)  = EyeAbertura2;
                                    AberLevS(3,ThisCarr,CurrentTest)  = EyeAbertura1;
                                    ValsLevS(1,ThisCarr,CurrentTest)  = ThisDeciLev3;
                                    ValsLevS(2,ThisCarr,CurrentTest)  = ThisDeciLev2;
                                    ValsLevS(3,ThisCarr,CurrentTest)  = ThisDeciLev1;

                                    %plot(t(NPPB/2),ThisDeciLev3,'b*');plot(t(NPPB/2),mean(LevHi3),'bv');plot(t(NPPB/2),mean(LevLo3),'b^');
                                    %plot(t(NPPB/2),ThisDeciLev2,'k*');plot(t(NPPB/2),mean(LevHi2),'kv');plot(t(NPPB/2),mean(LevLo2),'k^');
                                    %plot(t(NPPB/2),ThisDeciLev1,'g*');plot(t(NPPB/2),mean(LevHi1),'gv');plot(t(NPPB/2),mean(LevLo1),'g^');
                                    
                                    ThisDataSize = NPPB/2:NPPB:length(Ix);
                                    IxRecDeS = zeros(1,2*size(ThisDataSize,2));%Initialization of the vector that will store the income data
                                    ContBit2  = 1;
                                    for kk=1:length(ThisDataSize)%NPPB:length(Ix)                                       %The comparison process will be made for each symbol period
                                        %                 midaux = round(mean(SymLoc(1:round(end/2))));
                                        midaux = NPPB/2;%SymLoc(1);
                                        aux1 = Ix((kk-1)*NPPB+midaux+1);     %An small portion of the income signal is take for evaluation
                                        MeanRec = mean(aux1);                                      %Measuring the avarage value of the samples taken
                                        if MeanRec <= ThisDeciLev1                                   %If it is the lowest level the incoming data
                                            IxRecDeS(ContBit2) = 0;
                                            ContBit2  = 1 + ContBit2;
                                            IxRecDeS(ContBit2) = 0;
                                            ContBit2  = 1 + ContBit2;
                                            %IxRecDeS = [IxRecDeS 0 0];                             %is 01 (1)
                                        elseif (MeanRec <= ThisDeciLev2)&&(MeanRec > ThisDeciLev1)     %If it is the second level the incoming data
                                            IxRecDeS(ContBit2) = 0;
                                            ContBit2  = 1 + ContBit2;
                                            IxRecDeS(ContBit2) = 1;
                                            ContBit2  = 1 + ContBit2;
                                            %IxRecDeS = [IxRecDeS 0 1];                             %is 00 (0)
                                        elseif (MeanRec <= ThisDeciLev3)&&(MeanRec > ThisDeciLev2)     %If it is the tird level the incoming data
                                            IxRecDeS(ContBit2) = 1;
                                            ContBit2  = 1 + ContBit2;
                                            IxRecDeS(ContBit2) = 1;
                                            ContBit2  = 1 + ContBit2;
                                            %IxRecDeS = [IxRecDeS 1 1];                             %is 10 (2)
                                        elseif MeanRec > ThisDeciLev3                                %If it is the uper level the incoming data
                                            IxRecDeS(ContBit2) = 1;
                                            ContBit2  = 1 + ContBit2;
                                            IxRecDeS(ContBit2) = 0;
                                            ContBit2  = 1 + ContBit2;
                                            %IxRecDeS = [IxRecDeS 1 0];                             %is 11 (3)
                                        else                                                       %If for some misteriose reason neither of previous verification were sucedded
                                            IxRecDeS(ContBit2) = 0;
                                            ContBit2  = 1 + ContBit2;
                                            IxRecDeS(ContBit2) = 0;
                                            ContBit2  = 1 + ContBit2;
                                            %IxRecDeS = [IxRecDeS 0 0];                             %by default the current data is set to be 00 (0)
                                        end
                                    end
                                    
                                    BitErrS = sum(xor(TxData,IxRecDeS));%Comparison between the Transmited and received and counting the differences
                                    Ber4PAMS(CurrentTest,ThisCarr) = BitErrS/((Nb4Pam-(4*SyncPeriod))*CurTesSiz);%Calculating the ration of wrong bits and the total number of bits transmited
                               
                                end
                                
                                %%           Ploting for Qualitative Analizes
                                %if ThisCarr==126
                                %EyeToPlot(CurrentTest,1:length(Esync)) = Esync;
                                %save(['savingforplotingeye' num2str(VetSnr)],'EyeToPlot','VetSnr','IfftOrSum','UsedModula','T','NPPB');
                                %end
                                if PrintinEye==1
                                    PrintInfo(PlotingThis*36,Ix,T,NPPB);
                                    hold all;
                                    plot(t(NPPB/2),LevDec1,'b*');plot(t(NPPB/2),InterEF(SeqFinEF(LocMaxEF)-1),'bv');plot(t(NPPB/2),InterEF(round(SeqFinEF(LocMaxEF)-MaxValEF+1)),'b^');
                                    plot(t(NPPB/2),LevDec2,'gd');plot(t(NPPB/2),InterCD(SeqFinCD(LocMaxCD)-1),'gv');plot(t(NPPB/2),InterCD(round(SeqFinCD(LocMaxCD)-MaxValCD+1)),'g^');
                                    plot(t(NPPB/2),LevDec3,'kd');plot(t(NPPB/2),InterAB(SeqFinAB(LocMaxAB)-1),'kv');plot(t(NPPB/2),InterAB(round(SeqFinAB(LocMaxAB)-MaxValAB+1)),'k^');
                                  
                                    %figure;
                                    %hold all;
                                    %plotpos = zeros(1,length(IxAux));
                                    %plot(IxAux,plotpos,'o','color',[1 0.4 0]);
                                    %plot(EvmMatRef(ObsCarrUsed(ThisCarr),:),plotpos,'b*');
                                    set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                                    %a=a+1;
                                    %                     drawnow;
                                    %                     pause(0.1);
                                end
                                if PlotingThis
                                    figure;
                                    hold all;
                                    plot(TxData);
                                    plot(IxRec);
                                    set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                                    display(Ber4PAM);
                                end
                                close all;
                            else
                                Ber4PAM(CurrentTest,ThisCarr) = 1;
                            end
                        end
                    otherwise
                        %%              Receiver OOK
                        for ThisCarr=InitCarrUp:CarrPass:NumCarr                                           %For each carrier the same process of reception will be used.
                            if CarrUsedUp(ThisCarr)
                                %%  Reception Process: Handling the income optical field
                                %At the first step the income field will be selected from the
                                %optical FFT Output with the help of the maping vector
                                %(previously described).
                                Ix = EoutAux1(:,:,VetThisCarr==ObsCarrUsed(ThisCarr));
                                %%            Fiber Time Delay Compensation
                                switch Medium
                                    case 'Fiber'
                                        Ix = ifft(fft(Ix).*exp(1j*2*pi*f*(FiberDelay(ThisCarr)*Ta)));
                                    otherwise
                                end
                                
                                %The current incoming signal is them converted from the optical
                                %domain to the eletrical domain with the help of an photo
                                %detector.
                                Ix =Ix.*conj(Ix);
                                %%           Creating the Reception Filter
                                [BitFilt,~] = FiltroGaussiano(f,BWD,CenFeq,FiltOrd);           %Creating filter for selection of the received signal
                                BitFilt = fftshift(BitFilt);                                   %Shifting the filter for matching the received signal
                                %%
                                if RecFilBanPas==1
                                    Ix = ifft(fft(Ix).*BitFilt);
                                    IxMin = min(Ix);
                                    IxMin = repmat(IxMin,size(Ix,1),1);
                                    Ix = Ix - IxMin;%Removing the DC component from them Eletrical signal received
                                    IxMax = max(Ix);
                                    IxMax = repmat(IxMax,size(Ix,1),1);
                                    Ix = Ix./IxMax;
                                    %Ix = Ix - min(Ix);%Removing the DC component from them Eletrical signal received
                                    %Ix = Ix./max(abs(Ix));%Normalizing the eletrical signal (amplifying what is needed)
                                end
                                if ReceptorNoise==1
                                    [~,SigPower] = MeasPower(Ix);
                                    SigPower2 = SigPower-30;%10*log10(SigPower);
                                    Ix = awgn(Ix,SNR,SigPower2);
                                end
                                if ~RecFilBanPas
                                    Ix = ifft(fft(Ix).*BitFilt);
                                    IxMin = min(Ix);
                                    IxMin = repmat(IxMin,size(Ix,1),1);
                                    Ix = Ix - IxMin;%Removing the DC component from them Eletrical signal received
                                    IxMax = max(Ix);
                                    IxMax = repmat(IxMax,size(Ix,1),1);
                                    Ix = Ix./IxMax;
                                end
                                if CurrentModula==2
                                    EoutAuxF = 20*log10(abs(fftshift(fft(Ix)./length(Ix))));
                                    for jj=1:size(Ix,2)
                                        VetElecPowerUp(CurrentTest,ThisCarr,jj)= MeasPower(Ix(:,jj));
                                        [VetOptiPowerUp(CurrentTest,ThisCarr,jj),~]= findpeaks(EoutAuxF(:,jj),'SortStr','descend','NPeaks',1);
                                    end
                                else
                                %EoutAuxF = 20*log10(abs(fftshift(fft(Ix)./length(Ix))));
                                    %for jj=1:size(Ix,2)
                                        %VetElecPower(CurrentTest,ThisCarr,jj)= MeasPower(Ix(:,jj));
                                        %[VetOptiPower(CurrentTest,ThisCarr,jj),~]= findpeaks(EoutAuxF(:,jj),'SortStr','descend','NPeaks',1);
                                    %end
                                end
                                sn=((1/2)*(((1/2)^MaxNumStag)-1))/((1/2)-1);
                                AuxSyncCorr = round((sn/(2/2^IfftOrSum))*T/Ta);
                                Ix = [Ix(AuxSyncCorr+1:end,:);Ix(1:AuxSyncCorr,:)];            %Shift based on sampling sliding
                                %%                  Removing CP
                                %                 if ThisCarr==126
                                %                     EyeToPlot(CurrentTest,1:length(Ix)) = Ix;
                                %                     save(['savingforplotingeye' num2str(VetSnr)],'EyeToPlot','VetSnr','IfftOrSum','UsedModula','T','NPPB');
                                %                 end
                                if AddCP
                                    IxAux = Ix(1:end - StuffSampels,:);
                                    IxAux = reshape(IxAux,(2*NumAmosCP+NPPB),(Nb-NumBitDesc)*CurTesSiz);
                                    IxAux = IxAux(1+NumAmosCP:end-NumAmosCP,:);
                                    IxAux = reshape(IxAux,NPPB*(Nb-NumBitDesc),CurTesSiz);
                                    Ix    = IxAux;
                                end
                                %% Taking the sampling the EVM meassurement
                                %PosAuxEout = NPPB/2:NPPB:length(Ix);                       %Varriable respossible to take just the samples at the middle of the symbol
                                %IxAux      = Ix(PosAuxEout);                               %Normalizing the reference
                                %a=a+0;
                                %RxSymbAmos(CurrentTest,ObsCarrPos==ThisCarr,:) = IxAux;%RxSymbAmos = [];
                                %EvmMatRec(ObsCarrPos==ThisCarr,:) = IxAux;                       %Taking just the middle samples as references
                                %[EvmDB(CurrentTest,ThisCarr), EvmPer(CurrentTest,ThisCarr), EvmRms(CurrentTest,ThisCarr) ] = EvmCalc( EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux );
                                %[EvmDBJ(CurrentTest,ThisCarr),EvmPerJ(CurrentTest,ThisCarr),EvmRmsJ(CurrentTest,ThisCarr)] = evm1(2,'pam',EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux);
                                %% Ploting the result for qualitative analizes
                                
                                %##########################################################################
                                Ix = Ix(1+SyncPeriod*NPPB:end-SyncPeriod*NPPB,:);
                                Ix = reshape(Ix,1,size(Ix,1)*size(Ix,2));
                                if UsingMex==1
                                    [~,~,EyeOpen,EyeOpenHigh,EyeOpenLow] = Olho_mex((Ix),T,NPPB,0,1);
                                    [AberLev(CurrentTest,ThisCarr),ValsLev(CurrentTest,ThisCarr),Inter,SeqFin,LocMax,MaxVal,TxData,EoutCorrD,EoutCorr2,EoutCorr,BerOOK(CurrentTest,ThisCarr),LocMax2,SeqFin2,MaxVal2,Inter2] = RedDowOok_mex(Ix,TxDataMat,NPPB,ThisCarr,Nb,NumBitDesc,SyncPeriod,CurTesSiz,EyeOpenLow,EyeOpen);
                                else
                                    [~,~,EyeOpen,EyeOpenHigh,EyeOpenLow] = Olho((Ix),T,NPPB,0,1);
                                    [AberLev(CurrentTest,ThisCarr),ValsLev(CurrentTest,ThisCarr),Inter,SeqFin,LocMax,MaxVal,TxData,EoutCorrD,EoutCorr2,EoutCorr,BerOOK(CurrentTest,ThisCarr),LocMax2,SeqFin2,MaxVal2,Inter2] = RedDowOok(Ix,TxDataMat,NPPB,ThisCarr,Nb,NumBitDesc,SyncPeriod,CurTesSiz,EyeOpenLow,EyeOpen);
                                end
                                if CalcS==1
                                    AberLevS(CurrentTest,ThisCarr)= EyeOpen;
                                    ValsLevS(CurrentTest,ThisCarr)= EyeOpenLow + EyeOpen/2;
                                    
                                    ThisDataSize = NPPB/2:NPPB:length(Ix);
                                    EoutCorrS = zeros(1,length(ThisDataSize));                                                 %Initialization of the vector that will store the income data
                                    for kk=1:length(ThisDataSize)%length(Ix(ThisDataSize))                                       %The comparison process will be made for each symbol period
                                        %An small portion of the income signal is take for
                                        %evaluation by measuring the avarage value of the samples
                                        %taken
                                        %                 CalcMean = mean((Ix((kk-1)+SymLoc(1))));
                                        CalcMeanS = mean((Ix((kk-1)*NPPB+NPPB/2)));
                                        %Verifying the interval for each symbol received.
                                        if CalcMeanS >= EyeOpenLow + EyeOpen/2%EyeOpenLow+EyeOpen/2                        %If it is the uper level the incoming data
                                            EoutCorrS(kk) = 1;                               %is 1
                                        else                                                       %If it is the lowest level the incoming data
                                            EoutCorrS(kk) = 0;                               %is 0
                                        end
                                    end
                                    BitErrS = sum(xor(TxData,EoutCorrS));%Comparison between the Transmited and received and counting the differences
                                    BerOOKS(CurrentTest,ThisCarr) = BitErrS/(((Nb-NumBitDesc)-2*SyncPeriod)*CurTesSiz);%Calculating the ration of wrong bits and the total number of bits transmited
                                end
                                
                                %% Ploting the result for qualitative analizes
                                
                                %if ThisCarr==126
                                %EyeToPlot(CurrentTest,1:length(Ix)) = Ix;
                                %save(['savingforplotingeye' num2str(VetSnr)],'EyeToPlot','VetSnr','IfftOrSum','UsedModula','T','NPPB');
                                %end
                                if PrintinEye==1
                                    PrintInfo(PlotingThis*47,(Ix),T,NPPB);
                                    hold on;
                                    plot(t((NPPB)/2),Inter(round(SeqFin(LocMax)-MaxVal)),'k^');
                                    plot(t((NPPB)/2),Inter(SeqFin(LocMax)),'kv');
                                    plot(t((NPPB)/2),Inter(round(SeqFin(LocMax)-MaxVal/2)),'kd');
                                    plot(t((NPPB)/2),Inter2(round(SeqFin2(LocMax2)-MaxVal2)),'b^');
                                    plot(t((NPPB)/2),Inter2(SeqFin2(LocMax2)),'bv');
                                    plot(t((NPPB)/2),Inter2(round(SeqFin2(LocMax2)-MaxVal2/2)),'b*');
                                    if CalcS
                                        plot(t((NPPB)/2),EyeOpenLow,'mx');
                                        plot(t((NPPB)/2),EyeOpenHigh,'mo');
                                        plot(t((NPPB)/2),EyeOpenLow + EyeOpen/2,'md');
                                    end
                                    set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                                    %figure;
                                    %hold all;
                                    %plotpos = zeros(1,length(IxAux));
                                    %plot(IxAux,plotpos,'o','color',[1 0.4 0]);
                                    %plot(EvmMatRef(ObsCarrUsed(ThisCarr),:),plotpos,'b*');
                                    %set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                                    %a=a+0;
                                end
                                if PlotingThis
                                    figure;
                                    hold all;
                                    plot(TxData);
                                    plot(EoutCorr);
                                    set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
                                    display(BerOOK);
                                end
                                close all;
                            else
                                BerOOK(CurrentTest,ThisCarr) = 1;
                            end
                            %         end
                        end
                end
            end
        end
        %%              Saving
        if  ~exist('EvmMatRec','var')
            EvmMatRec = [];
            EvmMatRef = [];
            EvmRms = [];
            EvmDB = [];
            EvmPer = [];
            EvmRmsJ = [];
            EvmDBJ = [];
            EvmPerJ = [];
            EvmRmsUp = [];
            EvmDBUp = [];
            EvmPerUp = [];
        end
        if ~exist('CarrRecPowUp','var')
            CarrRecPowUp = 0;
        end
        if ~exist('CarrRecPowDo','var')
            CarrRecPowDo = 0;
        end
        if ~exist('VetOptiPower','var')
            VetOptiPower = 0;
            VetElecPower = 0;
        end
        if ~exist('VetElecPowerI','var')
            VetOptiPower = 0;
            VetElecPower = 0;
            VetElecPowerI = 0;
            VetElecPowerQ = 0;
        end
        if ~CalcS
            AberLevS  = [];
            ValsLevS  = [];
            BerOOKS   = [];
            Ber4PAMS  = [];
            BerDPSKS  = [];
            BerDQPSKS = [];
            AberLevIS = [];
            ValsLevIS = [];
            AberLevQS = [];
            ValsLevQS = [];
        end
        if ~SaveToScatter
            TxSigOfdm     = [];
            RxSigOfdm     = [];
            RxSigOfdmNoEq = [];
            TxSymbAmos    = [];
            RxSymbAmos    = [];
        end
        if CurrentModula == 2
            berpos1 = 1:2:size(Ber4PAM,2);
            berpos2 = 2:2:size(Ber4PAM,2);
            ber = [Ber4PAM(CurrentTest,berpos1);Ber4PAM(CurrentTest,berpos2)];
            b = Ber4PAM;
            
            ElapsedTime = toc - PastTime;
            PastTime = toc;
            if PastTime/60 > MinSavTim
                if CurrentModula == UpStModula
                    save(SavingAT,'NumbOfTest',...
                        'Ber4PAM','ValsLev','VetOptiPower','VetElecPower','AberLev',...
                        'EvmRms','EvmDB','EvmPer','DecLevDef1','DecLevDef2','DecLevDef3',...
                        'OSNRPC','OfcName','RefCarr', 'CarrUsed','SelecSetUP','FFTSplit',...
                        'CurrentMedium','UsedModula','NumCarr','PAM4Type','UpStModula',...
                        'IfftOrSum','NumCarr','lambdac','FiberLength','InitCarrUp', 'InitCarrDo'...
                        ,'CarrPass', 'CarrUsedDo','CarrUsedUp','ObsCarrUsed','CarrRecPowUp',...
                        'CarrRecPowDo','AddingNoiseF','AddingNoiseE',...
                        'AddingNoiseO','Ber4PAMS','AberLevS','ValsLevS','EvmMatRec','EvmMatRef',...
                        'EvmDB','EvmPer','EvmRms','EvmDBJ','EvmPerJ','EvmRmsJ','TxSymbAmos',...
                        'RxSymbAmos','CarSNR');
                else
                    save(SavingAT,'NumbOfTest',...
                        'Ber4PAM','ValsLev2','ValsLev','VetOptiPower','VetElecPower','AberLev',...
                        'EvmRms','EvmDB','EvmPer','DecLevDef1','DecLevDef2','DecLevDef3',...
                        'OSNRPC','OfcName','VetOptiPowerUp','VetElecPowerUp','ValsLevUp',...
                        'AberLevUp','EvmRmsUp','EvmDBUp','EvmPerUp','RefCarr', 'CarrUsed',...
                        'SelecSetUP','FFTSplit','CurrentMedium','UsedModula','NumCarr',...
                        'PAM4Type','UpStModula','IfftOrSum','NumCarr','lambdac','FiberLength',...
                        'InitCarrUp', 'InitCarrDo', 'CarrPass', 'CarrUsedDo','CarrUsedUp',...
                        'ObsCarrUsed','CarrRecPowUp', 'CarrRecPowDo',...
                        'AddingNoiseF','AddingNoiseE','AddingNoiseO','Ber4PAMS','AberLevS',...
                        'ValsLevS','EvmMatRec','EvmMatRef','EvmDB','EvmPer','EvmRms','EvmDBJ',...
                        'EvmPerJ','EvmRmsJ','TxSymbAmos','RxSymbAmos','CarSNR');
                end
                PastTime = 0;
                tic;
            end
        elseif CurrentModula == 1
            berpos1 = 1:2:size(BerOOK,2);
            berpos2 = 2:2:size(BerOOK,2);
            ber = [BerOOK(CurrentTest,berpos1);BerOOK(CurrentTest,berpos2)];
            b = BerOOK;
            
            ElapsedTime = toc - PastTime;
            PastTime = toc;
            if PastTime/60 > MinSavTim
                save(SavingAT,'NumbOfTest',...
                    'BerOOK','ValsLev','VetOptiPower','VetElecPower','AberLev',...
                    'EvmRms','EvmDB','EvmPer','DecLevDef1','DecLevDef2','DecLevDef3',...
                    'OSNRPC','OfcName','RefCarr', 'CarrUsed','SelecSetUP','FFTSplit',...
                    'CurrentMedium','UsedModula','NumCarr','PAM4Type','UpStModula',...
                    'IfftOrSum','NumCarr','lambdac','FiberLength','InitCarrUp', 'InitCarrDo'...
                    ,'CarrPass', 'CarrUsedDo','CarrUsedUp','ObsCarrUsed','CarrRecPowUp',...
                    'CarrRecPowDo','AddingNoiseF','AddingNoiseE',...
                    'AddingNoiseO','BerOOKS','AberLevS','ValsLevS','EvmMatRec','EvmMatRef',...
                    'EvmDB','EvmPer','EvmRms','EvmDBJ','EvmPerJ','EvmRmsJ','TxSymbAmos',...
                    'RxSymbAmos','CarSNR');
                PastTime = 0;
                tic;
            end
        elseif CurrentModula == 3
            berpos1 = 1:2:size(BerDQPSK,2);
            berpos2 = 2:2:size(BerDQPSK,2);
            ber=[BerDQPSK(CurrentTest,berpos1);BerDQPSK(CurrentTest,berpos2)];
            b = BerDQPSK;
            ElapsedTime = toc - PastTime;
            PastTime = toc;
            if PastTime/60 > MinSavTim
                save(SavingAT,'NumbOfTest','BerDQPSK','VetOptiPower',...
                    'AberLevQ','ValsLevQ','AberLevI','ValsLevI','DecLevDef1','DecLevDef2',...
                    'DecLevDef3','OSNRPC','OfcName','RefCarr','CarrUsed','SelecSetUP',...
                    'FFTSplit','CurrentMedium','UsedModula','NumCarr','PAM4Type',...
                    'UpStModula','IfftOrSum','NumCarr','lambdac','FiberLength',...
                    'InitCarrUp','InitCarrDo','CarrPass','CarrUsedDo','CarrUsedUp',...
                    'ObsCarrUsed','CarrRecPowUp','CarrRecPowDo',...
                    'AddingNoiseF','AddingNoiseE','AddingNoiseO','BerDQPSKS','AberLevIS',...
                    'ValsLevIS','AberLevQS','ValsLevQS','EvmMatRec','EvmMatRef','EvmDB',...
                    'EvmPer','EvmRms','EvmDBJ','EvmPerJ','EvmRmsJ','TxSymbAmos',...
                    'RxSymbAmos','CarSNR');
                PastTime = 0;
                tic;
            end
        elseif CurrentModula == 4
            berpos1 = 1:2:size(BerDPSK,2);
            berpos2 = 2:2:size(BerDPSK,2);
            ber = [BerDPSK(CurrentTest,berpos1);BerDPSK(CurrentTest,berpos2)];
            b = BerDPSK;
            
            ElapsedTime = toc - PastTime;
            PastTime = toc;
            if PastTime/60 > MinSavTim
                save(SavingAT,'NumbOfTest','BerDPSK','ValsLev',...
                    'VetOptiPower','VetElecPower','AberLev','DecLevDef1','DecLevDef2',...
                    'DecLevDef3','OSNRPC','OfcName','RefCarr', 'CarrUsed','SelecSetUP',...
                    'FFTSplit','CurrentMedium','UsedModula','NumCarr','PAM4Type',...
                    'UpStModula','IfftOrSum','NumCarr','lambdac','FiberLength','InitCarrUp',...
                    'InitCarrDo','CarrPass', 'CarrUsedDo','CarrUsedUp','ObsCarrUsed',...
                    'CarrRecPowUp','CarrRecPowDo','AddingNoiseF',...
                    'AddingNoiseE','AddingNoiseO','BerDPSKS','AberLevS','ValsLevS',...
                    'EvmMatRec','EvmMatRef','EvmDB','EvmPer','EvmRms','EvmDBJ','EvmPerJ',...
                    'EvmRmsJ','TxSymbAmos','RxSymbAmos','CarSNR');
                PastTime = 0;
                tic;
            end
        elseif CurrentModula == 5
            if (UseOfdmElect)&&(UsedModula==5)&&(NumCarr==2)
                berpos1 = 1:size(BerOFDM,2);
                berpos2 = 1:size(BerOFDM,2);
                ber = BerOFDM.';
                b = BerOFDM;
                ElapsedTime = toc - PastTime;
                PastTime = toc;
                if PastTime/60 > MinSavTim
                    save(SavingAT,'SelecGaus','OfdMod','SelModTp',...
                        'UsingHermitian','NuDmSuDi','DmtMinM','DmtPas','TapN',...
                        'ExtSam','NpEl','NFFT','OveT','M','ControlZt','CpLe','BW','DatSiz','Te',...
                        'Tu','Tg','Ts','g','Dtf','N','Ns','DifNsN','k','MZ_Input_File','NuSaTs',...
                        'NuSaTu','NuSaTg','OBw','Ofc','Ofs','OvSam','NuPart','NuSuSe','DmtMve',...
                        'nRb','nRbAgg','Zt','ZtC','NumFraPar','NumFra',...
                        'NPOFEX','NPPOF','NPSTUf','DmtKvep',...
                        'BerOFDM','EvmDB','EvmPer','EvmRms','EvmDBJ','EvmPerJ','EvmRmsJ',...
                        'SaveRxEq','SaveRxNotEq');
                end
            else
                berpos1 = 1:2:size(BerOFDM,2);
                berpos2 = 2:2:size(BerOFDM,2);
                ber = [BerOFDM(CurrentTest,berpos1);BerOFDM(CurrentTest,berpos2)];
                b = BerOFDM;
                
                ElapsedTime = toc - PastTime;
                PastTime = toc;
                if PastTime/60 > MinSavTim
                    save(SavingAT,'NumbOfTest','BerOFDM','OSNRPC','OfcName','DmtKvep',...
                        'RefCarr', 'CarrUsed','SelecSetUP','FFTSplit','CurrentMedium','OfdMod',...
                        'UsedModula','NumCarr','PAM4Type','UpStModula','IfftOrSum','NumCarr',...
                        'lambdac','FiberLength','InitCarrUp','InitCarrDo','CarrPass','CpLe',...
                        'CarrUsedDo','CarrUsedUp','ObsCarrUsed','AddingNoiseF','AddingNoiseE',...
                        'AddingNoiseO','EvmMatRec','EvmMatRef','EvmDB','EvmPer','EvmRms','M',...
                        'EvmDBJ','EvmPerJ','EvmRmsJ','TxSigOfdm','RxSigOfdm','RxSigOfdmNoEq',...
                        'NPSTUf','nRbAgg','nRb','DmtMve','NuSuSe','NuPart','OvSam','Ofs','Ofc',...
                        'OBw','MZ_Input_File','Tg','ControlZt','OveT','NpEl','ExtSam','TapN',...
                        'SelecGaus','NuDmSuDi','OfdMod','SelModTp','DmtPas','DmtMinM','CarSNR','BerToPlotOfdm');
                    PastTime = 0;
                    tic;
                end
            end
        end
        CurrentTest;
        VetSnr;
        baux = mean(b,1);
        bdis = [baux(berpos1);baux(berpos2)];
        display(bdis);
        display(sprintf(['Elapsed Time(s) = %f | Current Test = %d '...
            '| \nTested = %0.1f | Threshold = %0.2f | Decision = %0.2f '...
            '\n'],ElapsedTime,CurrentTest,VetSnr,((CurrentTest)*ceil(...
            0.8*Nb)),(5/min(mean(b(:,berpos2),1)))));
        a=a+1;
        if CurrentTest > MinTesNumb
            if (((CurrentTest)*ceil(0.8*Nb))>(5/min(mean(b(:,berpos2),1))))
                break;
            end
        end
    end
    if CurrentModula == 2
        if CurrentModula == UpStModula
            save(SavingAT,'NumbOfTest',...
                'Ber4PAM','ValsLev','VetOptiPower','VetElecPower','AberLev',...
                'EvmRms','EvmDB','EvmPer','DecLevDef1','DecLevDef2','DecLevDef3',...
                'OSNRPC','OfcName','RefCarr', 'CarrUsed','SelecSetUP','FFTSplit',...
                'CurrentMedium','UsedModula','NumCarr','PAM4Type','UpStModula',...
                'IfftOrSum','NumCarr','lambdac','FiberLength','InitCarrUp','InitCarrDo',...
                'CarrPass', 'CarrUsedDo','CarrUsedUp','ObsCarrUsed','CarrRecPowUp',...
                'CarrRecPowDo','AddingNoiseF','AddingNoiseE',...
                'AddingNoiseO','Ber4PAMS','AberLevS','ValsLevS','EvmMatRec','EvmMatRef',...
                'EvmDB','EvmPer','EvmRms','EvmDBJ','EvmPerJ','EvmRmsJ','TxSymbAmos',...
                'RxSymbAmos','CarSNR');
        else
            save(SavingAT,'NumbOfTest',...
                'Ber4PAM','ValsLev2','ValsLev','VetOptiPower','VetElecPower','AberLev',...
                'EvmRms','EvmDB','EvmPer','DecLevDef1','DecLevDef2','DecLevDef3',...
                'OSNRPC','OfcName','VetOptiPowerUp','VetElecPowerUp','ValsLevUp',...
                'AberLevUp','EvmRmsUp','EvmDBUp','EvmPerUp','RefCarr', 'CarrUsed',...
                'SelecSetUP','FFTSplit','CurrentMedium','UsedModula','NumCarr',...
                'PAM4Type','UpStModula','IfftOrSum','NumCarr','lambdac','FiberLength',...
                'InitCarrUp', 'InitCarrDo', 'CarrPass', 'CarrUsedDo','CarrUsedUp',...
                'ObsCarrUsed','CarrRecPowUp', 'CarrRecPowDo',...
                'AddingNoiseF','AddingNoiseE','AddingNoiseO','Ber4PAMS','AberLevS',...
                'ValsLevS','EvmMatRec','EvmMatRef','EvmDB','EvmPer','EvmRms','EvmDBJ',...
                'EvmPerJ','EvmRmsJ','TxSymbAmos','RxSymbAmos','CarSNR');
        end
    elseif CurrentModula == 1
        save(SavingAT,'NumbOfTest',...
            'BerOOK','ValsLev','VetOptiPower','VetElecPower','AberLev',...
            'EvmRms','EvmDB','EvmPer','DecLevDef1','DecLevDef2','DecLevDef3',...
            'OSNRPC','OfcName','RefCarr', 'CarrUsed','SelecSetUP','FFTSplit',...
            'CurrentMedium','UsedModula','NumCarr','PAM4Type','UpStModula',...
            'IfftOrSum','NumCarr','lambdac','FiberLength','InitCarrUp', 'InitCarrDo'...
            ,'CarrPass', 'CarrUsedDo','CarrUsedUp','ObsCarrUsed','CarrRecPowUp',...
            'CarrRecPowDo','AddingNoiseF','AddingNoiseE',...
            'AddingNoiseO','BerOOKS','AberLevS','ValsLevS','EvmMatRec','EvmMatRef',...
            'EvmDB','EvmPer','EvmRms','EvmDBJ','EvmPerJ','EvmRmsJ','TxSymbAmos',...
            'RxSymbAmos','CarSNR');
    elseif CurrentModula == 3
        save(SavingAT,'NumbOfTest','BerDQPSK','VetOptiPower',...
            'AberLevQ','ValsLevQ','AberLevI','ValsLevI','DecLevDef1','DecLevDef2',...
            'DecLevDef3','OSNRPC','OfcName','RefCarr','CarrUsed','SelecSetUP',...
            'FFTSplit','CurrentMedium','UsedModula','NumCarr','PAM4Type',...
            'UpStModula','IfftOrSum','NumCarr','lambdac','FiberLength',...
            'InitCarrUp','InitCarrDo','CarrPass','CarrUsedDo','CarrUsedUp',...
            'ObsCarrUsed','CarrRecPowUp','CarrRecPowDo',...
            'AddingNoiseF','AddingNoiseE','AddingNoiseO','BerDQPSKS','AberLevIS',...
            'ValsLevIS','AberLevQS','ValsLevQS','EvmMatRec','EvmMatRef','EvmDB',...
            'EvmPer','EvmRms','EvmDBJ','EvmPerJ','EvmRmsJ','TxSymbAmos',...
            'RxSymbAmos','CarSNR');
    elseif CurrentModula == 4
        save(SavingAT,'NumbOfTest','BerDPSK','ValsLev',...
            'VetOptiPower','VetElecPower','AberLev','DecLevDef1','DecLevDef2',...
            'DecLevDef3','OSNRPC','OfcName','RefCarr', 'CarrUsed','SelecSetUP',...
            'FFTSplit','CurrentMedium','UsedModula','NumCarr','PAM4Type',...
            'UpStModula','IfftOrSum','NumCarr','lambdac','FiberLength','InitCarrUp',...
            'InitCarrDo','CarrPass', 'CarrUsedDo','CarrUsedUp','ObsCarrUsed',...
            'CarrRecPowUp','CarrRecPowDo','AddingNoiseF',...
            'AddingNoiseE','AddingNoiseO','BerDPSKS','AberLevS','ValsLevS',...
            'EvmMatRec','EvmMatRef','EvmDB','EvmPer','EvmRms','EvmDBJ','EvmPerJ',...
            'EvmRmsJ','TxSymbAmos','RxSymbAmos','CarSNR');
    elseif CurrentModula == 5
        if (UseOfdmElect)&&(UsedModula==5)&&(NumCarr==2)
            save(SavingAT,'SelecGaus','OfdMod','SelModTp',...
                'UsingHermitian','NuDmSuDi','DmtMinM','DmtPas','TapN',...
                'ExtSam','NpEl','NFFT','OveT','M','ControlZt','CpLe','BW','DatSiz','Te',...
                'Tu','Tg','Ts','g','Dtf','N','Ns','DifNsN','k','MZ_Input_File','NuSaTs',...
                'NuSaTu','NuSaTg','OBw','Ofc','Ofs','OvSam','NuPart','NuSuSe','DmtMve',...
                'nRb','nRbAgg','Zt','ZtC','NumFraPar','NumFra','DmtKvep',...
                'NPOFEX','NPPOF','NPSTUf','BerOFDM','EvmDB','EvmPer',...
                'EvmRms','EvmDBJ','EvmPerJ','EvmRmsJ','SaveRxEq','SaveRxNotEq');
        else
            save(SavingAT,'NumbOfTest','BerOFDM','OSNRPC','OfcName','DmtKvep',...
                'RefCarr', 'CarrUsed','SelecSetUP','FFTSplit','CurrentMedium','OfdMod',...
                'UsedModula','NumCarr','PAM4Type','UpStModula','IfftOrSum','NumCarr',...
                'lambdac','FiberLength','InitCarrUp','InitCarrDo','CarrPass','CpLe',...
                'CarrUsedDo','CarrUsedUp','ObsCarrUsed','AddingNoiseF','AddingNoiseE',...
                'AddingNoiseO','EvmMatRec','EvmMatRef','EvmDB','EvmPer','EvmRms','M',...
                'EvmDBJ','EvmPerJ','EvmRmsJ','TxSigOfdm','RxSigOfdm','RxSigOfdmNoEq',...
                'NPSTUf','nRbAgg','nRb','DmtMve','NuSuSe','NuPart','OvSam','Ofs','Ofc',...
                'OBw','MZ_Input_File','Tg','ControlZt','OveT','NpEl','ExtSam','TapN',...
                'SelecGaus','NuDmSuDi','OfdMod','SelModTp','DmtPas','DmtMinM','CarSNR','BerToPlotOfdm');
        end
    end
    if ChanAwgn||ReceptorNoise
        PlotingNoise;
    end
end
if exist('CompEbn0','var')
    
    SaveBerEbn0 = [num2str(UsedModula) '_' num2str(ThisM) '_' VerThisMod];
    %     SaveBerEbn0Path = [pwd 'test0'];
    
    save([SaveBerEbn0PreFix SaveBerEbn0 SaveBerEbn0SuFix],'CompEbn0','BerTheo',...
        'ThisSimuBer','VerThisMod','LegendThis','ThisM','b');
    saveas(gcf,[SaveBerEbn0PreFix SaveBerEbn0 SaveBerEbn0FigSuFix]);
end
















