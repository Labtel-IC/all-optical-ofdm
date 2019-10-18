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
%c   ThisPath      : Storage local for the path of this current simulation
%c                   folder.
%c   UofcsLocal    : Store local for where the file to load the OFCS can be
%c                   found. It is the may file needed for the simulation.
%c                   If it doesn't exist it must be created.
%c   UofcsLocal    : It is the name of the main file to be load with its
%c                   full location path.
%c   SigAten       : Flag to indicate if the singal from the OFCS will be
%c                   attenuated or not.
%c   SnrRef        : Flag to indicate which signal will be used as
%c                   reference for the awg noise matlab function.
%c   FFTSplit      : Flag to indicate if at the receptor the income OFDM
%c                   signal will be splited by the OFFT or by Filters.
%c   CarSNR        : Actual sinal to noise ratio that will be applied.
%c   PastTime      : To mark the time elapsed by the tic toc function.
%c   Atenuation    : Actual attenuation applied to the income OFCS.
%c   NumbOfTest    : Number of repetitions for one single simulation
%c                   aiming to achive a target number of bits transmited.
%c   TestBundle    : Flag to set if there will be a bundle of tests or
%c                   just one per modulations.
%c   CurrentMedium : Flag to sellect if the channel will be SSMF our B2B.
%c   PAM4Type      : Variable used to set which type of PAM4 will be used.
%c   UsedModula    : This is a vector that store the sequency of all types
%c                   of modulations that will be used on this whole
%c                   simulation.
%c   ThisTest      : The actual number of testes passed to the for loop.
%c   ThisModula    : The actual modulation testes passed to the for loop.
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
%Checking if this main file was not called by another script
if ~exist('UsingAllCarriers','var')
    clear ; clc; close all;
end
%MainInputData;
%%                    OFC Generator Parameters
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

%%                 Control Variable
SelecSetUP    = 1;                                                         %Select the simulated setup 1 for original 0 for IFFT/OFFT at receptor
SendingDowStr = 1;                                                         %Set 1 to send data downstream 0 to not send.
SigAten       = 0;                                                         %Set 1 to atenuate the OFCS signal. Set to zero to turn off attenuation
SnrRef        = 0;                                                         %Set 1 to select the own signal for the noise reference at receptor.
%Set 0 for the pilote carrier as noise reference at receptor.
FFTSplit      = 1;                                                         %Set 1 to use the OFFT to split the income signal on receptor. Set 0 to
%use filter to split the income signal at receptor.
CarSNR        = VetSnr;                                                    %Set the snr to be applied at receptor.
Nsamp         = 2^9;%1;%                                                   %Number of samples per bit, it was necessary because the number of carriers here used
RecFilBanPas  = 1;                                                         %Filter after photodetector, 1 filter Before noise, 0 filter After noise.
SendingUp     = 1;                                                         %Select if the carrier upstream will be send. Set 1 to send, set 0 to not send.
SaveToScatter = 0;                                                         %Set 1 to save points to create scatter plot latter Set 0 to not save points (keep it zero to save data space)
TimeSys       = 1;                                                         %Select which time period will be accounted to the noise, set 1 to use T or set 0 to use t(2*NumAmosCP+NPPB).
%When CP is used Seting TimeSys to 1 gives higher SNR
PastTime      = 0;                                                         %Storage for the time flow on simulation
% Atenuation    = 10^(4.0);                                                %Set how much the OFCS will be attenuated
NumbOfTest    = 20;                                                      %The number of repetition that one test will be done
TestBundle    = 1;                                                         %Set if the scrip will be run once or many times
CurrentMedium = 1;                                                         %Sellection of channel 1 for Fiber other of B2B
PAM4Type      = [1 0];                                                     %Sellection the type of 4PAM that will be tested
GatingAtReception = 0;
UseOfdmElect  = 0;
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
PrintinEye    = 1;                                                         %Printing Eye Diagram of modulations

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
UsedModula = [2 []];
UpStModula = [2 []];
if UsedModula ~= 5
    UseOfdmElect = 0;
end
%The following verification will determinate which type of test will be on
%running, a test bundle or a single test. In the section the user may
%change with variables will be loaded for testing.
if TestBundle                                                              %Test Bundle selected
    ThisTest   = NumbOfTest;                                               %The number of tests per type of modulation
    ThisModula = UsedModula;                                               %The modulations under test
else
    ThisTest   = 1;                                                        %The number of tests per type of modulation
    ThisModula = UsedModula(1);                                            %The modulations under test
end
%%          Variables for the Medium Fiber
%This part of the simulation it is important to chose before hand some
%general parameters, that will affect the whole simulation.

if ~exist('UsingAllCarriers','var')                                        %Those variables will just be seted if any other outsider script haven't call them yet
    IfftOrSum    = 1;                                                      %Selection if the modulated signal will be assembled by just a adder or by the OIFFT
    %Set 1 for the OFFT or set 0 for the adder
    if (UsedModula==5)&&(UseOfdmElect==1)
        NumCarr  = 2;                                                      %Total number of carriers (DownStream and UpStream) - This variable must be even
    else
        NumCarr  = 2;                                                      %Total number of carriers (DownStream and UpStream) - This variable must be even
    end
    lambdac      = 1550e-9;                                                %Central wave length of the signal throughout the fiber
    FiberLength  = VetSnr;                                                     %Setting the length of the fiber in this simulation
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

a=0;

%if TestBundle                                                      %When many tests will be run, the OFCS will loaded from a file
    load([OfcName '.mat']);                                        %Loading a OFC from memory.
%else                                                               %When it will be just a general test it will ask if user wants a new OFCS
    %Check = input(['Do you want to load a new OFC?\n If not pre'...%Fristly the script check if the user want to
        %'ss Enter to continue. The OFC will be loaded from file.\n']...
        %,'s');%load an OFC from file or to generate a new one
    %if ~isempty(Check)
        %GenerateOFC;                                               %Generate a new OFC if one wishes so.
    %else
        %load([OfcName '.mat']);                                    %Otherwise load a OFC from memory.
    %end%
%end

if ~exist('freqGHz','var')
    tps        = t/1E-12;
    freqTHz    = time2freq_lamb_2(tps);
    freqGHz    = freqTHz*1e-3;			% Frequencia em GHz
    freqGHz    = -fftshift(freqGHz);
    freqGHz(1) = freqGHz(2);
end
% Ploting the OFC that will be used in this simulation
%opticalnoiseadded;
Ploting = 0;
PrintInfo(Ploting*1,f,Eout);
%xmax = max(20*log10(abs(fftshift(fft(Eout)./length(Eout)))));
%axis([-1e11 2e12 xmax-60 xmax+0.1]);

%% DatumElectricTransmission;%%       General Parameter for the simulation
ON          =  1;                                                          %DO NOT CHANGE THIS PARAMETER
OFF         =  0;                                                          %DO NOT CHANGE THIS PARAMETER
Ploting     = OFF;                                                         %Use On or Off to say if the figures should be plot or not
PlotingThis = ON;                                                          %Use On or Off to say if the figures should be plot or not
Selecting   = ON;                                                         %Use to filter or not the carriers to be transmited
ModAll      = OFF;                                                         %Set if all carriers will be modulate at once or not
AddCP       = ON;                                                          %Add the Cycle Prefix on the data stream
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
MaxNumCarrTx = 2^MaxNumStagTx;                                             %With the number of stages the actual nubmer of carriers can be calculated
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
    MzFilesGenerated = 1;
end
%%         Creating file for the DP-MZM
phi_0_steps  = 1;%Number of different phase between arms of the DP-MZM
V1bias_steps = 1;%Number of different Bias voltage at the Up arm
V2bias_steps = 1;%Number of different Bias voltage at the Lower arm
%
%IMPORTANT: The RF amplitude Vpp must respect the device limits!
phi_0_vet     = -2.05;                                                     %Variation of phase betwee arms.
V1bias        = 1.752;                                                     %Variation of bias voltage at the Up arm.
V2bias        = 1.052;                                                     %Variation of bias voltage at the Lower arm.

% Differently from the previouly Input Data file. This scrip does not have
% a previous configuration file for the MZM. Thus, it can create a new file
% to configure the MZM-IQ.
Local_Dp = [pwd '\input_files_dp\'];                                       %Seting the location where this file will be stored.

L          = 10;
U10        = 1.550;
U20        = 1.625;
U_pi1      = 3.1;
U_pi2      = 3.25;
eta1       = 89;% 89;
eta2       =-12.7;%-12.7;
V1pi0      = 0.350;
V2pi0      = 1.200;
Vphi0      = 0.100;
%
nopt       =  2.17;
nel        =  2.60;
alfa_ins   =  5.1;
phi_0      =  1.85;
alfa0      =  0.55;

C     = (eta1+eta2)/(eta1-eta2);       	%  Parametro de chirp
alfa0 = 10^(alfa0/20); 			% [1/cm.GHz^0.5]
if ~exist('MzFilesDpGenerated','var')
    if ~exist('L','var')
    MZ_Input_Data_DP;                                                      %Loading the basic data to be saved
end
S=1;                                                                       %Variabel that give the name for the input file
for Inc0=1:phi_0_steps 
    for Inc1=1:V1bias_steps                                                %Secundary loop for the possible combinations
        for Inc2=1:V2bias_steps                                            %Secundary loop for the possible combinations
            phi_0 = phi_0_vet(Inc0);
            U10               = V1bias(Inc1);                              %Variation of the Vbias
            U20               = V2bias(Inc2);                              %Variation of the Vbias
            [MZ_Input_Sdata] = Set_MZ_Input_Data_DP(S,L,U10,U20,U_pi1,...  %Functionresponsible to create the input files
            U_pi2,eta1,eta2,V1pi0,V2pi0,Vphi0,nopt,nel,alfa_ins,phi_0,...
                                                           alfa0,Local_Dp);
            S                = S+1;
        end
    end
end                                                    %Create the files for the MZM input data
    MzFilesDpGenerated = 1;
end
MZ_Input_File_Dp = 1;                                                      %Telling the script witch file will be used for modulation                                               %Telling the script witch file will be used for modulation
%%           Control Variables
%This variable will be used to select the type of modulation used for the
%datum Enconding. Type 'DQPSK' for DQPSK encoding, '4PAM' for 4PAM
%encoding. Type 'NoN' if an encoding is not necessary.
switch CurrentModula
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
if CurrentMedium
    if ~exist('FiberDelay','var')
        [FiberDelay,PowRef] = ProgationDelay(RefCarr,NumCarr,t,lambdac,T,...
            FiberLength,0,f,Ploting);
    end
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
        %% OfdmInputData;
        SelModFilt= 0;
        SelecGaus = 1;
        OfdMod    = 'qam';                                                %Chossing the type of modulation for the electrical subcarriers
        SelModTp  = 'BP';                                                  %Selecting with OFDM will be in Base Band "BP", Amplitude Modulation
        %with double side band "AM" or with single side band "AMSSB" when a
        %different option is choosen this script will implement phase modulation.
        UsingHermitian = 0;
        if ~UseOfdmElect
            UsingHermitian = 1;
        end
        ForcOvSam = 0;
        UsingParCar= 0;
        NuDmSuDi  = 1;
        DmtPas    = 2;
        DmtMinM   = 3;
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
                if UsingParCar
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
        
        % SNR = CarSNR + 10*log10(log2(max(DmtMve))) - 10*log10(1*Tu*OBw);
        SNR = CarSNR + 10*log10(log2(max(DmtMve))) - 10*log10(OvSam);
        
        %         NPOFEX = floor(NumFrE*OvSam*NFFT/NumFra);
        %         NPPOF  = floor((2^19)/((NPOFEX*NumFra)+(OvSam*NumFra*NFFT)));
        %         NPSTUf = 2^19 - NPPOF*NumFra*OvSam*NFFT;
        %Variables for the Selection of carriers - Wave Divided.
        fin           = RefCarr*fc;                                        %The frequency of the first carrier of the optical COMB
        FBWD          = 1.0*fc;                                            %The filter BandWidth to split each individual carrier
        Order         = 5;                                                 %The order of the filter used to split each carrier
        a=1;
    case 'DPSK' %Tested to fiberlength of 80 -> worked
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
        BWD2          = 1*fc;
        if TimeSys
            SNR = CarSNR + 10*log10(1) - 10*log10(1*T*BWD2);
        else
            SNR = CarSNR + 10*log10(1) - 10*log10(1*t(2*NumAmosCP+NPPB)*BWD2);
        end
    case 'DQPSK'%Tested to fiberlength of 40 -> work in progress
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
        BWD2          = 1*fc;
        if TimeSys
            SNR = CarSNR + 10*log10(2) - 10*log10(1*T*BWD2);
        else
            SNR = CarSNR + 10*log10(2) - 10*log10(1*t(2*NumAmosCP+NPPB)*BWD2);
        end
    case '4PAM'%Tested to fiberlength of 1 -> work in progress
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
        BWD2          = 1*fc;
        if TimeSys
            SNR = CarSNR + 10*log10(2) - 10*log10(1*T*BWD2);
        else
            SNR = CarSNR + 10*log10(2) - 10*log10(1*t(2*NumAmosCP+NPPB)*BWD2);
        end
    otherwise %Tested to fiberlength of 80 -> worked
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
        BWD           = fc;                                                %Badnwidth for the bit conformation
        CenFeq        = 0;                                                 %Center frequency for the bit conformation
        FiltOrd       = 1;                                                 %The order for the filter which will conform data
        MZ_Input_File = 5;                                                 %Telling the script witch file will be used for modulation
        %Variables for the Selection of carriers - Wave Divided.
        fin           = RefCarr*fc;                                        %The frequency of the first carrier of the optical COMB
        FBWD          = 1.0*fc;                                            %The filter BandWidth to split each individual carrier
        Order         = 5;                                                 %The order of the filter used to split each carrier
        BWD2          = 1*fc;
        if TimeSys
            SNR = CarSNR + 10*log10(1) - 10*log10(1*T*BWD2);
        else
            SNR = CarSNR + 10*log10(1) - 10*log10(1*t(2*NumAmosCP+NPPB)*BWD2);
        end
end
%%      Filter to Select the Carriers
[ SelecFilt ] = FiltroGaussiano(f,(NumCarr)*fc,((NumCarr+1)/2)*fc,11);
% [ SelecFilt ] = Filtro_Retangular( 34*fc,16*fc,f);
SelecFilt = fftshift(SelecFilt);
%%      Filter for the EVM measurement taking
EvmFilt = FiltroGaussiano(f,fc,0,1);
EvmFilt = fftshift(EvmFilt);
%%       Variables for Syncronism

%%    Control Variables of this program
ON            =  1; %DO NOT CHANGE THIS PARAMETER
OFF           =  0; %DO NOT CHANGE THIS PARAMETER
Ploting       = OFF; %It is used to select an ploting.
PlotingThis   = ON; %It is used to select an especific ploting among others
Elucidation   = OFF;%It is used to enable the OFFT example
UsePhotDiod   = OFF;%To select or not the usage of the photodiod
SendingUpStre = OFF; %To send the Data Upstream
% ReceptorNoise = ON;
BitFilBanWid  = 0.9*Rb;%Badnwidth for the bit conformation
BitFiltCenFre = 0;  %Center frequency for the bit conformation
BitFilOrd     = 3;  %The order for the filter which will conform data
DecMod        = 1;  %Selec how the decission level will be measured
% 1 for measuring the zeros on the eye diagram (best)
% 0 for taking the meddle of the zeros vector
InitCarr      = InitCarrDo;  %Selecting if the data processed will be down or up
%stream: 1 for downstream and 2 for upstream.
%%    Constants to controu the Optical FFT
% The optical FFT here implemented can take any number of carriers. The
% user must stay realistic because the number of actual devices used to
% implement this setup increases by the power of 2, thus higher orders of
% the FFT may not be feaseble to be implemented as the stability of each
% device may vary with temperature, aging and othter features.
MaxNumStag = nextpow2(NumCarr);                                        %Next number, power of 2, by the given number of carriers
MaxNumCarr = 2^MaxNumStag;                                             %With the number of stages the actual nubmer of carriers can be calculated
% end

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
% switch MaxNumStag
%     case 1
%         VetThisCarr = [1 2];
%     case 2
%         VetThisCarr = [1 3 2 4];
%     case 3
%         VetThisCarr = [1 5 3 7 2 6 4 8];
%     case 4
%         VetThisCarr = [1 9 5 13 3 11 7 15 2 10 6 14 4 12 8 16];
%     case 5
%         VetThisCarr = [1 17 9 25 5 21 13 29 3 19 11 27 7 23 15 31 2 18...
%                                   10 26 6 22 14 30 4 20 12 28 8 24 16 32];
%     otherwise
%         MaxNumStag = 1;
%         MaxNumCarr = 2^MaxNumStag;
%         VetThisCarr = [1 2];
% end
%At this point if for some misteriouse reason the Number of carriers or the
%Maxmum order of the OFFT were seted beyond the maximum point this
%verification will pointout a error while the OpticalOFFT function is not
%optmized for an 'n' order filter.
% if MaxNumStag>5
% error(['MaxNumStag allowed value for MaxNumStag was achieved. Please'...
%                               ' cange the number Maximum of carriers.']);
% end

Count = ones(MaxNumStag,1);                                                %Control Variable for adding correct variation
E0 = 0;                                                                    %Second input, usualy zero
ActStag = 1;                                                               %Starting point of Optical FFT
%%               Variables for syncronism
%Depending of each type of modulation, it is needed to define know symble
%for syncronism. As it deals with optcal devices and optical paths the
%signal takes time to travel, which leads to an sampling miss mach of the
%received and transmite data. Thus at the receiver the syncronism variable
%are created to approach this problem for each different type of
%modulation.
switch Modulation
    case 'OFDM'
        BWD           = 1.0*fc;                                        %Badnwidth for the bit conformation
        CenFeq        = 0;                                                 %Center frequency for the bit conformation
        FiltOrd       = 0;                                                 %The order for the filter which will conform data
    case 'DPSK'
        %%             Syncing DPSK signal
        
        BWD2          = 1*fc;                                          %Badnwidth for the bit conformation
        CenFeq2       = 0;                                                 %Center frequency for the bit conformation
        FiltOrd2      = 1;                                                 %The order for the filter which will conform data
        SyncSymb      = [zeros(1,10) ones(1,4) zeros(1,10)];               %Generating the sync-word at the beginning of the data stream
        IniSyncPos    = 1+5*(2*NumAmosCP+NPPB);
        SyncPos       = 19*(2*NumAmosCP+NPPB);                             %Seting the possition where the sync-symble is
        SyncSymbEnd   = SyncSymb;                                          %Generating the sync-word at the end of the data stream
        SyncPeriod    = length(SyncSymb);                                  %Seting how many symble periods will be within the sync-word
        AuxSyncLength = SyncPeriod*(2*NumAmosCP+NPPB);                     %auxiliar variable to set where is the sync-symble
        tsync         = linspace(0,SyncPeriod*T,(2*NumAmosCP+NPPB)*...
            SyncPeriod);%create a time vetor for the sync-word to help on the ploting of this signal
        fsync         = time2freq(tsync);                                  %creating the the frequency vector of the sync-word for hereafter filtering
        BWD           = 2*fc;                                          %Badnwidth for the bit conformation
        CenFeq        = 0;                                                 %Center frequency for the bit conformation
        FiltOrd       = 1;                                                 %The order for the filter which will conform data
        BitFilBanWid  = 10e9;                                              %Badnwidth for the bit conformation
        BitFiltCenFre = 0;                                                 %Center frequency for the bit conformation
        BitFilOrd     = 1;                                                 %The order for the filter which will conform data
        [BitFilt,~]   = FiltroGaussiano(fsync,BWD,CenFeq,FiltOrd);         %Creating filter for conformation of the input information
        BitFilt       = fftshift(BitFilt);
        
        SyncSymb      = rectpulse(SyncSymb,2*NumAmosCP+NPPB);
        SyncSymbEnd   = rectpulse(SyncSymbEnd,2*NumAmosCP+NPPB);
        SencondAdjust = OFF;                                               %Setting whether is needed to sync the information by the Syncronism Symble at the End
        LowSymPer     = 1;
        UpeSymPer     = 1;
    case 'DQPSK'
        %%             Syncing DQPSK signal
        
        BWD2          = 1.0*fc;                                        %Badnwidth for the bit conformation
        CenFeq2       = 0;                                                 %Center frequency for the bit conformation
        FiltOrd2      = 1;                                                 %The order for the filter which will conform data
        SyncSymb      = [zeros(1,10) ones(1,4) zeros(1,10)];               %Generating the sync-word at the beginning of the data stream
        IniSyncPos    = 1+5*(2*NumAmosCP+NPPB);
        SyncPos       = 19*(2*NumAmosCP+NPPB);                             %Seting the possition where the sync-symble is
        SyncSymbEnd   = SyncSymb;                                          %Generating the sync-word at the end of the data stream
        SyncPeriod    = length(SyncSymb);                                  %Seting how many symble periods will be within the sync-word
        AuxSyncLength = SyncPeriod*(2*NumAmosCP+NPPB);                     %auxiliar variable to set where is the sync-symble
        tsync         = linspace(0,SyncPeriod*T,(2*NumAmosCP+NPPB)*...
            SyncPeriod);%create a time vetor for the sync-word to help on the ploting of this signal
        fsync         = time2freq(tsync);                                  %creating the the frequency vector of the sync-word for hereafter filtering
        BWD           = 2*fc;                                          %Badnwidth for the bit conformation
        CenFeq        = 0;                                                 %Center frequency for the bit conformation
        FiltOrd       = 1;                                                 %The order for the filter which will conform data
        BitFilBanWid  = 10e9;                                              %Badnwidth for the bit conformation
        BitFiltCenFre = 0;                                                 %Center frequency for the bit conformation
        BitFilOrd     = 1;                                                 %The order for the filter which will conform data
        [BitFilt,~]   = FiltroGaussiano(fsync,BWD,CenFeq,FiltOrd);         %Creating filter for conformation of the input information
        BitFilt       = fftshift(BitFilt);
        
        SyncSymb      = rectpulse(SyncSymb,2*NumAmosCP+NPPB);
        SyncSymbEnd   = rectpulse(SyncSymbEnd,2*NumAmosCP+NPPB);
        SencondAdjust = ON;                                                %Setting whether is needed to sync the information by the Syncronism Symble at the End
        LowSymPer     = 1;
        UpeSymPer     = 1;
        if FiberLength>300
            ExtDel = 8;
        else
            ExtDel = 1;
        end
    case '4PAM'
        %%             Syncing 4PAM signal
        SyncSymb      = [zeros(1,10) ones(1,4) zeros(1,10)];               %Generating the sync-word at the beginning of the data stream
        IniSyncPos    = 1+5*(2*NumAmosCP+NPPB);
        SyncPos       = 19*(2*NumAmosCP+NPPB);                             %Seting the possition where the sync-symble is
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
    otherwise
        %%             Syncing OOK signal
        
        SyncSymb      = [zeros(1,10) ones(1,4) zeros(1,10)];               %Generating the sync-word at the beginning of the data stream
        SyncPos       = 19*(2*NumAmosCP+NPPB);                             %Seting the possition where the sync-symble is
        IniSyncPos    = 1+5*(2*NumAmosCP+NPPB);
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
        BitFilBanWid  = 10e9;                                              %Badnwidth for the bit conformation
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
end
%% DatumReceptionInputData;

%%       Variables for Syncronism

%%    Control Variables of this program
BitFilBanWid  = 0.9*Rb;%Badnwidth for the bit conformation
BitFiltCenFre = 0;  %Center frequency for the bit conformation
BitFilOrd     = 3;  %The order for the filter which will conform data
DecMod        = 1;  %Selec how the decission level will be measured
% 1 for measuring the zeros on the eye diagram (best)
% 0 for taking the meddle of the zeros vector

%stream: 1 for downstream and 2 for upstream.
%%    Constants to controu the Optical FFT
% The optical FFT here implemented can take any number of carriers. The
% user must stay realistic because the number of actual devices used to
% implement this setup increases by the power of 2, thus higher orders of
% the FFT may not be feaseble to be implemented as the stability of each
% device may vary with temperature, aging and othter features.
MaxNumStag = nextpow2(NumCarr);                                        %Next number, power of 2, by the given number of carriers
MaxNumCarr = 2^MaxNumStag;                                             %With the number of stages the actual nubmer of carriers can be calculated
% end

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
% switch MaxNumStag
%     case 1
%         VetThisCarr = [1 2];
%     case 2
%         VetThisCarr = [1 3 2 4];
%     case 3
%         VetThisCarr = [1 5 3 7 2 6 4 8];
%     case 4
%         VetThisCarr = [1 9 5 13 3 11 7 15 2 10 6 14 4 12 8 16];
%     case 5
%         VetThisCarr = [1 17 9 25 5 21 13 29 3 19 11 27 7 23 15 31 2 18...
%                                   10 26 6 22 14 30 4 20 12 28 8 24 16 32];
%     otherwise
%         MaxNumStag = 1;
%         MaxNumCarr = 2^MaxNumStag;
%         VetThisCarr = [1 2];
% end
%At this point if for some misteriouse reason the Number of carriers or the
%Maxmum order of the OFFT were seted beyond the maximum point this
%verification will pointout a error while the OpticalOFFT function is not
%optmized for an 'n' order filter.
% if MaxNumStag>5
% error(['MaxNumStag allowed value for MaxNumStag was achieved. Please'...
%                               ' cange the number Maximum of carriers.']);
% end

Count = ones(MaxNumStag,1);                                                %Control Variable for adding correct variation
E0 = 0;                                                                    %Second input, usualy zero
ActStag = 1;                                                               %Starting point of Optical FFT
%%               Variables for syncronism
%Depending of each type of modulation, it is needed to define know symble
%for syncronism. As it deals with optcal devices and optical paths the
%signal takes time to travel, which leads to an sampling miss mach of the
%received and transmite data. Thus at the receiver the syncronism variable
%are created to approach this problem for each different type of
%modulation.
switch Modulation
    case 'OFDM'
        BWD           = 1.0*fc;                                        %Badnwidth for the bit conformation
        CenFeq        = 0;                                                 %Center frequency for the bit conformation
        FiltOrd       = 0;                                                 %The order for the filter which will conform data
    case 'DPSK'
        %%             Syncing DPSK signal
        
        BWD2          = 1*fc;                                          %Badnwidth for the bit conformation
        CenFeq2       = 0;                                                 %Center frequency for the bit conformation
        FiltOrd2      = 1;                                                 %The order for the filter which will conform data
        SyncSymb      = [zeros(1,10) ones(1,4) zeros(1,10)];               %Generating the sync-word at the beginning of the data stream
        IniSyncPos    = 1+5*(2*NumAmosCP+NPPB);
        SyncPos       = 19*(2*NumAmosCP+NPPB);                             %Seting the possition where the sync-symble is
        SyncSymbEnd   = SyncSymb;                                          %Generating the sync-word at the end of the data stream
        SyncPeriod    = length(SyncSymb);                                  %Seting how many symble periods will be within the sync-word
        AuxSyncLength = SyncPeriod*(2*NumAmosCP+NPPB);                     %auxiliar variable to set where is the sync-symble
        tsync         = linspace(0,SyncPeriod*T,(2*NumAmosCP+NPPB)*...
            SyncPeriod);%create a time vetor for the sync-word to help on the ploting of this signal
        fsync         = time2freq(tsync);                                  %creating the the frequency vector of the sync-word for hereafter filtering
        BWD           = 2*fc;                                          %Badnwidth for the bit conformation
        CenFeq        = 0;                                                 %Center frequency for the bit conformation
        FiltOrd       = 1;                                                 %The order for the filter which will conform data
        BitFilBanWid  = 10e9;                                              %Badnwidth for the bit conformation
        BitFiltCenFre = 0;                                                 %Center frequency for the bit conformation
        BitFilOrd     = 1;                                                 %The order for the filter which will conform data
        [BitFilt,~]   = FiltroGaussiano(fsync,BWD,CenFeq,FiltOrd);         %Creating filter for conformation of the input information
        BitFilt       = fftshift(BitFilt);
        
        SyncSymb      = rectpulse(SyncSymb,2*NumAmosCP+NPPB);
        SyncSymbEnd   = rectpulse(SyncSymbEnd,2*NumAmosCP+NPPB);
        SencondAdjust = OFF;                                               %Setting whether is needed to sync the information by the Syncronism Symble at the End
        LowSymPer     = 1;
        UpeSymPer     = 1;
    case 'DQPSK'
        %%             Syncing DQPSK signal
        
        BWD2          = 1.0*fc;                                        %Badnwidth for the bit conformation
        CenFeq2       = 0;                                                 %Center frequency for the bit conformation
        FiltOrd2      = 1;                                                 %The order for the filter which will conform data
        SyncSymb      = [zeros(1,10) ones(1,4) zeros(1,10)];               %Generating the sync-word at the beginning of the data stream
        IniSyncPos    = 1+5*(2*NumAmosCP+NPPB);
        SyncPos       = 19*(2*NumAmosCP+NPPB);                             %Seting the possition where the sync-symble is
        SyncSymbEnd   = SyncSymb;                                          %Generating the sync-word at the end of the data stream
        SyncPeriod    = length(SyncSymb);                                  %Seting how many symble periods will be within the sync-word
        AuxSyncLength = SyncPeriod*(2*NumAmosCP+NPPB);                     %auxiliar variable to set where is the sync-symble
        tsync         = linspace(0,SyncPeriod*T,(2*NumAmosCP+NPPB)*...
            SyncPeriod);%create a time vetor for the sync-word to help on the ploting of this signal
        fsync         = time2freq(tsync);                                  %creating the the frequency vector of the sync-word for hereafter filtering
        BWD           = 2*fc;                                          %Badnwidth for the bit conformation
        CenFeq        = 0;                                                 %Center frequency for the bit conformation
        FiltOrd       = 1;                                                 %The order for the filter which will conform data
        BitFilBanWid  = 10e9;                                              %Badnwidth for the bit conformation
        BitFiltCenFre = 0;                                                 %Center frequency for the bit conformation
        BitFilOrd     = 1;                                                 %The order for the filter which will conform data
        [BitFilt,~]   = FiltroGaussiano(fsync,BWD,CenFeq,FiltOrd);         %Creating filter for conformation of the input information
        BitFilt       = fftshift(BitFilt);
        
        SyncSymb      = rectpulse(SyncSymb,2*NumAmosCP+NPPB);
        SyncSymbEnd   = rectpulse(SyncSymbEnd,2*NumAmosCP+NPPB);
        SencondAdjust = ON;                                                %Setting whether is needed to sync the information by the Syncronism Symble at the End
        LowSymPer     = 1;
        UpeSymPer     = 1;
        if FiberLength>300
            ExtDel = 8;
        else
            ExtDel = 1;
        end
    case '4PAM'
        %%             Syncing 4PAM signal
        SyncSymb      = [zeros(1,10) ones(1,4) zeros(1,10)];               %Generating the sync-word at the beginning of the data stream
        IniSyncPos    = 1+5*(2*NumAmosCP+NPPB);
        SyncPos       = 19*(2*NumAmosCP+NPPB);                             %Seting the possition where the sync-symble is
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
    otherwise
        %%             Syncing OOK signal
        
        SyncSymb      = [zeros(1,10) ones(1,4) zeros(1,10)];               %Generating the sync-word at the beginning of the data stream
        SyncPos       = 19*(2*NumAmosCP+NPPB);                             %Seting the possition where the sync-symble is
        IniSyncPos    = 1+5*(2*NumAmosCP+NPPB);
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
        BitFilBanWid  = 10e9;                                              %Badnwidth for the bit conformation
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
end

%% DatumUpstreamInputData;
%%    Simulation Parameters
NRZPolarity = 1;%Polarity of NRZ coding, 0 for one-polar, 1 for bi-polar
Which4PAM   = PAM4Type(1);%Which 4PAM will be used, electrical "0" or optical "1"
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
MZ_Input_Data;                                                             %Loading the basic data to be saved
%At this point the file with the MZM will be created. The first file will
%have the polarization point at Vpi/2 because it was used to generate the
%UfOCS. The second file is used by the 4PAM and by the OOK modulation,
%because they are polarized at Vpi/2 as well. The third file will be used
%by the DPSK modulation as it is polarized at Vpi. The last one will vary
%in value for testing, therefore it can assume any value.
Vbias_steps = 4;                                                           %Seting the number of variation of the Input file
Vbias = [ U_pi1/2 U_pi1/2 U_pi1 -1*U_pi1/2];                               %Creating the Vbias that will be writen on the Inpu Data file
Local = [pwd '\input_files\'];                                             %Seting the location where this file will be stored.
Make_MZ_Input_Files_Simp;                                                  %Create the files for the MZM input data

%%         Creating file for the DP-MZM
phi_0_steps  = 1;%Number of different phase between arms of the DP-MZM
V1bias_steps = 1;%Number of different Bias voltage at the Up arm
V2bias_steps = 1;%Number of different Bias voltage at the Lower arm
%
%IMPORTANT: The RF amplitude Vpp must respect the device limits!
phi_0_vet     = -2.05;                                                     %Variation of phase betwee arms.
V1bias        = 1.752;                                                     %Variation of bias voltage at the Up arm.
V2bias        = 1.052;                                                     %Variation of bias voltage at the Lower arm.

% Differently from the previouly Input Data file. This scrip does not have
% a previous configuration file for the MZM. Thus, it can create a new file
% to configure the MZM-IQ.
Local_Dp = [pwd '\input_files_dp\'];                                       %Seting the location where this file will be stored.
MZ_Input_Data_DP;                                                          %Loading the basic data to be saved
Make_MZ_Input_Files_DP;                                                    %Create the files for the MZM input data

MZ_Input_File_Dp = 1;                                                      %Telling the script witch file will be used for modulation
%%           Control Variables
%This variable will be used to select the medium in which the information
%will travel Type 'B2B' for a Back-To-Back transmission or 'Fiber' for an
%optical Fiber
% Medium = 'B2B';%'Fiber';%
switch CurrentMedium
    case 1
        Medium = 'Fiber';
    otherwise
        Medium = 'B2B';
end
%This variable will be used to select the type of modulation used for the
%datum Enconding. Type 'DQPSK' for DQPSK encoding, '4PAM' for 4PAM
%encoding. Type 'NoN' if an encoding is not necessary.
switch UpStModula
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

% %%          Variables for the Medium Fiber
% NumCarr      = 8;                                                          %Total number of carriers (DownStream and UpStream) - This variable must be even
% lambdac      = 1550e-9;                                                    %Central wave length of the signal throughout the fiber
% FiberLength  = 40;                                                         %Length of the optical fiber [km]
% %The pilot carrier must be close to the center of the OFDM frame for
% %probing the channel once the time delay varies from carrier to carrier.
% %The carries on the right side fo the pilot carrier will have a right shift
% %time delay whereas the carries on the left side of the pilot carrier will
% %have a left shit time delay, where the reference is the synchronization
% %symbol.
% % PilotCarrier = (NumCarr/2)*fc;
% %The main idea about the time delay compensation is that the receptor will
% %have previous knowledge about the length of the transmission link hence it
% %will posse the time delay of this link. Therefore, before any other
% %processing takes place first the income signal will have the time delay
% %adjusted by this information previously stored.
% if CurrentMedium
%     if ~exist('FiberDelay','var')
%         [FiberDelay,PowRef] = ProgationDelay(NumCarr,t,lambdac,T,...
%                                                   FiberLength,0,f,Ploting);
%     end
% end
%%          Variables for different types of Modulations
switch Modulation
    case 'OFDM'
        %         BW     = fc;
        %         M      = 4;
        %         Tu     = (2^10)/fc;
        %         Tg     = 0;
        %         Ts     = Tu + Tg;
        %         g      = Tu/Ts;
        %         Dtf    = 1/Tu;
        %         NFFT   = 2^nextpow2(ceil(BW/Dtf));
        %         N      = NFFT/2 - 1;
        %         Ns     = ceil(N/2);
        %         k      = log2(M);
        %         nRb    = (BW/(N+2))*((Ns*k)/(1+g));
        %         NumFra = Nb/NFFT;
        %         MZ_Input_File = 2;
        %Variables for the Selection of carriers - Wave Divided.
        %         fin           = RefCarr*fc;                                        %The frequency of the first carrier of the optical COMB
        %         FBWD          = 1.2*Rb;                                            %The filter BandWidth to split each individual carrier
        %         Order         = 1;                                                 %The order of the filter used to split each carrier
    case 'DPSK' %Tested to fiberlength of 80 -> worked
        %%      Variables for the DPSK Datum Modification
        %This variable allow the user to chose to send an random
        %information or just a high pulse within the midle of the datastram
        NbDPSK        = Nb - NumBitDesc;                                   %The total number of bits to be transmited
        DatGai        = 1;                                                 %Seting the maximum range for the eletrical signal
        JusVal        = [zeros(1,10) ones(1,4) zeros(1,10)];               %Value that will be load into the justificaiton space
        JusValEnd     = JusVal;                                            %Value that will be load into the justificaiton space
        JusPos        = length(JusVal);                                    %Length of the justification symbol
        BWD           = 2*fc;                                          %Badnwidth for the bit conformation
        CenFeq        = 0;                                                 %Center frequency for the bit conformation
        FiltOrd       = 1;                                                 %The order for the filter which will conform data
        MZ_Input_File = 3;                                                 %Telling the script witch file will be used for modulation
        %Variables for the Selection of carriers - Wave Divided.
        fin           = RefCarr*fc;                                            %The frequency of the first carrier of the optical COMB
        FBWD          = 1.2*Rb;                                            %The filter BandWidth to split each individual carrier
        Order         = 5;                                                 %The order of the filter used to split each carrier
    case 'DQPSK'%Tested to fiberlength of 40 -> work in progress
        %%      Variables for the DQPSK Datum Modification
        %This variable allow the user to chose to send an random
        %information or just a high pulse within the midle of the datastram
        NbDQPSK     = 2*(Nb - NumBitDesc);                                 %The total number of bits to be transmited
        JusVal      = [zeros(1,20) ones(1,8) zeros(1,20)];                 %Value that will be load into the justificaiton space
        JusValEnd   = JusVal;                                              %Value that will be load into the justificaiton space
        JusPos      = length(JusVal);                                      %Length of the justification symbol
        BWD         = 4*fc;                                            %Badnwidth for the bit conformation
        CenFeq      = 0;                                                   %Center frequency for the bit conformation
        FiltOrd     = 1;                                                   %The order for the filter which will conform data
        DatGai      = 1;                                                   %Gain to increase the modulatin signal amplitude
        V0          = 3.8;
        Vpi         = 3.8;
        %Variables for the Selection of carriers - Wave Divided.
        fin         = RefCarr*fc;                                              %The frequency of the first carrier of the optical COMB
        FBWD        = 1.2*Rb;                                              %The filter BandWidth to split each individual carrier
        Order       = 5;                                                   %The order of the filter used to split each carrier
    case '4PAM'%Tested to fiberlength of 1 -> work in progress
        %%      Variables for the 4PAM Datum Modification
        %This variable allow the user to chose to send an random
        %information or just a high pulse within the midle of the datastram
        Nb4Pam        = 2*(Nb - NumBitDesc);                               %The total number of bits to be transmited
        %This variable controls if the electrical signal will vary from
        %0[V] to max[V] (signal unipolar) or from -max[V] to max[V]
        %(signal bipolar).
        Polirized     = 1;
        MaxAmp4PAM    = -3.8/4;                                            %Seting the maximum range for the eletrical signal
        JusVal        = [zeros(1,20) 1 0 1 0 1 0 1 0 zeros(1,20)];         %Value that will be load into the justificaiton space
        JusValEnd     = JusVal;                                            %Value that will be load into the justificaiton space
        JusPos        = length(JusVal);                                    %Length of the justification symbol
        U_pi2         = 2.5;                                               %Biar voltage of the selected device
        %There are three PAM4 techniques implemented in this simulation,
        %Electical PAM4, Optical PAM4 with IQ-MZM and Optical PAM4 with
        %DD-MZM. The last one is the modulation with the best result for
        %the current setup configuration so far. It is selected when
        %Which4PAM set 1 and ModSchem set 0.
        if Which4PAM
            ModSchem      = PAM4Type(2);
            BWD           = 1.0*fc;                                    %Badnwidth for the bit conformation
            CenFeq        = 0;                                             %Center frequency for the bit conformation
            FiltOrd       = 5;                                             %The order for the filter which will conform data
            Vbias         = 0;                                             %Biar voltage of the selected device
            if ModSchem
                MZ_Input_File = 4;                                         %Telling the script witch file will be used for modulation
                Vmax          = 1;
                Vmin          = -1;
            else
                MZ_Input_File = 5;                                         %Telling the script witch file will be used for modulation
                Vmax          = 1.25;
                Vmin          = 0;
            end
        else
            BWD           = 1.0*fc;                                    %Badnwidth for the bit conformation
            CenFeq        = 0;                                             %Center frequency for the bit conformation
            FiltOrd       = 1;                                             %The order for the filter which will conform data
            VPI           = 3.8;                                           %Characteristic voltage of the MZM - Vbias at VPI is the point of minimum
            Vbias         = 1.9;                                           %Biar voltage of the selected device
            MZ_Input_File = 2;                                             %Telling the script witch file will be used for modulation
        end
        %Variables for the Selection of carriers - Wave Divided.
        fin           = RefCarr*fc;                                            %The frequency of the first carrier of the optical COMB
        FBWD          = 1.2*Rb;                                            %The filter BandWidth to split each individual carrier
        Order         = 5;                                                 %The order of the filter used to split each carrier
    otherwise %Tested to fiberlength of 80 -> worked
        %%      Variables for the OOK Datum Modification
        JusVal        = [zeros(1,10) ones(1,4) zeros(1,10)];               %Value that will be load into the justificaiton space
        JusValEnd     = [zeros(1,10) ones(1,4) zeros(1,10)];               %Value that will be load into the justificaiton space
        JusLen        = length(JusVal);                                    %Length for the information justification
        DatMin        = 0;                                                 %Minimum value that the Data can assume
        DatMax        = 1;                                                 %Maximum value that the DAta can assume
        NrzMin        = 0.96;%3.9;                                         %Minimum value for the NRZ bit
        NrzMax        = -0.96;%3.8;                                        %Maximum value for the NRZ bit
        IQAmp         = 1;                                                 %Amplitude of the signal IQ
        DatGai        = 1;                                                 %Gain to increase the modulatin signal amplitude
        BWD           = fc;                                            %Badnwidth for the bit conformation
        CenFeq        = 0;                                                 %Center frequency for the bit conformation
        FiltOrd       = 1;                                                 %The order for the filter which will conform data
        MZ_Input_File = 4;                                                 %Telling the script witch file will be used for modulation
        %Variables for the Selection of carriers - Wave Divided.
        fin           = RefCarr*fc;                                            %The frequency of the first carrier of the optical COMB
        FBWD          = 1.2*Rb;                                            %The filter BandWidth to split each individual carrier
        Order         = 5;                                                 %The order of the filter used to split each carrier
end
%%      Filter to Select the Carriers
[ SelecFilt ] = FiltroGaussiano(f,18*fc,10*fc,9);
% [ SelecFilt ] = Filtro_Retangular( 34*fc,16*fc,f);
SelecFilt = fftshift(SelecFilt);
%%      Filter for the EVM measurement taking
EvmFilt = FiltroGaussiano(f,fc,0,1);
EvmFilt = fftshift(EvmFilt);
clear TxData TxDataMat EoutAux1 VetThisCarr
%% DatumUpstrRecInputData;                                                    % Configuration File of DatumReception.


%%       General Parameter for the simulation
ON          =  1; %DO NOT CHANGE THIS PARAMETER
OFF         =  0; %DO NOT CHANGE THIS PARAMETER
Ploting     = OFF;%Use On or Off to say if the figures should be plot or not
PlotingThis = ON;%Use On or Off to say if the figures should be plot or not
Selecting   = OFF;%Use to filter or not the carriers to be transmited
AddCP       = ON; %Add the Cycle Prefix on the data stream
%%    Simulation Parameters
NRZPolarity = 1;%Polarity of NRZ coding, 0 for one-polar, 1 for bi-polar
Which4PAM   = PAM4Type(1);%Which 4PAM will be used, electrical "0" or optical "1"
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
MZ_Input_Data;                                                             %Loading the basic data to be saved
%At this point the file with the MZM will be created. The first file will
%have the polarization point at Vpi/2 because it was used to generate the
%UfOCS. The second file is used by the 4PAM and by the OOK modulation,
%because they are polarized at Vpi/2 as well. The third file will be used
%by the DPSK modulation as it is polarized at Vpi. The last one will vary
%in value for testing, therefore it can assume any value.
Vbias_steps = 4;                                                           %Seting the number of variation of the Input file
Vbias = [ U_pi1/2 U_pi1/2 U_pi1 -1*U_pi1/2];                               %Creating the Vbias that will be writen on the Inpu Data file
Local = [pwd '\input_files\'];                                             %Seting the location where this file will be stored.
Make_MZ_Input_Files_Simp;                                                  %Create the files for the MZM input data

%%         Creating file for the DP-MZM
phi_0_steps  = 1;%Number of different phase between arms of the DP-MZM
V1bias_steps = 1;%Number of different Bias voltage at the Up arm
V2bias_steps = 1;%Number of different Bias voltage at the Lower arm
%
%IMPORTANT: The RF amplitude Vpp must respect the device limits!
phi_0_vet     = -2.05;                                                     %Variation of phase betwee arms.
V1bias        = 1.752;                                                     %Variation of bias voltage at the Up arm.
V2bias        = 1.052;                                                     %Variation of bias voltage at the Lower arm.

% Differently from the previouly Input Data file. This scrip does not have
% a previous configuration file for the MZM. Thus, it can create a new file
% to configure the MZM-IQ.
Local_Dp = [pwd '\input_files_dp\'];                                       %Seting the location where this file will be stored.
MZ_Input_Data_DP;                                                          %Loading the basic data to be saved
Make_MZ_Input_Files_DP;                                                    %Create the files for the MZM input data

MZ_Input_File_Dp = 1;                                                      %Telling the script witch file will be used for modulation
%%           Control Variables
%This variable will be used to select the medium in which the information
%will travel Type 'B2B' for a Back-To-Back transmission or 'Fiber' for an
%optical Fiber
% Medium = 'B2B';%'Fiber';%
switch CurrentMedium
    case 1
        Medium = 'Fiber';
    otherwise
        Medium = 'B2B';
end
%This variable will be used to select the type of modulation used for the
%datum Enconding. Type 'DQPSK' for DQPSK encoding, '4PAM' for 4PAM
%encoding. Type 'NoN' if an encoding is not necessary.
switch UpStModula
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

% %%          Variables for the Medium Fiber
% NumCarr      = 8;                                                          %Total number of carriers (DownStream and UpStream) - This variable must be even
% lambdac      = 1550e-9;                                                    %Central wave length of the signal throughout the fiber
% FiberLength  = 40;                                                         %Length of the optical fiber [km]
% %The pilot carrier must be close to the center of the OFDM frame for
% %probing the channel once the time delay varies from carrier to carrier.
% %The carries on the right side fo the pilot carrier will have a right shift
% %time delay whereas the carries on the left side of the pilot carrier will
% %have a left shit time delay, where the reference is the synchronization
% %symbol.
% % PilotCarrier = (NumCarr/2)*fc;
% %The main idea about the time delay compensation is that the receptor will
% %have previous knowledge about the length of the transmission link hence it
% %will posse the time delay of this link. Therefore, before any other
% %processing takes place first the income signal will have the time delay
% %adjusted by this information previously stored.
% if CurrentMedium
%     if ~exist('FiberDelay','var')
%         [FiberDelay,PowRef] = ProgationDelay(NumCarr,t,lambdac,T,...
%                                                   FiberLength,0,f,Ploting);
%     end
% end
%%          Variables for different types of Modulations
switch Modulation
    case 'OFDM'
        %         BW     = fc;
        %         M      = 4;
        %         Tu     = (2^10)/fc;
        %         Tg     = 0;
        %         Ts     = Tu + Tg;
        %         g      = Tu/Ts;
        %         Dtf    = 1/Tu;
        %         NFFT   = 2^nextpow2(ceil(BW/Dtf));
        %         N      = NFFT/2 - 1;
        %         Ns     = ceil(N/2);
        %         k      = log2(M);
        %         nRb    = (BW/(N+2))*((Ns*k)/(1+g));
        %         NumFra = Nb/NFFT;
        %         MZ_Input_File = 2;
        %Variables for the Selection of carriers - Wave Divided.
        %         fin           = RefCarr*fc;                                        %The frequency of the first carrier of the optical COMB
        %         FBWD          = 1.2*Rb;                                            %The filter BandWidth to split each individual carrier
        %         Order         = 1;                                                 %The order of the filter used to split each carrier
    case 'DPSK' %Tested to fiberlength of 80 -> worked
        %%      Variables for the DPSK Datum Modification
        %This variable allow the user to chose to send an random
        %information or just a high pulse within the midle of the datastram
        NbDPSK        = Nb - NumBitDesc;                                   %The total number of bits to be transmited
        DatGai        = 1;                                                 %Seting the maximum range for the eletrical signal
        JusVal        = [zeros(1,10) ones(1,4) zeros(1,10)];               %Value that will be load into the justificaiton space
        JusValEnd     = JusVal;                                            %Value that will be load into the justificaiton space
        JusPos        = length(JusVal);                                    %Length of the justification symbol
        BWD           = 2*fc;                                          %Badnwidth for the bit conformation
        CenFeq        = 0;                                                 %Center frequency for the bit conformation
        FiltOrd       = 1;                                                 %The order for the filter which will conform data
        MZ_Input_File = 3;                                                 %Telling the script witch file will be used for modulation
        %Variables for the Selection of carriers - Wave Divided.
        fin           = RefCarr*fc;                                            %The frequency of the first carrier of the optical COMB
        FBWD          = 1.2*Rb;                                            %The filter BandWidth to split each individual carrier
        Order         = 5;                                                 %The order of the filter used to split each carrier
    case 'DQPSK'%Tested to fiberlength of 40 -> work in progress
        %%      Variables for the DQPSK Datum Modification
        %This variable allow the user to chose to send an random
        %information or just a high pulse within the midle of the datastram
        NbDQPSK     = 2*(Nb - NumBitDesc);                                 %The total number of bits to be transmited
        JusVal      = [zeros(1,20) ones(1,8) zeros(1,20)];                 %Value that will be load into the justificaiton space
        JusValEnd   = JusVal;                                              %Value that will be load into the justificaiton space
        JusPos      = length(JusVal);                                      %Length of the justification symbol
        BWD         = 4*fc;                                            %Badnwidth for the bit conformation
        CenFeq      = 0;                                                   %Center frequency for the bit conformation
        FiltOrd     = 1;                                                   %The order for the filter which will conform data
        DatGai      = 1;                                                   %Gain to increase the modulatin signal amplitude
        V0          = 3.8;
        Vpi         = 3.8;
        %Variables for the Selection of carriers - Wave Divided.
        fin         = RefCarr*fc;                                              %The frequency of the first carrier of the optical COMB
        FBWD        = 1.2*Rb;                                              %The filter BandWidth to split each individual carrier
        Order       = 5;                                                   %The order of the filter used to split each carrier
    case '4PAM'%Tested to fiberlength of 1 -> work in progress
        %%      Variables for the 4PAM Datum Modification
        %This variable allow the user to chose to send an random
        %information or just a high pulse within the midle of the datastram
        Nb4Pam        = 2*(Nb - NumBitDesc);                               %The total number of bits to be transmited
        %This variable controls if the electrical signal will vary from
        %0[V] to max[V] (signal unipolar) or from -max[V] to max[V]
        %(signal bipolar).
        Polirized     = 1;
        MaxAmp4PAM    = -3.8/4;                                            %Seting the maximum range for the eletrical signal
        JusVal        = [zeros(1,20) 1 0 1 0 1 0 1 0 zeros(1,20)];         %Value that will be load into the justificaiton space
        JusValEnd     = JusVal;                                            %Value that will be load into the justificaiton space
        JusPos        = length(JusVal);                                    %Length of the justification symbol
        U_pi2         = 2.5;                                               %Biar voltage of the selected device
        %There are three PAM4 techniques implemented in this simulation,
        %Electical PAM4, Optical PAM4 with IQ-MZM and Optical PAM4 with
        %DD-MZM. The last one is the modulation with the best result for
        %the current setup configuration so far. It is selected when
        %Which4PAM set 1 and ModSchem set 0.
        if Which4PAM
            ModSchem      = PAM4Type(2);
            BWD           = 1.0*fc;                                    %Badnwidth for the bit conformation
            CenFeq        = 0;                                             %Center frequency for the bit conformation
            FiltOrd       = 5;                                             %The order for the filter which will conform data
            Vbias         = 0;                                             %Biar voltage of the selected device
            if ModSchem
                MZ_Input_File = 4;                                         %Telling the script witch file will be used for modulation
                Vmax          = 1;
                Vmin          = -1;
            else
                MZ_Input_File = 5;                                         %Telling the script witch file will be used for modulation
                Vmax          = 1.25;
                Vmin          = 0;
            end
        else
            BWD           = 1.0*fc;                                    %Badnwidth for the bit conformation
            CenFeq        = 0;                                             %Center frequency for the bit conformation
            FiltOrd       = 1;                                             %The order for the filter which will conform data
            VPI           = 3.8;                                           %Characteristic voltage of the MZM - Vbias at VPI is the point of minimum
            Vbias         = 1.9;                                           %Biar voltage of the selected device
            MZ_Input_File = 2;                                             %Telling the script witch file will be used for modulation
        end
        %Variables for the Selection of carriers - Wave Divided.
        fin           = RefCarr*fc;                                            %The frequency of the first carrier of the optical COMB
        FBWD          = 1.2*Rb;                                            %The filter BandWidth to split each individual carrier
        Order         = 5;                                                 %The order of the filter used to split each carrier
    otherwise %Tested to fiberlength of 80 -> worked
        %%      Variables for the OOK Datum Modification
        JusVal        = [zeros(1,10) ones(1,4) zeros(1,10)];               %Value that will be load into the justificaiton space
        JusValEnd     = [zeros(1,10) ones(1,4) zeros(1,10)];               %Value that will be load into the justificaiton space
        JusLen        = length(JusVal);                                    %Length for the information justification
        DatMin        = 0;                                                 %Minimum value that the Data can assume
        DatMax        = 1;                                                 %Maximum value that the DAta can assume
        NrzMin        = 0.96;%3.9;                                         %Minimum value for the NRZ bit
        NrzMax        = -0.96;%3.8;                                        %Maximum value for the NRZ bit
        IQAmp         = 1;                                                 %Amplitude of the signal IQ
        DatGai        = 1;                                                 %Gain to increase the modulatin signal amplitude
        BWD           = fc;                                            %Badnwidth for the bit conformation
        CenFeq        = 0;                                                 %Center frequency for the bit conformation
        FiltOrd       = 1;                                                 %The order for the filter which will conform data
        MZ_Input_File = 4;                                                 %Telling the script witch file will be used for modulation
        %Variables for the Selection of carriers - Wave Divided.
        fin           = RefCarr*fc;                                            %The frequency of the first carrier of the optical COMB
        FBWD          = 1.2*Rb;                                            %The filter BandWidth to split each individual carrier
        Order         = 5;                                                 %The order of the filter used to split each carrier
end
%%      Filter to Select the Carriers
[ SelecFilt ] = FiltroGaussiano(f,18*fc,10*fc,9);
% [ SelecFilt ] = Filtro_Retangular( 34*fc,16*fc,f);
SelecFilt = fftshift(SelecFilt);
%%      Filter for the EVM measurement taking
EvmFilt = FiltroGaussiano(f,fc,0,1);
EvmFilt = fftshift(EvmFilt);




VetSnrIni     = 20;
VetSnrPass    = -10;
VerSnrEnd     = 20;
ThisPlotCont  = 1;
clear MzFilesGenerated MzFilesDpGenerated
tic;                                                                       %Starting timer to check the simulation time
for VetSnr=VetSnrIni:VetSnrPass:VerSnrEnd
    clear BerOOK BerPAM4 BerDQPSK BerDPSK BerOOKS BerPAM4S BerDQPSKS ...
        BerDPSKS FiberDelay BerOFDM
    
    
    SaveRxNotEq=[];
    SaveRxEq=[];
    %%              Loading the Ultra OFC Source
    % Initially this scrip will load the Ultra Optical Flat Comb Source, which
    % will be the carriers equaly spaced of our sistem. It is important to
    % remember that each carrier was originated from an optical ring in which a
    % single laser CW (Continuos Wave) was used. The main point of this scrip
    % is to demonstrate how it is possible to implement an All Optical FFT.
    % Therefore, the UOFCS needs to be precisely spaced according to the bit
    % ratio Rs.That is directly related to the Symbol ration of the sistem.
    
    %The two following loops are to control the test bundle. The first loop is
    %responsible to select the type of modulation that will be under test and
    %the second loop will determinate how many tests will be running.
    % load('DqpskCp122FullSysNoNoiseKm_100_1');
    for CurrentModula=ThisModula
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
                
                TxSig = 0;
                for CarDmtM=1:length(DmtMve)
                    M = DmtMve(CarDmtM);
                    TxData          = randi([0 M-1],NumFra,1);                %Generation random information
                    TxDataMat(1,1+NumFra*(CarDmtM-1):NumFra*CarDmtM) = TxData(:);                               %Saving data for future comparison
                    switch OfdMod                                              %Sellect which modulation was used
                        case 'qam'
                            TxSymbAux(1:NumFra,CarDmtM)  = qammod(TxData,M);                    %Modulating information by QAM
                        otherwise
                            TxSymbAux(1:NumFra,CarDmtM)  = dpskmod(TxData,M);                   %Modulating information DPSK as default
                    end
                end
                TxSymb      = TxSymbAux.';                                %Converting from serial to paralel
                EvmMatRef(1,:) = TxSymb(:);
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
                    TxSymbH = [zeros(1,NumFra);TxSymb;zeros(1,NumFra);TxSymbConj];%Hermitian Symetri
                    TxSymbC = zeros(NFFT,NumFra);
                    TxSymbC(1+Zt-ZtC:Zt-ZtC+size(TxSymbH,1)/2,:) = TxSymbH(1:size(TxSymbH,1)/2,:);
                    TxSymbC(end-Zt+ZtC-(size(TxSymbH,1)/2)+1:end-Zt+ZtC,:) = TxSymbH(1+size(TxSymbH,1)/2:end,:);
                else
                    TxSymbH = TxSymb;
                    %     TxSymbC = [zeros((NFFT-size(TxSymbH,1))/2,NumFra);TxSymbH;zeros((NFFT-size(TxSymbH,1))/2,NumFra)];
                    TxSymbC = [TxSymbH(1:size(TxSymbH,1)/2,:);zeros((NFFT-size(TxSymbH,1)),NumFra);TxSymbH(end+1-size(TxSymbH,1)/2:end,:)];
                end
                %Here was implemented the modulation through FFT process, it basica mix-up
                %all infomation following Fourier transform.
                TxSymbModA       = ifft(TxSymbC,NFFT);
                %Sellecting whether the OFDM signal will be transmited on
                %Base band or not
                switch SelModTp
                    case 'BP'                                              %Base Band transmission
                        TxSymbModA  = rectpulse(TxSymbModA,OvSam);         %Over sampling
                        %                         TxSymbMod = TxSymbModA;
                        tt = linspace(0,1*Te,size(TxSymbModA,1));
                        for jj=1:NumFra
                            tta(:,jj) = tt;
                        end
                        TxSymbMod   = TxSymbModA;
                        ff = time2freq(tt);
                        if SelecGaus
                            [ ModFilt ] = FiltroGaussiano(ff,OBw/2,0,1);
                        else
                            [ ModFilt ] = Filtro_Retangular(OBw,0,ff);
                        end
                        ModFilt = fftshift(ModFilt);
                        if  (Ploting)
                            figure;hold all;
                            plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                        end
                        for jj=1:NumFra
                            ModFilta(:,jj) = ModFilt;
                        end
                        if SelModFilt
                            TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
                        end
                        if  (Ploting)
                            %                             for jj=1:NumFra
                            %                                 ff2(:,jj) = ff;
                            %                             end
                            plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                            plot(ff,20*log10(abs(fftshift(ModFilt))));
                            set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                        end
                    case 'AM'                                              %Out of base band
                        TxSymbModA  = rectpulse(TxSymbModA,OvSam);         %Over sampling
                        tt = linspace(0,1*Te,size(TxSymbModA,1));
                        ff = time2freq(tt);
                        for jj=1:NumFra
                            tta(:,jj) = tt;
                        end
                        TxSymbMod   = TxSymbModA.*cos(2*pi*Ofc*tta);
                        if SelecGaus
                            [ ModFilt ] = FiltroGaussiano(ff,OBw,Ofc,1);
                        else
                            [ ModFilt ] = Filtro_Retangular( OBw,Ofc,ff);
                        end
                        ModFilt = fftshift(ModFilt);
                        if  (Ploting)
                            figure;hold all;
                            plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                        end
                        for jj=1:NumFra
                            ModFilta(:,jj) = ModFilt;
                        end
                        if SelModFilt
                            TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
                        end
                        if  (Ploting)
                            plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                            plot(ff,20*log10(abs(fftshift(ModFilt))));
                            set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                        end
                    case 'AMSSB'                                           %Out of base band
                        TxSymbModA  = rectpulse(TxSymbModA,OvSam);         %Over sampling
                        tt = linspace(0,1*Te,size(TxSymbModA,1));
                        ff = time2freq(tt);
                        for jj=1:NumFra
                            tta(:,jj) = tt;
                        end
                        TxSymbMod   = TxSymbModA.*exp(1j*2*pi*Ofc*tta);
                        if SelecGaus
                            [ ModFilt ] = FiltroGaussiano(ff,OBw,Ofc,1);
                        else
                            [ ModFilt ] = Filtro_Retangular( OBw,Ofc,ff);
                        end
                        ModFilt = fftshift(ModFilt);
                        if  (Ploting)
                            figure;hold all;
                            plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                        end
                        for jj=1:NumFra
                            ModFilta(:,jj) = ModFilt;
                        end
                        if SelModFilt
                            TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
                        end
                        if  (Ploting)
                            plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                            plot(ff,20*log10(abs(fftshift(ModFilt))));
                            set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                            a=1;
                        end
                    otherwise
                        TxSymbModA  = rectpulse(TxSymbModA,OvSam);         %Over sampling
                        [TxSymbMod,tt] = modulate(TxSymbModA,Ofc,Ofs,'pm');
                        ff = time2freq(tt);
                        if SelecGaus
                            [ ModFilt ] = FiltroGaussiano(ff,3*Ofc,0,1);
                        else
                            [ ModFilt ] = Filtro_Retangular( 6*Ofc,0,ff);
                        end
                        ModFilt = fftshift(ModFilt);
                        if  (Ploting)
                            figure;hold all;
                            plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                        end
                        for jj=1:NumFra
                            ModFilta(:,jj) = ModFilt;
                        end
                        TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
                        if  (Ploting)
                            plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                            plot(ff,20*log10(abs(fftshift(ModFilt))));
                            set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                        end
                end
                %
                %                 figure;plot(TxSymbMod(:,1));
                
                AuxTxSymbMod(1,:)  = TxSymbMod(:);
                TxSymbModA  = [TxSymbMod(1+end-NPOFEX:end,:);TxSymbMod];   %Adding Ciclyc prefix
                SigTx  = rectpulse(TxSymbModA,NPPOF);                      %Over sampling
                
                TxSigA = SigTx(:).';                                       %Serializing signal
                TxSig = TxSig + [TxSigA(1:NPSTUf) TxSigA];                         %Conforming vector sizes to the same length
                % CpPaFr = CpPaFr + 1;
                %Assigning the eletrical signal to one drive of the MZM - The aplitude of
                %the signal can be controlled by the variable DatGai, which can be
                %understood as an gain for the eletrical signal or an atenuation. The
                %second signal will be similar with the only difference a phase shift of pi
                NormFact = max(TxSig);
                TxSig = (TxSig./NormFact);                     %Normalizing the sinal to mach with the MZM espect to receive.
                
                % if  (PlotingThisThis)
                %     figure;hold all;grid on;
                %     plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
                %     axis([-25e9 37.5e9 -200 0]);
                %     set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                % end
                
                
                if ChanAwgn
                    [~,SigPower] = MeasPower(TxSig);
                    SigPower2 = SigPower-30;%20*log10(SigPowerI);
                    %     SNR = CarSNR + 10*log10(log2(max(DmtMve))) - 10*log10(1*Ts*OBw);
                    SNR = CarSNR + 10*log10(log2(max(DmtMve))) - 10*log10(OvSam);
                    %     SNR
                    EoutAux = awgn(TxSig,SNR,SigPower2);
                    %     EoutAux = TxSig;
                else
                    EoutAux = TxSig;
                end
                
                
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
                % ControlVars.SelecFilt     = SelecFilt;
                ControlVars.DifNsN        = Medium;
                ControlVars.DifNsN        = DifNsN;
                ControlVars.fin           = fin;
                ControlVars.FBWD          = FBWD;
                ControlVars.Order         = Order;
                ControlVars.FibLen        = FiberLength;
                ControlVars.Tu            = Tu;
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
                        if ~exist('ChanelEqualizer','var')
                            ChanelEqualizer = 0;
                            for kk=1:TapN
                                if CurrentMedium
                                    ChanelEqualizer = (ChanelEqualizer + FindChanelElectricalResponse(ControlVars,FiberDelay));
                                else
                                    ChanelEqualizer = (ChanelEqualizer + FindChanelElectricalResponse(ControlVars));
                                end
                            end
                            ChanelEqualizer = (ChanelEqualizer)./kk;
                        end
                    otherwise
                end
            else
                %% DatumTransmission;
                
                % The first step here is to configure all the parameters that will used
                % hereafter by the transmiter. Hence, implementing functionalities such as
                % generation of random data to fill the pay load, create the syncronism
                % symble, modulating the incoming data acordingly, transmiting this
                % information throughout a medium (fiber or B2B). All those functionalities
                % will need variables to be first initialyzed in this following Input Data.
                %% DatumTransmissionInputData;
                
                
                
                % The electrical OFDM when modulated with QAM the signal will need a
                % equalizer, because the system implyes on phase rotation. Therefore, this
                % system needs correction before comparing the TX and RX signals. It is
                % also important to mention that the DPSK modulation format don't need
                % equalizer as the information is stored at the phase difference on
                % subsequent symbols phase rotaion will affect adjacent symbol in a very
                % similar way hence information will have very little sencibility to phase
                % rotation.
                % ChanelEqualizerP = FindChanelResponseP(OfdMod,M,NFFT,Ns,NumFra,MZ_Input_File,SelecFilt,Medium,DifNsN,fin,FBWD,Order,ExtSam);
                if exist('OfdMod','var')
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
                            if ~exist('ChanelEqualizer','var')
                                ChanelEqualizer = 0;
                                for kk=1:TapN
                                    if CurrentMedium
                                        ChanelEqualizer = (ChanelEqualizer + FindChanelResponse(ControlVars,FiberDelay));
                                    else
                                        ChanelEqualizer = (ChanelEqualizer + FindChanelResponse(ControlVars));
                                    end
                                end
                                ChanelEqualizer = (ChanelEqualizer)./kk;
                            end
                        otherwise
                    end
                end
                
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
                
                if SendingDowStr
                    
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
                        if SelecSetUP                                                      %Vefify is the split process will be by filter...
                            EoutTx=SelectEachCarrier(Eout,NumCarr,f,fin,FBWD,Order,fc);    %This function is resposible to split each carrier
                            VetThisCarrTx = (RefCarr-1)+1:(RefCarr-1)+NumCarr;             %Keeping track of the right carrier sequence
                        else                                                               %... or by the OFFT
                            %As the OFFT has a periodic response the OFC needs to be constrained other
                            %whise carrier multiple carrier may interfir with other channels, This
                            %first selection was done with a broad pass-band filter. A higher OFFT
                            %order can also be used, although it may increase the computational time
                            %and my not be exactly feasible in the real world.
                            if Selecting
                                EoutT = ifft(fft(Eout).*SelecFilt);
                                %                 PrintInfo(Ploting*1,f,EoutT);
                                %                 a=0;
                            else
                                EoutT = Eout;
                            end
                            [EoutTx,~,VetThisCarrTx]=OpticalFFTN(t,T,MaxNumStagTx,EoutT);  %This function is resposible to split each carrier
                        end
                        %         PrintInfo(Ploting*51,EoutTx,f);                              %Printing for qualitative analizes.
                        %         axis([min(f) max(f) -400 0]);
                        %         a=0;
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
                            for kk=InitCarrDo:CarrPass:NumCarr                             %Generating different data for each carrier
                                TxSig = 0;
                                %                 CpPaFr = 1;
                                %for CarrOffSet=-1*(NumFraPar-1)/2:1:(NumFraPar-1)/2
                                for CarDmtM=1:length(DmtMve)
                                    M = DmtMve(CarDmtM);
                                    TxData          = randi([0 M-1],NumFra,1);                %Generation random information
                                    TxDataMat(kk,1+NumFra*(CarDmtM-1):NumFra*CarDmtM) = TxData(:);                               %Saving data for future comparison
                                    switch OfdMod                                              %Sellect which modulation was used
                                        case 'qam'
                                            TxSymbAux(1:NumFra,CarDmtM)  = qammod(TxData,M);                    %Modulating information by QAM
                                        otherwise
                                            TxSymbAux(1:NumFra,CarDmtM)  = dpskmod(TxData,M);                   %Modulating information DPSK as default
                                    end
                                end
                                TxSymb      = TxSymbAux.';                                %Converting from serial to paralel
                                TxSigOfdm(CurrentTest,ObsCarrPos==kk,:) = TxSymb(:);
                                EvmMatRef(ObsCarrPos==kk,:) = TxSymb(:);                   %Taking just the middle samples as references
                                %Construction the Hermitian Simetry
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
                                    TxSymbH = [zeros(1,NumFra);TxSymb;zeros(1,NumFra);TxSymbConj];%Hermitian Symetri
                                    TxSymbC = zeros(NFFT,NumFra);
                                    TxSymbC(1+Zt-ZtC:Zt-ZtC+size(TxSymbH,1)/2,:) = TxSymbH(1:size(TxSymbH,1)/2,:);
                                    TxSymbC(end-Zt+ZtC-(size(TxSymbH,1)/2)+1:end-Zt+ZtC,:) = TxSymbH(1+size(TxSymbH,1)/2:end,:);
                                else
                                    TxSymbH = TxSymb;
                                    %     TxSymbC = [zeros((NFFT-size(TxSymbH,1))/2,NumFra);TxSymbH;zeros((NFFT-size(TxSymbH,1))/2,NumFra)];
                                    TxSymbC = [TxSymbH(1:size(TxSymbH,1)/2,:);zeros((NFFT-size(TxSymbH,1)),NumFra);TxSymbH(end+1-size(TxSymbH,1)/2:end,:)];
                                end
                                %Here was implemented the modulation through FFT process, it basica mix-up
                                %all infomation following Fourier transform.
                                TxSymbModA       = ifft(TxSymbC,NFFT);
                                %Sellecting whether the OFDM signal will be transmited on
                                %Base band or not
                                switch SelModTp
                                    case 'BP'                                              %Base Band transmission
                                        TxSymbModA  = rectpulse(TxSymbModA,OvSam);         %Over sampling
                                        %                         TxSymbMod = TxSymbModA;
                                        tt = linspace(0,1*Te,size(TxSymbModA,1));
                                        for jj=1:NumFra
                                            tta(:,jj) = tt;
                                        end
                                        TxSymbMod   = TxSymbModA;%.*exp(1j*2*pi*CarrOffSet*FramSpac*tta);
                                        ff = time2freq(tt);
                                        if SelecGaus
                                            [ ModFilt ] = FiltroGaussiano(ff,OBw/2,0,1);%CarrOffSet*FramSpac,1);
                                        else
                                            [ ModFilt ] = Filtro_Retangular(OBw,0,ff);%CarrOffSet*FramSpac,ff);
                                        end
                                        ModFilt = fftshift(ModFilt);
                                        if (kk==1)&&(Ploting)
                                            figure;hold all;
                                            plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                                        end
                                        for jj=1:NumFra
                                            ModFilta(:,jj) = ModFilt;
                                        end
                                        if SelModFilt
                                            TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
                                        end
                                        if (kk==1)&&(Ploting)
                                            %                             for jj=1:NumFra
                                            %                                 ff2(:,jj) = ff;
                                            %                             end
                                            plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                                            plot(ff,20*log10(abs(fftshift(ModFilt))));
                                            set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                        end
                                    case 'AM'                                              %Out of base band
                                        TxSymbModA  = rectpulse(TxSymbModA,OvSam);         %Over sampling
                                        tt = linspace(0,1*Te,size(TxSymbModA,1));
                                        ff = time2freq(tt);
                                        for jj=1:NumFra
                                            tta(:,jj) = tt;
                                        end
                                        TxSymbMod   = TxSymbModA.*cos(2*pi*Ofc*tta);
                                        if SelecGaus
                                            [ ModFilt ] = FiltroGaussiano(ff,OBw,Ofc,1);
                                        else
                                            [ ModFilt ] = Filtro_Retangular( OBw,Ofc,ff);
                                        end
                                        ModFilt = fftshift(ModFilt);
                                        if (kk==1)&&(Ploting)
                                            figure;hold all;
                                            plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                                        end
                                        for jj=1:NumFra
                                            ModFilta(:,jj) = ModFilt;
                                        end
                                        if SelModFilt
                                            TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
                                        end
                                        if (kk==1)&&(Ploting)
                                            plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                                            plot(ff,20*log10(abs(fftshift(ModFilt))));
                                            set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                        end
                                    case 'AMSSB'                                           %Out of base band
                                        TxSymbModA  = rectpulse(TxSymbModA,OvSam);         %Over sampling
                                        tt = linspace(0,1*Te,size(TxSymbModA,1));
                                        ff = time2freq(tt);
                                        for jj=1:NumFra
                                            tta(:,jj) = tt;
                                        end
                                        TxSymbMod   = TxSymbModA.*exp(1j*2*pi*Ofc*tta);
                                        if SelecGaus
                                            [ ModFilt ] = FiltroGaussiano(ff,OBw,Ofc,1);
                                        else
                                            [ ModFilt ] = Filtro_Retangular( OBw,Ofc,ff);
                                        end
                                        ModFilt = fftshift(ModFilt);
                                        if (kk==1)&&(Ploting)
                                            figure;hold all;
                                            plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                                        end
                                        for jj=1:NumFra
                                            ModFilta(:,jj) = ModFilt;
                                        end
                                        if SelModFilt
                                            TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
                                        end
                                        if (kk==1)&&(Ploting)
                                            plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                                            plot(ff,20*log10(abs(fftshift(ModFilt))));
                                            set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                            a=1;
                                        end
                                    otherwise
                                        TxSymbModA  = rectpulse(TxSymbModA,OvSam);         %Over sampling
                                        [TxSymbMod,tt] = modulate(TxSymbModA,Ofc,Ofs,'pm');
                                        ff = time2freq(tt);
                                        if SelecGaus
                                            [ ModFilt ] = FiltroGaussiano(ff,3*Ofc,0,1);
                                        else
                                            [ ModFilt ] = Filtro_Retangular( 6*Ofc,0,ff);
                                        end
                                        ModFilt = fftshift(ModFilt);
                                        if (kk==1)&&(Ploting)
                                            figure;hold all;
                                            plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                                        end
                                        for jj=1:NumFra
                                            ModFilta(:,jj) = ModFilt;
                                        end
                                        TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
                                        if (kk==1)&&(Ploting)
                                            plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                                            plot(ff,20*log10(abs(fftshift(ModFilt))));
                                            set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                        end
                                end
                                %
                                %                 figure;plot(TxSymbMod(:,1));
                                
                                AuxTxSymbMod(1,:)  = TxSymbMod(:);
                                TxSymbModA  = [TxSymbMod(1+end-NPOFEX:end,:);TxSymbMod];   %Adding Ciclyc prefix
                                SigTx  = rectpulse(TxSymbModA,NPPOF);                      %Over sampling
                                
                                TxSigA = SigTx(:).';                                       %Serializing signal
                                TxSig = TxSig + [TxSigA(1:NPSTUf) TxSigA];                         %Conforming vector sizes to the same length
                                %CpPaFr = CpPaFr + 1;
                                %end
                                %Assigning the eletrical signal to one drive of the MZM - The aplitude of
                                %the signal can be controlled by the variable DatGai, which can be
                                %understood as an gain for the eletrical signal or an atenuation. The
                                %second signal will be similar with the only difference a phase shift of pi
                                NormFact = max((TxSig));
                                TxSig = 0.95*(TxSig./NormFact);                     %Normalizing the sinal to mach with the MZM espect to receive.
                                U.U1t = TxSig;
                                U.U2t = exp(-1j*pi).*TxSig;
                                %As both signals will have the mostrly the same characteristics with the
                                %only difference the phase shift of 180 degress. The MZM-I will be working
                                %on the Push-Pull configuration. It is necessary to reduce the Chirp noise
                                %to zero.
                                [EoutModAux,~]=Mach_Zehnder_Modulator_simplificado_2(freqGHz,...
                                    EoutTx(VetThisCarrTx==(RefCarr-1+kk),:),U,MZ_Input_File);%Modulating individual carriers
                                EoutModTem(ObsCarrPos(kk),1:length(EoutModAux)) = ...
                                    EoutModAux;%Adding the current Modulated output to the final OutPut
                                if (kk==1)&&(Ploting)
                                    figure;hold all;grid on;
                                    plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
                                    axis([-25e9 37.5e9 -200 0]);
                                    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                end
                                %This section adds the unmodulated upstream carrier on the optical OFDM
                                %signal
                                if (~mod(CarrPass,2))&&(SendingUp)
                                    clear EoutModAux;
                                    EoutModAux = EoutTx(VetThisCarrTx==(RefCarr-1+kk...
                                        +1),:);
                                    if AddingNoiseO                                            %Adding optical noise, if needed
                                        EsN0 = 10*log10(snr);
                                        [~,sigpower] = MeasPower(EoutModAux);
                                        EoutModAux = awgn(EoutModAux,osnr,sigpower-30);
                                        [~,sigpower] = MeasPower(EoutModAux);
                                    end
                                    EoutModTem(ObsCarrPos(kk+1),1:length(EoutModAux)...
                                        ) = EoutModAux;%Adding the current Unmodulated output to the final OutPut
                                end
                            end
                            %At this section individual and parralel carriers will be assemble to
                            %compose the optical OFDM signal. It can be done in two way, with a simple
                            %adder, which gives a small ICI clearance. Or by the optical IFFT, which
                            %gives a better ICI clearance.
                            if size(EoutModTem,1)<NumCarr
                                EoutModTem(end+1,:) = 0;
                            end
                            if IfftOrSum                                                   %Use the OIFFT
                                if size(EoutModTem,1)>1                                    %Just if the signal has more than 1 dimension
                                    [EoutMod,~,~] = OpticalIFFTN(t,T,MaxNumStagT,EoutModTem);
                                else
                                    EoutMod = EoutModTem;                                  %otherwise just let the signal pass
                                end
                            else                                                           %Use a simple adder which don't display a better result
                                if size(EoutModTem,1)>1
                                    EoutMod = sum(EoutModTem);
                                else
                                    EoutMod = EoutModTem;
                                end
                            end
                        case 'DPSK'
                            %%            Generate the data DPSK
                            for kk=InitCarrDo:CarrPass:NumCarr                             %Generating different data for each carrier
                                if CarrUsedDo(kk)
                                    if SendingData                                         % Frist it is chosen to transmite a random information...
                                        TxData = (randi(2,1,NbDPSK)-1);                    %Creating the data stream to be transmited
                                        TxData(1:JusPos) = JusVal;                         %Adding the Justification simble at the begining of the stream to sincronize received data frame
                                        TxData(end-(JusPos-1):end) = JusValEnd;            %Adding the Justification simble at the end of the stream to sincronize received data frame
                                    else                                                   %... or just one high pulse for testing the channel
                                        TxData        = zeros(1,NbDPSK);                   %Creating the base fo the data stream
                                        TxData(end/2) = 1;                                 %Adding the code in the data stream to produce the highest level possible
                                    end
                                    TxDataMat(kk,:)   = TxData;                            %Storring the transmited information for latter evaluation
                                    %Once that we have the data to be transmited it needs to be converted to
                                    %the eletrical signal that hereafter will modulate our carrier through the
                                    %MZM-I
                                    SigTx           = DpskEncodEq( TxData);
                                    SigTx(SigTx==1) =  1.9;
                                    SigTx(SigTx==0) = -1.9;
                                    SigTx           = rectpulse(SigTx,NPPB);
                                    % Adding CP to the data
                                    if AddCP
                                        TxAux = reshape(SigTx,NPPB,NbDPSK);
                                        if SetCpSampZer
                                            TxAux = [zeros(NumAmosCP,size(TxAux,2));...
                                                TxAux;zeros(NumAmosCP,size(TxAux,2))];
                                        else
                                            TxAux = [TxAux(1:NumAmosCP,:);TxAux;TxAux(...
                                                end-(NumAmosCP-1):end,:)];
                                        end
                                        TxAux = reshape(TxAux,1,(2*NumAmosCP+NPPB)*NbDPSK);
                                        SigTx = [TxAux TxAux(end-(StuffSampels-1):end)];
                                    end
                                    
                                    %Because of reasons, the real world does not have components capable
                                    %instantly change the voltage level. Thus, to emulate this behaviour the
                                    %eletrical PAM signal pass through a gaussian filter that will remove the
                                    %frequencies of higher order, which will result in small slope in the level
                                    %variation
                                    
                                    [BitFilt,~] = FiltroGaussiano(f,BWD,CenFeq,FiltOrd);   %Creating filter for conformation of the input information
                                    BitFilt     = fftshift(BitFilt);                       %Doing a shift on the vector for matching the transmited data
                                    %TxSig      = ifft(fft(SigTx).*BitFilt);               %Conforming the information and Creating the modulation signal
                                    TxSig       = SigTx;                                   %Adding a gain to the eletrical signal
                                    
                                    %                     if (kk==1)&&(Ploting)%Conforming the information and Creating the modulation signal
                                    %                     figure;hold all;grid on;
                                    %                     plot(t,TxSig)
                                    %                     end
                                    if AddingNoiseE
                                        EsN0 = 10*log10(snr);
                                        [~,sigpower] = MeasPower(TxSig);
                                        TxSig = awgn(TxSig,snr,sigpower-30);
                                        [~,sigpower] = MeasPower(TxSig);
                                    end
                                    %                     if (kk==1)&&(Ploting)
                                    %                       plot(t,TxSig)
                                    %                       PrintInfo(Ploting*2,TxSig(kk,1:end-...
                                    %                       StuffSampels),t(2*NumAmosCP+NPPB),2*NumAmosCP+...
                                    %                                                                    NPPB);
                                    %                       a=4;
                                    %                     end
                                    
                                    %Assigning the eletrical signal to one drive of the MZM - The aplitude of
                                    %the signal can be controlled by the variable DatGai, which can be
                                    %understood as an gain for the eletrical signal or an atenuation. The
                                    %second signal will be similar with the only difference a phase shift of pi
                                    U.U1t = TxSig;
                                    U.U2t = exp(-1j*pi).*TxSig;
                                    %As both signals will have the mostrly the same characteristics with the
                                    %only difference the phase shift of 180 degress. The MZM-I will be working
                                    %on the Push-Pull configuration. It is necessary to reduce the Chirp noise
                                    %to zero.
                                    [EoutModAux,~]=Mach_Zehnder_Modulator_simplificado_2(freqGHz...
                                        ,EoutTx(VetThisCarrTx==(RefCarr-1+kk),:),U,...
                                        MZ_Input_File);%Modulating individual carriers
                                    
                                    %                     if (kk==1)&&(Ploting)
                                    %                         figure;hold all;grid on;
                                    %                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
                                    %                     end
                                    if AddingNoiseO
                                        EsN0 = 10*log10(snr);
                                        [~,sigpower] = MeasPower(EoutModAux);
                                        EoutModAux = awgn(EoutModAux,osnr,sigpower-30);
                                        [~,sigpower] = MeasPower(EoutModAux);
                                    end
                                    %                     if (kk==1)&&(Ploting)
                                    %                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
                                    %                         PrintInfo(Ploting*4,abs(EoutModAux(kk,1:end-...
                                    %                         StuffSampels)).^2./max(abs(EoutModAux(kk,1:end-...
                                    %                         StuffSampels)).^2),t(2*NumAmosCP+NPPB),2*...
                                    %                         NumAmosCP+NPPB,abs(EoutModAux(1:end-StuffSampels...
                                    %                         )).^2./max(abs(EoutModAux(kk,1:end-StuffSampels)...
                                    %                                                                     ).^2));
                                    %                         a=4;
                                    %                     end
                                    EoutModTem(ObsCarrPos(kk),1:length(EoutModAux)) = ...
                                        EoutModAux;%Adding the current Modulated output to the final OutPut
                                    %% Taking the sampling the EVM meassurement
                                    PhaDel          = 0;
                                    TimDel          = T;
                                    [ESync1,ESync2] = DelayInterf(t,TimDel,PhaDel,EoutModAux); %Coverting phase shift to amplitude variation
                                    ESync1 = ESync1.*conj(ESync1);
                                    ESync2 = ESync2.*conj(ESync2);
                                    Esync  = ESync2-ESync1;
                                    PosAuxEout      = (2*NumAmosCP+NPPB)/2:(2*NumAmosCP+...
                                        NPPB):(length(Esync)-StuffSampels);%Varriable respossible to take just the samples at the middle of the symbol
                                    Esync           = ifft(fft(Esync).*EvmFilt);           %Removing higher signal generated by the receiving process
                                    %                     Esync  = Esync-min(Esync);
                                    Esync  = Esync-mean(Esync);
                                    Esync  = Esync./max(abs(Esync));
                                    IxAux  = Esync(PosAuxEout);
                                    
                                    TxSymbAmos(CurrentTest,ObsCarrPos==kk,:) = IxAux;%RxSymbAmosDo = [];
                                    a=0;
                                    EvmMatRef(ObsCarrPos==kk,:) = IxAux;                   %Taking just the middle samples as references
                                    %                     if (kk==1)&&(Ploting)
                                    %                         figure;hold all;grid on;
                                    %                         plot(IxAux,zeros(1,length(IxAux)),'rx');
                                    %                         set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                    %                         a=0;
                                    %                     end
                                    %%
                                    if (~mod(CarrPass,2))&&(SendingUp)
                                        clear EoutModAux;
                                        EoutModAux = EoutTx(VetThisCarrTx==(RefCarr-1+kk...
                                            +1),:);
                                        
                                        %                     if (kk==1)&&(Ploting)
                                        %                         figure;hold all;grid on;
                                        %                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
                                        %                     end
                                        if AddingNoiseO
                                            EsN0 = 10*log10(snr);
                                            [~,sigpower] = MeasPower(EoutModAux);
                                            EoutModAux = awgn(EoutModAux,osnr,sigpower-30);
                                            [~,sigpower] = MeasPower(EoutModAux);
                                        end
                                        %                     if (kk==1)&&(Ploting)
                                        % %                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
                                        %                         PrintInfo(Ploting*4,abs(EoutModAux(kk,1:end-...
                                        %                         StuffSampels)).^2./max(abs(EoutModAux(kk,1:end-...
                                        %                         StuffSampels)).^2),t(2*NumAmosCP+NPPB),2*...
                                        %                         NumAmosCP+NPPB,abs(EoutModAux(1:end-StuffSampels...
                                        %                         )).^2./max(abs(EoutModAux(kk,1:end-StuffSampels)...
                                        %                                                                     ).^2));
                                        %                         a=4;
                                        %                     end
                                        EoutModTem(ObsCarrPos(kk+1),1:length(EoutModAux)...
                                            ) = EoutModAux;%Adding the current Modulated output to the final OutPut
                                    end
                                else
                                    EoutModAux = EoutTx(VetThisCarrTx==(RefCarr-1+kk),:);
                                    %EoutMod = EoutMod + EoutModAux + EoutTx(kk+1,:);      %Adding the current Modulated output to the final OutPut
                                    %                     if (kk==1)&&(Ploting)
                                    %                         figure;hold all;grid on;
                                    %                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
                                    %                     end
                                    if AddingNoiseO
                                        EsN0 = 10*log10(snr);
                                        [~,sigpower] = MeasPower(EoutModAux);
                                        EoutModAux = awgn(EoutModAux,osnr,sigpower-30);
                                        [~,sigpower] = MeasPower(EoutModAux);
                                    end
                                    %                     if (kk==1)&&(Ploting)
                                    %                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
                                    %                         PrintInfo(Ploting*4,abs(EoutModAux(kk,1:end-...
                                    %                         StuffSampels)).^2./max(abs(EoutModAux(kk,1:end-...
                                    %                         StuffSampels)).^2),t(2*NumAmosCP+NPPB),2*...
                                    %                         NumAmosCP+NPPB,abs(EoutModAux(1:end-StuffSampels...
                                    %                         )).^2./max(abs(EoutModAux(kk,1:end-StuffSampels)...
                                    %                                                                     ).^2));
                                    %                         a=4;
                                    %                     end
                                    EoutModTem(ObsCarrPos(kk),1:length(EoutModAux)) = ...
                                        EoutModAux;%Adding the current Modulated output to the final OutPut
                                    if ~mod(CarrPass,2)
                                        clear EoutModAux;
                                        EoutModAux = EoutTx(VetThisCarrTx==(RefCarr-1+kk...
                                            +1),:);
                                        %                     if (kk==1)&&(Ploting)
                                        %                         figure;hold all;grid on;
                                        %                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
                                        %                     end
                                        if AddingNoiseO
                                            EsN0 = 10*log10(snr);
                                            [~,sigpower] = MeasPower(EoutModAux);
                                            EoutModAux = awgn(EoutModAux,osnr,sigpower-30);
                                            [~,sigpower] = MeasPower(EoutModAux);
                                        end
                                        %                     if (kk==1)&&(Ploting)
                                        %                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
                                        %                         PrintInfo(Ploting*4,abs(EoutModAux(kk,1:end-...
                                        %                         StuffSampels)).^2./max(abs(EoutModAux(kk,1:end-...
                                        %                         StuffSampels)).^2),t(2*NumAmosCP+NPPB),2*...
                                        %                         NumAmosCP+NPPB,abs(EoutModAux(1:end-StuffSampels...
                                        %                         )).^2./max(abs(EoutModAux(kk,1:end-StuffSampels)...
                                        %                                                                     ).^2));
                                        %                         a=4;
                                        %                     end
                                        EoutModTem(ObsCarrPos(kk+1),1:length(EoutModAux)...
                                            ) = EoutModAux;%Adding the current Modulated output to the final OutPut
                                    end
                                end
                            end
                            if size(EoutModTem,1)<NumCarr
                                EoutModTem(end+1,:) = 0;
                            end
                            if IfftOrSum
                                %                 Ekaux=sum(EoutModTem);
                                %                 figure;hold all;
                                %                 plot(f,20*log10(abs(fftshift(fft(Ekaux)./length(Ekaux)))));
                                %                 set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                if size(EoutModTem,1)>1
                                    [EoutMod,~,~] = OpticalIFFTN(t,T,MaxNumStagT,EoutModTem);
                                else
                                    EoutMod = EoutModTem;
                                end
                                %                 plot(f,20*log10(abs(fftshift(fft(EoutMod)./length(...
                                %                                                               EoutMod)))));
                                %                 axis([0.9*fin 1.2*NumCarr*fc+fin -140 -20]);
                                %                 a=1;
                            else
                                if size(EoutModTem,1)>1
                                    EoutMod = sum(EoutModTem);
                                else
                                    EoutMod = EoutModTem;
                                end
                            end
                        case 'DQPSK'
                            %%            Generate the data DQPSK
                            EoutMod    = 0;
                            EoutModTem = 0;
                            for kk=InitCarrDo:CarrPass:NumCarr                             %Generating different data for each carrier
                                if CarrUsedDo(kk)
                                    if SendingData == 1                                    % Frist it is chosen to transmite a random information...
                                        TxData = (randi(2,1,NbDQPSK)-1);                   %Creating the data stream to be transmited
                                        TxData(1:JusPos) = JusVal;                         %Adding the Justification simble at the begining of the stream to sincronize received data frame
                                        TxData(end-(JusPos-1):end) = JusValEnd;            %Adding the Justification simble at the end of the stream to sincronize received data frame
                                        PreBits  = [0 0];                                  %Seting the start point for the DQPSK maping
                                    elseif SendingData == 2
                                        %For previous test it was needed a known bit sequence. This
                                        %section creates it. The [I Q] output will be:
                                        %[1 1 0 0 0 1 0 0 0 0 1 0 0 1 0 1 1 1 1 1 0 1 1 0 1 1 1 0 1
                                        % 0 0 0]
                                        %This sequency covers all possible combinations of I and Q
                                        %which is:
                                        %I = 1 0 0 0 0 1 0 0 1 1 0 1 1 1 1 0
                                        %Q = 1 0 1 0 0 0 1 1 1 1 1 0 1 0 0 0
                                        TxDataSeg = [0 0; %01
                                            0 0; %02
                                            0 1; %03
                                            1 0; %04
                                            1 1; %05
                                            1 0; %06
                                            0 0; %07
                                            1 1; %08
                                            0 1; %09
                                            1 1; %10
                                            1 0; %11
                                            0 0; %12
                                            1 0; %13
                                            0 1; %14
                                            1 1; %15
                                            0 1];%16
                                        TxData   = [];
                                        PreBits  = [0 0];
                                        countseg = 1;
                                        for jj=1:Nb
                                            TxData   = [TxData TxDataSeg(countseg,:)];
                                            countseg = countseg + 1;
                                            if countseg > 16
                                                countseg = 1;
                                            end
                                        end
                                    else                                                   %... or just one high pulse for testing the channel
                                        TxData = zeros(1,NbDQPSK);                         %Creating the base fo the data stream
                                        TxData((end/2)-1:end/2) = [1 1];                   %Adding the code in the data stream to produce the highest level possible
                                        PreBits  = [0 0];
                                    end
                                    TxDataPos       = linspace(1,NbDQPSK,NbDQPSK);         %Auxiliar variable to split the TxData
                                    [DataI,DataQ]   = DqpskEncodEq(TxData,PreBits);        %Maping the income data to the DQPSK format
                                    %Converting the I and Q components to the polirezed NRZ format
                                    DataI(DataI==1) =  1.9;
                                    DataI(DataI==0) = -1.9;
                                    DataQ(DataQ==1) =  1.9;
                                    DataQ(DataQ==0) = -1.9;
                                    TxOdd  = TxData(logical(mod(TxDataPos,2)));            %Spliting the information of odd positions
                                    TxEven = TxData(~(mod(TxDataPos,2)));                  %Spliting the information of even positions
                                    TxDataMat(kk,:)         = TxOdd;                       %Storring the transmited information for latter evaluation
                                    TxDataMat(kk+NumCarr,:) = TxEven;                      %Storring the transmited information for latter evaluation
                                    %Once that we have the data to be transmited it needs to be converted to
                                    %the eletrical signal that hereafter will modulate our carrier through the
                                    %MZM-I
                                    SigI = rectpulse(DataI,NPPB);
                                    SigQ = rectpulse(DataQ,NPPB);
                                    
                                    % Adding CP to the data
                                    if AddCP
                                        TxAux1 = reshape(SigI,NPPB,NbDQPSK/2);
                                        if SetCpSampZer
                                            TxAux1 = [zeros(NumAmosCP,size(TxAux1,2));...
                                                TxAux1;zeros(NumAmosCP,size(TxAux1,2))];
                                        else
                                            TxAux1 = [TxAux1(1:NumAmosCP,:);TxAux1;...
                                                TxAux1(end-(NumAmosCP-1):end,:)];
                                        end
                                        TxAux1 = reshape(TxAux1,1,(2*NumAmosCP+NPPB)*...
                                            NbDQPSK/2);
                                        SigI = [TxAux1 TxAux1(end-(StuffSampels-1):end)];
                                        
                                        TxAux2 = reshape(SigQ,NPPB,NbDQPSK/2);
                                        if SetCpSampZer
                                            TxAux2 = [zeros(NumAmosCP,size(TxAux2,2));...
                                                TxAux2;zeros(NumAmosCP,size(TxAux2,2))];
                                        else
                                            TxAux2 = [TxAux2(1:NumAmosCP,:);TxAux2;...
                                                TxAux2(end-(NumAmosCP-1):end,:)];
                                        end
                                        TxAux2 = reshape(TxAux2,1,(2*NumAmosCP+NPPB)*...
                                            NbDQPSK/2);
                                        SigQ = [TxAux2 TxAux2(end-(StuffSampels-1):end)];
                                    end
                                    
                                    %Because of reasons, the real world does not have components capable
                                    %instantly change the voltage level. Thus, to emulate this behaviour the
                                    %eletrical PAM signal pass through a gaussian filter that will remove the
                                    %frequencies of higher order, which will result in small slope in the level
                                    %variation
                                    [BitFilt,~] = FiltroGaussiano(f,BWD,CenFeq,FiltOrd);   %Creating filter for conformation of the input information
                                    BitFilt = fftshift(BitFilt);                           %Doing a shift on the vector for matching the transmited data
                                    TxSigI = SigI;                                         %Conforming the information and Creating the modulation signal
                                    TxSigQ = SigQ;                                         %Conforming the information and Creating the modulation signal
                                    %TxSigI = ifft(fft(SigI).*BitFilt);                    %Conforming the information and Creating the modulation signal
                                    %TxSigQ = ifft(fft(SigQ).*BitFilt);                    %Conforming the information and Creating the modulation signal
                                    
                                    %                     if (kk==1)&&(Ploting)
                                    %                       figure;hold all;grid on;
                                    %                       plot(t,TxSig1,t,TxSig2)
                                    %                     end
                                    if AddingNoiseE
                                        EsN0 = 10*log10(snr);
                                        [~,sigpower] = MeasPower(TxSigI);
                                        [~,sigpower2] = MeasPower(TxSigQ);
                                        TxSigI = awgn(TxSigI,snr,sigpower-30);
                                        TxSigQ = awgn(TxSigQ,snr,sigpower2-30);
                                        [~,sigpower] = MeasPower(TxSigI);
                                        [~,sigpower2] = MeasPower(TxSigQ);
                                    end
                                    %                     if (kk==1)&&(Ploting)
                                    %                       plot(t,TxSig1,t,TxSig2)
                                    %                       PrintInfo(Ploting*3,TxSigI(1:end-StuffSampels),2*t(2*NumAmosCP+NPPB),2*NumAmosCP+NPPB,TxSigQ(1:end-StuffSampels));
                                    %                       a=4;
                                    %                     end
                                    [EoutModAux] = IqMod(EoutTx(VetThisCarrTx==(RefCarr-...
                                        1+kk),:),TxSigI,TxSigQ,Vpi,V0);
                                    %                     if (kk==1)&&(Ploting)
                                    %                         figure;hold all;grid on;
                                    %                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
                                    %                     end
                                    if AddingNoiseO
                                        EsN0 = 10*log10(snr);
                                        [~,sigpower] = MeasPower(EoutModAux);
                                        EoutModAux = awgn(EoutModAux,osnr,sigpower-30);
                                        [~,sigpower] = MeasPower(EoutModAux);
                                    end
                                    %                     if (kk==1)&&(Ploting)
                                    % %                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
                                    %                         PrintInfo(Ploting*4,abs(EoutModAux(kk,1:end-...
                                    %                         StuffSampels)).^2./max(abs(EoutModAux(kk,1:end-...
                                    %                         StuffSampels)).^2),t(2*NumAmosCP+NPPB),2*...
                                    %                         NumAmosCP+NPPB,abs(EoutModAux(1:end-StuffSampels...
                                    %                         )).^2./max(abs(EoutModAux(kk,1:end-StuffSampels)...
                                    %                                                                     ).^2));
                                    %                         a=4;
                                    %                     end
                                    EoutModTem(ObsCarrPos(kk),1:length(EoutModAux)) = ...
                                        EoutModAux;%Adding the current Modulated output to the final OutPut
                                    %% Taking the sampling the EVM meassurement
                                    PhaDel = 1*pi/4;
                                    TimDel = T;
                                    [EoutA,EoutB] = DelayInterf(t,TimDel,PhaDel,EoutModAux);%Coverting phase shift to amplitude variation
                                    PhaDel = -1*pi/4;
                                    TimDel = T;
                                    [EoutC,EoutD] = DelayInterf(t,TimDel,PhaDel,EoutModAux);%Coverting phase shift to amplitude variation
                                    EoutA = EoutA.*conj(EoutA);
                                    EoutB = EoutB.*conj(EoutB);
                                    EoutC = EoutC.*conj(EoutC);
                                    EoutD = EoutD.*conj(EoutD);                   %Shifting the filter for matching the received signal
                                    EoutA = ifft(fft(EoutA).*EvmFilt);
                                    EoutB = ifft(fft(EoutB).*EvmFilt);
                                    EoutC = ifft(fft(EoutC).*EvmFilt);
                                    EoutD = ifft(fft(EoutD).*EvmFilt);
                                    EoutI = (EoutB - EoutA);
                                    EoutQ = (EoutD - EoutC);
                                    EoutI = EoutI-mean(EoutI);
                                    EoutQ = EoutQ-mean(EoutQ);
                                    EoutI = EoutI./max(abs(EoutI));                                %Normalizing the signal
                                    EoutQ = EoutQ./max(abs(EoutQ));                                %Normalizing the signal
                                    PosAuxEout      = (2*NumAmosCP+NPPB)/2:(2*NumAmosCP+...
                                        NPPB):(length(EoutModAux)-StuffSampels);%Varriable respossible to take just the samples at the middle of the symbol
                                    IxAux = EoutI(PosAuxEout) + 1j.*EoutQ(PosAuxEout);
                                    
                                    TxSymbAmos(CurrentTest,ObsCarrPos==kk,:) = IxAux;%RxSymbAmosDo = [];
                                    
                                    EvmMatRef(ObsCarrPos==kk,:) = IxAux;       %Taking just the middle samples as references
                                    %                     if (kk==1)&&(Ploting)
                                    %                         figure;hold all;grid on;
                                    %                         plot(real(IxAux),imag(IxAux),'rx');
                                    %                         set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                    %                         a=0;
                                    %                     end
                                    %%
                                    if (~mod(CarrPass,2))&&(SendingUp)
                                        clear EoutModAux;
                                        EoutModAux = EoutTx(VetThisCarrTx==(RefCarr-1+kk+1),:);
                                        %                     if (kk==1)&&(Ploting)
                                        %                         figure;hold all;grid on;
                                        %                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
                                        %                     end
                                        if AddingNoiseO
                                            EsN0 = 10*log10(snr);
                                            [~,sigpower] = MeasPower(EoutModAux);
                                            EoutModAux = awgn(EoutModAux,osnr,sigpower-30);
                                            [~,sigpower] = MeasPower(EoutModAux);
                                        end
                                        %                     if (kk==1)&&(Ploting)
                                        %                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
                                        %                         PrintInfo(Ploting*4,abs(EoutModAux(kk,1:end-...
                                        %                         StuffSampels)).^2./max(abs(EoutModAux(kk,1:end-...
                                        %                         StuffSampels)).^2),2*t(2*NumAmosCP+NPPB),2*...
                                        %                         NumAmosCP+NPPB,abs(EoutModAux(1:end-StuffSampels...
                                        %                         )).^2./max(abs(EoutModAux(kk,1:end-StuffSampels)...
                                        %                                                                     ).^2));
                                        %                         a=4;
                                        %                     end
                                        EoutModTem(ObsCarrPos(kk+1),1:length(EoutModAux)...
                                            ) = EoutModAux;%Adding the current Modulated output to the final OutPut
                                    end
                                else
                                    EoutModAux = EoutTx(VetThisCarrTx==(RefCarr-1+kk),:);
                                    
                                    %                     if (kk==1)&&(Ploting)
                                    %                         figure;hold all;grid on;
                                    %                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
                                    %                     end
                                    if AddingNoiseO
                                        EsN0 = 10*log10(snr);
                                        [~,sigpower] = MeasPower(EoutModAux);
                                        EoutModAux = awgn(EoutModAux,osnr,sigpower-30);
                                        [~,sigpower] = MeasPower(EoutModAux);
                                    end
                                    %                     if (kk==1)&&(Ploting)
                                    % %                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
                                    %                         PrintInfo(Ploting*4,abs(EoutModAux(kk,1:end-...
                                    %                         StuffSampels)).^2./max(abs(EoutModAux(kk,1:end-...
                                    %                         StuffSampels)).^2),t(2*NumAmosCP+NPPB),2*...
                                    %                         NumAmosCP+NPPB,abs(EoutModAux(1:end-StuffSampels...
                                    %                         )).^2./max(abs(EoutModAux(kk,1:end-StuffSampels)...
                                    %                                                                     ).^2));
                                    %                         a=4;
                                    %                     end
                                    EoutModTem(ObsCarrPos(kk),1:length(EoutModAux)) = ...
                                        EoutModAux;%Adding the current Modulated output to the final OutPut
                                    if ~mod(CarrPass,2)
                                        clear EoutModAux;
                                        EoutModAux = EoutTx(VetThisCarrTx==(RefCarr-1+kk...
                                            +1),:);
                                        %                     if (kk==1)&&(Ploting)
                                        %                         figure;hold all;grid on;
                                        %                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
                                        %                     end
                                        if AddingNoiseO
                                            EsN0 = 10*log10(snr);
                                            [~,sigpower] = MeasPower(EoutModAux);
                                            EoutModAux = awgn(EoutModAux,osnr,sigpower-30);
                                            [~,sigpower] = MeasPower(EoutModAux);
                                        end
                                        %                     if (kk==1)&&(Ploting)
                                        %                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
                                        %                         PrintInfo(Ploting*4,abs(EoutModAux(kk,1:end-...
                                        %                         StuffSampels)).^2./max(abs(EoutModAux(kk,1:end-...
                                        %                         StuffSampels)).^2),t(2*NumAmosCP+NPPB),2*...
                                        %                         NumAmosCP+NPPB,abs(EoutModAux(1:end-StuffSampels...
                                        %                         )).^2./max(abs(EoutModAux(kk,1:end-StuffSampels)...
                                        %                                                                     ).^2));
                                        %                         a=4;
                                        %                     end
                                        EoutModTem(ObsCarrPos(kk+1),1:length(EoutModAux)...
                                            ) = EoutModAux;%Adding the current Modulated output to the final OutPut
                                    end
                                end
                            end
                            if size(EoutModTem,1)<NumCarr
                                EoutModTem(end+1,:) = 0;
                            end
                            if IfftOrSum
                                %                 Ekaux=sum(EoutModTem);
                                %                 figure;hold all;
                                %                 plot(f,20*log10(abs(fftshift(fft(Ekaux)./length(Ekaux)))));
                                %                 set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                if size(EoutModTem,1)>1
                                    [EoutMod,~,~] = OpticalIFFTN(t,T,MaxNumStagT,EoutModTem);
                                else
                                    EoutMod = EoutModTem;
                                end
                                %                 plot(f,20*log10(abs(fftshift(fft(EoutMod)./length(...
                                %                                                               EoutMod)))));
                                %                 axis([0.9*fin 1.2*NumCarr*fc+fin -140 -20]);
                                %                 a=1;
                            else
                                if size(EoutModTem,1)>1
                                    EoutMod = sum(EoutModTem);
                                else
                                    EoutMod = EoutModTem;
                                end
                            end
                        case '4PAM'
                            %%        Generate the data 4PAM
                            EoutMod    = 0;
                            EoutModTem = 0;
                            for kk=InitCarrDo:CarrPass:NumCarr                             %Generating different data for each carrier
                                if CarrUsedDo(kk)                                          % Frist it is chosen to transmite a random information...
                                    if SendingData
                                        TxData = (randi(2,1,Nb4Pam)-1);                    %Creating the data stream to be transmited
                                        TxData(1:JusPos) = JusVal;                         %Adding the Justification simble at the begining of the stream to sincronize received data frame
                                        TxData(end-(JusPos-1):end) = JusValEnd;            %Adding the Justification simble at the end of the stream to sincronize received data frame
                                    else                                                   %... or just one high pulse for testing the channel
                                        TxData = zeros(1,Nb4Pam);                          %Creating the base fo the data stream
                                        TxDataPos = linspace(1,Nb4Pam,Nb4Pam);             %Creating a vector to auxilliate the addresing process
                                        TxData(~mod(TxDataPos,2))=1;                       %Adding the code in the data stream to produce the lowest level possible
                                        TxData(end/2 - 1) = 1;                             %Adding the code in the data stream to produce the highest level possible
                                    end
                                    
                                    TxDataMat(kk,:) = TxData;                              %Storring the transmited information for latter evaluation
                                    %Once that we have the data to be transmited it needs to be
                                    %converted to the eletrical signal that hereafter will modulate
                                    %our carrier through the MZM-I
                                    if Which4PAM                                           %Selecting if the PAM will be done on Electrical or Optical domain
                                        %  Good for others kms
                                        [Phi1,Phi2] = Maping4PamIq(TxData,Vmin,Vmax,...
                                            ModSchem,FiberLength,SetCpSampZer);%Generating the eletrical signal for the optical PAM4
                                    else
                                        [Phi1,Phi2] = Maping4Pam(TxData,VPI,Polirized,...
                                            MaxAmp4PAM);%Generating the eletrical signal for the electrical PAM4
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
                                    if AddCP
                                        TxAux1 = reshape(TxSig1,NPPB,Nb4Pam/2);
                                        if SetCpSampZer
                                            TxAux1 = [zeros(NumAmosCP,size(TxAux1,2));...
                                                TxAux1;zeros(NumAmosCP,size(TxAux1,2))];
                                        else
                                            TxAux1 = [TxAux1(1:NumAmosCP,:);TxAux1;...
                                                TxAux1(end-(NumAmosCP-1):end,:)];
                                        end
                                        TxAux1 = reshape(TxAux1,1,(2*NumAmosCP+NPPB)*...
                                            Nb4Pam/2);
                                        TxSig1 = [TxAux1 TxAux1(end-(StuffSampels-1):end)];
                                        TxAux2 = reshape(TxSig2,NPPB,Nb4Pam/2);
                                        if SetCpSampZer
                                            TxAux2 = [zeros(NumAmosCP,size(TxAux2,2));...
                                                TxAux2;zeros(NumAmosCP,size(TxAux2,2))];
                                        else
                                            TxAux2 = [TxAux2(1:NumAmosCP,:);TxAux2;...
                                                TxAux2(end-(NumAmosCP-1):end,:)];
                                        end
                                        TxAux2 = reshape(TxAux2,1,(2*NumAmosCP+NPPB)*...
                                            Nb4Pam/2);
                                        TxSig2 = [TxAux2 TxAux2(end-(StuffSampels-1):end)];
                                    end
                                    
                                    %                     if (kk==1)&&(Ploting)
                                    %                       figure;hold all;grid on;
                                    %                       plot(t,TxSig1,t,TxSig2)
                                    %                     end
                                    if AddingNoiseE
                                        EsN0 = 10*log10(snr);
                                        [~,sigpower] = MeasPower(TxSig1);
                                        [~,sigpower2] = MeasPower(TxSig2);
                                        TxSig1 = awgn(TxSig1,snr,sigpower-30);
                                        TxSig2 = awgn(TxSig2,snr,sigpower2-30);
                                        [~,sigpower] = MeasPower(TxSig1);
                                        [~,sigpower2] = MeasPower(TxSig2);
                                    end
                                    %                     if (kk==1)&&(Ploting)
                                    %                       plot(t,TxSig1,t,TxSig2)
                                    %                       PrintInfo(Ploting*4,TxSig1(1:end-StuffSampels),2*t(2*NumAmosCP+NPPB),2*NumAmosCP+NPPB,TxSig2(1:end-StuffSampels));
                                    %                     end
                                    %Because of reasons, the real world does not have components capable
                                    %instantly change the voltage level. Thus, to emulate this behaviour the
                                    %eletrical PAM signal pass through a gaussian filter that will remove the
                                    %frequencies of higher order, which will result in small slope in the level
                                    %variation
                                    [BitFilt,~] = FiltroGaussiano(f,BWD,CenFeq,FiltOrd);   %Creating filter for conformation of the input information
                                    BitFilt     = fftshift(BitFilt);                       %Doing a shift on the vector for matching the transmited data
                                    %TxSig1     = ifft(fft(TxSig1).*BitFilt);              %Conforming the information and Creating the modulation signal
                                    %TxSig2     = ifft(fft(TxSig2).*BitFilt);              %Conforming the information and Creating the modulation signal
                                    %                     if (kk==1)&&(Ploting)
                                    %                         PrintInfo(Ploting*4,TxSig1(kk,1:end-...
                                    %                         StuffSampels),t(2*NumAmosCP+NPPB),2*NumAmosCP+...
                                    %                                           NPPB,TxSig2(1:end-StuffSampels));
                                    %                         a=4;
                                    %                     end
                                    
                                    %If the user choses to modulate all carriers at the same time there is no
                                    %need of this script to generate data for each individual carrier.
                                    %Therefore, this For loop can be halted
                                    U.U1t = TxSig2;                                        %Assigning the electrical signal to one drive of the MZM
                                    U.U2t = TxSig1;                                        %Assigning the electrical signal to another drive of the MZM
                                    if Which4PAM
                                        if ModSchem
                                            [EoutModAux] = ...
                                                Mach_Zehnder_Modulator_simplificado_2(freqGHz,EoutTx...
                                                (VetThisCarrTx==(RefCarr-1+kk),:),U,...
                                                MZ_Input_File);
                                        else
                                            [EoutModAux] = ...
                                                Mach_Zehnder_Modulator_simplificado_2(freqGHz,EoutTx...
                                                (VetThisCarrTx==(RefCarr-1+kk),:),U,...
                                                MZ_Input_File);
                                        end
                                    else
                                        [EoutModAux] = ...
                                            Mach_Zehnder_Modulator_simplificado_2(freqGHz,EoutTx(...
                                            VetThisCarrTx==(RefCarr-1+kk),:),U,MZ_Input_File);
                                    end
                                    
                                    if AddingNoiseO
                                        EsN0 = 10*log10(snr);
                                        [~,sigpower] = MeasPower(EoutModAux);
                                        EoutModAux = awgn(EoutModAux,osnr,sigpower-30);
                                        [~,sigpower] = MeasPower(EoutModAux);
                                    end
                                    %                     if (kk==1)&&(Ploting)
                                    % %                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
                                    %                         PrintInfo(Ploting*4,abs(EoutModAux(kk,1:end-...
                                    %                         StuffSampels)).^2./max(abs(EoutModAux(kk,1:end-...
                                    %                         StuffSampels)).^2),2*t(2*NumAmosCP+NPPB),2*...
                                    %                         NumAmosCP+NPPB,abs(EoutModAux(1:end-StuffSampels...
                                    %                         )).^2./max(abs(EoutModAux(kk,1:end-StuffSampels)...
                                    %                                                                     ).^2));
                                    %                         a=4;
                                    %                     end
                                    EoutModTem(ObsCarrPos(kk),1:length(EoutModAux)) = ...
                                        EoutModAux;%Adding the current Modulated output to the final OutPut
                                    %% Taking the sampling the EVM meassurement
                                    PosAuxEout = (2*NumAmosCP+NPPB)/2:(2*NumAmosCP+NPPB)...
                                        :(length(EoutModAux)-StuffSampels);%Varriable respossible to take just the samples at the middle of the symbol
                                    Ix         = EoutModAux.*conj(EoutModAux);             %Recovering the signal that will be transmited
                                    Ix         = ifft(fft(Ix).*EvmFilt);                   %Removing higher signal generated by the receiving process
                                    Ix         = Ix - min(Ix);
                                    Ix         = Ix./max(abs(Ix));                         %Normalizing the reference
                                    IxAux      = Ix(PosAuxEout);
                                    a=0;
                                    
                                    TxSymbAmos(CurrentTest,ObsCarrPos==kk,:) = IxAux;%RxSymbAmosDo = [];
                                    EvmMatRef(ObsCarrPos==kk,:) = IxAux;                   %Taking just the middle samples as references
                                    %                     if (kk==1)&&(Ploting)
                                    %                         figure;hold all;grid on;
                                    %                         plot(t,Ix);plot(t(PosAuxEout),IxAux,'x');
                                    %                         set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                    %                         figure;hold all;grid on;
                                    %                         plot(IxAux,zeros(1,length(IxAux)),'rx');
                                    %                         set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                    %                         a=0;
                                    %                     end
                                    %%
                                    if (~mod(CarrPass,2))&&(SendingUp)
                                        clear EoutModAux;
                                        EoutModAux = EoutTx(VetThisCarrTx==(RefCarr-1+kk...
                                            +1),:);
                                        %                     if (kk==1)&&(Ploting)
                                        %                         figure;hold all;grid on;
                                        %                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
                                        %                     end
                                        if AddingNoiseO
                                            EsN0 = 10*log10(snr);
                                            [~,sigpower] = MeasPower(EoutModAux);
                                            EoutModAux = awgn(EoutModAux,osnr,sigpower-30);
                                            [~,sigpower] = MeasPower(EoutModAux);
                                        end
                                        %                     if (kk==1)&&(Ploting)
                                        %                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
                                        %                         PrintInfo(Ploting*4,abs(EoutModAux(kk,1:end-...
                                        %                         StuffSampels)).^2./max(abs(EoutModAux(kk,1:end-...
                                        %                         StuffSampels)).^2),2*t(2*NumAmosCP+NPPB),2*...
                                        %                         NumAmosCP+NPPB,abs(EoutModAux(1:end-StuffSampels...
                                        %                         )).^2./max(abs(EoutModAux(kk,1:end-StuffSampels)...
                                        %                                                                     ).^2));
                                        %                         a=4;
                                        %                     end
                                        EoutModTem(ObsCarrPos(kk+1),1:length(EoutModAux)...
                                            ) = EoutModAux;%Adding the current Modulated output to the final OutPut
                                    end
                                    %                     if (kk==1)&&(Ploting)
                                    %                         figure;hold all;grid on;
                                    %                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
                                    %                     end
                                    if AddingNoiseO
                                        EsN0 = 10*log10(snr);
                                        [~,sigpower] = MeasPower(EoutModAux);
                                        EoutModAux = awgn(EoutModAux,osnr,sigpower-30);
                                        [~,sigpower] = MeasPower(EoutModAux);
                                    end
                                    %                     if (kk==1)&&(Ploting)
                                    %                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
                                    %                         PrintInfo(Ploting*4,abs(EoutModAux(kk,1:end-...
                                    %                         StuffSampels)).^2./max(abs(EoutModAux(kk,1:end-...
                                    %                         StuffSampels)).^2),2*t(2*NumAmosCP+NPPB),2*...
                                    %                         NumAmosCP+NPPB,abs(EoutModAux(1:end-StuffSampels...
                                    %                         )).^2./max(abs(EoutModAux(kk,1:end-StuffSampels)...
                                    %                                                                     ).^2));
                                    %                         a=4;
                                    %                     end
                                else
                                    EoutModAux = EoutTx(VetThisCarrTx==(RefCarr-1+kk),:);
                                    %                     if (kk==1)&&(Ploting)
                                    %                         figure;hold all;grid on;
                                    %                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
                                    %                     end
                                    if AddingNoiseO
                                        EsN0 = 10*log10(snr);
                                        [~,sigpower] = MeasPower(EoutModAux);
                                        EoutModAux = awgn(EoutModAux,osnr,sigpower-30);
                                        [~,sigpower] = MeasPower(EoutModAux);
                                    end
                                    %                     if (kk==1)&&(Ploting)
                                    % %                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
                                    %                         PrintInfo(Ploting*4,abs(EoutModAux(kk,1:end-...
                                    %                         StuffSampels)).^2./max(abs(EoutModAux(kk,1:end-...
                                    %                         StuffSampels)).^2),t(2*NumAmosCP+NPPB),2*...
                                    %                         NumAmosCP+NPPB,abs(EoutModAux(1:end-StuffSampels...
                                    %                         )).^2./max(abs(EoutModAux(kk,1:end-StuffSampels)...
                                    %                                                                     ).^2));
                                    %                         a=4;
                                    %                     end
                                    EoutModTem(ObsCarrPos(kk),1:length(EoutModAux)) = ...
                                        EoutModAux;%Adding the current Modulated output to the final OutPut
                                    if ~mod(CarrPass,2)
                                        clear EoutModAux;
                                        EoutModAux = EoutTx(VetThisCarrTx==(RefCarr-1+kk...
                                            +1),:);
                                        %                     if (kk==1)&&(Ploting)
                                        %                         figure;hold all;grid on;
                                        %                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
                                        %                     end
                                        if AddingNoiseO
                                            EsN0 = 10*log10(snr);
                                            [~,sigpower] = MeasPower(EoutModAux);
                                            EoutModAux = awgn(EoutModAux,osnr,sigpower-30);
                                            [~,sigpower] = MeasPower(EoutModAux);
                                        end
                                        %                     if (kk==1)&&(Ploting)
                                        %                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
                                        %                         PrintInfo(Ploting*4,abs(EoutModAux(kk,1:end-...
                                        %                         StuffSampels)).^2./max(abs(EoutModAux(kk,1:end-...
                                        %                         StuffSampels)).^2),t(2*NumAmosCP+NPPB),2*...
                                        %                         NumAmosCP+NPPB,abs(EoutModAux(1:end-StuffSampels...
                                        %                         )).^2./max(abs(EoutModAux(kk,1:end-StuffSampels)...
                                        %                                                                     ).^2));
                                        %                         a=4;
                                        %                     end
                                        EoutModTem(ObsCarrPos(kk+1),1:length(EoutModAux)...
                                            ) = EoutModAux;%Adding the current Modulated output to the final OutPut
                                    end
                                end
                            end
                            if size(EoutModTem,1)<NumCarr
                                EoutModTem(end+1,:) = 0;
                            end
                            if IfftOrSum
                                %                 Ekaux=sum(EoutModTem);
                                %                 figure;hold all;
                                %                 plot(f,20*log10(abs(fftshift(fft(Ekaux)./length(Ekaux)))));
                                %                 set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                if size(EoutModTem,1)>1
                                    [EoutMod,~,~] = OpticalIFFTN(t,T,MaxNumStagT,EoutModTem);
                                else
                                    EoutMod = EoutModTem;
                                end
                                %                 plot(f,20*log10(abs(fftshift(fft(EoutMod)./length(...
                                %                                                               EoutMod)))));
                                %                 axis([0.9*fin 1.2*NumCarr*fc+fin -140 -20]);
                                %                 a=1;
                            else
                                if size(EoutModTem,1)>1
                                    EoutMod = sum(EoutModTem);
                                else
                                    EoutMod = EoutModTem;
                                end
                            end
                        otherwise
                            %%        Generate the data OOK
                            EoutMod    = 0;
                            EoutModTem = 0;
                            for kk=InitCarrDo:CarrPass:NumCarr                             %Generating different data for each carrier
                                if CarrUsedDo(kk)                                          % Frist it is chosen to transmite just one high pulse for
                                    if AdjusData                                           %testing the channel...
                                        TxData        = zeros(1,Nb - NumBitDesc);          %creating vector of zerros for correcting the delay caused by the transmission
                                        TxData(end/2) = 1;                                 %Adding an strategical pulse for measuring the result.
                                        TxDataMat(kk,:) = TxData;                          %Storring the transmited information for latter evaluation
                                        if NRZPolarity
                                            TxData(TxData==0) = NrzMin;
                                            TxData(TxData==1) = NrzMax;
                                        end
                                    else                                                   %... or just a random information
                                        
                                        TxData = (randi(2,1,Nb - NumBitDesc)-1);           %Creating Random Information that will be loaded in each individual subcarrier
                                        TxData(1:JusLen)           = JusVal;               %Making the First 4 bits equal to zero
                                        TxData(end-(JusLen-1):end) = JusVal;               %Making the Last 4 bits equal to zero
                                        TxDataMat(kk,:) = TxData;                          %Storring the transmited information for latter evaluation
                                        if NRZPolarity
                                            TxData(TxData==0) = NrzMin;
                                            TxData(TxData==1) = NrzMax;
                                        end
                                    end
                                    
                                    %The signal generated are not yet with the same number of samples as the
                                    %OFCS loaded. These nexte lines do oversampling
                                    TxDataRes = rectpulse(TxData,NPPB);                    %Changing the length of the Data acordingly with the time vector
                                    %Thus, if it would be required to add cycle prefix the number of samples
                                    %per symbol needs to change as well as some adjustments needs to be done
                                    %for the new signal match in size with the size of the vector time. This
                                    %problem just exist on simulation, at practice the main point is the
                                    %syncronism of the signals.
                                    % Adding CP to the data
                                    if AddCP
                                        TxAux = reshape(TxDataRes,NPPB,Nb-NumBitDesc);
                                        if SetCpSampZer
                                            TxAux = [zeros(NumAmosCP,size(TxAux,2));...
                                                TxAux;zeros(NumAmosCP,size(TxAux,2))];
                                        else
                                            TxAux = [TxAux(1:NumAmosCP,:);TxAux;TxAux(...
                                                end-(NumAmosCP-1):end,:)];
                                        end
                                        TxAux = reshape(TxAux,1,(2*NumAmosCP+NPPB)*(Nb-...
                                            NumBitDesc));
                                        TxDataRes = [TxAux TxAux(end-(StuffSampels-1):...
                                            end)];
                                    end
                                    
                                    
                                    %Because of reasons, the real world does not have components capable
                                    %instantly change the voltage level. Thus, to emulate this behaviour the
                                    %eletrical PAM signal pass through a gaussian filter that will remove the
                                    %frequencies of higher order, which will result in small slope in the level
                                    %variation
                                    
                                    [BitFilt,~] = FiltroGaussiano(f,BWD,CenFeq,FiltOrd);   %Creating filter for conformation of the input information
                                    BitFilt = fftshift(BitFilt);                           %Doing a shift on the vector for matching the transmited data
                                    %                   TxSig = ifft(fft(TxDataRes).*BitFilt);                 %Conforming the information and Creating the modulation signal
                                    TxSig = TxDataRes;
                                    %                     if (kk==1)&&(Ploting)%Conforming the information and Creating the modulation signal
                                    %                     figure;hold all;grid on;
                                    %                     plot(t,TxSig)
                                    %                     end
                                    if AddingNoiseE
                                        EsN0 = 10*log10(snr);
                                        [~,sigpower] = MeasPower(TxSig);
                                        TxSig = awgn(TxSig,snr,sigpower-30);
                                        [~,sigpower] = MeasPower(TxSig);
                                    end
                                    %                     if (kk==1)&&(Ploting)
                                    %                     plot(t,TxSig)
                                    %                     PrintInfo(Ploting*4,TxSig(1:end-StuffSampels),2*t(2*NumAmosCP+NPPB),2*NumAmosCP+NPPB,TxSig(1:end-StuffSampels));
                                    %                     end
                                    %                   TxSigMat(kk,:) = TxSig;                                %Storring the transmited information for latter evaluation
                                    %                   if (kk==1)&&(Ploting)
                                    %                       PrintInfo(Ploting*5,t,TxSigMat(kk,:),TxDataRes);
                                    %                   end
                                    %If the user choses to modulate all carriers at the same time there is no
                                    %need of this script to generate data for each individual carrier.
                                    %Therefore, this For loop can be halted
                                    
                                    
                                    %Assigning the eletrical signal to one drive of the MZM - The aplitude of
                                    %the signal can be controlled by the variable DatGai, which can be
                                    %understood as an gain for the eletrical signal or an atenuation. The
                                    %second signal will be similar with the only difference a phase shift of pi
                                    U.U1t = TxSig;
                                    U.U2t = exp(-1j*pi).*TxSig;
                                    %As both signals will have the mostrly the same characteristics with the
                                    %only difference the phase shift of 180 degress. The MZM-I will be working
                                    %on the Push-Pull configuration. It is necessary to reduce the Chirp noise
                                    %to zero.
                                    [EoutModAux,~]=Mach_Zehnder_Modulator_simplificado_2(freqGHz...
                                        ,EoutTx(VetThisCarrTx==(RefCarr-1+kk),:),U,...
                                        MZ_Input_File);%Modulating individual carriers
                                    %                     if (kk==1)&&(Ploting)
                                    %                         figure;hold all;grid on;
                                    %                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
                                    %                     end
                                    if AddingNoiseO
                                        EsN0 = 10*log10(snr);
                                        [~,sigpower] = MeasPower(EoutModAux);
                                        EoutModAux = awgn(EoutModAux,osnr,sigpower-30);
                                        [~,sigpower] = MeasPower(EoutModAux);
                                    end
                                    %                     if (kk==1)&&(Ploting)
                                    % %                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
                                    %                         PrintInfo(Ploting*4,abs(EoutModAux(kk,1:end-...
                                    %                         StuffSampels)).^2./max(abs(EoutModAux(kk,1:end-...
                                    %                         StuffSampels)).^2),t(2*NumAmosCP+NPPB),2*...
                                    %                         NumAmosCP+NPPB,abs(EoutModAux(1:end-StuffSampels...
                                    %                         )).^2./max(abs(EoutModAux(kk,1:end-StuffSampels)...
                                    %                                                                     ).^2));
                                    %                         a=4;
                                    %                     end
                                    
                                    EoutModTem(ObsCarrPos(kk),1:length(EoutModAux)) = ...
                                        EoutModAux;%Adding the current Modulated output to the final OutPut
                                    %% Taking the sampling the EVM meassurement
                                    PosAuxEout = (2*NumAmosCP+NPPB)/2:(2*NumAmosCP+NPPB)...
                                        :(length(EoutModAux)-StuffSampels);%Varriable respossible to take just the samples at the middle of the symbol
                                    Ix         = EoutModAux.*conj(EoutModAux);             %Recovering the signal that will be transmited
                                    Ix         = ifft(fft(Ix).*EvmFilt);                   %Removing higher signal generated by the receiving process
                                    Ix         = Ix - min(Ix);
                                    Ix         = Ix./max(abs(Ix));                         %Normalizing the reference
                                    IxAux      = Ix(PosAuxEout);                           %Normalizing the reference
                                    a=0;
                                    
                                    TxSymbAmos(CurrentTest,ObsCarrPos==kk,:) = IxAux;%RxSymbAmosDo = [];
                                    EvmMatRef(ObsCarrPos==kk,:) = IxAux;                   %Taking just the middle samples as references
                                    %                     if (kk==1)&&(Ploting)
                                    %                         figure;hold all;grid on;
                                    %                         plot(IxAux,zeros(1,length(IxAux)),'rx');
                                    %                         set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                    %                         a=0;
                                    %                     end
                                    %%
                                    %                     if kk~=1
                                    if (~mod(CarrPass,2))&&(SendingUp)
                                        clear EoutModAux;
                                        EoutModAux = EoutTx(VetThisCarrTx==(RefCarr-1+kk+1),:);
                                        %                         if (kk==1)&&(Ploting)
                                        %                         figure;hold all;grid on;
                                        %                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
                                        %                         end
                                        if AddingNoiseO
                                            EsN0 = 10*log10(snr);
                                            [~,sigpower] = MeasPower(EoutModAux);
                                            EoutModAux = awgn(EoutModAux,osnr,sigpower-30);
                                            [~,sigpower] = MeasPower(EoutModAux);
                                        end
                                        %                         if (kk==1)&&(Ploting)
                                        %                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
                                        %                         PrintInfo(Ploting*4,abs(EoutModAux(kk,1:end-...
                                        %                         StuffSampels)).^2./max(abs(EoutModAux(kk,1:end-...
                                        %                         StuffSampels)).^2),t(2*NumAmosCP+NPPB),2*...
                                        %                         NumAmosCP+NPPB,abs(EoutModAux(1:end-StuffSampels...
                                        %                         )).^2./max(abs(EoutModAux(kk,1:end-StuffSampels)...
                                        %                                                                     ).^2));
                                        %                         end
                                        EoutModTem(ObsCarrPos(kk+1),1:length(EoutModAux)...
                                            ) = EoutModAux;%Adding the current Modulated output to the final OutPut
                                    end
                                    %                     a=1;
                                    %                     else
                                    %                     if ~mod(CarrPass,2)
                                    %                         EoutModTem(ObsCarrPos(kk+1),1:length(EoutModAux)) = 0;%Adding the current Modulated output to the final OutPut
                                    %                     end
                                    %                     end
                                else
                                    EoutModAux = EoutTx(VetThisCarrTx==(RefCarr-1+kk),:);
                                    %                     if (kk==1)&&(Ploting)
                                    %                         figure;hold all;grid on;
                                    %                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
                                    %                     end
                                    if AddingNoiseO
                                        EsN0 = 10*log10(snr);
                                        [~,sigpower] = MeasPower(EoutModAux);
                                        EoutModAux = awgn(EoutModAux,osnr,sigpower-30);
                                        [~,sigpower] = MeasPower(EoutModAux);
                                    end
                                    %                     if (kk==1)&&(Ploting)
                                    %                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
                                    %                         PrintInfo(Ploting*4,abs(EoutModAux(kk,1:end-...
                                    %                         StuffSampels)).^2./max(abs(EoutModAux(kk,1:end-...
                                    %                         StuffSampels)).^2),t(2*NumAmosCP+NPPB),2*...
                                    %                         NumAmosCP+NPPB,abs(EoutModAux(1:end-StuffSampels...
                                    %                         )).^2./max(abs(EoutModAux(kk,1:end-StuffSampels)...
                                    %                                                                     ).^2));
                                    %                         a=4;
                                    %                     end
                                    
                                    EoutModTem(ObsCarrPos(kk),1:length(EoutModAux)) = ...
                                        EoutModAux;%Adding the current Modulated output to the final OutPut
                                    if ~mod(CarrPass,2)
                                        clear EoutModAux;
                                        EoutModAux = EoutTx(VetThisCarrTx==(RefCarr-1+kk+1),:);
                                        %                         if (kk==1)&&(Ploting)
                                        %                         figure;hold all;grid on;
                                        %                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
                                        %                         end
                                        if AddingNoiseO
                                            EsN0 = 10*log10(snr);
                                            [~,sigpower] = MeasPower(EoutModAux);
                                            EoutModAux = awgn(EoutModAux,osnr,sigpower-30);
                                            [~,sigpower] = MeasPower(EoutModAux);
                                        end
                                        %                         if (kk==1)&&(Ploting)
                                        %                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
                                        %                         PrintInfo(Ploting*4,abs(EoutModAux(kk,1:end-...
                                        %                         StuffSampels)).^2./max(abs(EoutModAux(kk,1:end-...
                                        %                         StuffSampels)).^2),t(2*NumAmosCP+NPPB),2*...
                                        %                         NumAmosCP+NPPB,abs(EoutModAux(1:end-StuffSampels...
                                        %                         )).^2./max(abs(EoutModAux(kk,1:end-StuffSampels)...
                                        %                                                                     ).^2));
                                        %                         end
                                        EoutModTem(ObsCarrPos(kk+1),1:length(EoutModAux)...
                                            ) = EoutModAux;%Adding the current Modulated output to the final OutPut
                                    end
                                    %                     EoutModTem(ObsCarrPos(kk),1:length(EoutModAux)) = 0;%Adding the current Modulated output to the final OutPut
                                    %                     if ~mod(CarrPass,2)
                                    %                         EoutModTem(ObsCarrPos(kk+1),1:length(EoutModAux)) = 0;%Adding the current Modulated output to the final OutPut
                                    %                     end
                                end
                            end
                            if size(EoutModTem,1)<NumCarr
                                EoutModTem(end+1,:) = 0;
                            end
                            if IfftOrSum
                                %                 Ekaux=sum(EoutModTem);
                                %                 figure;hold all;
                                %                 plot(f,20*log10(abs(fftshift(fft(Ekaux)./length(Ekaux)))));
                                %                 set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                if size(EoutModTem,1)>1
                                    [EoutMod,~,~] = OpticalIFFTN(t,T,MaxNumStagT,EoutModTem);
                                else
                                    EoutMod = EoutModTem;
                                end
                                %                 plot(f,20*log10(abs(fftshift(fft(EoutMod)./length(...
                                %                                                               EoutMod)))));
                                %                 axis([0.9*fin 1.2*NumCarr*fc+fin -140 -20]);
                                %                 a=1;
                            else
                                if size(EoutModTem,1)>1
                                    EoutMod = sum(EoutModTem);
                                else
                                    EoutMod = EoutModTem;
                                end
                            end
                    end
                else
                    %There are two ways, implemented, to split the incoming OFCS. One by
                    %filters e another by the optical FFT. But this step just need to be done
                    %once in the simulation.
                    if ~exist('EoutTx','var')                                              %Verify is this step was previouly done
                        if SelecSetUP                                                      %Vefify is the split process will be by filter...
                            EoutTx=SelectEachCarrier(EoutT,NumCarr,f,fin,FBWD,Order,fc);   %This function is resposible to split each carrier
                            VetThisCarrTx = (RefCarr-1)+1:(RefCarr-1)+NumCarr;
                        else                                                               %... or by the OFFT
                            [EoutTx,~,VetThisCarrTx]=OpticalFFTN(t,T,MaxNumStagTx,EoutT);  %This function is resposible to split each carrier
                        end
                        %         PrintInfo(Ploting*51,EoutTx,f);                                %Printing for qualitative analizes.
                    end
                    EoutMod = 0;
                    if IfftOrSum
                        [EoutMod,~,~] = OpticalIFFTN(t,T,MaxNumStagT,EoutTx);
                    else
                        EoutMod = sum(EoutTx);
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
                
                if AddingNoiseF
                    EsN0 = 10*log10(snr);
                    [~,sigpower] = MeasPower(EoutRec);
                    EoutRec = awgn(EoutRec,osnrf,sigpower-30);
                    [~,sigpower] = MeasPower(EoutRec);
                end
                
                clear TxSigMat TxSig TxSigI TxSigQ EoutAux1 VetThisCarr EoutMod Ix...
                    EoutModAux Esync Esync1 Esync2 EoutA EoutB EoutC EoutD EoutModTem
                % clearvars -except fc Rb NPPB fsample T Ta NumberOf_T FinalTime Nb NbDPSK...
                % TotalSamples t f Eout EoutMod EoutRec TxDataMat NumCarr NbDQPSK ValsLev...
                % Modulation Polirized MaxAmp4PAM NumbOfTest UsedModula TestBundle Nb4Pam ...
                % ThisModula CurrentTest CurrentModula ThisTest ThisModula OfcName BerOOK ...
                % Ber4PAM BerDQPSK BerDPSK FiberDelay Medium CurrentMedium OSNRPC AddCP V0...
                % OcsToTest CurrentOCS PulseResp NumAmosCP NumBitDesc StuffSampels JusPos ...
                % contcp Ber4PAM PAM4Type Vpi MZ_Input_File ModAll JusValEnd NRZPolarity ...
                % VetElecPower AberLev fin FFTSplit PastTime PowRef SnrRef CarSNR ValsLev2...
                % AberLevI ValsLevI AberLevQ ValsLevQ VetElecPowerQ VetElecPowerI SigAten ...
                %                                      lambdac Atenuation VetOptiPower JusVal
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
                %                 DatumElectricReception;
                %%
                
                Ix = EoutAux.*NormFact;
                
                if Ploting
                    figure;hold all;grid on;plot(f,20*log10(abs(fftshift(fft(Ix)./length(Ix)))));axis([-25e9 37.5e9 -200 0]);set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                end
                Ix = Ix(1+NPSTUf:end);
                Ix = reshape(Ix,length(Ix)/NumFra,NumFra);
                Ix = intdump(Ix,NPPOF);                                        %Downsampling the income signal
                SigRecepA = Ix(1+NPOFEX:end,:);
                
                switch SelModTp
                    case 'BP'
                        %                     SigRecep  = intdump(SigRecepA,OvSam);
                        [BitFiltEle,~] = FiltroGaussiano(ff,OBw/2,0,1);        %Creating filter for selection of the received signal
                        %                     [ BitFiltEle ] = Filtro_Retangular(OBw,0,ff);
                        BitFiltEle = fftshift(BitFiltEle);
                        for jj=1:NumFra
                            ModFilta(:,jj) = BitFiltEle;
                        end
                        SigRecepB   = SigRecepA;
                        if Ploting
                            figure;
                            hold all;
                            plot(ff,20*log10(abs(fftshift(fft(SigRecepB)./length(SigRecepB)))));
                            plot(ff,20*log10(abs(fftshift(BitFiltEle))));
                        end
                        if SelModFilt
                            SigRecepB = ifft(fft(SigRecepB).*ModFilta);
                        end
                        SigRecep  = intdump(SigRecepB,OvSam);
                        if Ploting
                            plot(ff,20*log10(abs(fftshift(fft(SigRecepB)./length(SigRecepB)))));
                            set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                        end
                    case 'AM'
                        [BitFiltEle,~] = FiltroGaussiano(ff,OBw/2,0,1);           %Creating filter for selection of the received signal
                        BitFiltEle = fftshift(BitFiltEle);
                        for jj=1:NumFra
                            ModFilta(:,jj) = BitFiltEle;
                        end
                        SigRecepB   = SigRecepA.*cos(-2*pi*Ofc*tta);
                        if Ploting
                            figure;
                            hold all;
                            plot(ff,20*log10(abs(fftshift(fft(SigRecepB)./length(SigRecepB)))));
                            plot(ff,20*log10(abs(fftshift(BitFiltEle))));
                        end
                        if SelModFilt
                            SigRecepB = ifft(fft(SigRecepB).*ModFilta);
                        end
                        SigRecep  = intdump(SigRecepB,OvSam);
                        if Ploting
                            plot(ff,20*log10(abs(fftshift(fft(SigRecepB)./length(SigRecepB)))));
                            set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                        end
                        a=1;
                    case 'AMSSB'
                        [BitFiltEle,~] = FiltroGaussiano(ff,OBw/2,0,1);           %Creating filter for selection of the received signal
                        BitFiltEle = fftshift(BitFiltEle);
                        for jj=1:NumFra
                            ModFilta(:,jj) = BitFiltEle;
                        end
                        SigRecepB   = SigRecepA.*exp(-1j*2*pi*Ofc*tta);
                        if Ploting
                            figure;
                            hold all;
                            plot(ff,20*log10(abs(fftshift(fft(SigRecepB)./length(SigRecepB)))));
                            plot(ff,20*log10(abs(fftshift(BitFiltEle))));
                        end
                        if SelModFilt
                            SigRecepB = ifft(fft(SigRecepB).*ModFilta);
                        end
                        SigRecep  = intdump(SigRecepB,OvSam);              %Over sampling
                        if Ploting
                            plot(ff,20*log10(abs(fftshift(fft(SigRecepB)./length(SigRecepB)))));
                            set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                        end
                        a=1;
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
                if Ploting
                    switch OfdMod                                                  %Modulation Tx signal for future comparison
                        case 'qam'
                            TxDataA = TxDataMat(1,:);
                            TxDataA = reshape(TxDataA,NumFra,length(TxDataA)/NumFra);
                            TxDataA = TxDataA.';
                            for CarDmtM=1:length(DmtMve)
                                M = DmtMve(CarDmtM);
                                TxSigToPlot(CarDmtM,:)  = qammod(TxDataA(CarDmtM,:),M);            %Modulating information
                            end
                        otherwise
                            TxDataA = TxDataMat(1,:);
                            TxDataA = reshape(TxDataA,NumFra,length(TxDataA)/NumFra);
                            TxDataA = TxDataA.';
                            for CarDmtM=1:length(DmtMve)
                                M = DmtMve(CarDmtM);
                                TxSigToPlot(CarDmtM,:)  = dpskmod(TxDataA(CarDmtM,:),M);            %Modulating information
                            end
                            %                     SigRecepA = reshape(TxDataMat(ThisCarr,:),NumFra,NFFT);
                            %                     TxSigToPlot = dpskmod(TxDataA,M);                  %Modulating information
                            %                     TxSigToPlot = TxSigToPlot(:);
                    end
                    %         TxSigToPlot = dpskmod(TxDataMat(ThisCarr,2:end),M);
                    if Ploting
                        figure;
                        txcolor = [0.2 0 1];
                        rxcolor = [1 0.4 0];
                        hold all;
                        plot(TxSigToPlot(:),'*','LineWidth',2,'color',txcolor);
                        plot(SigRecep3(:),'o','LineWidth',2,'color',rxcolor);
                        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                    end
                end
                SaveRxNotEq = [SaveRxNotEq SigRecep3(:)];
                
                switch SelModTp
                    case 'BP'
                        switch OfdMod
                            case 'qam'
                                RxSigOfdmNoEq(CurrentTest,1,:) = SigRecep3(:).';
                                clear SigRecep4;
                                for CarDmtM=1:length(DmtMve)
                                    M = DmtMve(CarDmtM);
                                    SigRecep4(CarDmtM,:)  = qamdemod(SigRecep3(CarDmtM,:),M);            %Modulating information
                                end
                            otherwise
                                clear SigRecep4;
                                for CarDmtM=1:length(DmtMve)
                                    M = DmtMve(CarDmtM);
                                    SigRecep4(CarDmtM,:)  = dpskdemod(SigRecep3(CarDmtM,:),M);            %Modulating information
                                end
                        end
                    case 'AMSSB'
                        switch OfdMod
                            case 'qam'
                                RxSigOfdmNoEq(CurrentTest,1,:) = SigRecep3(:).';
                                equa = ChanelEqualizer(1,:);
                                equ = reshape(equa,length(equa)/NumFra,NumFra);
                                if ~ChanAwgn
                                    SigRecep3 = ((SigRecep3)./equ);%Fixing phase rotation and amplitude variation with the equalizer
                                end
                                clear SigRecep4;
                                for CarDmtM=1:length(DmtMve)
                                    M = DmtMve(CarDmtM);
                                    SigRecep4(CarDmtM,:)  = qamdemod(SigRecep3(CarDmtM,:),M);            %Modulating information
                                end
                            otherwise
                                clear SigRecep4;
                                for CarDmtM=1:length(DmtMve)
                                    M = DmtMve(CarDmtM);
                                    SigRecep4(CarDmtM,:)  = dpskdemod(SigRecep3(CarDmtM,:),M);            %Modulating information
                                end
                        end
                    otherwise
                        switch OfdMod
                            case 'qam'
                                RxSigOfdmNoEq(CurrentTest,1,:) = SigRecep3(:).';
                                equa = ChanelEqualizer(1,:);
                                equ = reshape(equa,length(equa)/NumFra,NumFra);
                                SigRecep3 = ((SigRecep3)./equ);%Fixing phase rotation and amplitude variation with the equalizer
                                clear SigRecep4;
                                for CarDmtM=1:length(DmtMve)
                                    M = DmtMve(CarDmtM);
                                    SigRecep4(CarDmtM,:)  = qamdemod(SigRecep3(CarDmtM,:),M);            %Modulating information
                                end
                            otherwise
                                clear SigRecep4;
                                for CarDmtM=1:length(DmtMve)
                                    M = DmtMve(CarDmtM);
                                    SigRecep4(CarDmtM,:)  = dpskdemod(SigRecep3(CarDmtM,:),M);            %Modulating information
                                end
                        end
                end
                EvmMatRec(1,:) = SigRecep3(:);%Taking just the middle samples as references
                [EvmDB(CurrentTest,1), EvmPer(CurrentTest,1), EvmRms(CurrentTest,1) ] = ...
                    EvmCalc( EvmMatRef(1,:),SigRecep3(:).' );
                [EvmDBJ(CurrentTest,1),EvmPerJ(CurrentTest,1),EvmRmsJ(CurrentTest,1)] = ...
                    evm1(M,OfdMod,EvmMatRef(1,:),SigRecep3(:).');
                %         SigRecep4 = dpskdemod(SigRecep3,M);
                SigRecep4 = SigRecep4.';
                SigRecep4 = SigRecep4(:).';
                RxSigOfdm(CurrentTest,1,:) = SigRecep4;
                SaveRxEq = [SaveRxEq SigRecep4];
                if Ploting
                    figure(111);
                    txcolor = [0.2 0 1];
                    rxcolor = [1 0.4 0];
                    hold all;
                    plot(TxSigToPlot(:),'*','LineWidth',2,'color',txcolor);
                    plot(SigRecep3(:),'o','LineWidth',2,'color',rxcolor);
                    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                    figure;
                    hold all;
                    plot(TxDataMat(1,:));
                    plot(SigRecep4);
                    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                end
                SymbErr = ~(TxDataMat(1,:)==SigRecep4);%Measuring system bit error ration
                DmtMvep = repmat(DmtMve.',1,NumFra);
                DmtKvep = log2(DmtMvep(:).');
                %             BerOFDM(CurrentTest,ThisCarr) = sum(SymbErr)/length(SymbErr);
                BerOFDM(CurrentTest,1) = sum(SymbErr.*DmtKvep)/sum(DmtKvep);
                % BerOFDM
                a=0;
                % close all;
            else
                if SendingDowStr
                    %                     if GatingAtReception
                    %                         DatumReceptionGating;
                    %                     else
                    %% DatumReception;
                    
                    %%                Main of All optical FFT
                    % The principal importance of the script is to load the All Optical FFT.
                    % Also, it is responsible to receive the transmited data.
                    
                    %% DatumReceptionInputData;                                                   % Configuration File of DatumReception.
                    
                    
                    
                    
                    
                    
                    
                    %%            Recovering the signal
                    % At this point the FFT will be inplemented. The received signal need to be
                    % given as an parameter. This following function recursively implement the
                    % FFT Operation. One may need to investigate how would work a receptor with
                    % filters insted of a OFFT. The second implementation is an array of
                    % filters.
                    if FFTSplit
                        [EoutAux1,~,VetThisCarr]=OpticalFFTN(t,T,MaxNumStag,EoutRec);
                        %     [EoutAux1,~,VetThisCarr]=OpticalFFTN(t,T,1,EoutRec);
                        %     for kk=1:size(EoutAux1,1)
                        %         [RecFilt,~]   = FiltroGaussiano(f,2*fc,(kk-1)*fc+fin,5);%Creating filter for conformation of the input information
                        %         RecFilt       = fftshift(RecFilt);
                        %         EoutAux1(kk,:) = ifft(fft(EoutAux1(kk,:)).*RecFilt);
                        %     end
                    else
                        VetThisCarr = ObsCarrPos;
                        if SendingUp%&&(NumCarr>2)
                            EoutAux1    = SelectEachCarrier(EoutRec,NumCarr,f,fin,1.0*fc,11,fc);
                        else
                            EoutAux1(1,:)    = EoutRec;
                            EoutAux1(2,:)    = 0;
                        end
                    end
                    for kk=1:size(EoutAux1,1)
                        [~,CarrRecPowDo(CurrentTest,kk)] = MeasPower(EoutAux1(VetThisCarr==...
                            ObsCarrUsed(kk),:),t);
                    end
                    %%  Ploting some results for qualitative analizes
                    % figure; plot(f,db(abs(fftshift(fft(EoutRec)./length(EoutRec))))); set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                    % PrintInfo(Ploting*12,f,db(abs((fftshift(fft(EoutRec))))/length(EoutRec)),EoutAux1,RefCarr*fc);
                    %%      Recovering the Data
                    %This is basicaly the final step, so far, in this process. Here, the
                    %transmited signal will be received and processed to recover the actual
                    %data within it.
                    switch Modulation
                        case 'OFDM'
                            for ThisCarr=InitCarr:CarrPass:NumCarr
                                EoutAux = EoutAux1(VetThisCarr==ObsCarrUsed(ThisCarr),:);
                                switch Medium
                                    case 'Fiber'
                                        EoutAux = ifft(fft(EoutAux).*exp(1j*2*pi*f*(...
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
                                Ix = EoutAux.*conj(EoutAux);
                                %             Ix = Ix.*NormFact;
                                [BitFilt,~] = FiltroGaussiano(f,BWD/2,CenFeq,FiltOrd);           %Creating filter for selection of the received signal
                                BitFilt = fftshift(BitFilt);                                   %Shifting the filter for matching the received signal
                                Ix = ifft(fft(Ix).*BitFilt);                                   %Filtering the signal for better detection
                                Ix = Ix - min(Ix);                                             %Removing any off-set that may exist
                                Ix = Ix./max(abs(Ix));                                         %Normalizing the eletrical signal (amplifying what is needed)
                                
                                if ReceptorNoise                                               %Verify whether the noise is to be added or not
                                    [~,SigPower] = MeasPower(Ix);
                                    SigPower2 = SigPower-30;%20*log10(SigPowerI);
                                    Ix = awgn(Ix,SNR,SigPower2);
                                end
                                if Ploting
                                    figure;hold all;grid on;plot(f,20*log10(abs(fftshift(fft(Ix)./length(Ix)))));axis([-25e9 37.5e9 -200 0]);set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                end
                                Ix = Ix(1+NPSTUf:end);
                                Ix = reshape(Ix,length(Ix)/NumFra,NumFra);
                                Ix = intdump(Ix,NPPOF);                                        %Downsampling the income signal
                                SigRecepA = Ix(1+NPOFEX:end,:);
                                %CpPaFr = 1;
                                %for CarrOffSet=-1*(NumFraPar-1)/2:1:(NumFraPar-1)/2
                                switch SelModTp
                                    case 'BP'
                                        %                     SigRecep  = intdump(SigRecepA,OvSam);
                                        if SelecGaus
                                            [ BitFiltEle ] = FiltroGaussiano(ff,OBw/2,0,1);%CarrOffSet*FramSpac,1);
                                        else
                                            [ BitFiltEle ] = Filtro_Retangular(OBw,0,ff);%CarrOffSet*FramSpac,ff);
                                        end
                                        BitFiltEle = fftshift(BitFiltEle);
                                        for jj=1:NumFra
                                            ModFilta(:,jj) = BitFiltEle;
                                        end
                                        SigRecepB   = SigRecepA;%.*exp(-1j*2*pi*CarrOffSet*FramSpac*tta);
                                        if Ploting
                                            figure;
                                            hold all;
                                            plot(ff,20*log10(abs(fftshift(fft(SigRecepB)./length(SigRecepB)))));
                                            plot(ff,20*log10(abs(fftshift(BitFiltEle))));
                                        end
                                        if SelModFilt
                                            SigRecepB = ifft(fft(SigRecepB).*ModFilta);
                                        end
                                        SigRecep  = intdump(SigRecepB,OvSam);
                                        if Ploting
                                            plot(ff,20*log10(abs(fftshift(fft(SigRecepB)./length(SigRecepB)))));
                                            set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                        end
                                    case 'AM'
                                        %                     tt = linspace(0,1*Te,size(SigRecepA,1));
                                        %                     ff = time2freq(tt);
                                        %                     for jj=1:NumFra
                                        %                         tta(:,jj) = tt;
                                        %                     end
                                        [BitFiltEle,~] = FiltroGaussiano(ff,OBw/2,0,1);           %Creating filter for selection of the received signal
                                        BitFiltEle = fftshift(BitFiltEle);
                                        for jj=1:NumFra
                                            ModFilta(:,jj) = BitFiltEle;
                                        end
                                        SigRecepB   = SigRecepA.*cos(2*pi*Ofc*tta);
                                        if Ploting
                                            figure;
                                            hold all;
                                            plot(ff,20*log10(abs(fftshift(fft(SigRecepB)./length(SigRecepB)))));
                                            plot(ff,20*log10(abs(fftshift(BitFiltEle))));
                                        end
                                        if SelModFilt
                                            SigRecepB = ifft(fft(SigRecepB).*ModFilta);
                                        end
                                        SigRecep  = intdump(SigRecepB,OvSam);
                                        if Ploting
                                            plot(ff,20*log10(abs(fftshift(fft(SigRecepB)./length(SigRecepB)))));
                                            set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                        end
                                        %                         x_ofdm = x_ofdm_ssb;
                                    case 'AMSSB'
                                        %                     tt = linspace(0,1*Te,size(SigRecepA,1));
                                        %                     ff = time2freq(tt);
                                        %                     for jj=1:NumFra
                                        %                         tta(:,jj) = tt;
                                        %                     end
                                        [BitFiltEle,~] = FiltroGaussiano(ff,OBw/2,0,1);           %Creating filter for selection of the received signal
                                        BitFiltEle = fftshift(BitFiltEle);
                                        for jj=1:NumFra
                                            ModFilta(:,jj) = BitFiltEle;
                                        end
                                        SigRecepB   = SigRecepA.*exp(1j*2*pi*Ofc*tta);
                                        if Ploting
                                            figure;
                                            hold all;
                                            plot(ff,20*log10(abs(fftshift(fft(SigRecepB)./length(SigRecepB)))));
                                            plot(ff,20*log10(abs(fftshift(BitFiltEle))));
                                        end
                                        if SelModFilt
                                            SigRecepB = ifft(fft(SigRecepB).*ModFilta);
                                        end
                                        SigRecep  = intdump(SigRecepB,OvSam);              %Over sampling
                                        if Ploting
                                            plot(ff,20*log10(abs(fftshift(fft(SigRecepB)./length(SigRecepB)))));
                                            set(gcf,'units','normalized','outerposition',[0 0 1 1]);
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
                                if Ploting
                                    switch OfdMod                                                  %Modulation Tx signal for future comparison
                                        case 'qam'
                                            TxDataA = TxDataMat(ThisCarr,:);
                                            TxDataA = reshape(TxDataA,NumFra,length(TxDataA)/NumFra);
                                            TxDataA = TxDataA.';
                                            for CarDmtM=1:length(DmtMve)
                                                M = DmtMve(CarDmtM);
                                                TxSigToPlot(CarDmtM,:)  = qammod(TxDataA(CarDmtM,:),M);            %Modulating information
                                            end
                                        otherwise
                                            TxDataA = TxDataMat(ThisCarr,:);
                                            TxDataA = reshape(TxDataA,NumFra,length(TxDataA)/NumFra);
                                            TxDataA = TxDataA.';
                                            for CarDmtM=1:length(DmtMve)
                                                M = DmtMve(CarDmtM);
                                                TxSigToPlot(CarDmtM,:)  = dpskmod(TxDataA(CarDmtM,:),M);            %Modulating information
                                            end
                                            %                     SigRecepA = reshape(TxDataMat(ThisCarr,:),NumFra,NFFT);
                                            %                     TxSigToPlot = dpskmod(TxDataA,M);                  %Modulating information
                                            %                     TxSigToPlot = TxSigToPlot(:);
                                    end
                                    %         TxSigToPlot = dpskmod(TxDataMat(ThisCarr,2:end),M);
                                    if Ploting
                                        figure;
                                        txcolor = [0.2 0 1];
                                        rxcolor = [1 0.4 0];
                                        hold all;
                                        plot(TxSigToPlot(:),'*','LineWidth',2,'color',txcolor);
                                        plot(SigRecep3(:),'o','LineWidth',2,'color',rxcolor);
                                        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                    end
                                end
                                
                                switch OfdMod
                                    case 'qam'
                                        RxSigOfdmNoEq(CurrentTest,ThisCarr,:) = SigRecep3(:).';
                                        equa = ChanelEqualizer(ObsCarrUsed(ThisCarr),:);
                                        equ = reshape(equa,length(equa)/NumFra,NumFra);
                                        SigRecep3 = ((SigRecep3)./equ);%Fixing phase rotation and amplitude variation with the equalizer
                                        clear SigRecep4;
                                        for CarDmtM=1:length(DmtMve)
                                            M = DmtMve(CarDmtM);
                                            SigRecep4(CarDmtM,:)  = qamdemod(SigRecep3(CarDmtM,:),M);            %Modulating information
                                        end
                                        %                     SigRecep4  = qamdemod(SigRecep3,M);                    %Modulating information
                                    otherwise
                                        clear SigRecep4;
                                        for CarDmtM=1:length(DmtMve)
                                            M = DmtMve(CarDmtM);
                                            SigRecep4(CarDmtM,:)  = dpskdemod(SigRecep3(CarDmtM,:),M);            %Modulating information
                                        end
                                        %                     SigRecep3 = reshape(SigRecep3,NumFra,NFFT);            %Reshaping the signal if multiple carriers were transmited
                                        %                     SigRecep4 = dpskdemod(SigRecep3,M);                    %Modulating information
                                        %                     SigRecep4 = SigRecep4(:);                              %Serializing the signal
                                end
                                EvmMatRec(ObsCarrPos==ThisCarr,:) = ...
                                    SigRecep3(:);%Taking just the middle samples as references
                                [EvmDB(CurrentTest,ThisCarr), ...
                                    EvmPer(CurrentTest,ThisCarr), ...
                                    EvmRms(CurrentTest,ThisCarr) ] = ...
                                    EvmCalc( EvmMatRef(ObsCarrPos==ThisCarr,:),SigRecep3(:).' );
                                [EvmDBJ(CurrentTest,ThisCarr),...
                                    EvmPerJ(CurrentTest,ThisCarr),...
                                    EvmRmsJ(CurrentTest,ThisCarr)] = ...
                                    evm1(M,OfdMod,EvmMatRef(ObsCarrPos==ThisCarr,:),SigRecep3(:).');
                                %         SigRecep4 = dpskdemod(SigRecep3,M);
                                SigRecep4 = SigRecep4.';
                                SigRecep4 = SigRecep4(:).';
                                RxSigOfdm(CurrentTest,ThisCarr,:) = SigRecep4;
                                if Ploting
                                    figure(111);
                                    txcolor = [0.2 0 1];
                                    rxcolor = [1 0.4 0];
                                    hold all;
                                    plot(TxSigToPlot(:),'*','LineWidth',2,'color',txcolor);
                                    plot(SigRecep3(:),'o','LineWidth',2,'color',rxcolor);
                                    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                end
                                if Ploting
                                    figure;
                                    hold all;
                                    plot(TxDataMat(ThisCarr,:));
                                    plot(SigRecep4);
                                    SampPos = 1:length(SigRecep4);
                                    plot(SampPos(~(TxDataMat(ThisCarr,:)==SigRecep4)),...
                                        SigRecep4(~(TxDataMat(ThisCarr,:)==SigRecep4)),'o');
                                    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                end
                                SymbErr = ~(TxDataMat(ThisCarr,:)==SigRecep4);%Measuring system bit error ration
                                %CpPaFr = CpPaFr + 1;
                                %end
                                DmtMvep = repmat(DmtMve,NumFra,1);
                                DmtKvep = log2(DmtMvep(:).');
                                BerToPlotOfdm(CurrentTest,ThisCarr,:) = SymbErr.*DmtKvep;
                                %             BerOFDM(CurrentTest,ThisCarr) = sum(SymbErr)/length(SymbErr);
                                BerOFDM(CurrentTest,ThisCarr) = sum(SymbErr.*DmtKvep)/sum(DmtKvep);
                                close all;
                            end
                        case 'DPSK'
                            %%                   Receiver DPSK
                            for ThisCarr=InitCarr:CarrPass:NumCarr                                    %For each carrier the same process of reception will be used.
                                if CarrUsedDo(ThisCarr)
                                    %%  Reception Process: Handling the income optical field
                                    %At the first step the income field will be selected from the
                                    %optical FFT Output with the help of the maping vector
                                    %(previously described).
                                    EoutAux = EoutAux1(VetThisCarr==ObsCarrUsed(ThisCarr),:);
                                    if AddingNoiseP
                                        [~,sigpower] = MeasPower(EoutAux);
                                        EoutAux      = awgn(EoutAux,osnrp,sigpower-30);
                                        [~,sigpower] = MeasPower(EoutAux);
                                    end
                                    %%            Fiber Time Delay Compensation
                                    switch Medium
                                        case 'Fiber'
                                            EoutAux = ifft(fft(EoutAux).*exp(1j*2*pi*f*(...
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
                                    [ESync1,ESync2] = DelayInterf(t,TimDel,PhaDel,EoutAux);        %Coverting phase shift to amplitude variation
                                    
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
                                    
                                    
                                    [BitFilt,~] = FiltroGaussiano(f,BWD2,CenFeq2,FiltOrd2);        %Creating filter for selection of the received signal
                                    BitFilt = fftshift(BitFilt);                                   %Shifting the filter for matching the received signal
                                    
                                    %%
                                    Esync  = ESync2-ESync1;
                                    
                                    if RecFilBanPas
                                        Esync  = ifft(fft(Esync).*BitFilt);                             %Filter is used to remove higher order components
                                        Esync  = Esync-mean(Esync);
                                        Esync  = Esync./max(abs(Esync));
                                    end
                                    %                 ESync1 = ESync1./max(abs(ESync1));                             %Normalizing the signal
                                    %                 ESync2 = ESync2./max(abs(ESync2));                             %Normalizing the signal
                                    %                 Esync  = Esync./max(abs(Esync));                               %Normalizing the signal
                                    
                                    %%         Adding Noise to Received Signal
                                    %Just after the photo diod is where the noise must be added as
                                    %it is the main constrain of the system. Once the signal is
                                    %converted from optical to electrical and it digital values
                                    %were recevered, there is no point to add noise, because the
                                    %information was already received, it just needs decoding.
                                    %The noise here added will be gaussian basically it represents
                                    %the reception sensibility. In another words, how much the
                                    %signal must be about the noise for an acceptable recovering.
                                    %                 if ReceptorNoise                                               %Verify whether the noise is to be added or not
                                    %                     [BitFiltN,~] = FiltroGaussiano(f,2*BWD2,CenFeq2,FiltOrd2);        %Creating filter for selection of the received signal
                                    %                     BitFiltN = fftshift(BitFiltN);                                   %Shifting the filter for matching the received signal
                                    %                     [~,SigPower] = MeasPower(Esync);                           %it own received signal or
                                    %                     SigPower2    = SigPower-30;%10*log10(SigPower);
                                    %                     SigAux       = linspace(SigPower2,SigPower2,length(t));
                                    %                     IrecNoise = awgn(SigAux,CarSNR,'measured');
                                    % %                     IrecNoise = ifft(fft(IrecNoise).*BitFiltN);                             %Filter is used to remove higher order components
                                    %                     Esync     = Esync + IrecNoise;
                                    %                 end
                                    if ReceptorNoise                                               %Verify whether the noise is to be added or not
                                        [~,SigPower] = MeasPower(Esync);
                                        SigPower2 = SigPower-30;%20*log10(SigPowerI);
                                        SNR2 = CarSNR + 10*log10(1) - 10*log10(Nsamp) + 10*log10(10^0.36);
                                        %                     if TimeSys
                                        %                         SNR = CarSNR + 10*log10(1) - 10*log10(1*T*BWD2);
                                        %                     else
                                        %                         SNR = CarSNR + 10*log10(1) - 10*log10(1*t(2*NumAmosCP+NPPB)*BWD2);
                                        %                     end
                                        Esync = awgn(Esync,SNR,SigPower2);
                                        %                     Esync = ifft(fft(Esync).*BitFiltN);                             %Filter is used to remove higher order components
                                    end
                                    if ~RecFilBanPas
                                        Esync  = ifft(fft(Esync).*BitFilt);                             %Filter is used to remove higher order components
                                        Esync  = Esync-mean(Esync);
                                        Esync  = Esync./max(abs(Esync));
                                    end
                                    
                                    VetElecPower(CurrentTest,ThisCarr)= MeasPower(Esync);
                                    EoutAuxF = 20*log10(abs(fftshift(fft(EoutAux)./length(...
                                        EoutAux))));
                                    [VetOptiPower(CurrentTest,ThisCarr),~]= findpeaks(EoutAuxF,...
                                        'SortStr','descend','NPeaks',1);
                                    
                                    SyncAux   = Esync(IniSyncPos:SyncPos);                         %Selecting just the symbol to synchronize
                                    SyncedAux = SyncSymb(IniSyncPos:SyncPos);                      %Selecting just the symbol to synchronize
                                    
                                    %%                   Plot for Qualitative analizes
                                    %             PrintInfo(Ploting*15,t(IniSyncPos:SyncPos),Esync(IniSyncPos:...
                                    %             SyncPos),SyncSymb(IniSyncPos:SyncPos),ESync1(IniSyncPos:...
                                    %                                       SyncPos),ESync2(IniSyncPos:SyncPos));
                                    
                                    %%                   Synchronizing
                                    %This synchronization process is based on the mean expected
                                    %value. That means, the information mean value should be within
                                    %one period of symbol. Thus, the mean value of the received
                                    %signal is acquired and compare of the known sync-word to
                                    %verify if this mean value is at the right possition. Which is
                                    %the midel point (peak) of the highest level at the sync period
                                    SyncAux(SyncAux<0)              = 0;                           %To keep the mean value above zero anything under is neglected
                                    SyncAux(SyncAux>=mean(SyncAux)) = 1;                           %Adding a flag to the first sample of the received mean value
                                    SyncAux(SyncAux<mean(SyncAux))  = -1;                          %All the others samples at set to the lowest level
                                    
                                    PosToSyn  = find(ismember(SyncAux,1));                         %Finding where is the location of the samples to synchronize
                                    PosSynced = find(ismember(SyncedAux,1));                       %Finding where is the location of the samples to synchronize
                                    
                                    DiffPos = PosToSyn(round(end/2)) - PosSynced(round(end/2));    %Accounting the peak (midel point) displacement
                                    sn=((1/2)*(((1/2)^MaxNumStag)-1))/((1/2)-1);
                                    DiffPos = round((sn/(3/2.55^IfftOrSum))*T/Ta);
                                    if DiffPos~=0
                                        if DiffPos>0%If the difference is positive, left-shift...
                                            %                     Esync = ifft(fft(Esync).*exp(1j*2*pi*f*(DiffPos*Ta))); %Shift based on time change
                                            Esync = [Esync(DiffPos+1:end) Esync(1:DiffPos)];       %Shift based on sampling sliding
                                        else%... but if the difference is negative, right-shift
                                            %                     Esync = ifft(fft(Esync).*exp(1j*2*pi*f*(DiffPos*Ta))); %Shift based on time change
                                            Esync = [Esync(end+DiffPos+1:end) Esync(1:end + ...
                                                DiffPos)];%Shift based on sampling sliding
                                        end
                                    end
                                    
                                    %Because of reasons, sometimes it may be required to make a
                                    %synchronization process with the  end of the data stream as
                                    %well. This following verification check if the user set (or
                                    %not) a second synchronization process to be done.
                                    
                                    %%          Ploting the result for qualitative analizes
                                    %             PrintInfo(Ploting*16,t(end-SyncPos+1:end-IniSyncPos+1),Esync...
                                    %             (end-SyncPos+1:end-IniSyncPos+1),SyncSymb(end-SyncPos+1:end-...
                                    %             IniSyncPos+1),ESync1(end-SyncPos+1:end-IniSyncPos+1),ESync2(...
                                    %                                           end-SyncPos+1:end-IniSyncPos+1));
                                    if SencondAdjust
                                        %%                   Synchronizing
                                        %This synchronization process is based on the mean expected
                                        %value. That means, the information mean value should be
                                        %within one period of symbol. Thus, the mean value of the
                                        %received signal is acquired and compare of the known
                                        %sync-word to verify if this mean value is at the right
                                        %possition.
                                        
                                        SyncAuxEnd   = Esync(end-SyncPos+1:end-IniSyncPos+1);
                                        SyncedAuxEnd = SyncSymbEnd(end-SyncPos+1:end-IniSyncPos+1);
                                        SyncAuxEnd(SyncAuxEnd<0)                 = 0;              %To keep the mean value above zero anything under is neglected
                                        SyncAuxEnd(SyncAuxEnd>=mean(SyncAuxEnd)) = 1;              %Adding a flag to the first sample of the received mean value
                                        SyncAuxEnd(SyncAuxEnd<mean(SyncAuxEnd))  = -1;             %All the others samples at set to the lowest level
                                        
                                        PosToSynEnd  = find(ismember(SyncAuxEnd,1));               %Finding where is the location of the first sample to synchronize
                                        PosSyncedEnd = find(ismember(SyncedAuxEnd,1));             %Finding where is the location of the first sample to synchronize
                                        
                                        DiffPosEnd = PosToSynEnd(round(end/2))-PosSyncedEnd(...
                                            round(end/2));
                                        if DiffPosEnd~=0
                                            if DiffPosEnd>0%If positive difference, left-shift...
                                                %                         Esync = ifft(fft(Esync).*exp(1j*2*pi*f*(...
                                                %                                                           DiffPosEnd*Ta)));%Shift based on time change
                                                Esync = [Esync(DiffPosEnd+1:end) Esync(1:...
                                                    DiffPosEnd)];%Shift based on sampling sliding
                                            else%... but if the difference is negative, right-shift
                                                %                         Esync = ifft(fft(Esync).*exp(1j*2*pi*f*(...
                                                %                                                           DiffPosEnd*Ta)));%Shift based on time change
                                                Esync = [Esync(end-DiffPosEnd+1:end) Esync(...
                                                    1:end-DiffPosEnd)];%Shift based on sampling sliding
                                            end
                                        end
                                    end
                                    
                                    %% Removing CP
                                    if AddCP
                                        IxAux  = Esync(1:end - StuffSampels);
                                        IxAux  = reshape(IxAux,(2*NumAmosCP+NPPB),NbDPSK);
                                        IxAux  = IxAux(1+NumAmosCP:end-NumAmosCP,:);
                                        IxAux  = reshape(IxAux,1,NPPB*NbDPSK);
                                        Esync  = IxAux;
                                    end
                                    %                 PassGate = round(Nb*NPPB/round(t(end)*ThisCarr*fc*4));
                                    %                 VectGate = 1:PassGate:(Nb*NPPB);
                                    %                 RecSig   = Esync(VectGate);
                                    % %                 RecSig   = RecSig./max(abs(RecSig));
                                    %                 RecF     = f(VectGate);
                                    %                 figure; hold all;
                                    %                 plot(f,20*log10(abs(fftshift(fft(Esync)./length(Esync)))));
                                    %                 plot(RecF,20*log10(abs(fftshift(fft(RecF)./length(RecF)))));
                                    %                  set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                    %                 figure; hold all;
                                    %                 plot(t,Esync);
                                    %                 plot(t(VectGate),RecSig);
                                    %                  set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                    %% Taking the sampling the EVM meassurement
                                    clear IxAux;
                                    PosAuxEout1 = NPPB/2:NPPB:length(Esync);                   %Varriable respossible to take just the samples at the middle of the symbol
                                    PosAuxEout2 = ((NPPB/2)+(NPPB/16)):NPPB:length(Esync);      %Varriable respossible to take just the samples at the middle of the symbol
                                    PosAuxEout3 = ((NPPB/2)-(NPPB/16)):NPPB:length(Esync);      %Varriable respossible to take just the samples at the middle of the symbol
                                    IxAux1      = Esync(PosAuxEout1);                          %Normalizing the reference
                                    IxAux2      = Esync(PosAuxEout2);    %Normalizing the reference
                                    IxAux3      = Esync(PosAuxEout3);    %Normalizing the reference
                                    a=0;
                                    EvmMatRecA(1,:) = IxAux1;                       %Taking just the middle samples as references
                                    EvmMatRecA(2,:) = IxAux2;                       %Taking just the middle samples as references
                                    EvmMatRecA(3,:) = IxAux3;                       %Taking just the middle samples as references
                                    EvmMatRecA(4,:) = IxAux1;                       %Taking just the middle samples as references
                                    [EvmDBA(1), EvmPerA(1), EvmRmsA(1) ] = EvmCalc( EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux1 );
                                    [EvmDBA(2), EvmPerA(2), EvmRmsA(2) ] = EvmCalc( EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux2 );
                                    [EvmDBA(3), EvmPerA(3), EvmRmsA(3) ] = EvmCalc( EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux3 );
                                    [EvmDBA(4), EvmPerA(4), EvmRmsA(4) ] = EvmCalc( EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux1 );
                                    [EvmDBJA(1),EvmPerJA(1),EvmRmsJA(1)] = evm1(2,'dpsk',EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux1);
                                    [EvmDBJA(2),EvmPerJA(2),EvmRmsJA(2)] = evm1(2,'dpsk',EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux2);
                                    [EvmDBJA(3),EvmPerJA(3),EvmRmsJA(3)] = evm1(2,'dpsk',EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux3);
                                    [EvmDBJA(4),EvmPerJA(4),EvmRmsJA(4)] = evm1(2,'dpsk',EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux1);
                                    %                 if ThisCarr==1
                                    %                     figure;hold all;grid on;
                                    %                     plot(real(IxAux1),zeros(1,length(IxAux1)),'rx');
                                    %                     plot(real(IxAux2),zeros(1,length(IxAux2)),'yx');
                                    %                     plot(real(IxAux3),zeros(1,length(IxAux3)),'g^');
                                    %                     set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                    %                     a=0;
                                    %                 end
                                    %##########################################################################
                                    %                 Esync = Esync(1+SyncPeriod:end-SyncPeriod);
                                    PosIx = NPPB/2:NPPB:length(Esync);                                %Possition of the central samples - but the number of samples per symbol is even ????
                                    IxAux = Esync(PosIx);                                             %From the main received signal just a few samples are taken for further evaluation
                                    n     = 100;                                                   %The number of boxes to be filled up on the histogram process
                                    Perce = 0.7;
                                    %At this process each eye will be evaluated separetelly. It is
                                    %important to remember that the PAM4 has 4 levels, which means
                                    %three level of decissions that are exaclty the center of the
                                    %eye diagram.
                                    IxAuxAB = Esync(PosIx);       %Taking just those values relative to the uper eye
                                    InterAB = linspace(Perce*(min(Esync)),Perce*(max(Esync)),n);                     %Building the histogram boxes
                                    EyeAB = hist(IxAuxAB,InterAB);                                 %filling up the boxes with samples that fit on them.
                                    
                                    %The same process described for the uper level will be done at
                                    %the middle and lower eyes levels.
                                    IxAuxCD = Esync(PosIx-NPPB/16);
                                    InterCD = linspace(Perce*(min(Esync)),Perce*(max(Esync)),n);%NPPB*2^n);
                                    EyeCD = hist(IxAuxCD,InterCD);
                                    
                                    IxAuxEF = Esync(PosIx+NPPB/16);
                                    InterEF = linspace(Perce*(min(Esync)),Perce*(max(Esync)),n);%NPPB*2^n);
                                    EyeEF = hist(IxAuxEF,InterEF);
                                    
                                    %What we are looking for here are not where exist the
                                    %occurrences of values. The eye is where there is not samples,
                                    %this means, the decission level is a reagion between actual
                                    %levels thus this reagion does not contain any signal.
                                    EyeAB = ~EyeAB;                                                %Changing zeros to one - Zeros compose the eye region
                                    EyeCD = ~EyeCD;
                                    EyeEF = ~EyeEF;
                                    %Starting variables that will be used to identify the eye
                                    %diagram.
                                    CountAB   = 1;
                                    CountCD   = 1;
                                    CountEF   = 1;
                                    
                                    SeqOnesAB = 0;
                                    SeqOnesCD = 0;
                                    SeqOnesEF = 0;
                                    
                                    SeqFinAB  = 0;
                                    SeqFinCD  = 0;
                                    SeqFinEF  = 0;
                                    
                                    SeqIniAB  = 1;
                                    SeqIniCD  = 1;
                                    SeqIniEF  = 1;
                                    %The for loop will take account of every box with ones. It is
                                    %important to take note that the not operator was used in this
                                    %vector, therefore ones means zeros (the eye diagram -
                                    %possibly) and zeros means values abouve zeroa (not the eye).
                                    for kk=1:length(EyeAB)                                         %For every box
                                        if EyeAB(kk)                                               %if it contains "1"
                                            SeqOnesAB(SeqIniAB)=CountAB;                           %count this element as part of a consecutive sequency
                                            CountAB = CountAB + 1;                                 %adds one to the counter of consecutive elements "1"
                                            if kk==length(EyeAB)                                   %if the current box is the last box we got to an end
                                                SeqFinAB(SeqIniAB) = kk;                           %The final sequency element is equal to its possition (kk)
                                            end
                                        else                                                       %else if the current box contains "0"
                                            SeqFinAB(SeqIniAB) = kk-1;                             %the previous element was the last element of a consecutive sequency of ones
                                            SeqIniAB = SeqIniAB + 1;                               %adds one to the counter of consecutive number (take in account how many sequencies there are)
                                            CountAB = 1;                                           %reset the counter of consecutive elements to mark the start of a new consecutive sequence
                                        end
                                        
                                        if EyeCD(kk)
                                            SeqOnesCD(SeqIniCD)=CountCD;
                                            CountCD = CountCD + 1;
                                            if kk==length(EyeCD)
                                                SeqFinCD(SeqIniCD) = kk;
                                            end
                                        else
                                            SeqFinCD(SeqIniCD) = kk-1;
                                            SeqIniCD = SeqIniCD + 1;
                                            CountCD = 1;
                                        end
                                        
                                        if EyeEF(kk)
                                            SeqOnesEF(SeqIniEF)=CountEF;
                                            CountEF = CountEF + 1;
                                            if kk==length(EyeEF)
                                                SeqFinEF(SeqIniEF) = kk;
                                            end
                                        else
                                            SeqFinEF(SeqIniEF) = kk-1;
                                            SeqIniEF = SeqIniEF + 1;
                                            CountEF = 1;
                                        end
                                    end
                                    
                                    
                                    %If the eye is open, which means there is a clear difference
                                    %between adjacent levels, the eye diagram will be the longest
                                    %sequence of ones.
                                    [MaxValAB,LocMaxAB]=max(SeqOnesAB);
                                    if LocMaxAB<2 || MaxValAB<2                                    %if any sequency was found or there is just one sequency it is a error thus
                                        LevDec3 = 0.0;                                            %the decission level will be by default 0.67. Also, the other variables will
                                        LocMaxAB = 1;                                              %will be set with values to not cause errors in the future
                                        SeqFinAB(1)=2;
                                        MaxValAB = 0;
                                        InterAB(1)=LevDec3;
                                    else                                                           %if a sequency was found the middle of the eye will be the middle of the sequency
                                        if (SeqFinAB(LocMaxAB)-MaxValAB/2)<1                       %if for some reason the final element of the sequency minus the half of its
                                            LevDec3 = 0.0;                                         %length results in a negative value, something went very wrong, and by
                                        else                                                       %default it will be set to 0.7
                                            LevDec3 = InterAB(round(SeqFinAB(LocMaxAB)-MaxValAB/...
                                                2));%Otherwise, the decission level is the middle of the sequency
                                        end
                                    end
                                    
                                    [MaxValCD,LocMaxCD]=max(SeqOnesCD);
                                    if LocMaxCD<2 || MaxValCD<2
                                        LevDec2 = 0.0;
                                        LocMaxCD = 1;
                                        SeqFinCD(1)=2;
                                        MaxValCD = 0;
                                        InterCD(1)=LevDec2;
                                    else
                                        if (SeqFinCD(LocMaxCD)-MaxValCD/2)<1
                                            LevDec2 = 0.0;
                                        else
                                            LevDec2 = InterCD(round(SeqFinCD(LocMaxCD)-MaxValCD/2));
                                        end
                                    end
                                    
                                    [MaxValEF,LocMaxEF]=max(SeqOnesEF);
                                    if LocMaxEF<2 || MaxValEF<2
                                        LevDec1 = 0.0;
                                        LocMaxEF = 1;
                                        SeqFinEF(1)=2;
                                        MaxValEF = 0;
                                        InterEF(1)=LevDec1;
                                    else
                                        if (SeqFinEF(LocMaxEF)-MaxValEF/2)<1
                                            LevDec1 = 0.0;
                                        else
                                            LevDec1 = InterEF(round(SeqFinEF(LocMaxEF)-MaxValEF/2));
                                        end
                                        
                                    end
                                    
                                    %another way to measure the eye opening is the get all the
                                    %boxes and find all peaks on it, that will be a plato created
                                    %by the sequences of ones (which was zeros). From thos peaks,
                                    %the eye diagram will be the longer of them hence it will take
                                    %the most part of the vector that store the findpeaks result.
                                    %Thus, the middle of the eye will be basically the middle of
                                    %the peaks vector.
                                    LocAB = find(EyeAB);
                                    if isempty(LocAB)                                              %if for some reason there are no peaks, something went wrong.
                                        LocAB = 1;
                                        LevelDec3 = 0.0;%mean(Levels(3:4));                       %by default the decission level will be set to 0.65
                                    else
                                        LevelDec3 = LevDec3;
                                    end
                                    
                                    LocCD = find(EyeCD);
                                    if isempty(LocCD)
                                        LocCD = 1;
                                        LevelDec2 = 0.0;%mean(Levels(2:3));
                                    else
                                        LevelDec2 = LevDec2;
                                    end
                                    
                                    LocEF = find(EyeEF);
                                    if isempty(LocEF)
                                        LocEF = 1;
                                        LevelDec1 = 0.0;%mean(Levels(1:2));
                                    else
                                        LevelDec1 = LevDec1;
                                    end
                                    %##########################################################################
                                    if CalcS
                                        [~,~,EyeOpenI,EyeOpenHighI,EyeOpenLowI] = Olho(Esync,T,NPPB,0);
                                        AberLevS(CurrentTest,ThisCarr)  = EyeOpenI;
                                        ValsLevS(CurrentTest,ThisCarr)  = EyeOpenLowI + EyeOpenI/2;          %Comparison between the Transmited and received and counting the differences
                                    end
                                    %                 [~,~,EyeOpenQ,EyeOpenHighQ,EyeOpenLowQ] = Olho(EoutQ,T,NPPB,0);
                                    AberLevAuxI(1)= InterAB(SeqFinAB(LocMaxAB))-InterAB(round(SeqFinAB(LocMaxAB)-MaxValAB));
                                    AberLevAuxI(3)= InterCD(SeqFinCD(LocMaxCD))-InterCD(round(SeqFinCD(LocMaxCD)-MaxValCD));
                                    AberLevAuxI(2)= InterEF(SeqFinEF(LocMaxEF))-InterEF(round(SeqFinEF(LocMaxEF)-MaxValEF));
                                    ValsLevAuxI(1)= LevDec3;
                                    ValsLevAuxI(3)= LevDec2;
                                    ValsLevAuxI(2)= LevDec1;
                                    if ThisCarr==126
                                        EyeToPlot(CurrentTest,1:length(Esync)) = Esync;
                                        save(['savingforplotingeye' num2str(VetSnr)],'EyeToPlot','VetSnr','IfftOrSum','UsedModula','T','NPPB');
                                    end
                                    %% Ploting the result for qualitative analizes
                                    if PrintinEye
                                        PrintInfo(Ploting*17,Esync,T,NPPB);
                                        hold on;
                                        plot(t((NPPB)/2),InterAB(round(SeqFinAB(LocMaxAB)-MaxValAB)),'k^');
                                        plot(t((NPPB)/2),InterAB(SeqFinAB(LocMaxAB)),'kv');
                                        plot(t((NPPB)/2),InterAB(round(SeqFinAB(LocMaxAB)-MaxValAB/2)),'k*');
                                        plot(t(((NPPB)/2)-(NPPB)/16),InterCD(round(SeqFinCD(LocMaxCD)-MaxValCD)),'g^');
                                        plot(t(((NPPB)/2)-(NPPB)/16),InterCD(SeqFinCD(LocMaxCD)),'gv');
                                        plot(t(((NPPB)/2)-(NPPB)/16),InterCD(round(SeqFinCD(LocMaxCD)-MaxValCD/2)),'g*');
                                        plot(t(((NPPB)/2)+(NPPB)/16),InterEF(round(SeqFinEF(LocMaxEF)-MaxValEF)),'b^');
                                        plot(t(((NPPB)/2)+(NPPB)/16),InterEF(SeqFinEF(LocMaxEF)),'bv');
                                        plot(t(((NPPB)/2)+(NPPB)/16),InterEF(round(SeqFinEF(LocMaxEF)-MaxValEF/2)),'b*');
                                        if CalcS
                                            plot(t((NPPB)/2),EyeOpenLowI,'mx');
                                            plot(t((NPPB)/2),EyeOpenHighI,'mo');
                                            plot(t((NPPB)/2),EyeOpenLowI + EyeOpenI/2,'md');
                                        end
                                        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                        figure;
                                        hold all;
                                        plotpos = zeros(1,length(IxAux1));
                                        plot(IxAux1,plotpos,'o','color',[1 0.4 0]);
                                        plot(EvmMatRef(ObsCarrUsed(ThisCarr),:),plotpos,'b*');
                                        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                        a=0;
                                    end
                                    %                 [~,~,EyeOpen,EyeOpenHigh,EyeOpenLow] = Olho(Esync,T,NPPB,0);
                                    
                                    %                 AberLev(CurrentTest,ThisCarr)= Inter(SeqFin(LocMax))-Inter(round(SeqFin(LocMax)-MaxVal));
                                    %                 ValsLev(CurrentTest,ThisCarr)= LevDec;
                                    %                 AberLev(CurrentTest,ThisCarr)  = EyeOpen;
                                    %                 ValsLev(CurrentTest,ThisCarr)  = EyeOpenLow + EyeOpen/2;
                                    
                                    %                 Txaux1 = rectpulse(TxDataMat(ThisCarr,:),NPPB);
                                    %                 Txaux1(Txaux1==0) = -1;
                                    %% Ploting the result for qualitative analizes
                                    %                 hold on;
                                    %                 plot(t((NPPB)/2),Inter(round(SeqFin(LocMax)-MaxVal)),'b^');
                                    %                 plot(t((NPPB)/2),Inter(SeqFin(LocMax)),'bv');
                                    %                 plot(t((NPPB)/2),Inter(round(SeqFin(LocMax)-MaxVal/2)),'b*');
                                    %                 plot(t((NPPB)/2),EyeOpenLow,'cx');
                                    %                 plot(t((NPPB)/2),EyeOpenHigh,'co');
                                    % %                 plot(t((NPPB)/2),EyeOpenLow + EyeOpen/2,'cd');
                                    %                 a=0;
                                    %                 PrintInfo(Ploting*18,t(1:length(Esync)),Txaux1,Esync);
                                    %%               Recovering the Information
                                    %Once the Income signal was synchronized it is possible to
                                    %recover the signal.
                                    %
                                    %After passing the optical signal to the eletrical domain, for
                                    %actually detect the data withing the signal the following
                                    %steps are needed.
                                    %
                                    %Finding Decission Levels:
                                    %The process for decoding the income signal will be based on
                                    %eletronic comparators. Inasmuch as the right decission level
                                    %must be acquired for accurately decide, within a symbol
                                    %periode, what that current leavel means (ones or zeros).
                                    %
                                    %The process hereafter of chosing the  decission levels is not
                                    %deterministic rather it is a statistic process. The main idea
                                    %is to take the decission level from the histogram generated
                                    %from the income signal stream.
                                    %
                                    %This process is realized inside the function Olho.
                                    %
                                    %Basicaly the decission level will be the minimal value of the
                                    %currente eye under evaluation plus the half of the its eye
                                    %opening.The following ilustration better describe this process
                                    %
                                    %Eye Limits:   + (1/2 * Eye Opening:)  =    Comparisson Limit:
                                    %
                                    %UperLevel  ______     ______________     _____
                                    %                 \   /      |       \   /
                                    %                  \ /       |        \ /
                                    %                   \ Half Eye Opening /     Decission Level
                                    %                  / \       |        / \
                                    %LowerLevel ______/   \______|_______/   \_____
                                    %
                                    %
                                    %Once the signal was processed the next step is through a
                                    %comparator decide the actual information received.
                                    %
                                    %For accurately read the income signal is necessary more than
                                    %just have the decision level, it is also needed to know where
                                    %it occurs. For instance, if we decide to use just to take the
                                    %decision level and measure the mean value across and a portion
                                    %of the time period, which portion should we take? It would be
                                    %logical to take the central portion, but the bit may be
                                    %deformed in such way that the information gets concentrated
                                    %out of the middle part, which means the symbol is not
                                    %symmetric. The symmetry information can be acquired from the
                                    %eye diagram by measuring the longitudinal opening. The
                                    %following sketch better describes this process:
                                    %
                                    %                   Point of Symmetry
                                    %         _____     _______|_______     _____
                                    %              \   /               \   /
                                    %               \ /  Longitudinal   \ /
                                    %                \ ----------------- /
                                    %               / \    Opening      / \
                                    %         _____/   \_______________/   \_____
                                    %
                                    %With those two pieces of information, decision level and point
                                    %of symmetry, we have the X, Y coordinates for the centre of
                                    %the Eye Diagram. Therefore, as long as there is an opening on
                                    %it it will be possible to recover the transmitted information
                                    %without error... theoretically.
                                    %
                                    %As this process is also statistical, first we reshape the
                                    %income vector to analyze all periods at the same time.
                                    %                 EyeSymMat = reshape(Esync(1+SyncPeriod*NPPB:end-SyncPeriod*...
                                    %                     NPPB),NPPB,NbDPSK-2*SyncPeriod);
                                    %Then we take the values that compose the decision level
                                    %because they will mark the point of symmetry.
                                    %
                                    %Firstly it was set the interval in which the histogram will be
                                    %build. It is based on the number of samples per bit period.
                                    %                 Interval = linspace(min(Esync(1+SyncPeriod*NPPB:end-...
                                    %                     SyncPeriod*NPPB)),max(Esync(1+SyncPeriod*NPPB:end-SyncPeriod...
                                    %                     *NPPB)),2*NPPB);
                                    %Therefore, the MATLAB hist function returns the number of
                                    %occurrence of each interval.
                                    %                 EyeMax = hist(Esync,Interval);
                                    %                 EyeMaxaux = [0 EyeMax 0];                                      %Zeros are added at the EyeMax to auxiliate the finding peaks process
                                    %                 [~,EyeLoc] = findpeaks(EyeMaxaux,'MinPeakDistance',NPPB/4,...
                                    %                     'SortStr','descend','NPeaks',2);%The peaks on the Eye profile will be the levels at the Eyes limit
                                    
                                    %From the location of the max values that occured, which means
                                    %the uper and lower level of the eye diagram it needs to take
                                    %the actual value that those occurences represent that is
                                    %withing the the Interval variable.
                                    %                 ValToSeek = Interval(EyeLoc-1);
                                    %The number of ocurrences is a statical measure therefore one
                                    %does not have control which interval will have the highest
                                    %peak, thus it is important to ordenate the values to be seek
                                    %from the lower part of the eye diagram to the uper part of the
                                    %eye diagram.
                                    %                 ValToSeek = sort(ValToSeek,'ascend');
                                    %                 OccuCount = zeros(1,size(EyeSymMat,1));                        %Auxiliar Variable for accounting.
                                    %                 for kk=1:size(EyeSymMat,1)                                     %For every sample within a symbol period
                                    
                                    %                     OccuCount(kk) = OccuCount(kk)+sum((EyeSymMat(kk,:)>=min(...
                                    %                         Esync))&(EyeSymMat(kk,:)<=UpeSymPer*EyeOpenLow)); %Account all occurencies of the value 1
                                    %                     OccuCount(kk) = OccuCount(kk)+sum((EyeSymMat(kk,:)>=...
                                    %                         LowSymPer*EyeOpenHigh)&(EyeSymMat(kk,:)<=max(Esync))); %Account all occurencies of the value 2
                                    %                 end
                                    %The point of symmetry of the eye diagram will be where the
                                    %maximum number of occurrences were measured inasmuch as those
                                    %are the points where all the bits go to the center of the
                                    %symbol.  From the maximum number of occurrences, it can happen
                                    %for more than one sample within one symbol period, in the best
                                    %case, all samples would have the same accounting as it is
                                    %shown the ilustration above hence the symmetry will be at the
                                    %middle sample of this group of maximum occurrences. This value
                                    %can be found by the mean of the samples positions within the
                                    %symbol period. The problem with this approach is that the
                                    %signal must be synchronized with the maximum displacement of
                                    %a symbol period minus 25% of the eye Longitudinal opening if
                                    %the displacement is higher than that the point of symmetry
                                    %will be wrongly measured.
                                    %                 [SymLoc] = round(mean(find(ismember(OccuCount,max(OccuCount)...
                                    %                     ))));
                                    %             figure;findpeaks(OccuCount,'SortStr','descend');
                                    
                                    %% Actualy Receiving Data:
                                    Data = [];                                                  %Initialization of the vector that will store the income data
                                    DataM = [];                                                  %Initialization of the vector that will store the income data
                                    DataL = [];                                                  %Initialization of the vector that will store the income data
                                    DataU = [];                                                  %Initialization of the vector that will store the income data
                                    DataS = [];                                                  %Initialization of the vector that will store the income data
                                    for kk=1:NPPB:length(Esync)                                    %The comparison process will be made for each symbol period
                                        %                 MeanOfData = mean(EoutI((kk-1)+SymLocI(1)));
                                        MeanOfData = mean(Esync((kk-1)+NPPB/2));
                                        if MeanOfData > LevelDec3%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                                            Data = [Data 1];                                 %is 1
                                        else                                                       %If it is the lowest level the incoming data
                                            Data = [Data 0];                                 %is 0
                                        end
                                        MeanOfData = mean(Esync((kk-1)+(NPPB/2)+(+NPPB/16)));
                                        if MeanOfData > LevelDec1%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                                            DataM = [DataM 1];                                 %is 1
                                        else                                                       %If it is the lowest level the incoming data
                                            DataM = [DataM 0];                                 %is 0
                                        end
                                        MeanOfData = mean(Esync((kk-1)+(NPPB/2)-(+NPPB/16)));
                                        if MeanOfData > LevelDec2%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                                            DataL = [DataL 1];                                 %is 1
                                        else                                                       %If it is the lowest level the incoming data
                                            DataL = [DataL 0];                                 %is 0
                                        end
                                        MeanOfData = mean(Esync((kk-1)+(NPPB/2)));
                                        if MeanOfData > LevDefDpqsk%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                                            DataU = [DataU 1];                                 %is 1
                                        else                                                       %If it is the lowest level the incoming data
                                            DataU = [DataU 0];                                 %is 0
                                        end
                                        if CalcS
                                            MeanOfData = mean(Esync((kk-1)+(NPPB/2)));
                                            if MeanOfData > EyeOpenLowI + EyeOpenI/2%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                                                DataS = [DataS 1];                                 %is 1
                                            else                                                       %If it is the lowest level the incoming data
                                                DataS = [DataS 0];                                 %is 0
                                            end
                                        end
                                    end
                                    %%       Calculating the Bit Error Ratio (BER)
                                    %The final process here is to count the number of wrongdoings
                                    %of this whole process upon the transmited data for
                                    %quantitative analizes
                                    TxData  = TxDataMat(ThisCarr,1+SyncPeriod:end-SyncPeriod);
                                    
                                    if CalcS
                                        DataS  = DataS(1+SyncPeriod:end-SyncPeriod);
                                        BitErrS    = sum(xor(TxData,DataS));                       %Comparison between the Transmited and received and counting the differences
                                        BerDPSKS(CurrentTest,ThisCarr) = BitErrS/(NbDPSK-(2*SyncPeriod));                  %Calculating the ration of wrong bits and the total number of bits transmited
                                    end
                                    Data  = Data(1+SyncPeriod:end-SyncPeriod);
                                    DataM  = DataM(1+SyncPeriod:end-SyncPeriod);
                                    DataL  = DataL(1+SyncPeriod:end-SyncPeriod);
                                    DataU  = DataU(1+SyncPeriod:end-SyncPeriod);
                                    AberLevAuxI(4) = 0;
                                    ValsLevAuxI(4) = 0.01;
                                    
                                    BitErr(1)  = sum(xor(TxData,Data));                        %Comparison between the Transmited and received and counting the differences
                                    BitErr(2)  = sum(xor(TxData,DataM));                       %Comparison between the Transmited and received and counting the differences
                                    BitErr(3)  = sum(xor(TxData,DataL));                       %Comparison between the Transmited and received and counting the differences
                                    BitErr(4)  = sum(xor(TxData,DataU));                       %Comparison between the Transmited and received and counting the differences
                                    BerDPSK(CurrentTest,ThisCarr) = (min(BitErr))/(NbDPSK-(2*SyncPeriod));%Calculating the ration of wrong bits and the total number of bits transmited
                                    AberIAux = max(AberLevAuxI(BitErr<=min(BitErr)));
                                    AberLev(CurrentTest,ThisCarr)  = AberIAux;
                                    ValsLev(CurrentTest,ThisCarr)  = max(ValsLevAuxI(AberLevAuxI==AberIAux));          %Comparison between the Transmited and received and counting the differences
                                    if BerDPSK(CurrentTest,ThisCarr)==0
                                        LevDefDpqsk = max(ValsLevAuxI(AberLevAuxI==AberIAux));
                                        if LevDefDpqsk==0
                                            LevDefDpqsk = 0.01;
                                        end
                                    end
                                    [~,LocAux] = max(AberLevAuxI==AberIAux);
                                    EvmDB(CurrentTest,ThisCarr)   = EvmDBA(LocAux);
                                    EvmPer(CurrentTest,ThisCarr)  = EvmPerA(LocAux);
                                    EvmRms(CurrentTest,ThisCarr)  = EvmRmsA(LocAux);
                                    EvmDBJ(CurrentTest,ThisCarr)  = EvmDBJA(LocAux);
                                    EvmPerJ(CurrentTest,ThisCarr) = EvmPerJA(LocAux);
                                    EvmRmsJ(CurrentTest,ThisCarr) = EvmRmsJA(LocAux);
                                    
                                    RxSymbAmos(CurrentTest,ObsCarrPos==ThisCarr,:) = EvmMatRecA(LocAux,:);%RxSymbAmosUp = [];
                                    EvmMatRec(ObsCarrPos==ThisCarr,:) = EvmMatRecA(LocAux,:);
                                    %% Ploting the result for qualitative analizes
                                    %PrintInfo(Ploting*19,TxData,Data);
                                    %%
                                    %             berpos = 1:2:size(BerDPSK,2);
                                    %             BerDPSK(size(BerDPSK,1),berpos)
                                    a=6;
                                    %                         close all;
                                    %                         BerDPSK(CurrentTest,ThisCarr+1) = 1;
                                end
                            end
                        case 'DQPSK'
                            
                            %%                   Receiver DQPSK
                            for ThisCarr=InitCarr:CarrPass:NumCarr                                           %For each carrier the same process of reception will be used.
                                if CarrUsedDo(ThisCarr)
                                    %%  Reception Process: Handling the income optical field
                                    %At the first step the income field will be selected from the
                                    %optical FFT Output with the help of the maping vector
                                    %(previously described).
                                    EoutAux = EoutAux1(VetThisCarr==ObsCarrUsed(ThisCarr),:);
                                    if AddingNoiseP
                                        [~,sigpower] = MeasPower(EoutAux);
                                        EoutAux      = awgn(EoutAux,osnrp,sigpower-30);
                                        [~,sigpower] = MeasPower(EoutAux);
                                    end
                                    %%            Fiber Time Delay Compensation
                                    switch Medium
                                        case 'Fiber'
                                            EoutAux = ifft(fft(EoutAux).*exp(1j*2*pi*f*(...
                                                FiberDelay(ThisCarr)*Ta)));
                                        otherwise
                                    end
                                    %             EoutAux = EoutRec;
                                    %%  Plot for Qualitative analizes
                                    %PrintInfo(Ploting*20,t,abs(EoutAux));
                                    %             PrintInfo(Ploting*21,f,20*log10(abs(fftshift(fft(EoutAux)./...
                                    %                                                        length(EoutAux)))));
                                    
                                    %%                   synchronizing
                                    %Differently from the previuos process at the DQPSK
                                    %demodulation it is needed to make the synchronism before
                                    %interferometric process. Because the recovery of the
                                    %transmited data is extremally dependent of the time delay.
                                    %Thus, for the reception process work properly, it is needed to
                                    %sincronized the recieved signal or the sampling process will
                                    %not work.
                                    
                                    %At this first moment a small part of the income signal needs
                                    %be analized for the sincronism process. As we are not looking
                                    %for the actual data rather for the synchronis symbol, which is
                                    %the bigest phase shift of the income signal. The next delay
                                    %interferometer convert this phase shift to a amplitude
                                    %variation.
                                    
                                    %             %The phase delay is not needed because the optical field will
                                    %             %be analized as a whole.
                                    %             PhaDel          = 0;
                                    %             %Remember that at DQPSK modulation the information is stored at
                                    %             %the difference of phase between the current and previous
                                    %             %symbol hence this time delay is needed.
                                    %             TimDel          = T;
                                    %             [ESync1,ESync2] = DelayInterf(t,TimDel,PhaDel,EoutAux);        %Coverting phase shift to amplitude variation
                                    %
                                    %             %The second moment is to transfer this information from the
                                    %             %optical domain to the eletrical domain for an eletronic
                                    %             %processing. It is done with the help of an photo diode. The
                                    %             %configuration here used is an balanced receiver as the output
                                    %             %of the Delay Interferometer has two signals resulting from
                                    %             %the constructive and destructive signal interaction.
                                    %             ESync1 = ESync1.*conj(ESync1);
                                    %             ESync2 = ESync2.*conj(ESync2);
                                    %             Esync  = ESync2-ESync1;
                                    %
                                    %             ESync1 = ESync1./max(abs(ESync1));                             %Normalizing the signal
                                    %             ESync2 = ESync2./max(abs(ESync2));                             %Normalizing the signal
                                    %             Esync  = Esync./max(abs(Esync));                               %Normalizing the signal
                                    %
                                    %             SyncAux   = Esync(IniSyncPos:SyncPos);                         %Selecting just the symbol to synchronize
                                    %             SyncedAux = SyncSymb(IniSyncPos:SyncPos);                      %Selecting just the symbol to synchronize
                                    %             %%                   Plot for Qualitative analizes
                                    % %              PrintInfo(Ploting*22,t(IniSyncPos:SyncPos),Esync(IniSyncPos...
                                    % %              :SyncPos),SyncSymb(IniSyncPos:SyncPos),ESync1(IniSyncPos:...
                                    % %                                       SyncPos),ESync2(IniSyncPos:SyncPos));
                                    %             %%                   Synchronizing
                                    %             %This synchronization process is based on the mean expected
                                    %             %value. That means, the information mean value should be within
                                    %             %one period of symbol. Thus, the mean value of the received
                                    %             %signal is acquired and compare of the known sync-word to
                                    %             %verify if this mean value is at the right possition. Which is
                                    %             %the midel point (peak) of the highest level at the sync period
                                    %             SyncAux(SyncAux<0)              = 0;                           %To keep the mean value above zero anything under is neglected
                                    %             SyncAux(SyncAux>=mean(SyncAux)) = 1;                           %Adding a flag to the first sample of the received mean value
                                    %             SyncAux(SyncAux<mean(SyncAux))  = -1;                          %All the others samples at set to the lowest level
                                    %
                                    %             PosToSyn  = find(ismember(SyncAux,1));                         %Finding where is the location of the samples to synchronize
                                    %             PosSynced = find(ismember(SyncedAux,1));                       %Finding where is the location of the samples to synchronize
                                    %
                                    %             DiffPos = PosToSyn(round(end/2)) - PosSynced(round(end/2));    %Accounting the peak (midel point) displacement
                                    %             if DiffPos>=0%If the difference is positive, left-shift...
                                    %                 EoutAux = ifft(fft(EoutAux).*exp(1j*2*pi*f*(DiffPos*Ta))); %Shift based on time change
                                    % %                 EoutAux = [EoutAux(DiffPos+1:end) EoutAux(1:DiffPos)];   %Shift based on sampling sliding
                                    %             else%... but if the difference is negative, right-shift
                                    %                 EoutAux = ifft(fft(EoutAux).*exp(1j*2*pi*f*(DiffPos*Ta))); %Shift based on time change
                                    % %                 EoutAux = [EoutAux(end+DiffPos+1:end) EoutAux(1:end + ...
                                    % %                                                                 DiffPos)];%Shift based on sampling sliding
                                    %             end
                                    %
                                    %                 %                   Synchronizing
                                    %             %Because of reasons, sometimes it may be required to make a
                                    %             %synchronization process with the  end of the data stream as
                                    %             %well. This following verification check if the user set (or
                                    %             %not) a second synchronization process to be done.
                                    %             if SencondAdjust
                                    %                 PhaDel          = 0;
                                    %                 TimDel          = T;
                                    %                 [ESync1,ESync2] = DelayInterf(t,TimDel,PhaDel,EoutAux);
                                    %
                                    %                 ESync1 = ESync1.*conj(ESync1);
                                    %                 ESync2 = ESync2.*conj(ESync2);
                                    %                 Esync  = ESync2-ESync1;
                                    %
                                    %                 ESync1 = ESync1./max(abs(ESync1));
                                    %                 ESync2 = ESync2./max(abs(ESync2));
                                    %                 Esync  = Esync./max(abs(Esync));
                                    %
                                    %
                                    %                 %This synchronization process is based on the mean expected
                                    %                 %value. That means, the information mean value should be
                                    %                 %within one period of symbol. Thus, the mean value of the
                                    %                 %received signal is acquired and compare of the known
                                    %                 %sync-word to verify if this mean value is at the right
                                    %                 %possition.
                                    %
                                    %                 SyncAuxEnd   = Esync(end-SyncPos+1:end-IniSyncPos+1);
                                    %                 SyncedAuxEnd = SyncSymbEnd(end-SyncPos+1:end-IniSyncPos+1);
                                    %                 SyncAuxEnd(SyncAuxEnd<0)                 = 0;              %To keep the mean value above zero anything under is neglected
                                    %                 SyncAuxEnd(SyncAuxEnd>=mean(SyncAuxEnd)) = 1;              %Adding a flag to the first sample of the received mean value
                                    %                 SyncAuxEnd(SyncAuxEnd<mean(SyncAuxEnd))  = -1;             %All the others samples at set to the lowest level
                                    %
                                    %                 PosToSynEnd  = find(ismember(SyncAuxEnd,1));               %Finding where is the location of the first sample to synchronize
                                    %                 PosSyncedEnd = find(ismember(SyncedAuxEnd,1));             %Finding where is the location of the first sample to synchronize
                                    %
                                    %                 DiffPosEnd = PosToSynEnd(round(end/2))-PosSyncedEnd(...
                                    %                                                              round(end/2));
                                    %                 if DiffPosEnd>=0%If positive difference, left-shift...
                                    %                     EoutAux = ifft(fft(EoutAux).*exp(1j*2*pi*f*(...
                                    %                                                           DiffPosEnd*Ta)));%Shift based on time change
                                    % %                     EoutAux = [EoutAux(DiffPosEnd+1:end) EoutAux(1:...
                                    % %                                                             DiffPosEnd)];%Shift based on sampling sliding
                                    %                 else%... but if the difference is negative, right-shift
                                    %                     EoutAux = ifft(fft(EoutAux).*exp(1j*2*pi*f*(...
                                    %                                                           DiffPosEnd*Ta)));%Shift based on time change
                                    % %                     EoutAux = [EoutAux(end+DiffPosEnd+1:end) EoutAux(...
                                    % %                                                       1:end+DiffPosEnd)];%Shift based on sampling sliding
                                    %                 end
                                    %             end
                                    
                                    %%                   Plot for Qualitative analizes
                                    %             PrintInfo(Ploting*23,t(IniSyncPos:SyncPos),Esync(IniSyncPos:...
                                    %                 SyncPos),SyncSymb(IniSyncPos:SyncPos),ESync1(IniSyncPos:...
                                    %                                       SyncPos),ESync2(IniSyncPos:SyncPos));
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
                                    
                                    E_rec3 = EoutAux./max(abs(EoutAux));                           %Normalizing income signal
                                    %For the interferometric process  take in account just the real
                                    %component it is needed a phase delay of 45° degrees;
                                    PhaDel = 1*pi/4;
                                    %Remember that at DQPSK modulation the information is stored at
                                    %the difference of phase between the current and previous
                                    %symbol hence this time delay of one symbol period is needed.
                                    TimDel = T;
                                    [EoutA,EoutB] = DelayInterf(t,TimDel,PhaDel,EoutAux);          %Coverting phase shift to amplitude variation
                                    %For the interferometric process  take in account just the real
                                    %component it is needed a phase delay of -45° degrees;
                                    PhaDel = -1*pi/4;
                                    %Remember that at DQPSK modulation the information is stored at
                                    %the difference of phase between the current and previous
                                    %symbol hence this time delay of one symbol period is needed.
                                    TimDel = T;
                                    [EoutC,EoutD] = DelayInterf(t,TimDel,PhaDel,EoutAux);          %Coverting phase shift to amplitude variation
                                    
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
                                    
                                    if ReceptorNoise                                               %Verify whether the noise is to be added or not
                                        [~,SigPowerA] = MeasPower(EoutA);                          %it own received signal or
                                        [~,SigPowerB] = MeasPower(EoutB);
                                        [~,SigPowerC] = MeasPower(EoutC);                          %it own received signal or
                                        [~,SigPowerD] = MeasPower(EoutD);
                                        SNR2 = CarSNR + 10*log10(2) - 10*log10(Nsamp) + 10*log10(10^0.3);
                                        
                                        %                     if TimeSys
                                        %                         SNR = CarSNR + 10*log10(2) - 10*log10(T*BWD2);
                                        %                     else
                                        %                         SNR = CarSNR + 10*log10(2) - 10*log10(1*t(2*NumAmosCP+NPPB)*BWD2);
                                        %                     end
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
                                    
                                    EoutI = (EoutB - EoutA);
                                    EoutQ = (EoutD - EoutC);
                                    EoutI = EoutI-mean(EoutI);
                                    EoutQ = EoutQ-mean(EoutQ);
                                    EoutI = EoutI./max(abs(EoutI));                                %Normalizing the signal
                                    EoutQ = EoutQ./max(abs(EoutQ));                                %Normalizing the signal
                                    VetElecPowerI(CurrentTest,ThisCarr)= MeasPower(EoutI);
                                    VetElecPowerQ(CurrentTest,ThisCarr)= MeasPower(EoutQ);
                                    EoutAuxF = 20*log10(abs(fftshift(fft(EoutAux)./length(...
                                        EoutAux))));
                                    [VetOptiPower(CurrentTest,ThisCarr),~]= findpeaks(EoutAuxF,...
                                        'SortStr','descend','NPeaks',1);
                                    
                                    
                                    SyncAuxI   = EoutI(IniSyncPos:SyncPos);                        %Selecting just the symbol to synchronize
                                    SyncAuxQ   = EoutQ(IniSyncPos:SyncPos);                        %Selecting just the symbol to synchronize
                                    %                 PrintInfo(Ploting*22,t(IniSyncPos:SyncPos),EoutI(IniSyncPos:SyncPos),SyncSymb(IniSyncPos:SyncPos),EoutA(IniSyncPos:SyncPos),EoutB(IniSyncPos:SyncPos));
                                    
                                    SyncedAux  = SyncSymb(IniSyncPos:SyncPos);                     %Selecting just the symbol to synchronize
                                    %%                   Synchronizing
                                    %This synchronization process is based on the mean expected
                                    %value. That means, the information mean value should be within
                                    %one period of symbol. Thus, the mean value of the received
                                    %signal is acquired and compare of the known sync-word to
                                    %verify if this mean value is at the right possition. Which is
                                    %the midel point (peak) of the highest level at the sync period
                                    SyncAuxI(EoutI(IniSyncPos:SyncPos)<0.5*max((EoutI(IniSyncPos:SyncPos))))               = 0;   %To keep the mean value above zero anything under is neglected
                                    SyncAuxQ(EoutQ(IniSyncPos:SyncPos)<0.5*max((EoutQ(IniSyncPos:SyncPos))))               = 0;   %To keep the mean value above zero anything under is neglected
                                    SyncAuxI(SyncAuxI>=mean(SyncAuxI)) = 1;   %Adding a flag to the first sample of the received mean value
                                    SyncAuxQ(SyncAuxQ>=mean(SyncAuxQ)) = 1;   %Adding a flag to the first sample of the received mean value
                                    SyncAuxI(SyncAuxI<mean(SyncAuxI))  = -1;  %All the others samples at set to the lowest level
                                    SyncAuxQ(SyncAuxQ<mean(SyncAuxQ))  = -1;  %All the others samples at set to the lowest level
                                    
                                    PosToSynI  = find(ismember(SyncAuxI,1));                       %Finding where is the location of the samples to synchronize
                                    PosToSynQ  = find(ismember(SyncAuxQ,1));                       %Finding where is the location of the samples to synchronize
                                    PosSynced = find(ismember(SyncedAux,1));                       %Finding where is the location of the samples to synchronize
                                    
                                    DiffPosI = ExtDel*(PosToSynI(round(end/2)) - PosSynced(round(end/2)));  %Accounting the peak (midel point) displacement
                                    DiffPosQ = ExtDel*(PosToSynQ(round(end/2)) - PosSynced(round(end/2)));  %Accounting the peak (midel point) displacement
                                    sn=((1/2)*(((1/2)^MaxNumStag)-1))/((1/2)-1);
                                    DiffPosI = round((sn/(3/2.55^IfftOrSum))*T/Ta);
                                    DiffPosQ = round((sn/(3/2.55^IfftOrSum))*T/Ta);
                                    if DiffPosI>=0%If the difference is positive, left-shift...
                                        %                 EoutI = ifft(fft(EoutI).*exp(1j*2*pi*f*(DiffPosI*Ta))); %Shift based on time change
                                        EoutI = [EoutI(DiffPosI+1:end) EoutI(1:DiffPosI)];   %Shift based on sampling sliding
                                        EoutQ = [EoutQ(DiffPosI+1:end) EoutQ(1:DiffPosI)];   %Shift based on sampling sliding
                                    else%... but if the difference is negative, right-shift
                                        %                 EoutI = ifft(fft(EoutI).*exp(-1j*2*pi*f*(DiffPosI*Ta))); %Shift based on time change
                                        EoutI = [EoutI(end+DiffPosI+1:end) EoutI(1:end+DiffPosI)]; %Shift based on sampling sliding
                                        EoutQ = [EoutQ(end+DiffPosI+1:end) EoutQ(1:end+DiffPosI)]; %Shift based on sampling sliding
                                    end
                                    %             if DiffPosQ>=0%If the difference is positive, left-shift...
                                    % %                 EoutQ = (ifft(fft(EoutQ).*exp(1j*2*pi*f*(DiffPosQ*Ta)))); %Shift based on time change
                                    %                 EoutQ = [EoutQ(DiffPosQ+1:end) EoutQ(1:DiffPosQ)];   %Shift based on sampling sliding
                                    %             else%... but if the difference is negative, right-shift
                                    % %                 EoutQ = (ifft(fft(EoutQ).*exp(-1j*2*pi*f*(DiffPosQ*Ta)))); %Shift based on time change
                                    %                 EoutQ = [EoutQ(end+DiffPosQ+1:end) EoutQ(1:end+DiffPosQ)]; %Shift based on sampling sliding
                                    %             end
                                    %%                   Plot for Qualitative analizes
                                    %                              PrintInfo(Ploting*22,t(IniSyncPos:SyncPos),EoutI(IniSyncPos...
                                    %                              :SyncPos),SyncSymb(IniSyncPos:SyncPos),EoutA(IniSyncPos:...
                                    %                                                       SyncPos),EoutB(IniSyncPos:SyncPos));
                                    %                              PrintInfo(Ploting*22,t(IniSyncPos:SyncPos),EoutQ(IniSyncPos...
                                    %                              :SyncPos),SyncSymb(IniSyncPos:SyncPos),EoutC(IniSyncPos:...
                                    %                                                       SyncPos),EoutD(IniSyncPos:SyncPos));
                                    %% Removing CP
                                    if AddCP
                                        IxAux1  = EoutI(1:end - StuffSampels);
                                        IxAux1  = reshape(IxAux1,(2*NumAmosCP+NPPB),NbDQPSK/2);
                                        IxAux1  = IxAux1(1+NumAmosCP:end-NumAmosCP,:);
                                        IxAux1  = reshape(IxAux1,1,NPPB*NbDQPSK/2);
                                        EoutI  = IxAux1;
                                        
                                        IxAux2  = EoutQ(1:end - StuffSampels);
                                        IxAux2  = reshape(IxAux2,(2*NumAmosCP+NPPB),NbDQPSK/2);
                                        IxAux2  = IxAux2(1+NumAmosCP:end-NumAmosCP,:);
                                        IxAux2  = reshape(IxAux2,1,NPPB*NbDQPSK/2);
                                        EoutQ  = IxAux2;
                                        
                                        IxAux3  = E_rec3(1:end - StuffSampels);
                                        IxAux3  = reshape(IxAux3,(2*NumAmosCP+NPPB),NbDQPSK/2);
                                        IxAux3  = IxAux3(1+NumAmosCP:end-NumAmosCP,:);
                                        IxAux3  = reshape(IxAux3,1,NPPB*NbDQPSK/2);
                                        E_rec3  = IxAux3;
                                    end
                                    %% Taking the sampling the EVM meassurement
                                    clear IxAux;
                                    PosAuxEout1 = NPPB/2:NPPB:length(EoutI);%Varriable respossible to take just the samples at the middle of the symbol
                                    PosAuxEout2 = ((NPPB/2)+(NPPB/4)):NPPB:length(EoutI);%Varriable respossible to take just the samples at the middle of the symbol
                                    PosAuxEout3 = ((NPPB/2)-(NPPB/4)):NPPB:length(EoutI);%Varriable respossible to take just the samples at the middle of the symbol
                                    IxAux1      = EoutI(PosAuxEout1) + 1j.*EoutQ(PosAuxEout1);    %Normalizing the reference
                                    IxAux2      = EoutI(PosAuxEout2) + 1j.*EoutQ(PosAuxEout2);    %Normalizing the reference
                                    IxAux3      = EoutI(PosAuxEout3) + 1j.*EoutQ(PosAuxEout3);    %Normalizing the reference
                                    a=0;
                                    EvmMatRecA(1,:) = IxAux1;                       %Taking just the middle samples as references
                                    EvmMatRecA(2,:) = IxAux2;                       %Taking just the middle samples as references
                                    EvmMatRecA(3,:) = IxAux3;                       %Taking just the middle samples as references
                                    EvmMatRecA(4,:) = IxAux1;                       %Taking just the middle samples as references
                                    [EvmDBA(1), EvmPerA(1), EvmRmsA(1) ] = EvmCalc( EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux1 );
                                    [EvmDBA(2), EvmPerA(2), EvmRmsA(2) ] = EvmCalc( EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux2 );
                                    [EvmDBA(3), EvmPerA(3), EvmRmsA(3) ] = EvmCalc( EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux3 );
                                    [EvmDBA(4), EvmPerA(4), EvmRmsA(4) ] = EvmCalc( EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux1 );
                                    [EvmDBJA(1),EvmPerJA(1),EvmRmsJA(1)] = evm1(4,'dpsk',EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux1);
                                    [EvmDBJA(2),EvmPerJA(2),EvmRmsJA(2)] = evm1(4,'dpsk',EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux2);
                                    [EvmDBJA(3),EvmPerJA(3),EvmRmsJA(3)] = evm1(4,'dpsk',EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux3);
                                    [EvmDBJA(4),EvmPerJA(4),EvmRmsJA(4)] = evm1(4,'dpsk',EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux1);
                                    %                 if ThisCarr==1
                                    %                     figure;hold all;grid on;
                                    %                     plot(real(IxAux1),imag(IxAux1),'rx');
                                    %                     plot(real(IxAux2),imag(IxAux2),'yx');
                                    %                     plot(real(IxAux3),imag(IxAux3),'g^');
                                    %                     set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                    %                     a=0;
                                    %                 end
                                    %%
                                    %                 figure;plot(t(1:length(E_rec3)),abs(E_rec3),t(1:length(EoutQ)),EoutQ,t(1:length(EoutI)),EoutI);set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                    %One important step on this project was the confirmation of our
                                    %results from a closed and theoretical equation that relates
                                    %the income optical field witht its respective data (as
                                    %descrived above in the receiving process). The result from
                                    %this equation was further compared with the result from the
                                    %MZ-Interferometer as an proof of concept. This equation can be
                                    %found at the book of Optical Fiber Telecommunications V B,
                                    %which one of the authors is Ivan P. Kaminow at the page 144.
                                    taux = t(1:length(E_rec3));
                                    faux = time2freq(taux);
                                    Ui = real(exp(1j*pi/4).*E_rec3.*conj(ifft(fft(E_rec3).*exp(-...
                                        1j*2*pi*faux*T))));%The data at odd position
                                    Ui = Ui./max(abs(Ui));                                         %Normalizing the signal
                                    Uq = imag(exp(1j*pi/4).*E_rec3.*conj(ifft(fft(E_rec3).*exp(-...
                                        1j*2*pi*faux*T))));%The data at the even position
                                    Uq = Uq./max(abs(Uq));                                         %Normalizing the signal
                                    %##########################################################################
                                    PosIx = NPPB/2:NPPB:length(EoutI);                                %Possition of the central samples - but the number of samples per symbol is even ????
                                    IxAux = EoutI(PosIx);                                             %From the main received signal just a few samples are taken for further evaluation
                                    n     = 100;                                                   %The number of boxes to be filled up on the histogram process
                                    Perce = 0.7;
                                    %At this process each eye will be evaluated separetelly. It is
                                    %important to remember that the PAM4 has 4 levels, which means
                                    %three level of decissions that are exaclty the center of the
                                    %eye diagram.
                                    IxAuxAB = EoutI(PosIx);       %Taking just those values relative to the uper eye
                                    InterAB = linspace(Perce*(min(EoutI)),Perce*(max(EoutI)),n);                     %Building the histogram boxes
                                    EyeAB = hist(IxAuxAB,InterAB);                                 %filling up the boxes with samples that fit on them.
                                    
                                    QIxAuxAB = EoutQ(PosIx);       %Taking just those values relative to the uper eye
                                    QInterAB = linspace(Perce*(min(EoutQ)),Perce*(max(EoutQ)),n);                     %Building the histogram boxes
                                    QEyeAB = hist(QIxAuxAB,QInterAB);                                 %filling up the boxes with samples that fit on them.
                                    %The same process described for the uper level will be done at
                                    %the middle and lower eyes levels.
                                    IxAuxCD = EoutI(PosIx-NPPB/4);
                                    InterCD = linspace(Perce*(min(EoutI)),Perce*(max(EoutI)),n);%NPPB*2^n);
                                    EyeCD = hist(IxAuxCD,InterCD);
                                    
                                    QIxAuxCD = EoutQ(PosIx-NPPB/4);
                                    QInterCD = linspace(Perce*(min(EoutQ)),Perce*(max(EoutQ)),n);%NPPB*2^n);
                                    QEyeCD = hist(QIxAuxCD,QInterCD);
                                    
                                    IxAuxEF = EoutI(PosIx+NPPB/4);
                                    InterEF = linspace(Perce*(min(EoutI)),Perce*(max(EoutI)),n);%NPPB*2^n);
                                    EyeEF = hist(IxAuxEF,InterEF);
                                    
                                    QIxAuxEF = EoutQ(PosIx+NPPB/4);
                                    QInterEF = linspace(Perce*(min(EoutQ)),Perce*(max(EoutQ)),n);%NPPB*2^n);
                                    QEyeEF = hist(QIxAuxEF,QInterEF);
                                    
                                    %What we are looking for here are not where exist the
                                    %occurrences of values. The eye is where there is not samples,
                                    %this means, the decission level is a reagion between actual
                                    %levels thus this reagion does not contain any signal.
                                    EyeAB = ~EyeAB;                                                %Changing zeros to one - Zeros compose the eye region
                                    EyeCD = ~EyeCD;
                                    EyeEF = ~EyeEF;
                                    
                                    QEyeAB = ~QEyeAB;                                                %Changing zeros to one - Zeros compose the eye region
                                    QEyeCD = ~QEyeCD;
                                    QEyeEF = ~QEyeEF;
                                    %Starting variables that will be used to identify the eye
                                    %diagram.
                                    CountAB   = 1;
                                    CountCD   = 1;
                                    CountEF   = 1;
                                    
                                    QCountAB   = 1;
                                    QCountCD   = 1;
                                    QCountEF   = 1;
                                    
                                    SeqOnesAB = 0;
                                    SeqOnesCD = 0;
                                    SeqOnesEF = 0;
                                    
                                    QSeqOnesAB = 0;
                                    QSeqOnesCD = 0;
                                    QSeqOnesEF = 0;
                                    
                                    SeqFinAB  = 0;
                                    SeqFinCD  = 0;
                                    SeqFinEF  = 0;
                                    
                                    QSeqFinAB  = 0;
                                    QSeqFinCD  = 0;
                                    QSeqFinEF  = 0;
                                    
                                    SeqIniAB  = 1;
                                    SeqIniCD  = 1;
                                    SeqIniEF  = 1;
                                    
                                    QSeqIniAB  = 1;
                                    QSeqIniCD  = 1;
                                    QSeqIniEF  = 1;
                                    %The for loop will take account of every box with ones. It is
                                    %important to take note that the not operator was used in this
                                    %vector, therefore ones means zeros (the eye diagram -
                                    %possibly) and zeros means values abouve zeroa (not the eye).
                                    for kk=1:length(EyeAB)                                         %For every box
                                        if EyeAB(kk)                                               %if it contains "1"
                                            SeqOnesAB(SeqIniAB)=CountAB;                           %count this element as part of a consecutive sequency
                                            CountAB = CountAB + 1;                                 %adds one to the counter of consecutive elements "1"
                                            if kk==length(EyeAB)                                   %if the current box is the last box we got to an end
                                                SeqFinAB(SeqIniAB) = kk;                           %The final sequency element is equal to its possition (kk)
                                            end
                                        else                                                       %else if the current box contains "0"
                                            SeqFinAB(SeqIniAB) = kk-1;                             %the previous element was the last element of a consecutive sequency of ones
                                            SeqIniAB = SeqIniAB + 1;                               %adds one to the counter of consecutive number (take in account how many sequencies there are)
                                            CountAB = 1;                                           %reset the counter of consecutive elements to mark the start of a new consecutive sequence
                                        end
                                        
                                        if QEyeAB(kk)                                               %if it contains "1"
                                            QSeqOnesAB(QSeqIniAB)=QCountAB;                           %count this element as part of a consecutive sequency
                                            QCountAB = QCountAB + 1;                                 %adds one to the counter of consecutive elements "1"
                                            if kk==length(QEyeAB)                                   %if the current box is the last box we got to an end
                                                QSeqFinAB(QSeqIniAB) = kk;                           %The final sequency element is equal to its possition (kk)
                                            end
                                        else                                                       %else if the current box contains "0"
                                            QSeqFinAB(QSeqIniAB) = kk-1;                             %the previous element was the last element of a consecutive sequency of ones
                                            QSeqIniAB = QSeqIniAB + 1;                               %adds one to the counter of consecutive number (take in account how many sequencies there are)
                                            QCountAB = 1;                                           %reset the counter of consecutive elements to mark the start of a new consecutive sequence
                                        end
                                        
                                        if EyeCD(kk)
                                            SeqOnesCD(SeqIniCD)=CountCD;
                                            CountCD = CountCD + 1;
                                            if kk==length(EyeCD)
                                                SeqFinCD(SeqIniCD) = kk;
                                            end
                                        else
                                            SeqFinCD(SeqIniCD) = kk-1;
                                            SeqIniCD = SeqIniCD + 1;
                                            CountCD = 1;
                                        end
                                        
                                        if QEyeCD(kk)
                                            QSeqOnesCD(QSeqIniCD)=QCountCD;
                                            QCountCD = QCountCD + 1;
                                            if kk==length(QEyeCD)
                                                QSeqFinCD(QSeqIniCD) = kk;
                                            end
                                        else
                                            QSeqFinCD(QSeqIniCD) = kk-1;
                                            QSeqIniCD = QSeqIniCD + 1;
                                            QCountCD = 1;
                                        end
                                        
                                        if EyeEF(kk)
                                            SeqOnesEF(SeqIniEF)=CountEF;
                                            CountEF = CountEF + 1;
                                            if kk==length(EyeEF)
                                                SeqFinEF(SeqIniEF) = kk;
                                            end
                                        else
                                            SeqFinEF(SeqIniEF) = kk-1;
                                            SeqIniEF = SeqIniEF + 1;
                                            CountEF = 1;
                                        end
                                        
                                        if QEyeEF(kk)
                                            QSeqOnesEF(QSeqIniEF)=QCountEF;
                                            QCountEF = QCountEF + 1;
                                            if kk==length(QEyeEF)
                                                QSeqFinEF(QSeqIniEF) = kk;
                                            end
                                        else
                                            QSeqFinEF(QSeqIniEF) = kk-1;
                                            QSeqIniEF = QSeqIniEF + 1;
                                            QCountEF = 1;
                                        end
                                    end
                                    
                                    
                                    %If the eye is open, which means there is a clear difference
                                    %between adjacent levels, the eye diagram will be the longest
                                    %sequence of ones.
                                    [MaxValAB,LocMaxAB]=max(SeqOnesAB);
                                    if LocMaxAB<2 || MaxValAB<2                                    %if any sequency was found or there is just one sequency it is a error thus
                                        LevDec3 = 0.1;                                            %the decission level will be by default 0.67. Also, the other variables will
                                        LocMaxAB = 1;                                              %will be set with values to not cause errors in the future
                                        SeqFinAB(1)=2;
                                        MaxValAB = 0;
                                        InterAB(1)=LevDec3;
                                    else                                                           %if a sequency was found the middle of the eye will be the middle of the sequency
                                        if (SeqFinAB(LocMaxAB)-MaxValAB/2)<1                       %if for some reason the final element of the sequency minus the half of its
                                            LevDec3 = 0.1;                                         %length results in a negative value, something went very wrong, and by
                                        else                                                       %default it will be set to 0.7
                                            LevDec3 = InterAB(round(SeqFinAB(LocMaxAB)-MaxValAB/...
                                                2));%Otherwise, the decission level is the middle of the sequency
                                        end
                                    end
                                    [QMaxValAB,QLocMaxAB]=max(QSeqOnesAB);
                                    if QLocMaxAB<2 || QMaxValAB<2                                    %if any sequency was found or there is just one sequency it is a error thus
                                        QLevDec3 = 0.1;                                            %the decission level will be by default 0.67. Also, the other variables will
                                        QLocMaxAB = 1;                                              %will be set with values to not cause errors in the future
                                        QSeqFinAB(1)=2;
                                        QMaxValAB = 0;
                                        QInterAB(1)=QLevDec3;
                                    else                                                           %if a sequency was found the middle of the eye will be the middle of the sequency
                                        if (QSeqFinAB(QLocMaxAB)-QMaxValAB/2)<1                       %if for some reason the final element of the sequency minus the half of its
                                            QLevDec3 = 0.1;                                         %length results in a negative value, something went very wrong, and by
                                        else                                                       %default it will be set to 0.7
                                            QLevDec3 = QInterAB(round(QSeqFinAB(QLocMaxAB)-QMaxValAB/...
                                                2));%Otherwise, the decission level is the middle of the sequency
                                        end
                                    end
                                    
                                    [MaxValCD,LocMaxCD]=max(SeqOnesCD);
                                    if LocMaxCD<2 || MaxValCD<2
                                        LevDec2 = 0.1;
                                        LocMaxCD = 1;
                                        SeqFinCD(1)=2;
                                        MaxValCD = 0;
                                        InterCD(1)=LevDec2;
                                    else
                                        if (SeqFinCD(LocMaxCD)-MaxValCD/2)<1
                                            LevDec2 = 0.1;
                                        else
                                            LevDec2 = InterCD(round(SeqFinCD(LocMaxCD)-MaxValCD/2));
                                        end
                                    end
                                    [QMaxValCD,QLocMaxCD]=max(QSeqOnesCD);
                                    if QLocMaxCD<2 || QMaxValCD<2
                                        QLevDec2 = 0.1;
                                        QLocMaxCD = 1;
                                        QSeqFinCD(1)=2;
                                        QMaxValCD = 0;
                                        QInterCD(1)=LevDec2;
                                    else
                                        if (QSeqFinCD(QLocMaxCD)-QMaxValCD/2)<1
                                            QLevDec2 = 0.1;
                                        else
                                            QLevDec2 = QInterCD(round(QSeqFinCD(QLocMaxCD)-QMaxValCD/2));
                                        end
                                    end
                                    
                                    [MaxValEF,LocMaxEF]=max(SeqOnesEF);
                                    if LocMaxEF<2 || MaxValEF<2
                                        LevDec1 = 0.1;
                                        LocMaxEF = 1;
                                        SeqFinEF(1)=2;
                                        MaxValEF = 0;
                                        InterEF(1)=LevDec1;
                                    else
                                        if (SeqFinEF(LocMaxEF)-MaxValEF/2)<1
                                            LevDec1 = 0.1;
                                        else
                                            LevDec1 = InterEF(round(SeqFinEF(LocMaxEF)-MaxValEF/2));
                                        end
                                        
                                    end
                                    [QMaxValEF,QLocMaxEF]=max(QSeqOnesEF);
                                    if QLocMaxEF<2 || QMaxValEF<2
                                        QLevDec1 = 0.1;
                                        QLocMaxEF = 1;
                                        QSeqFinEF(1)=2;
                                        QMaxValEF = 0;
                                        QInterEF(1)=QLevDec1;
                                    else
                                        if (QSeqFinEF(QLocMaxEF)-QMaxValEF/2)<1
                                            QLevDec1 = 0.1;
                                        else
                                            QLevDec1 = QInterEF(round(QSeqFinEF(QLocMaxEF)-QMaxValEF/2));
                                        end
                                        
                                    end
                                    
                                    %another way to measure the eye opening is the get all the
                                    %boxes and find all peaks on it, that will be a plato created
                                    %by the sequences of ones (which was zeros). From thos peaks,
                                    %the eye diagram will be the longer of them hence it will take
                                    %the most part of the vector that store the findpeaks result.
                                    %Thus, the middle of the eye will be basically the middle of
                                    %the peaks vector.
                                    LocAB = find(EyeAB);
                                    if isempty(LocAB)                                              %if for some reason there are no peaks, something went wrong.
                                        LocAB = 1;
                                        LevelDec3 = 0.1;%mean(Levels(3:4));                       %by default the decission level will be set to 0.65
                                    else
                                        LevelDec3 = LevDec3;
                                    end
                                    QLocAB = find(QEyeAB);
                                    if isempty(QLocAB)                                              %if for some reason there are no peaks, something went wrong.
                                        QLocAB = 1;
                                        QLevelDec3 = 0.1;%mean(Levels(3:4));                       %by default the decission level will be set to 0.65
                                    else
                                        QLevelDec3 = QLevDec3;
                                    end
                                    
                                    LocCD = find(EyeCD);
                                    if isempty(LocCD)
                                        LocCD = 1;
                                        LevelDec2 = 0.1;%mean(Levels(2:3));
                                    else
                                        LevelDec2 = LevDec2;
                                    end
                                    QLocCD = find(QEyeCD);
                                    if isempty(QLocCD)
                                        QLocCD = 1;
                                        QLevelDec2 = 0.1;%mean(Levels(2:3));
                                    else
                                        QLevelDec2 = QLevDec2;
                                    end
                                    
                                    LocEF = find(EyeEF);
                                    if isempty(LocEF)
                                        LocEF = 1;
                                        LevelDec1 = 0.1;%mean(Levels(1:2));
                                    else
                                        LevelDec1 = LevDec1;
                                    end
                                    QLocEF = find(QEyeEF);
                                    if isempty(QLocEF)
                                        QLocEF = 1;
                                        QLevelDec1 = 0.1;%mean(Levels(1:2));
                                    else
                                        QLevelDec1 = QLevDec1;
                                    end
                                    %##########################################################################
                                    if CalcS
                                        [~,~,EyeOpenI,EyeOpenHighI,EyeOpenLowI] = Olho(EoutI,T,NPPB,0);
                                        [~,~,EyeOpenQ,EyeOpenHighQ,EyeOpenLowQ] = Olho(EoutQ,T,NPPB,0);
                                        AberLevIS(CurrentTest,ThisCarr)  = EyeOpenI;
                                        ValsLevIS(CurrentTest,ThisCarr)  = EyeOpenLowI + EyeOpenI/2;
                                        AberLevQS(CurrentTest,ThisCarr)  = EyeOpenQ;
                                        ValsLevQS(CurrentTest,ThisCarr)  = EyeOpenLowQ + EyeOpenQ/2;
                                    end
                                    
                                    AberLevAuxI(1)= InterAB(SeqFinAB(LocMaxAB))-InterAB(round(SeqFinAB(LocMaxAB)-MaxValAB));
                                    AberLevAuxI(3)= InterCD(SeqFinCD(LocMaxCD))-InterCD(round(SeqFinCD(LocMaxCD)-MaxValCD));
                                    AberLevAuxI(2)= InterEF(SeqFinEF(LocMaxEF))-InterEF(round(SeqFinEF(LocMaxEF)-MaxValEF));
                                    ValsLevAuxI(1)= LevDec3;
                                    ValsLevAuxI(3)= LevDec2;
                                    ValsLevAuxI(2)= LevDec1;
                                    AberLevAuxQ(1)= QInterAB(QSeqFinAB(QLocMaxAB))-QInterAB(round(QSeqFinAB(QLocMaxAB)-QMaxValAB));
                                    AberLevAuxQ(3)= QInterCD(QSeqFinCD(QLocMaxCD))-QInterCD(round(QSeqFinCD(QLocMaxCD)-QMaxValCD));
                                    AberLevAuxQ(2)= QInterEF(QSeqFinEF(QLocMaxEF))-QInterEF(round(QSeqFinEF(QLocMaxEF)-QMaxValEF));
                                    ValsLevAuxQ(1)= QLevDec3;
                                    ValsLevAuxQ(3)= QLevDec2;
                                    ValsLevAuxQ(2)= QLevDec1;
                                    
                                    %% Ploting the result for qualitative analizes
                                    
                                    if ThisCarr==126
                                        EyeToPlot(CurrentTest,1:length([EoutI EoutQ])) = [EoutI EoutQ];
                                        save(['savingforplotingeye' num2str(VetSnr)],'EyeToPlot','VetSnr','IfftOrSum','UsedModula','T','NPPB');
                                    end
                                    if PrintinEye
                                        PrintInfo(Ploting*24,EoutI,T,NPPB);
                                        hold on;
                                        plot(t((NPPB)/2),InterAB(round(SeqFinAB(LocMaxAB)-MaxValAB)),'k^');
                                        plot(t((NPPB)/2),InterAB(SeqFinAB(LocMaxAB)),'kv');
                                        plot(t((NPPB)/2),InterAB(round(SeqFinAB(LocMaxAB)-MaxValAB/2)),'kd');
                                        plot(t((NPPB)/4),InterCD(round(SeqFinCD(LocMaxCD)-MaxValCD)),'g^');
                                        plot(t((NPPB)/4),InterCD(SeqFinCD(LocMaxCD)),'gv');
                                        plot(t((NPPB)/4),InterCD(round(SeqFinCD(LocMaxCD)-MaxValCD/2)),'gd');
                                        plot(t(3*(NPPB)/4),InterEF(round(SeqFinEF(LocMaxEF)-MaxValEF)),'b^');
                                        plot(t(3*(NPPB)/4),InterEF(SeqFinEF(LocMaxEF)),'bv');
                                        plot(t(3*(NPPB)/4),InterEF(round(SeqFinEF(LocMaxEF)-MaxValEF/2)),'b*');
                                        %PrintInfo(Ploting*25,Ui,T,NPPB);
                                        PrintInfo(Ploting*26,EoutQ,T,NPPB);
                                        hold on;
                                        plot(t((NPPB)/2),QInterAB(round(QSeqFinAB(QLocMaxAB)-QMaxValAB)),'k^');
                                        plot(t((NPPB)/2),QInterAB(QSeqFinAB(QLocMaxAB)),'kv');
                                        plot(t((NPPB)/2),QInterAB(round(QSeqFinAB(QLocMaxAB)-QMaxValAB/2)),'kd');
                                        plot(t((NPPB)/4),QInterCD(round(QSeqFinCD(QLocMaxCD)-QMaxValCD)),'g^');
                                        plot(t((NPPB)/4),QInterCD(QSeqFinCD(QLocMaxCD)),'gv');
                                        plot(t((NPPB)/4),QInterCD(round(QSeqFinCD(QLocMaxCD)-QMaxValCD/2)),'gd');
                                        plot(t(3*(NPPB)/4),QInterEF(round(QSeqFinEF(QLocMaxEF)-QMaxValEF)),'b^');
                                        plot(t(3*(NPPB)/4),QInterEF(QSeqFinEF(QLocMaxEF)),'bv');
                                        plot(t(3*(NPPB)/4),QInterEF(round(QSeqFinEF(QLocMaxEF)-QMaxValEF/2)),'b*');
                                        if CalcS
                                            plot(t((NPPB)/2),EyeOpenLowI,'mx');
                                            plot(t((NPPB)/2),EyeOpenHighI,'mo');
                                            plot(t((NPPB)/2),EyeOpenLowI + EyeOpenI/2,'md');
                                            plot(t((NPPB)/2),EyeOpenLowQ,'cx');
                                            plot(t((NPPB)/2),EyeOpenHighQ,'co');
                                            plot(t((NPPB)/2),EyeOpenLowQ + EyeOpenQ/2,'cd');
                                        end
                                        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                        figure;
                                        hold all;
                                        plotpos = zeros(1,length(IxAux1));
                                        plot(IxAux1,'o','color',[1 0.4 0]);
                                        plot(EvmMatRef(ObsCarrUsed(ThisCarr),:),'b*');
                                        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                        a=0;
                                    end
                                    %PrintInfo(Ploting*27,Uq,T,NPPB);
                                    %%
                                    
                                    %                 Txaux1 = rectpulse(TxDataMat(ThisCarr,:),NPPB);
                                    %                 Txaux2 = rectpulse(TxDataMat(ThisCarr+NumCarr,:),NPPB);
                                    %                 Txaux1(Txaux1==0) = -1;
                                    %                 Txaux2(Txaux2==0) = -1;
                                    %%  Ploting some results for qualitative analizes
                                    %             PrintInfo(Ploting*28,t(1:length(EoutI)),Txaux1,Txaux2,EoutI,...
                                    %                           EoutA,EoutB,EoutQ,EoutC,EoutD,real(Ui),real(Uq));
                                    %                       a=1;
                                    %%               Recovering the Information
                                    %After passing the optical signal to the eletrical domain, for
                                    %actually detect the data withing the signal the following
                                    %steps are needed.
                                    %
                                    %Finding Decission Levels:
                                    %The process for decoding the income signal will be based on
                                    %eletronic comparators. Inasmuch as the right decission level
                                    %must be acquired for accurately decide, within a symbol
                                    %periode, what that current leavel means (ones or zeros).
                                    %
                                    %The process hereafter of chosing the  decission levels is not
                                    %deterministic rather it is a statistic process. The main idea
                                    %is to take the decission level from the histogram generated
                                    %from the income signal stream.
                                    %
                                    %This process is realized inside the function Olho.
                                    %
                                    %Basicaly the decission level will be the minimal value of the
                                    %currente eye under evaluation plus the half of the its eye
                                    %opening.The following ilustration better describe this process
                                    %
                                    %Eye Limits:   + (1/2 * Eye Opening:)  =    Comparisson Limit:
                                    %
                                    %UperLevel  ______     ______________     _____
                                    %                 \   /      |       \   /
                                    %                  \ /       |        \ /
                                    %                   \ Half Eye Opening /     Decission Level
                                    %                  / \       |        / \
                                    %LowerLevel ______/   \______|_______/   \_____
                                    %
                                    %
                                    %Actualy Receiving Data:
                                    %Once the signal was processed the next step is through a
                                    %comparator decide the actual information received.
                                    %
                                    %
                                    %For accurately read the income signal is necessary more than
                                    %just have the decision level, it is also needed to know where
                                    %it occurs. For instance, if we decide to use just to take the
                                    %decision level and measure the mean value across and a portion
                                    %of the time period, which portion should we take? It would be
                                    %logical to take the central portion, but the bit may be
                                    %deformed in such way that the information gets concentrated
                                    %out of the middle part, which means the symbol is not
                                    %symmetric. The symmetry information can be acquired from the
                                    %eye diagram by measuring the longitudinal opening. The
                                    %following sketch better describes this process:
                                    %
                                    %                   Point of Symmetry
                                    %         _____     _______|_______     _____
                                    %              \   /               \   /
                                    %               \ /  Longitudinal   \ /
                                    %                \ ----------------- /
                                    %               / \    Opening      / \
                                    %         _____/   \_______________/   \_____
                                    %
                                    %With those two pieces of information, decision level and point
                                    %of symmetry, we have the X, Y coordinates for the centre of
                                    %the Eye Diagram. Therefore, as long as there is an opening on
                                    %it it will be possible to recover the transmitted information
                                    %without error... theoretically.
                                    %
                                    %As this process is also statistical, first we reshape the
                                    %income vector to analyze all periods at the same time.
                                    %                 EyeSymMatI = reshape(EoutI(1+SyncPeriod*NPPB:end-SyncPeriod*...
                                    %                     NPPB),NPPB,NbDQPSK/2-2*SyncPeriod);
                                    %                 %Then we take the values that compose the decision level
                                    %                 %because they will mark the point of symmetry.
                                    %                 %
                                    %                 %Firstly it was set the interval in which the histogram will be
                                    %                 %build. It is based on the number of samples per bit period.
                                    %                 IntervalI = linspace(min(EoutI(1+SyncPeriod*NPPB:end-...
                                    %                     SyncPeriod*NPPB)),max(EoutI(1+SyncPeriod*NPPB:end-...
                                    %                     SyncPeriod*NPPB)),2*NPPB);
                                    %                 %Therefore, the MATLAB hist function returns the number of
                                    %                 %occurrence of each interval.
                                    %                 EyeMaxI = hist(EoutI,IntervalI);
                                    %                 EyeMaxauxI = [0 EyeMaxI 0];                                    %Zeros are added at the EyeMax to auxiliate the finding peaks process
                                    %                 [~,EyeLocI] = findpeaks(EyeMaxauxI,'MinPeakDistance',NPPB/4,...
                                    %                     'SortStr','descend','NPeaks',2);%The peaks on the Eye profile will be the levels at the Eyes limit
                                    %                 %From the location of the max values that occured, which means
                                    %                 %the uper and lower level of the eye diagram it needs to take
                                    %                 %the actual value that those occurences represent that is
                                    %                 %withing the the Interval variable.
                                    %                 ValToSeekI = IntervalI(EyeLocI-1);
                                    %                 ValToSeekI = sort(ValToSeekI,'ascend');
                                    %                 %The number of ocurrences is a statical measure therefore one
                                    %                 %does not have control which interval will have the highest
                                    %                 %peak, thus it is important to ordenate the values to be seek
                                    %                 %from the lower part of the eye diagram to the uper part of the
                                    %                 %eye diagram.
                                    %                 OccuCountI = zeros(1,size(EyeSymMatI,1));                      %Auxiliar Variable for accounting.
                                    %                 for kk=1:size(EyeSymMatI,1)                                    %For every sample within a symbol period
                                    %
                                    %                     OccuCountI(kk) = OccuCountI(kk)+sum((EyeSymMatI(kk,:)>=...
                                    %                         min(EoutI))&(EyeSymMatI(kk,:)<=UpeSymPer*EyeOpenLowI)); %Account all occurencies of the valeu 1
                                    %                     OccuCountI(kk) = OccuCountI(kk)+sum((EyeSymMatI(kk,:)>=...
                                    %                         LowSymPer*EyeOpenHighI)&(EyeSymMatI(kk,:)<=max(EoutI))); %Account all occurencies of the valeu 2
                                    %                 end
                                    %The point of symmetry of the eye diagram will be where the
                                    %maximum number of occurrences were measured inasmuch as those
                                    %are the points where all the bits go to the center of the
                                    %symbol.  From the maximum number of occurrences, it can happen
                                    %for more than one sample within one symbol period, in the best
                                    %case, all samples would have the same accounting as it is
                                    %shown the ilustration above hence the symmetry will be at the
                                    %middle sample of this group of maximum occurrences. This value
                                    %can be found by the mean of the samples positions within the
                                    %symbol period. The problem with this approach is that the
                                    %signal must be synchronized with the maximum displacement of
                                    %a symbol period minus 25% of the eye Longitudinal opening if
                                    %the displacement is higher than that the point of symmetry
                                    %will be wrongly measured.
                                    %                 [SymLocI] = round(mean(find(ismember(OccuCountI,max(...
                                    %                     OccuCountI)))));
                                    %             [~,SymLocI] = findpeaks(OccuCountI,'SortStr','descend');       %The peak on the Eye profile will be the Symmetry level
                                    DataOdd = [];                                                  %Initialization of the vector that will store the income data
                                    DataOddM = [];                                                  %Initialization of the vector that will store the income data
                                    DataOddL = [];                                                  %Initialization of the vector that will store the income data
                                    DataOddU = [];                                                  %Initialization of the vector that will store the income data
                                    DataOddU2 = [];                                                  %Initialization of the vector that will store the income data
                                    DataOddU3 = [];                                                  %Initialization of the vector that will store the income data
                                    DataOddS = [];                                                  %Initialization of the vector that will store the income data
                                    for kk=1:NPPB:length(EoutI)                                    %The comparison process will be made for each symbol period
                                        %                 MeanOfData = mean(EoutI((kk-1)+SymLocI(1)));
                                        MeanOfData = mean(EoutI((kk-1)+NPPB/2));
                                        if MeanOfData > LevelDec3%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                                            DataOdd = [DataOdd 1];                                 %is 1
                                        else                                                       %If it is the lowest level the incoming data
                                            DataOdd = [DataOdd 0];                                 %is 0
                                        end
                                        MeanOfData = mean(EoutI((kk-1)+3*NPPB/4));
                                        if MeanOfData > LevelDec1%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                                            DataOddM = [DataOddM 1];                                 %is 1
                                        else                                                       %If it is the lowest level the incoming data
                                            DataOddM = [DataOddM 0];                                 %is 0
                                        end
                                        MeanOfData = mean(EoutI((kk-1)+NPPB/4));
                                        if MeanOfData > LevelDec2%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                                            DataOddL = [DataOddL 1];                                 %is 1
                                        else                                                       %If it is the lowest level the incoming data
                                            DataOddL = [DataOddL 0];                                 %is 0
                                        end
                                        MeanOfData = mean(EoutI((kk-1)+(NPPB/2)));
                                        if MeanOfData > LevDefDpqsk%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                                            DataOddU = [DataOddU 1];                                 %is 1
                                        else                                                       %If it is the lowest level the incoming data
                                            DataOddU = [DataOddU 0];                                 %is 0
                                        end
                                        MeanOfData = mean(EoutI((kk-1)+(NPPB/4)));
                                        if MeanOfData > LevDefDpqsk%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                                            DataOddU2 = [DataOddU2 1];                                 %is 1
                                        else                                                       %If it is the lowest level the incoming data
                                            DataOddU2 = [DataOddU2 0];                                 %is 0
                                        end
                                        MeanOfData = mean(EoutI((kk-1)+(3*NPPB/4)));
                                        if MeanOfData > LevDefDpqsk%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                                            DataOddU3 = [DataOddU3 1];                                 %is 1
                                        else                                                       %If it is the lowest level the incoming data
                                            DataOddU3 = [DataOddU3 0];                                 %is 0
                                        end
                                        if CalcS
                                            MeanOfData = mean(EoutI((kk-1)+NPPB/2));
                                            if MeanOfData > EyeOpenLowI + EyeOpenI/2%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                                                DataOddS = [DataOddS 1];                                 %is 1
                                            else                                                       %If it is the lowest level the incoming data
                                                DataOddS = [DataOddS 0];                                 %is 0
                                            end
                                        end
                                    end
                                    %The identical process just described above will also be used
                                    %to recover the data at the even positions.
                                    %
                                    %As this process is also statistical, first we reshape the
                                    %income vector to analyze all periods at the same time.
                                    %                 EyeSymMatQ = reshape(EoutQ(1+SyncPeriod*NPPB:end-SyncPeriod*...
                                    %                     NPPB),NPPB,NbDQPSK/2-2*SyncPeriod);
                                    %                 %Then we take the values that compose the decision level
                                    %                 %because they will mark the point of symmetry.
                                    %                 %
                                    %                 %Firstly it was set the interval in which the histogram will be
                                    %                 %build. It is based on the number of samples per bit period.
                                    %                 IntervalQ = linspace(min(EoutQ(1+SyncPeriod*NPPB:end-...
                                    %                     SyncPeriod*NPPB)),max(EoutQ(1+SyncPeriod*NPPB:end-...
                                    %                     SyncPeriod*NPPB)),2*NPPB);
                                    %                 %Therefore, the MATLAB hist function returns the number of
                                    %                 %occurrence of each interval.
                                    %                 EyeMaxQ = hist(EoutQ,IntervalQ);
                                    %                 EyeMaxauxQ = [0 EyeMaxI 0];                                    %Zeros are added at the EyeMax to auxiliate the finding peaks process
                                    %                 [~,EyeLocQ] = findpeaks(EyeMaxauxQ,'MinPeakDistance',NPPB/4,...
                                    %                     'SortStr','descend','NPeaks',2);%The peaks on the Eye profile will be the levels at the Eyes limit
                                    %                 ValToSeekQ = IntervalQ(EyeLocQ-1);
                                    %                 ValToSeekQ = sort(ValToSeekQ,'ascend');
                                    %                 OccuCountQ = zeros(1,size(EyeSymMatQ,1));                      %Auxiliar Variable for accounting.
                                    %                 for kk=1:size(EyeSymMatQ,1)                                    %For every sample within a symbol period
                                    %                     OccuCountQ(kk) = OccuCountQ(kk)+sum((EyeSymMatI(kk,:)>=...
                                    %                         LowSymPer*ValToSeekQ(1))&(EyeSymMatQ(kk,:)<=UpeSymPer*...
                                    %                         ValToSeekQ(1)));%Account all occurencies of the valeu 1
                                    %                     OccuCountQ(kk) = OccuCountQ(kk)+sum((EyeSymMatQ(kk,:)>=...
                                    %                         LowSymPer*ValToSeekQ(2))&(EyeSymMatQ(kk,:)<=UpeSymPer*...
                                    %                         ValToSeekQ(2)));%Account all occurencies of the valeu 2
                                    %                 end
                                    %                 %             [~,SymLocQ] = findpeaks(OccuCountQ,'SortStr','descend');       %The peak on the Eye profile will be the Symmetry level
                                    %                 [SymLocQ] = round(mean(find(ismember(OccuCountQ,max(...
                                    %                     OccuCountQ)))));
                                    %##############################################################
                                    %######################Important###############################
                                    %The ber results for the Data Even is not as good as the
                                    %results of the Data Odd. One possible reason is the decision
                                    %point of symmetry hence, for testing, we change SymLocQ to
                                    %SymLocI for evaluation of improvement. It is expected as bouth
                                    %signal will the same point of symetry will perform in an equal
                                    %way. If it is confirmed the creation of SymLocQ will be
                                    %erased.
                                    %##############################################################
                                    %##############################################################
                                    DataEven = [];                                                 %Initialization of the vector that will store the income data
                                    DataEvenM = [];                                                 %Initialization of the vector that will store the income data
                                    DataEvenL = [];                                                 %Initialization of the vector that will store the income data
                                    DataEvenU = [];                                                 %Initialization of the vector that will store the income data
                                    DataEvenU2 = [];                                                 %Initialization of the vector that will store the income data
                                    DataEvenU3 = [];                                                 %Initialization of the vector that will store the income data
                                    DataEvenS = [];                                                 %Initialization of the vector that will store the income data
                                    for kk=1:NPPB:length(EoutQ)                                    %The comparison process will be made for each symbol period
                                        %                 MeanOfData = mean(EoutQ((kk-1)+SymLocI(1)));
                                        MeanOfData = mean(EoutQ((kk-1)+NPPB/2));
                                        if MeanOfData > QLevelDec3%EyeOpenLowQ+EyeOpenQ/2                     %If it is the uper level the incoming data
                                            DataEven = [DataEven 1];                               %is 1
                                        else                                                       %If it is the lowest level the incoming data
                                            DataEven = [DataEven 0];                               %is 0
                                        end
                                        MeanOfData = mean(EoutQ((kk-1)+3*NPPB/4));
                                        if MeanOfData > QLevelDec1%EyeOpenLowQ+EyeOpenQ/2                     %If it is the uper level the incoming data
                                            DataEvenM = [DataEvenM 1];                               %is 1
                                        else                                                       %If it is the lowest level the incoming data
                                            DataEvenM = [DataEvenM 0];                               %is 0
                                        end
                                        MeanOfData = mean(EoutQ((kk-1)+NPPB/4));
                                        if MeanOfData > QLevelDec3%EyeOpenLowQ+EyeOpenQ/2                     %If it is the uper level the incoming data
                                            DataEvenL = [DataEvenL 1];                               %is 1
                                        else                                                       %If it is the lowest level the incoming data
                                            DataEvenL = [DataEvenL 0];                               %is 0
                                        end
                                        MeanOfData = mean(EoutQ((kk-1)+(NPPB/2)));
                                        if MeanOfData > LevDefDpqsk%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                                            DataEvenU = [DataEvenU 1];                                 %is 1
                                        else                                                       %If it is the lowest level the incoming data
                                            DataEvenU = [DataEvenU 0];                                 %is 0
                                        end
                                        MeanOfData = mean(EoutQ((kk-1)+(NPPB/4)));
                                        if MeanOfData > LevDefDpqsk%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                                            DataEvenU2 = [DataEvenU2 1];                                 %is 1
                                        else                                                       %If it is the lowest level the incoming data
                                            DataEvenU2 = [DataEvenU2 0];                                 %is 0
                                        end
                                        MeanOfData = mean(EoutQ((kk-1)+(3*NPPB/4)));
                                        if MeanOfData > LevDefDpqsk%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                                            DataEvenU3 = [DataEvenU3 1];                                 %is 1
                                        else                                                       %If it is the lowest level the incoming data
                                            DataEvenU3 = [DataEvenU3 0];                                 %is 0
                                        end
                                        if CalcS
                                            MeanOfData = mean(EoutI((kk-1)+NPPB/2));
                                            if MeanOfData > EyeOpenLowQ + EyeOpenQ/2%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                                                DataEvenS = [DataEvenS 1];                                 %is 1
                                            else                                                       %If it is the lowest level the incoming data
                                                DataEvenS = [DataEvenS 0];                                 %is 0
                                            end
                                        end
                                    end
                                    %%       Calculating the Bit Error Ratio (BER)
                                    %The final process here is to count the number of wrongdoings
                                    %of this whole process upon the transmited data for
                                    %quantitative analizes
                                    TxDataOdd  = TxDataMat(ThisCarr,1+SyncPeriod:end-SyncPeriod);
                                    TxDataEven = TxDataMat(ThisCarr+NumCarr,1+SyncPeriod:end-...
                                        SyncPeriod);
                                    if CalcS
                                        DataOddS  = DataOddS(1+SyncPeriod:end-SyncPeriod);
                                        DataEvenS = DataEvenS(1+SyncPeriod:end-SyncPeriod);
                                        BitErrOddS    = sum(xor(TxDataOdd,DataOddS));                      %Comparison between the Transmited and received and counting the differences
                                        BitErrEvenS   = sum(xor(TxDataEven,DataEvenS));                    %Comparison between the Transmited and received and counting the differences
                                        BerDQPSKS(CurrentTest,ThisCarr) = (BitErrOddS+BitErrEvenS)/...
                                            ((NbDQPSK)-(2*SyncPeriod));%Calculating the ration of wrong bits and the total number of bits transmited
                                    end
                                    DataOdd  = DataOdd(1+SyncPeriod:end-SyncPeriod);
                                    DataEven = DataEven(1+SyncPeriod:end-SyncPeriod);
                                    DataOddM  = DataOddM(1+SyncPeriod:end-SyncPeriod);
                                    DataEvenM = DataEvenM(1+SyncPeriod:end-SyncPeriod);
                                    DataOddL  = DataOddL(1+SyncPeriod:end-SyncPeriod);
                                    DataEvenL = DataEvenL(1+SyncPeriod:end-SyncPeriod);
                                    DataOddU  = DataOddU(1+SyncPeriod:end-SyncPeriod);
                                    DataEvenU = DataEvenU(1+SyncPeriod:end-SyncPeriod);
                                    DataOddU2  = DataOddU2(1+SyncPeriod:end-SyncPeriod);
                                    DataEvenU2 = DataEvenU2(1+SyncPeriod:end-SyncPeriod);
                                    DataOddU3  = DataOddU3(1+SyncPeriod:end-SyncPeriod);
                                    DataEvenU3 = DataEvenU3(1+SyncPeriod:end-SyncPeriod);
                                    AberLevAuxI(4) = 0;
                                    ValsLevAuxI(4) = 0.01;
                                    AberLevAuxQ(4) = 0;
                                    ValsLevAuxQ(4) = 0.01;
                                    AberLevAuxI(5) = 0;
                                    ValsLevAuxI(5) = 0.01;
                                    AberLevAuxQ(5) = 0;
                                    ValsLevAuxQ(5) = 0.01;
                                    AberLevAuxI(6) = 0;
                                    ValsLevAuxI(6) = 0.01;
                                    AberLevAuxQ(6) = 0;
                                    ValsLevAuxQ(6) = 0.01;
                                    
                                    BitErrOdd(1)  = sum(xor(TxDataOdd,DataOdd));                      %Comparison between the Transmited and received and counting the differences
                                    BitErrEven(1) = sum(xor(TxDataEven,DataEven));                    %Comparison between the Transmited and received and counting the differences
                                    BitErrOdd(2)  = sum(xor(TxDataOdd,DataOddM));                      %Comparison between the Transmited and received and counting the differences
                                    BitErrEven(2) = sum(xor(TxDataEven,DataEvenM));                    %Comparison between the Transmited and received and counting the differences
                                    BitErrOdd(3)  = sum(xor(TxDataOdd,DataOddL));                      %Comparison between the Transmited and received and counting the differences
                                    BitErrEven(3) = sum(xor(TxDataEven,DataEvenL));                    %Comparison between the Transmited and received and counting the differences
                                    BitErrOdd(4)  = sum(xor(TxDataOdd,DataOddU));                      %Comparison between the Transmited and received and counting the differences
                                    BitErrEven(4) = sum(xor(TxDataEven,DataEvenU));                    %Comparison between the Transmited and received and counting the differences
                                    BitErrOdd(5)  = sum(xor(TxDataOdd,DataOddU2));                      %Comparison between the Transmited and received and counting the differences
                                    BitErrEven(5) = sum(xor(TxDataEven,DataEvenU2));                    %Comparison between the Transmited and received and counting the differences
                                    BitErrOdd(6)  = sum(xor(TxDataOdd,DataOddU3));                      %Comparison between the Transmited and received and counting the differences
                                    BitErrEven(6) = sum(xor(TxDataEven,DataEvenU3));                    %Comparison between the Transmited and received and counting the differences
                                    BerDQPSK(CurrentTest,ThisCarr) = (min(BitErrOdd)+min(BitErrEven))/...
                                        ((NbDQPSK)-(2*SyncPeriod));%Calculating the ration of wrong bits and the total number of bits transmited
                                    
                                    AberIAux = max(AberLevAuxI(BitErrOdd<=min(BitErrOdd)));
                                    AberQAux = max(AberLevAuxQ(BitErrEven<=min(BitErrEven)));
                                    AberLevI(CurrentTest,ThisCarr)  = AberIAux;
                                    ValsLevI(CurrentTest,ThisCarr)  = max(ValsLevAuxI(AberLevAuxI==AberIAux));
                                    
                                    AberLevQ(CurrentTest,ThisCarr)  = AberQAux;
                                    ValsLevQ(CurrentTest,ThisCarr)  = max(ValsLevAuxQ(AberLevAuxQ==AberQAux));
                                    
                                    [~,LocAux] = max(AberLevAuxI==AberIAux);
                                    EvmDB(CurrentTest,ThisCarr)   = EvmDBA(LocAux);
                                    EvmPer(CurrentTest,ThisCarr)  = EvmPerA(LocAux);
                                    EvmRms(CurrentTest,ThisCarr)  = EvmRmsA(LocAux);
                                    EvmDBJ(CurrentTest,ThisCarr)  = EvmDBJA(LocAux);
                                    EvmPerJ(CurrentTest,ThisCarr) = EvmPerJA(LocAux);
                                    EvmRmsJ(CurrentTest,ThisCarr) = EvmRmsJA(LocAux);
                                    
                                    RxSymbAmos(CurrentTest,ObsCarrPos==ThisCarr,:) = EvmMatRecA(LocAux,:);%RxSymbAmosUp = [];
                                    EvmMatRec(ObsCarrPos==ThisCarr,:) = EvmMatRecA(LocAux,:);                       %Taking just the middle samples as references
                                    %% Ploting the result for qualitative analizes
                                    %PrintInfo(Ploting*29,TxDataOdd,DataOdd);
                                    %PrintInfo(Ploting*30,TxDataEven,DataEven);
                                    %%
                                    %             berpos = 1:2:size(BerDQPSK,2);
                                    %             BerDQPSK(size(BerDQPSK,1),berpos)
                                    %             a=6;
                                    close all;
                                    %                 BerDQPSK(CurrentTest,ThisCarr+1) = 1;
                                end
                            end
                        case '4PAM'
                            
                            %%              Receiver 4PAM
                            for ThisCarr=InitCarr:CarrPass:NumCarr                                           %For each carrier the same process of reception will be used.
                                if CarrUsedDo(ThisCarr)
                                    %%  Reception Process: Handling the income optical field
                                    %At the first step the income field will be selected from the
                                    %optical FFT Output with the help of the maping vector
                                    %(previously described).
                                    EoutAux = EoutAux1(VetThisCarr==ObsCarrUsed(ThisCarr),:);
                                    if AddingNoiseP
                                        %                 figure;plot(f,20*log10(abs(fftshift(fft(EoutAux)./length(EoutAux)))));
                                        %                 [~,sigpower] = MeasPower(EoutAux);
                                        [p,l] = findpeaks(20*log10(abs(fftshift(fft(EoutAux)./length(EoutAux)))),'SortStr','descend','NPeaks',1);
                                        EoutAuxA = fftshift(fft(EoutAux));
                                        EoutAuxA(f==f(l)) = EoutAuxA(f==f(l-3));
                                        %                 [testfilt,~] = FiltroGaussiano(f,12.5e9,f(l),5);           %Creating filter for selection of the received signal
                                        %                 testfilt = fftshift(testfilt);                                   %Shifting the filter for matching the received signal
                                        %                 EoutAux = ifft(fftshift(EoutAuxA.*testfilt));
                                        EoutAuxA = ifft(fftshift(EoutAuxA));
                                        %                 figure;plot(f,20*log10(abs(fftshift(fft(EoutAuxA)./length(EoutAuxA)))));
                                        %                 [~,sigpower] = MeasPower(EoutAuxA);
                                        [~,sigpower] = MeasPower(EoutAuxA);
                                        EoutAux      = awgn(EoutAux,osnrp,sigpower-30);
                                        [~,sigpower] = MeasPower(EoutAux);
                                    end
                                    %%            Fiber Time Delay Compensation
                                    switch Medium
                                        case 'Fiber'
                                            EoutAux = ifft(fft(EoutAux).*exp(1j*2*pi*f*(...
                                                FiberDelay(ThisCarr)*Ta)));
                                        otherwise
                                    end
                                    %The current incoming signal is them converted from the optical
                                    %domain to the eletrical domain with the help of an photo
                                    %detector.
                                    Ix = EoutAux.*conj(EoutAux);
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
                                    if RecFilBanPas
                                        Ix = ifft(fft(Ix).*BitFilt);
                                        Ix(1:4*(2*NumAmosCP+NPPB))         = Ix(6*(2*NumAmosCP+NPPB));
                                        Ix(end+1-4*(2*NumAmosCP+NPPB):end) = Ix(6*(2*NumAmosCP+NPPB));
                                        Ix = Ix - min(Ix);                                         %Removing the DC component from them Eletrical signal received
                                        Ix = Ix./max(abs(Ix));                                     %Normalizing the eletrical signal (amplifying what is needed)
                                    end
                                    %                 figure; plot(f,db(abs(fftshift(fft(Ix)./length(Ix))))); set(gcf,'units','normalized','outerposition',[0 0 1 1]);
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
                                    if ReceptorNoise                                               %Verify whether the noise is to be added or not
                                        [~,SigPower] = MeasPower(Ix);                              %it own received signal or
                                        SigPower2 = SigPower-30;%10*log10(SigPower);
                                        %SNR = CarSNR - 10*log10(2) + 10*log10(Nsamp) + 10*log10(1/2);
                                        SNR2 = CarSNR + 10*log10(2) - 10*log10(Nsamp) - 10*log10(1/2) + 10*log10(10^0.6);
                                        
                                        %                     if TimeSys
                                        %                         SNR = CarSNR + 10*log10(2) - 10*log10(0.25*T*BWD2);
                                        %                     else
                                        %                         SNR = CarSNR + 10*log10(2) - 10*log10(0.25*t(2*NumAmosCP+NPPB)*BWD2);
                                        %                     end
                                        
                                        Ix = awgn(Ix,SNR,'measured');
                                    end
                                    
                                    if ~RecFilBanPas
                                        Ix = ifft(fft(Ix).*BitFilt);
                                        Ix(1:4*(2*NumAmosCP+NPPB))         = Ix(6*(2*NumAmosCP+NPPB));
                                        Ix(end+1-4*(2*NumAmosCP+NPPB):end) = Ix(6*(2*NumAmosCP+NPPB));
                                        Ix = Ix - min(Ix);                                         %Removing the DC component from them Eletrical signal received
                                        
                                        Ix = Ix./max(abs(Ix));                                     %Normalizing the eletrical signal (amplifying what is needed)
                                    end
                                    
                                    %                 [BitFilt,~] = FiltroGaussiano(f,200*BWD2,CenFeq2,FiltOrd2);        %Creating filter for selection of the received signal
                                    %                 BitFilt = fftshift(BitFilt);                                   %Shifting the filter for matching the received signal
                                    %                 Ix = ifft(fft(Ix).*BitFilt);
                                    
                                    VetElecPower(CurrentTest,ThisCarr)= MeasPower(Ix);
                                    EoutAuxF = 20*log10(abs(fftshift(fft(EoutAux)./length(...
                                        EoutAux))));
                                    [VetOptiPower(CurrentTest,ThisCarr),~]= findpeaks(EoutAuxF,...
                                        'SortStr','descend','NPeaks',1);
                                    %% Ploting the result for qualitative analizes
                                    %             PrintInfo(Ploting*31,f,20.*log10(abs(fftshift(fft(Ix)./...
                                    %                                                             length(Ix)))));
                                    %             PrintInfo(Ploting*32,t,Ix);
                                    %             PrintInfo(Ploting*33,Ix,T,NPPB);
                                    %             PrintInfo(Ploting*34,t(IniSyncPos:SyncPos),Ix(IniSyncPos:...
                                    %                                     SyncPos),SyncSymb(IniSyncPos:SyncPos));
                                    %%        synchronizing
                                    %For the reception process work properly, it is needed to
                                    %sincronized the recieved signal or the sampling process will
                                    %not work.
                                    AuxSync = (Ix(IniSyncPos:SyncPos));                            %Selecting the sync-word within the received signal
                                    AuxSync1 = AuxSync;
                                    %This synchronization process is based on the mean expected
                                    %value. That means, the information mean value should be within
                                    %one period of symbol. Thus, the mean value of the received
                                    %signal is acquired and compare of the known sync-word to
                                    %verify if this mean value is at the right possition
                                    AuxSync1(AuxSync1<0)               = 0;                        %To keep the mean value above zero anything under is neglected
                                    AuxSync1(AuxSync1>=mean(AuxSync1)) = 1;                        %Adding a flag to the first sample of the received mean value
                                    AuxSync1(AuxSync1<mean(AuxSync1))  = -1;                       %All the others samples at set to the lowest level
                                    AuxSync2                           = SyncSymb(IniSyncPos:...
                                        SyncPos);%Selecting the sync-word within the known signal
                                    AuxSync2(AuxSync2>=mean(AuxSync2)) = 1;                        %Adding a flag to the first sample of the known mean value
                                    AuxSync2(AuxSync2<mean(AuxSync2))  = -1;                       %All the others samples at set to the lowest level
                                    
                                    PosToSyn  = find(ismember(AuxSync1,1));                        %Finding where is the location of the first sample to synchronize
                                    PosSyn = find(ismember(AuxSync2,1));                           %Finding where is the location of the first sample to synchronize
                                    
                                    %AuxSyncCorr = PosToSyn(round(end/2,0))-PosSyn(round(end/2,0));
                                    sn=((1/2)*(((1/2)^MaxNumStag)-1))/((1/2)-1);
                                    AuxSyncCorr = round((sn/(2/2^IfftOrSum))*T/Ta);
                                    
                                    %The difference between the PossitionTosynchronize and
                                    %Possitionsynchronized will be used to correct the time
                                    %shifting on the transmition and reception process.
                                    if AuxSyncCorr>=0%If the difference is positive, left-shift...
                                        %                 Ix = ifft(fft(Ix).*exp(1j*2*pi*f*AuxSyncCorr));          %Shift based on time change
                                        Ix = [Ix(AuxSyncCorr+1:end) Ix(1:AuxSyncCorr)];            %Shift based on sampling sliding
                                    else%... but if the difference is negative, right-shift
                                        %                 Ix = ifft(fft(Ix).*exp(-1j*2*pi*f*AuxSyncCorr));         %Shift based on time change
                                        Ix = [Ix(end+AuxSyncCorr+1:end) Ix(1:end + AuxSyncCorr)];  %Shift based on sampling sliding
                                    end
                                    if SencondAdjust
                                        %For some reason that we could not understand sometimes the
                                        %time (sampling) sliding of the signal is not equal
                                        %throught the data stream. Thus, the second part of the
                                        %synchronism process will be turn ON or OFF according to
                                        %the user's will.
                                        AuxSyncEnd     = (Ix(end-SyncPos+1:end-IniSyncPos-1));     %Selecting the sync-word within the received signal
                                        SyncSymbEndAux = SyncSymbEnd(end-SyncPos+1:end-...
                                            IniSyncPos-1);%Selecting the sync-word within the known signal
                                        AuxSyncEnd1 = AuxSyncEnd;
                                        AuxSyncEnd1(AuxSyncEnd1<0) = 0;                            %To keep the mean value above zero anything under is neglected
                                        AuxSyncEnd1(AuxSyncEnd1>=mean(AuxSyncEnd1)) = 1;           %Adding a flag to the first sample of the received mean value
                                        AuxSyncEnd1(AuxSyncEnd1<mean(AuxSyncEnd1)) = -1;           %All the others samples at set to the lowest level
                                        AuxSyncEnd2 = SyncSymbEndAux;
                                        AuxSyncEnd2(AuxSyncEnd2>=mean(AuxSyncEnd2)) = 1;           %Adding a flag to the first sample of the known mean value
                                        AuxSyncEnd2(AuxSyncEnd2<mean(AuxSyncEnd2)) = -1;           %All the others samples at set to the lowest level
                                        
                                        
                                        PosToSynEnd  = find(ismember(AuxSyncEnd1,1));              %Finding where is the location of the first sample to synchronize
                                        PosSynEnd = find(ismember(AuxSyncEnd2,1));                 %Finding where is the location of the first sample to synchronize
                                        
                                        AuxSyncCorrEnd = PosToSynEnd(round(end/2))-PosSynEnd(...
                                            round(end/2));
                                        
                                        %The difference between the PossitionTosynchronize and
                                        %Possitionsynchronized will be used to correct the time
                                        %shifting on the transmition and reception process.
                                        if AuxSyncCorrEnd>=0%If possitive difference, left-shift...
                                            %                 Ix = ifft(fft(Ix).*exp(1j*2*pi*f*AuxSyncCorrEnd));   %Shift based on time change
                                            Ix = [Ix(AuxSyncCorrEnd+1:end) Ix(1:AuxSyncCorrEnd)];  %Shift based on sampling sliding
                                        else%... but if the difference is negative, right-shift
                                            %                 Ix = ifft(fft(Ix).*exp(-1j*2*pi*f*AuxSyncCorrEnd));  %Shift based on time change
                                            Ix = [Ix(end-AuxSyncCorrEnd+1:end) Ix(1:end - ...
                                                AuxSyncCorrEnd)];%Shift based on sampling sliding
                                        end
                                    end
                                    
                                    %%          Ploting the result for qualitative analizes
                                    %                             PrintInfo(Ploting*35,t(end-SyncPos+1:end-IniSyncPos+1),Ix(...
                                    %                                 end-SyncPos+1:end-IniSyncPos+1),SyncSymb(end-SyncPos+1:...
                                    %                                                                         end-IniSyncPos+1));
                                    %% Removing CP
                                    if AddCP
                                        IxAux = Ix(1:end - StuffSampels);
                                        IxAux = reshape(IxAux,(2*NumAmosCP+NPPB),Nb4Pam/2);
                                        IxAux = IxAux(1+NumAmosCP:end-NumAmosCP,:);
                                        IxAux = reshape(IxAux,1,NPPB*Nb4Pam/2);
                                        Ix    = IxAux;
                                    end
                                    %% Taking the sampling the EVM meassurement
                                    PosAuxEout = NPPB/2:NPPB:length(Ix);%Varriable respossible to take just the samples at the middle of the symbol
                                    IxAux      = Ix(PosAuxEout);                           %Normalizing the reference
                                    RxSymbAmos(CurrentTest,ObsCarrPos==ThisCarr,:) = IxAux;%RxSymbAmosUp = [];
                                    EvmMatRec(ObsCarrPos==ThisCarr,:) = IxAux;                       %Taking just the middle samples as references
                                    [EvmDB(CurrentTest,ThisCarr), EvmPer(CurrentTest,ThisCarr), EvmRms(CurrentTest,ThisCarr) ] = EvmCalc( EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux );
                                    [EvmDBJ(CurrentTest,ThisCarr),EvmPerJ(CurrentTest,ThisCarr),EvmRmsJ(CurrentTest,ThisCarr)] = evm1(4,'pam',EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux);
                                    a=0;
                                    %                 if ThisCarr==1
                                    %                     figure;hold all;grid on;
                                    %                     plot(IxAux,zeros(1,length(IxAux)),'rx');
                                    %                     set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                    %                     a=0;
                                    %                 end
                                    %                 %%  Measuring the EVM
                                    %                 IxPosAux = NPPB/2:NPPB:length(Ix);
                                    %                 IxRec    = Ix(IxPosAux);
                                    %                 [EvmRms(CurrentTest,ThisCarr),EvmDB(CurrentTest,ThisCarr),...
                                    %                     EvmPer(CurrentTest,ThisCarr)] = EvmCalc( EvmMatRef(ThisCarr,...
                                    %                     :),IxRec );
                                    %%          Ploting the result for qualitative analizes
                                    %%         Finding Decission Levels
                                    %The process for decoding the income signal will be based on
                                    %eletronic comparators. Inasmuch as the right decission level
                                    %must be acquired for accurately decide, within a symbol
                                    %periode, what that current leavel means (ones or zeros).
                                    %
                                    %The process hereafter of chosing the  decission levels is not
                                    %deterministic rather it is a statistic process. The main idea
                                    %is to take the decission level from the histogram generated
                                    %from the income signal stream.
                                    %
                                    %Firstly it was set the interval in which the histogram will be
                                    %build. It is based on the number of samples per bit period.
                                    %             IxToSeek = reshape(Ix(1+SyncPeriod*NPPB:end-SyncPeriod*NPPB...
                                    %                                        ),NPPB,(Nb4Pam/2) - (2*SyncPeriod));
                                    Interval = linspace(min(Ix(1+SyncPeriod:end-SyncPeriod)),max...
                                        (Ix(1+SyncPeriod:end-SyncPeriod)),IntervalStep);
                                    %Therefore, the MATLAB hist function returns the number of
                                    %occurrence of each interval.
                                    EyeMax = hist(Ix,Interval);
                                    EyeMaxaux = [0 EyeMax 0];                                      %Zeros are added at the EyeMax to auxiliate the finding peaks process
                                    [~,EyeLoc] = findpeaks(EyeMaxaux,'MinPeakDistance',MinDist,...
                                        'SortStr','descend','NPeaks',4);%,'MinPeakHeight',mean(EyeMax));%The peaks on the Eye profile will be the levels at the Eyes limit
                                    if length(EyeLoc)<4
                                        [~,EyeLoc] = findpeaks(EyeMaxaux,'MinPeakDistance',...
                                            MinDist*0.8,'SortStr','descend','NPeaks',4);%,'MinPeakHeight',mean(EyeMax)/4);%The peaks on the Eye profile will be the levels at the Eyes limit
                                    end
                                    %This variable will brings the eye profile of the input signal
                                    %that will be used to generate the decision level. Basicaly the
                                    %decission level will be the minimal value of the currente eye
                                    %under evaluation plus the half of the its eye opening. The
                                    %following ilustration better describe this process.
                                    %
                                    %Eye Limits:   + (1/2 * Eye Opening:)  =    Comparisson Limits:
                                    %
                                    %UperLevel  ______     ______________     _____
                                    %                 \   /      |       \   /
                                    %                  \ /       |        \ /
                                    %                   \ Half Eye Opening /     Decission Level 3
                                    %                  / \       |        / \
                                    %LowerLevel 3_____/   \______|_______/   \_____
                                    %                 \   /      |       \   /
                                    %                  \ /       |        \ /
                                    %                   \ Half Eye Opening /     Decission Level 2
                                    %                  / \       |        / \
                                    %LowerLevel 2_____/   \______|_______/   \_____
                                    %                 \   /      |       \   /
                                    %                  \ /       |        \ /
                                    %                   \ Half Eye Opening /     Decission Level 1
                                    %                  / \       |        / \
                                    %LowerLevel 1_____/   \______|_______/   \_____
                                    %
                                    %%           Ploting for Qualitative Analizes
                                    %             PrintInfo(Ploting*37,Interval,EyeMax);
                                    %%         Finding Decission Levels
                                    %It is not always possible to totaly recover the signal.
                                    %Depending on the configuration of the transmition and
                                    %reception system the eye diagram may be nonexistent. Which
                                    %means, there will not be a profile to be found therefore the
                                    %EyeLoc will not return the correct location. Inasmuch as the
                                    %detection process to works limts will be set accordling to the
                                    %amplitude of the received signal.
                                    if length(EyeLoc)<4%If it was not able to find the eye profile.
                                        EyeLoc = [2 3 4 5];
                                        Levels = [0 0 0.35 0.55 0.85];
                                        %                 Levels = sort(Levels);
                                        %                 LevelDec1 = mean(Ix) - 2*mean(Ix)/3 ;
                                        %                 LevelDec2 = mean(Ix) ;
                                        %                 LevelDec3 = mean(Ix) + 2*mean(Ix)/3 ;
                                    else%Whereas, if there is an profile the decission can be found
                                        Levels = [Interval(EyeLoc(1)-1) Interval(EyeLoc(2)-1) ...
                                            Interval(EyeLoc(3)-1) Interval(EyeLoc(4)-1)];
                                        Levels = sort(Levels);
                                        %                 LevelDec1 = Levels(1) + (Levels(2)-Levels(1))/2 ;
                                        %                 LevelDec2 = Levels(2) + (Levels(3)-Levels(2))/2 ;
                                        %                 LevelDec3 = Levels(3) + (Levels(4)-Levels(3))/2 ;
                                    end
                                    %##############################################################
                                    if CalcS
                                        limiar3 = (4/6)*abs(max(Ix));%(1/2)*abs(max(Ix)-min(Ix))/3;
                                        limiar2 = (2.5/6)*abs(max(Ix));%((1/2)*abs(max(Ix)-min(Ix))/3) + abs(max(Ix)-min(Ix))/3;
                                        limiar1 = (1/6)*abs(max(Ix));%((1/2)*abs(max(Ix)-min(Ix))/3) + 2*(abs(max(Ix)-min(Ix))/3);
                                        limiarPos1 = Levels(1);%abs(min(Ix));%Interval>limiar1;
                                        limiarPos2 = Levels(2);%0.33*(abs(max(Ix))-abs(min(Ix)));%(Interval<=limiar1)&(Interval>limiar2);
                                        limiarPos3 = Levels(3);%0.63*(abs(max(Ix))-abs(min(Ix)));%(Interval<=limiar2)&(Interval>limiar3);
                                        limiarPos4 = Levels(4);%abs(max(Ix));%Interval<=limiar3;
                                        
                                        EyeSymMa1 = reshape(Ix(1+SyncPeriod*NPPB:end-SyncPeriod*NPPB...
                                            ),NPPB,(Nb4Pam/2) - (2*SyncPeriod));
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
                                    end
                                    
                                    %                             for ii = 1:N1
                                    %                               media(ii) = mean(olho(ii,:));
                                    %                               nhi       = find(olho(ii,:) > media(ii));
                                    %                               nlow      = find(olho(ii,:) < media(ii));
                                    %                               hi(ii)    = mean(olho(ii,nhi));
                                    %                               low(ii)   =  mean(olho(ii,nlow));
                                    %                             %
                                    %                               Lhi(ii)   = mean(hi)  - std(olho(ii,nhi));
                                    %                               Llow(ii)  = mean(low) + std(olho(ii,nlow));
                                    %                               Ab(ii)    = Lhi(ii) - Llow(ii);
                                    %                             end
                                    %                             Abertura   = mean(Ab);
                                    a=1;
                                    %##############################################################
                                    %For accurately read the income signal is necessary more than
                                    %just have the decision level, it is also needed to know where
                                    %it occurs. For instance, if we decide to use just to take the
                                    %decision level and measure the mean value across and a portion
                                    %of the time period, which portion should we take? It would be
                                    %logical to take the central portion, but the bit may be
                                    %deformed in such way that the information gets concentrated
                                    %out of the middle part, which means the symbol is not
                                    %symmetric. The symmetry information can be acquired from the
                                    %eye diagram by measuring the longitudinal opening. The
                                    %following sketch better describes this process:
                                    %
                                    %                   Point of Symmetry
                                    %         _____     _______|_______     _____
                                    %              \   /               \   /
                                    %               \ /  Longitudinal   \ /
                                    %                \ ----------------- /
                                    %               / \    Opening      / \
                                    %         _____/   \_______________/   \_____
                                    %
                                    %With those two pieces of information, decision level and point
                                    %of symmetry, we have the X, Y coordinates for the centre of
                                    %the Eye Diagram. Therefore, as long as there is an opening on
                                    %it it will be possible to recover the transmitted information
                                    %without error, theoretically.
                                    
                                    %As this process is also statistical, first we reshape the
                                    %income vector to analyze all periods at the same time.
                                    %             EyeSymMat = reshape(Ix(1+SyncPeriod*NPPB:end-SyncPeriod*NPPB...
                                    %                 ),NPPB,(Nb4Pam/2) - (2*SyncPeriod));
                                    %             EyeSymMat = EyeSymMat(1+NPPB/4:end-NPPB/4,:);
                                    %Then we take the values that compose the decision level
                                    %because they will mark the point of symmetry.
                                    %From the location of the max values that occured, which means
                                    %the uper and lower level of the eye diagram it needs to take
                                    %the actual value that those occurences represent that is
                                    %withing the the Interval variable.
                                    %             ValToSeek = Interval(EyeLoc-1);
                                    %The number of ocurrences is a statical measure therefore one
                                    %does not have control which interval will have the highest
                                    %peak, thus it is important to ordenate the values to be seek
                                    %from the lower part of the eye diagram to the uper part of the
                                    %eye diagram.
                                    %             ValToSeek = sort(ValToSeek,'ascend');
                                    %             OccuCount = zeros(1,size(EyeSymMat,1));                        %Auxiliar Variable for accounting.
                                    %             for kk=1:size(EyeSymMat,1)                                     %For every sample within a symbol period
                                    %                 if (kk>=NPPB/4)&&(kk<=NPPB-NPPB/4)
                                    %                     OccuCount(kk) = OccuCount(kk) + sum((EyeSymMat(kk,:)...
                                    %                         >=LowSymPer*ValToSeek(1))&(EyeSymMat(kk,:)<=...
                                    %                         UpeSymPer*ValToSeek(1)));%Account all occurencies of the valeu 1
                                    %                     OccuCount(kk) = OccuCount(kk) + sum((EyeSymMat(kk,:)...
                                    %                         >=LowSymPer*ValToSeek(2))&(EyeSymMat(kk,:)<=...
                                    %                         UpeSymPer*ValToSeek(2)));%Account all occurencies of the valeu 2
                                    %                     OccuCount(kk) = OccuCount(kk) + sum((EyeSymMat(kk,:)...
                                    %                         >=LowSymPer*ValToSeek(3))&(EyeSymMat(kk,:)<=...
                                    %                         UpeSymPer*ValToSeek(3)));%Account all occurencies of the valeu 3
                                    %                     OccuCount(kk) = OccuCount(kk) + sum((EyeSymMat(kk,:)...
                                    %                         >=LowSymPer*ValToSeek(4))&(EyeSymMat(kk,:)<=...
                                    %                         UpeSymPer*ValToSeek(4)));%Account all occurencies of the valeu 4
                                    %                 end
                                    %             end
                                    
                                    %The point of symmetry of the eye diagram will be where the
                                    %maximum number of occurrences were measured inasmuch as those
                                    %are the points where all the bits go to the center of the
                                    %symbol.  From the maximum number of occurrences, it can happen
                                    %for more than one sample within one symbol period, in the best
                                    %case, all samples would have the same accounting as it is
                                    %shown the ilustration above hence the symmetry will be at the
                                    %middle sample of this group of maximum occurrences. This value
                                    %can be found at the highest peak founded.
                                    %             [~,SymLoc] = findpeaks(OccuCount,'SortStr','descend');       %The peak on the Eye profile will be the Symmetry level
                                    
                                    %             PosIx = SymLoc(1):NPPB:length(Ix);
                                    
                                    
                                    
                                    
                                    
                                    %To find the Eye opening this script will concider just the
                                    %samples at the center of a symbol. Thus, it is based on the
                                    %fact that the syncronism process was done succefully.
                                    %##############################################################
                                    PosIx = NPPB/2:NPPB:length(Ix);                                %Possition of the central samples - but the number of samples per symbol is even ????
                                    IxAux = Ix(PosIx);                                             %From the main received signal just a few samples are taken for further evaluation
                                    n     = 100;                                                   %The number of boxes to be filled up on the histogram process
                                    %At this process each eye will be evaluated separetelly. It is
                                    %important to remember that the PAM4 has 4 levels, which means
                                    %three level of decissions that are exaclty the center of the
                                    %eye diagram.
                                    IxAuxAB = IxAux((IxAux<=Levels(4))&(IxAux>=Levels(3)));        %Taking just those values relative to the uper eye
                                    InterAB = linspace(Levels(3),Levels(4),n);                     %Building the histogram boxes
                                    EyeAB = hist(IxAuxAB,InterAB);                                 %filling up the boxes with samples that fit on them.
                                    %The same process described for the uper level will be done at
                                    %the middle and lower eyes levels.
                                    IxAuxCD = IxAux((IxAux<=Levels(3))&(IxAux>=Levels(2)));
                                    InterCD = linspace(Levels(2),Levels(3),n);%NPPB*2^n);
                                    EyeCD = hist(IxAuxCD,InterCD);
                                    
                                    IxAuxEF = IxAux((IxAux<=Levels(2))&(IxAux>=Levels(1)));
                                    InterEF = linspace(Levels(1),Levels(2),n);%NPPB*2^n);
                                    EyeEF = hist(IxAuxEF,InterEF);
                                    
                                    %What we are looking for here are not where exist the
                                    %occurrences of values. The eye is where there is not samples,
                                    %this means, the decission level is a reagion between actual
                                    %levels thus this reagion does not contain any signal.
                                    EyeAB = ~EyeAB;                                                %Changing zeros to one - Zeros compose the eye region
                                    EyeCD = ~EyeCD;
                                    EyeEF = ~EyeEF;
                                    %Starting variables that will be used to identify the eye
                                    %diagram.
                                    CountAB   = 1;
                                    CountCD   = 1;
                                    CountEF   = 1;
                                    
                                    SeqOnesAB = 0;
                                    SeqOnesCD = 0;
                                    SeqOnesEF = 0;
                                    
                                    SeqFinAB  = 0;
                                    SeqFinCD  = 0;
                                    SeqFinEF  = 0;
                                    
                                    SeqIniAB  = 1;
                                    SeqIniCD  = 1;
                                    SeqIniEF  = 1;
                                    %The for loop will take account of every box with ones. It is
                                    %important to take note that the not operator was used in this
                                    %vector, therefore ones means zeros (the eye diagram -
                                    %possibly) and zeros means values abouve zeroa (not the eye).
                                    for kk=1:length(EyeAB)                                         %For every box
                                        if EyeAB(kk)                                               %if it contains "1"
                                            SeqOnesAB(SeqIniAB)=CountAB;                           %count this element as part of a consecutive sequency
                                            CountAB = CountAB + 1;                                 %adds one to the counter of consecutive elements "1"
                                            if kk==length(EyeAB)                                   %if the current box is the last box we got to an end
                                                SeqFinAB(SeqIniAB) = kk;                           %The final sequency element is equal to its possition (kk)
                                            end
                                        else                                                       %else if the current box contains "0"
                                            SeqFinAB(SeqIniAB) = kk-1;                             %the previous element was the last element of a consecutive sequency of ones
                                            SeqIniAB = SeqIniAB + 1;                               %adds one to the counter of consecutive number (take in account how many sequencies there are)
                                            CountAB = 1;                                           %reset the counter of consecutive elements to mark the start of a new consecutive sequence
                                        end
                                        
                                        if EyeCD(kk)
                                            SeqOnesCD(SeqIniCD)=CountCD;
                                            CountCD = CountCD + 1;
                                            if kk==length(EyeCD)
                                                SeqFinCD(SeqIniCD) = kk;
                                            end
                                        else
                                            SeqFinCD(SeqIniCD) = kk-1;
                                            SeqIniCD = SeqIniCD + 1;
                                            CountCD = 1;
                                        end
                                        
                                        if EyeEF(kk)
                                            SeqOnesEF(SeqIniEF)=CountEF;
                                            CountEF = CountEF + 1;
                                            if kk==length(EyeEF)
                                                SeqFinEF(SeqIniEF) = kk;
                                            end
                                        else
                                            SeqFinEF(SeqIniEF) = kk-1;
                                            SeqIniEF = SeqIniEF + 1;
                                            CountEF = 1;
                                        end
                                    end
                                    
                                    
                                    %If the eye is open, which means there is a clear difference
                                    %between adjacent levels, the eye diagram will be the longest
                                    %sequence of ones.
                                    [MaxValAB,LocMaxAB]=max(SeqOnesAB);
                                    if LocMaxAB<1 || MaxValAB<2                                    %if any sequency was found or there is just one sequency it is a error thus
                                        LevDec3 = 0.6953;                                            %the decission level will be by default 0.67. Also, the other variables will
                                        LocMaxAB = 1;                                              %will be set with values to not cause errors in the future
                                        SeqFinAB(1)=2;
                                        MaxValAB = 0;
                                        InterAB(1)=LevDec3;
                                    else                                                           %if a sequency was found the middle of the eye will be the middle of the sequency
                                        if (SeqFinAB(LocMaxAB)-MaxValAB/2)<1                       %if for some reason the final element of the sequency minus the half of its
                                            LevDec3 = 0.6953;                                         %length results in a negative value, something went very wrong, and by
                                        else                                                       %default it will be set to 0.7
                                            LevDec3 = InterAB(round(SeqFinAB(LocMaxAB)-MaxValAB/...
                                                2));%Otherwise, the decission level is the middle of the sequency
                                        end
                                    end
                                    [MaxValCD,LocMaxCD]=max(SeqOnesCD);
                                    if LocMaxCD<1 || MaxValCD<2
                                        LevDec2 = 0.4013;
                                        LocMaxCD = 1;
                                        SeqFinCD(1)=2;
                                        MaxValCD = 0;
                                        InterCD(1)=LevDec2;
                                    else
                                        if (SeqFinCD(LocMaxCD)-MaxValCD/2)<1
                                            LevDec2 = 0.4013;
                                        else
                                            LevDec2 = InterCD(round(SeqFinCD(LocMaxCD)-MaxValCD/2));
                                        end
                                    end
                                    
                                    [MaxValEF,LocMaxEF]=max(SeqOnesEF);
                                    if LocMaxEF<1 || MaxValEF<2
                                        LevDec1 = 0.1877;
                                        LocMaxEF = 1;
                                        SeqFinEF(1)=2;
                                        MaxValEF = 0;
                                        InterEF(1)=LevDec1;
                                    else
                                        if (SeqFinEF(LocMaxEF)-MaxValEF/2)<1
                                            LevDec1 = 0.1877;
                                        else
                                            LevDec1 = InterEF(round(SeqFinEF(LocMaxEF)-MaxValEF/2));
                                        end
                                        
                                    end
                                    
                                    %another way to measure the eye opening is the get all the
                                    %boxes and find all peaks on it, that will be a plato created
                                    %by the sequences of ones (which was zeros). From thos peaks,
                                    %the eye diagram will be the longer of them hence it will take
                                    %the most part of the vector that store the findpeaks result.
                                    %Thus, the middle of the eye will be basically the middle of
                                    %the peaks vector.
                                    LocAB = find(EyeAB);
                                    if isempty(LocAB)                                              %if for some reason there are no peaks, something went wrong.
                                        LocAB = 1;
                                        LevelDec3 = 0.6953;%mean(Levels(3:4));                       %by default the decission level will be set to 0.65
                                    else
                                        if DecMod%InterAB(LocAB(round(end/2)))<=LevDec3
                                            LevelDec3 = LevDec3;
                                        else
                                            LevelDec3 = InterAB(LocAB(round(end/2)));
                                        end
                                    end
                                    
                                    LocCD = find(EyeCD);
                                    if isempty(LocCD)
                                        LocCD = 1;
                                        LevelDec2 = 0.4013;%mean(Levels(2:3));
                                    else
                                        if DecMod%InterCD(LocCD(round(end/2)))>=LevDec2
                                            LevelDec2 = LevDec2;
                                        else
                                            LevelDec2 = InterCD(LocCD(round(end/2)));
                                        end
                                    end
                                    
                                    LocEF = find(EyeEF);
                                    if isempty(LocEF)
                                        LocEF = 1;
                                        LevelDec1 = 0.1877;%mean(Levels(1:2));
                                    else
                                        if DecMod%InterEF(LocEF(round(end/2)))<=LevDec1
                                            LevelDec1 = LevDec1;
                                        else
                                            LevelDec1 = InterEF(LocEF(round(end/2)));
                                        end
                                    end
                                    
                                    AberLev(1,ThisCarr,CurrentTest)  = abs(InterAB(SeqFinAB(LocMaxAB)-1) - InterAB(SeqFinAB(LocMaxAB)-MaxValAB+1));
                                    AberLev(2,ThisCarr,CurrentTest)  = abs(InterCD(SeqFinCD(LocMaxCD)-1) - InterCD(SeqFinCD(LocMaxCD)-MaxValCD+1));
                                    AberLev(3,ThisCarr,CurrentTest)  = abs(InterEF(SeqFinEF(LocMaxEF)-1) - InterEF(SeqFinEF(LocMaxEF)-MaxValEF+1));
                                    ValsLev(1,ThisCarr,CurrentTest)  = LevDec3;
                                    ValsLev(2,ThisCarr,CurrentTest)  = LevDec2;
                                    ValsLev(3,ThisCarr,CurrentTest)  = LevDec1;
                                    ValsLev2(1,ThisCarr,CurrentTest) = InterAB(LocAB(round(size(LocAB,2)/2)));
                                    ValsLev2(2,ThisCarr,CurrentTest) = InterCD(LocCD(round(size(LocCD,2)/2)));
                                    ValsLev2(3,ThisCarr,CurrentTest) = InterEF(LocEF(round(size(LocEF,2)/2)));
                                    %%           Ploting for Qualitative Analizes
                                    if ThisCarr==126
                                        EyeToPlot(CurrentTest,1:length(Esync)) = Esync;
                                        save(['savingforplotingeye' num2str(VetSnr)],'EyeToPlot','VetSnr','IfftOrSum','UsedModula','T','NPPB');
                                    end
                                    if PrintinEye
                                        PrintInfo(PlotingThis*36,Ix,T,NPPB);
                                        hold all;
                                        plot(t(NPPB/2),LevDec1,'b*');plot(t(NPPB/2),InterEF(SeqFinEF(LocMaxEF)-1),'bv');plot(t(NPPB/2),InterEF(round(SeqFinEF(LocMaxEF)-MaxValEF+1)),'b^');
                                        plot(t(NPPB/2),LevDec2,'gd');plot(t(NPPB/2),InterCD(SeqFinCD(LocMaxCD)-1),'gv');plot(t(NPPB/2),InterCD(round(SeqFinCD(LocMaxCD)-MaxValCD+1)),'g^');
                                        plot(t(NPPB/2),LevDec3,'kd');plot(t(NPPB/2),InterAB(SeqFinAB(LocMaxAB)-1),'kv');plot(t(NPPB/2),InterAB(round(SeqFinAB(LocMaxAB)-MaxValAB+1)),'k^');
                                        if CalcS
                                            plot(t(NPPB/2),ThisDeciLev3,'b*');plot(t(NPPB/2),mean(LevHi3),'bv');plot(t(NPPB/2),mean(LevLo3),'b^');
                                            plot(t(NPPB/2),ThisDeciLev2,'k*');plot(t(NPPB/2),mean(LevHi2),'kv');plot(t(NPPB/2),mean(LevLo2),'k^');
                                            plot(t(NPPB/2),ThisDeciLev1,'g*');plot(t(NPPB/2),mean(LevHi1),'gv');plot(t(NPPB/2),mean(LevLo1),'g^');
                                        end
                                        figure;
                                        hold all;
                                        plotpos = zeros(1,length(IxAux));
                                        plot(IxAux,plotpos,'o','color',[1 0.4 0]);
                                        plot(EvmMatRef(ObsCarrUsed(ThisCarr),:),plotpos,'b*');
                                        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                        a=1;
                                        %                     drawnow;
                                        %                     pause(0.1);
                                    end
                                    %                             plot(t(NPPB/2),InterAB(LocAB(round(end/2))),'k^');
                                    %                             plot(t(NPPB/2),InterCD(LocCD(round(end/2))),'g^');
                                    %                             plot(t(NPPB/2),InterEF(LocEF(round(end/2))),'b^');
                                    %                             plot([t(SymLoc(1)) t(SymLoc(1))],[0 1]);
                                    %
                                    %             drawnow;
                                    %             PrintInfo(Ploting*38,EyeMax,EyeLoc-1,Interval);
                                    % %                         figure;plot(SeqOnesEF);figure;plot(SeqOnesCD);figure;plot(SeqOnesAB)
                                    %             LevDec3
                                    %             a1=InterAB(LocAB(round(end/2)))
                                    %             b1=ThisDeciLev1
                                    %             LevDec2
                                    %             a2=InterCD(LocCD(round(end/2)))
                                    %             b2=ThisDeciLev2
                                    %             LevDec1
                                    %             c1=InterEF(LocEF(round(end/2)))
                                    %             c2=ThisDeciLev3
                                    %             PrintInfo(Ploting*39,NPPB,OccuCount);
                                    %             PrintInfo(Ploting*40,OccuCount,SymLoc);
                                    %%      Actualy Receiving Data
                                    %Once the signal was processed the next step is through a
                                    %comparator decide the actual information received.
                                    IxRec    = [];                                                 %Initialization of the vector that will store the income data
                                    IxRecDef = [];                                                 %Initialization of the vector that will store the income data
                                    IxRecDeS = [];                                                 %Initialization of the vector that will store the income data
                                    for kk=1:NPPB:length(Ix)                                       %The comparison process will be made for each symbol period
                                        %                 midaux = round(mean(SymLoc(1:round(end/2))));
                                        midaux = NPPB/2;%SymLoc(1);
                                        aux1 = Ix((kk-1)+midaux);     %An small portion of the income signal is take for evaluation
                                        MeanRec = mean(aux1);                                      %Measuring the avarage value of the samples taken
                                        %Verifying the interval for each symbol received.
                                        if MeanRec <= LevelDec1                                    %If it is the lowest level the incoming data
                                            IxRec = [IxRec 0 0];                                   %is 01 (1)
                                        elseif (MeanRec <= LevelDec2)&&(MeanRec > LevelDec1)       %If it is the second level the incoming data
                                            IxRec = [IxRec 0 1];                                   %is 00 (0)
                                        elseif (MeanRec <= LevelDec3)&&(MeanRec > LevelDec2)       %If it is the tird level the incoming data
                                            IxRec = [IxRec 1 1];                                   %is 10 (2)
                                        elseif MeanRec > LevelDec3                                 %If it is the uper level the incoming data
                                            IxRec = [IxRec 1 0];                                   %is 11 (3)
                                        else                                                       %If for some misteriose reason neither of previous verification were sucedded
                                            IxRec = [IxRec 0 0];                                   %by default the current data is set to be 00 (0)
                                        end
                                        
                                        %Verifying the interval for each symbol received.
                                        if MeanRec <= DecLevDef1                                   %If it is the lowest level the incoming data
                                            IxRecDef = [IxRecDef 0 0];                             %is 01 (1)
                                        elseif (MeanRec <= DecLevDef2)&&(MeanRec > DecLevDef1)     %If it is the second level the incoming data
                                            IxRecDef = [IxRecDef 0 1];                             %is 00 (0)
                                        elseif (MeanRec <= DecLevDef3)&&(MeanRec > DecLevDef2)     %If it is the tird level the incoming data
                                            IxRecDef = [IxRecDef 1 1];                             %is 10 (2)
                                        elseif MeanRec > DecLevDef3                                %If it is the uper level the incoming data
                                            IxRecDef = [IxRecDef 1 0];                             %is 11 (3)
                                        else                                                       %If for some misteriose reason neither of previous verification were sucedded
                                            IxRecDef = [IxRecDef 0 0];                             %by default the current data is set to be 00 (0)
                                        end
                                        
                                        if CalcS
                                            if MeanRec <= ThisDeciLev1                                   %If it is the lowest level the incoming data
                                                IxRecDeS = [IxRecDeS 0 0];                             %is 01 (1)
                                            elseif (MeanRec <= ThisDeciLev2)&&(MeanRec > ThisDeciLev1)     %If it is the second level the incoming data
                                                IxRecDeS = [IxRecDeS 0 1];                             %is 00 (0)
                                            elseif (MeanRec <= ThisDeciLev3)&&(MeanRec > ThisDeciLev2)     %If it is the tird level the incoming data
                                                IxRecDeS = [IxRecDeS 1 1];                             %is 10 (2)
                                            elseif MeanRec > ThisDeciLev3                                %If it is the uper level the incoming data
                                                IxRecDeS = [IxRecDeS 1 0];                             %is 11 (3)
                                            else                                                       %If for some misteriose reason neither of previous verification were sucedded
                                                IxRecDeS = [IxRecDeS 0 0];                             %by default the current data is set to be 00 (0)
                                            end
                                        end
                                    end
                                    %%           Ploting for Qualitative Analizes
                                    %             PrintInfo(Ploting*41,t(length(TxDataMat(ThisCarr,:))),Nb4Pam...
                                    %                                            /2,TxDataMat(ThisCarr,:),IxRec);
                                    %%       Calculating the Bit Error Ratio (BER)
                                    %The final process here is to count the number of wrongdoings
                                    %of this whole process upon the transmited data for
                                    %quantitative analizes
                                    BitErr = sum(xor(TxDataMat(ThisCarr,1+2*SyncPeriod:end-2*...
                                        SyncPeriod),IxRec(1+2*SyncPeriod:end-2*SyncPeriod)));%Comparison between the Transmited and received and counting the differences
                                    BitErrAux1 = BitErr;
                                    BitErrAux2 = sum(xor(TxDataMat(ThisCarr,1+2*SyncPeriod:end-...
                                        2*SyncPeriod),IxRecDef(1+2*SyncPeriod:end-2*SyncPeriod)));%Comparison between the Transmited and received and counting the differences
                                    if BitErr ~= 0
                                        if BitErrAux2<BitErrAux1
                                            BitErr = BitErrAux2;
                                        end
                                    else
                                        DecLevDef1 = LevelDec1;
                                        DecLevDef2 = LevelDec2;
                                        DecLevDef3 = LevelDec3;
                                    end
                                    if CalcS
                                        BitErrS = sum(xor(TxDataMat(ThisCarr,1+2*SyncPeriod:end-2*...
                                            SyncPeriod),IxRecDeS(1+2*SyncPeriod:end-2*SyncPeriod)));%Comparison between the Transmited and received and counting the differences
                                        Ber4PAMS(CurrentTest,ThisCarr) = BitErrS/(Nb4Pam-(2*SyncPeriod));%Calculating the ration of wrong bits and the total number of bits transmited
                                    end
                                    if BitErrS<BitErr
                                        Ber4PAM(CurrentTest,ThisCarr) = BitErrS/(Nb4Pam-(2*SyncPeriod));%Calculating the ration of wrong bits and the total number of bits transmited
                                    else
                                        Ber4PAM(CurrentTest,ThisCarr) = BitErr/(Nb4Pam-(2*SyncPeriod));%Calculating the ration of wrong bits and the total number of bits transmited
                                    end
                                    %             if BitErr>0
                                    %                 berpos = 1:2:size(Ber4PAM,2);
                                    %                 Ber4PAM(size(Ber4PAM,1),berpos)
                                    %                 drawnow;%pause(1);
                                    %                 a=2;
                                    %             end
                                    %             berpos = 1:2:size(Ber4PAM,2);
                                    %             Ber4PAM(size(Ber4PAM,1),berpos)
                                    %             drawnow;%pause(1);
                                    %             a=3;
                                    %             drawnow;
                                    %             pause(0.5);
                                    close all;
                                    %                             Ber4PAM(CurrentTest,ThisCarr+1) = 1;
                                    
                                end
                            end
                        otherwise
                            %%              Receiver OOK
                            %         AberLev(CurrentTest,:) = zeros(1,NumCarr);
                            for ThisCarr=InitCarr:CarrPass:NumCarr                                           %For each carrier the same process of reception will be used.
                                if CarrUsedDo(ThisCarr)
                                    %%  Reception Process: Handling the income optical field
                                    %At the first step the income field will be selected from the
                                    %optical FFT Output with the help of the maping vector
                                    %(previously described).
                                    EoutAux = EoutAux1(VetThisCarr==ObsCarrUsed(ThisCarr),:);
                                    if AddingNoiseP
                                        %                 figure;plot(f,20*log10(abs(fftshift(fft(EoutAux)./length(EoutAux)))));
                                        %                 [~,sigpower] = MeasPower(EoutAux);
                                        [p,l] = findpeaks(20*log10(abs(fftshift(fft(EoutAux)./length(EoutAux)))),'SortStr','descend','NPeaks',1);
                                        EoutAuxA = fftshift(fft(EoutAux));
                                        EoutAuxA(f==f(l)) = EoutAuxA(f==f(l-3));
                                        %                 [testfilt,~] = FiltroGaussiano(f,12.5e9,f(l),5);           %Creating filter for selection of the received signal
                                        %                 testfilt = fftshift(testfilt);                                   %Shifting the filter for matching the received signal
                                        %                 EoutAux = ifft(fftshift(EoutAuxA.*testfilt));
                                        EoutAuxA = ifft(fftshift(EoutAuxA));
                                        %                 figure;plot(f,20*log10(abs(fftshift(fft(EoutAuxA)./length(EoutAuxA)))));
                                        %                 [~,sigpower] = MeasPower(EoutAux);
                                        [~,sigpower] = MeasPower(EoutAuxA);
                                        EoutAux      = awgn(EoutAux,osnrp,sigpower-30);
                                        [~,sigpower] = MeasPower(EoutAux);
                                    end
                                    %%            Fiber Time Delay Compensation
                                    switch Medium
                                        case 'Fiber'
                                            EoutAux = ifft(fft(EoutAux).*exp(1j*2*pi*f*(...
                                                FiberDelay((ThisCarr))*Ta)));
                                        otherwise
                                    end
                                    
                                    %The current incoming signal is them converted from the optical
                                    %domain to the eletrical domain with the help of an photo
                                    %detector.
                                    Ix =EoutAux.*conj(EoutAux);
                                    %%           Creating the Reception Filter
                                    [BitFilt,~] = FiltroGaussiano(f,BWD,CenFeq,FiltOrd);           %Creating filter for selection of the received signal
                                    BitFilt = fftshift(BitFilt);                                   %Shifting the filter for matching the received signal
                                    %%
                                    if RecFilBanPas
                                        Ix = ifft(fft(Ix).*BitFilt);
                                        Ix = Ix - min(Ix);                                         %Removing the DC component from them Eletrical signal received
                                        Ix = Ix./max(abs(Ix));                                     %Normalizing the eletrical signal (amplifying what is needed)
                                    end
                                    
                                    if ReceptorNoise
                                        [~,SigPower] = MeasPower(Ix);
                                        SigPower2 = SigPower-30;%10*log10(SigPower);
                                        SNR2 = CarSNR + 10*log10(1) - 10*log10(Nsamp) - 10*log10(1/2) + 10*log10(10^0.3);
                                        
                                        %                     if TimeSys
                                        %                         SNR = CarSNR + 10*log10(1) - 10*log10(0.25*T*BWD);
                                        %                     else
                                        %                         SNR = CarSNR + 10*log10(1) - 10*log10(0.25*t(2*NumAmosCP+NPPB)*BWD);
                                        %                     end
                                        
                                        Ix = awgn(Ix,SNR,SigPower2);
                                    end
                                    if ~RecFilBanPas
                                        Ix = ifft(fft(Ix).*BitFilt);
                                        Ix = Ix - min(Ix);                                         %Removing the DC component from them Eletrical signal received
                                        Ix = Ix./max(abs(Ix));                                     %Normalizing the eletrical signal (amplifying what is needed)
                                    end
                                    
                                    VetElecPower(CurrentTest,ThisCarr)= MeasPower(Ix);
                                    EoutAuxF = 20*log10(abs(fftshift(fft(EoutAux)./length(...
                                        EoutAux))));
                                    [VetOptiPower(CurrentTest,ThisCarr),~]= findpeaks(EoutAuxF,...
                                        'SortStr','descend','NPeaks',1);
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
                                    
                                    AuxSync = (Ix(IniSyncPos:SyncPos));                            %Selecting the sync-word within the received signal
                                    AuxSync1 = AuxSync;
                                    %This synchronization process is based on the mean expected
                                    %value. That means, the information mean value should be within
                                    %one period of symbol. Thus, the mean value of the received
                                    %signal is acquired and compare of the known sync-word to
                                    %verify if this mean value is at the right possition
                                    AuxSync1(AuxSync1<0) = 0;                                      %To keep the mean value above zero anything under is neglected
                                    AuxSync1(AuxSync1>=mean(AuxSync1)) = 1;                        %Adding a flag to the first sample of the received mean value
                                    AuxSync1(AuxSync1<mean(AuxSync1)) = -1;                        %All the others samples at set to the lowest level
                                    AuxSync2 = SyncSymb(IniSyncPos:SyncPos);                       %Selecting the sync-word within the known signal
                                    AuxSync2(AuxSync2>=mean(AuxSync2)) = 1;                        %Adding a flag to the first sample of the known mean value
                                    AuxSync2(AuxSync2<mean(AuxSync2)) = -1;                        %All the others samples at set to the lowest level
                                    
                                    PosToSyn  = find(ismember(AuxSync1,1));
                                    PosSyn    = find(ismember(AuxSync2,1));
                                    
                                    AuxSyncCorr = PosToSyn(round(end/2))-PosSyn(round(end/2));
                                    sn=((1/2)*(((1/2)^MaxNumStag)-1))/((1/2)-1);
                                    AuxSyncCorr = round((sn/(2/2^IfftOrSum))*T/Ta);
                                    %The difference between the PossitionTosynchronize and
                                    %Possitionsynchronized will be used to correct the time
                                    %shifting on the transmition and reception process.
                                    
                                    if AuxSyncCorr>=0%If the difference is positive, left-shift...
                                        %                 Ix = ifft(fft(Ix).*exp(1j*2*pi*f*(AuxSyncCorr*Ta)));     %Shift based on time change
                                        Ix = [Ix(AuxSyncCorr+1:end) Ix(1:AuxSyncCorr)];            %Shift based on sampling sliding
                                    else%... but if the difference is negative, right-shift
                                        %                 Ix = ifft(fft(Ix).*exp(-1j*2*pi*f*(AuxSyncCorr*Ta)));    %Shift based on time change
                                        Ix = [Ix(end+AuxSyncCorr+1:end) Ix(1:end + AuxSyncCorr)];  %Shift based on sampling sliding
                                    end
                                    
                                    
                                    %%        synchronizing
                                    if SencondAdjust
                                        %For some reason that we could not understand sometimes the
                                        %time (sampling) sliding of the signal is not equal
                                        %throught the data stream. Thus, the second part of the
                                        %synchronism process will be turn ON or OFF according to the
                                        %user's will.
                                        AuxSyncEnd = (Ix(end-(SyncPos+1*abs(AuxSyncCorr))+1:end-...
                                            (2*NPPB+1*abs(AuxSyncCorr))));%Selecting the sync-word within the received signal
                                        SyncSymbEndAux = SyncSymbEnd(end-SyncPos+1:end-2*NPPB);    %Selecting the sync-word within the known signal
                                        AuxSyncEnd1 = AuxSyncEnd;
                                        AuxSyncEnd1(AuxSyncEnd1<0) = 0;                            %To keep the mean value above zero anything under is neglected
                                        AuxSyncEnd1(AuxSyncEnd1>=mean(AuxSyncEnd1)) = 1;           %Adding a flag to the first sample of the received mean value
                                        AuxSyncEnd1(AuxSyncEnd1<mean(AuxSyncEnd1)) = -1;           %All the others samples at set to the lowest level
                                        AuxSyncEnd2 = SyncSymbEndAux;
                                        AuxSyncEnd2(AuxSyncEnd2>=mean(AuxSyncEnd2)) = 1;           %Adding a flag to the first sample of the known mean value
                                        AuxSyncEnd2(AuxSyncEnd2<mean(AuxSyncEnd2)) = -1;           %All the others samples at set to the lowest level
                                        
                                        PosToSynEnd  = find(ismember(AuxSyncEnd1,1));
                                        PosSynEnd = find(ismember(AuxSyncEnd2,1));
                                        
                                        AuxSyncCorrEnd = PosToSynEnd(round(end/2))-PosSynEnd(...
                                            round(end/2));
                                        
                                        %The difference between the PossitionTosynchronize and
                                        %Possitionsynchronized will be used to correct the time
                                        %shifting on the transmition and reception process.
                                        
                                        if AuxSyncCorrEnd>=0%If possitive difference, left-shift...
                                            %                     Ix = ifft(fft(Ix).*exp(1j*2*pi*f*(AuxSyncCorrEnd*...
                                            %                                                                    Ta)));%Shift based on time change
                                            Ix = [Ix(AuxSyncCorrEnd+1:end) Ix(1:AuxSyncCorrEnd)];
                                        else%... but if the difference is negative, right-shift
                                            %                     Ix = ifft(fft(Ix).*exp(-1j*2*pi*f*(AuxSyncCorrEnd*...
                                            %                                                                    Ta)));%Shift based on time change
                                            Ix = [Ix(end-AuxSyncCorrEnd+1:end) Ix(1:end - ...
                                                AuxSyncCorrEnd)];%Shift based on sampling sliding
                                        end
                                    end
                                    
                                    %% Ploting the result for qualitative analizes
                                    
                                    %                             PrintInfo(Ploting*45,t(IniSyncPos:SyncPos),Ix(IniSyncPos:SyncPos),SyncSymb(IniSyncPos:SyncPos));
                                    %                             PrintInfo(Ploting*46,t(1:length(Ix)),Ix);
                                    %%                  Removing CP
                                    if AddCP
                                        IxAux = Ix(1:end - StuffSampels);
                                        IxAux = reshape(IxAux,(2*NumAmosCP+NPPB),Nb-NumBitDesc);
                                        IxAux = IxAux(1+NumAmosCP:end-NumAmosCP,:);
                                        IxAux = reshape(IxAux,1,NPPB*(Nb-NumBitDesc));
                                        Ix    = IxAux;
                                    end
                                    %% Taking the sampling the EVM meassurement
                                    %                 evm = comm.EVM('MaximumEVMOutputPort',true,...
                                    %                 'XPercentileEVMOutputPort',true, 'XPercentileValue',90,...
                                    %                 'SymbolCountOutputPort',true);
                                    PosAuxEout = NPPB/2:NPPB:length(Ix);                       %Varriable respossible to take just the samples at the middle of the symbol
                                    IxAux      = Ix(PosAuxEout);                               %Normalizing the reference
                                    a=0;
                                    RxSymbAmos(CurrentTest,ObsCarrPos==ThisCarr,:) = IxAux;%RxSymbAmosUp = [];
                                    EvmMatRec(ObsCarrPos==ThisCarr,:) = IxAux;                       %Taking just the middle samples as references
                                    [EvmDB(CurrentTest,ThisCarr), EvmPer(CurrentTest,ThisCarr), EvmRms(CurrentTest,ThisCarr) ] = EvmCalc( EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux );
                                    [EvmDBJ(CurrentTest,ThisCarr),EvmPerJ(CurrentTest,ThisCarr),EvmRmsJ(CurrentTest,ThisCarr)] = evm1(2,'pam',EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux);
                                    %                 IxRef = EvmMatRef(ObsCarrPos==ThisCarr,:);
                                    %                 [MatEvem,maxEVM,pctEVM,numSym] = evm(abs(IxRef),abs(IxAux));
                                    %                 EvmRmsM(CurrentTest,ThisCarr) = MatEvem;
                                    %                 EvmPerM(CurrentTest,ThisCarr) = MatEvem*100;
                                    %                 EvmDBM(CurrentTest,ThisCarr) = 20*log10(MatEvem);
                                    %                 if ThisCarr==1
                                    %                     figure;hold all;grid on;
                                    %                     plot(IxAux,zeros(1,length(IxAux)),'rx');
                                    %                     set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                    %                     a=0;
                                    %                 end
                                    %%  Measuring the EVM
                                    %                 IxPosAux = NPPB/2:NPPB:length(Ix);
                                    %                 IxRec    = Ix(IxPosAux);
                                    %                 [EvmRms(CurrentTest,ThisCarr),EvmDB(CurrentTest,ThisCarr),...
                                    %                     EvmPer(CurrentTest,ThisCarr)] = EvmCalc( EvmMatRef(ThisCarr,...
                                    %                     :),IxRec );
                                    %##########################################################################
                                    PosIx = NPPB/2:NPPB:length(Ix);                                %Possition of the central samples - but the number of samples per symbol is even ????
                                    IxAux = Ix(PosIx);                                             %From the main received signal just a few samples are taken for further evaluation
                                    n     = 100;                                                   %The number of boxes to be filled up on the histogram process
                                    MinDist = n/3;
                                    Perce = 0.5;
                                    %At this process each eye will be evaluated separetelly. It is
                                    %important to remember that the PAM4 has 4 levels, which means
                                    %three level of decissions that are exaclty the center of the
                                    %eye diagram.
                                    Interval = linspace(Perce*(min(Ix)),Perce*(max(Ix)),n);                     %Building the histogram boxes
                                    
                                    %Therefore, the MATLAB hist function returns the number of
                                    %occurrence of each interval.
                                    EyeMax = hist(Ix,Interval);
                                    EyeMaxaux = [0 EyeMax 0];                                      %Zeros are added at the EyeMax to auxiliate the finding peaks process
                                    [~,EyeLoc] = findpeaks(EyeMaxaux,'MinPeakDistance',MinDist,'SortStr','descend','NPeaks',2);%The peaks on the Eye profile will be the levels at the Eyes limit
                                    
                                    if length(EyeLoc)<2%If it was not able to find the eye profile.
                                        [~,EyeLoc] = findpeaks(EyeMaxaux,'MinPeakDistance',MinDist/2,'SortStr','descend','NPeaks',2);%The peaks on the Eye profile will be the levels at the Eyes limit
                                    end
                                    
                                    if length(EyeLoc)<2%If it was not able to find the eye profile.
                                        EyeLoc = [2 3];
                                        Levels = [0 -0.2 0.35];
                                    else%Whereas, if there is an profile the decission can be found
                                        Levels = [Interval(EyeLoc(1)-1) Interval(EyeLoc(2)-1)];
                                        Levels = sort(Levels);
                                        %                 LevelDec1 = Levels(1) + (Levels(2)-Levels(1))/2 ;
                                        %                 LevelDec2 = Levels(2) + (Levels(3)-Levels(2))/2 ;
                                        %                 LevelDec3 = Levels(3) + (Levels(4)-Levels(3))/2 ;
                                    end
                                    IxAuxAB = IxAux((IxAux<=Levels(2))&(IxAux>=Levels(1)));
                                    Inter = linspace(Levels(1),Levels(2),n);                     %Building the histogram boxes
                                    Inter2 = linspace(min(IxAux),max(IxAux),n);                     %Building the histogram boxes
                                    EyeAB = hist(IxAuxAB,Inter);                                 %filling up the boxes with samples that fit on them.
                                    EyeCD = hist(IxAux,Inter2);                                 %filling up the boxes with samples that fit on them.
                                    
                                    %What we are looking for here are not where exist the
                                    %occurrences of values. The eye is where there is not samples,
                                    %this means, the decission level is a reagion between actual
                                    %levels thus this reagion does not contain any signal.
                                    EyeAB = ~EyeAB;                                                %Changing zeros to one - Zeros compose the eye region
                                    EyeCD = ~EyeCD;                                                %Changing zeros to one - Zeros compose the eye region
                                    %Starting variables that will be used to identify the eye
                                    %diagram.
                                    Count   = 1;
                                    Count2   = 1;
                                    
                                    SeqOnes = 0;
                                    SeqOnes2 = 0;
                                    
                                    SeqIni  = 1;
                                    SeqIni2  = 1;
                                    
                                    SeqFin  = 0;
                                    SeqFin2  = 0;
                                    %The for loop will take account of every box with ones. It is
                                    %important to take note that the not operator was used in this
                                    %vector, therefore ones means zeros (the eye diagram -
                                    %possibly) and zeros means values abouve zeroa (not the eye).
                                    for kk=1:length(EyeAB)                                         %For every box
                                        if EyeAB(kk)                                               %if it contains "1"
                                            SeqOnes(SeqIni)=Count;                           %count this element as part of a consecutive sequency
                                            Count = Count + 1;                                 %adds one to the counter of consecutive elements "1"
                                            if kk==length(EyeAB)                                   %if the current box is the last box we got to an end
                                                SeqFin(SeqIni) = kk;                           %The final sequency element is equal to its possition (kk)
                                            end
                                        else                                                       %else if the current box contains "0"
                                            SeqFin(SeqIni) = kk-1;                             %the previous element was the last element of a consecutive sequency of ones
                                            SeqIni = SeqIni + 1;                               %adds one to the counter of consecutive number (take in account how many sequencies there are)
                                            Count = 1;                                           %reset the counter of consecutive elements to mark the start of a new consecutive sequence
                                        end
                                        
                                        if EyeCD(kk)                                               %if it contains "1"
                                            SeqOnes2(SeqIni2)=Count2;                           %count this element as part of a consecutive sequency
                                            Count2 = Count2 + 1;                                 %adds one to the counter of consecutive elements "1"
                                            if kk==length(EyeCD)                                   %if the current box is the last box we got to an end
                                                SeqFin2(SeqIni2) = kk;                           %The final sequency element is equal to its possition (kk)
                                            end
                                        else                                                       %else if the current box contains "0"
                                            SeqFin2(SeqIni2) = kk-1;                             %the previous element was the last element of a consecutive sequency of ones
                                            SeqIni2 = SeqIni2 + 1;                               %adds one to the counter of consecutive number (take in account how many sequencies there are)
                                            Count2 = 1;                                           %reset the counter of consecutive elements to mark the start of a new consecutive sequence
                                        end
                                    end
                                    
                                    
                                    %If the eye is open, which means there is a clear difference
                                    %between adjacent levels, the eye diagram will be the longest
                                    %sequence of ones.
                                    [MaxVal,LocMax]=max(SeqOnes);
                                    if LocMax<2 || MaxVal<2                                    %if any sequency was found or there is just one sequency it is a error thus
                                        LevDec = 0.0;                                            %the decission level will be by default 0.67. Also, the other variables will
                                        LocMax = 1;                                              %will be set with values to not cause errors in the future
                                        SeqFin(1)=2;
                                        MaxVal = 0;
                                        Inter(1)=LevDec;
                                    else                                                           %if a sequency was found the middle of the eye will be the middle of the sequency
                                        if (SeqFin(LocMax)-MaxVal/2)<1                       %if for some reason the final element of the sequency minus the half of its
                                            LevDec = 0.0;                                         %length results in a negative value, something went very wrong, and by
                                        else                                                       %default it will be set to 0.7
                                            LevDec = Inter(round(SeqFin(LocMax)-MaxVal/2));%Otherwise, the decission level is the middle of the sequency
                                        end
                                    end
                                    [MaxVal2,LocMax2]=max(SeqOnes2);
                                    if LocMax2<2 || MaxVal2<2                                    %if any sequency was found or there is just one sequency it is a error thus
                                        LevDec2 = 0.0;                                            %the decission level will be by default 0.67. Also, the other variables will
                                        LocMax2 = 1;                                              %will be set with values to not cause errors in the future
                                        SeqFin2(1)=2;
                                        MaxVal2 = 0;
                                        Inter2(1)=LevDec2;
                                    else                                                           %if a sequency was found the middle of the eye will be the middle of the sequency
                                        if (SeqFin2(LocMax2)-MaxVal2/2)<1                       %if for some reason the final element of the sequency minus the half of its
                                            LevDec2 = 0.0;                                         %length results in a negative value, something went very wrong, and by
                                        else                                                       %default it will be set to 0.7
                                            LevDec2 = Inter2(round(SeqFin2(LocMax2)-MaxVal2/2));%Otherwise, the decission level is the middle of the sequency
                                        end
                                    end
                                    
                                    %another way to measure the eye opening is the get all the
                                    %boxes and find all peaks on it, that will be a plato created
                                    %by the sequences of ones (which was zeros). From thos peaks,
                                    %the eye diagram will be the longer of them hence it will take
                                    %the most part of the vector that store the findpeaks result.
                                    %Thus, the middle of the eye will be basically the middle of
                                    %the peaks vector.
                                    Loc = find(EyeAB);
                                    if isempty(Loc)                                              %if for some reason there are no peaks, something went wrong.
                                        Loc = 1;
                                        LevelDec = 0.0;%mean(Levels(3:4));                       %by default the decission level will be set to 0.65
                                    else
                                        LevelDec = LevDec;
                                    end
                                    Loc2 = find(EyeCD);
                                    if isempty(Loc2)                                              %if for some reason there are no peaks, something went wrong.
                                        Loc2 = 1;
                                        LevelDec2 = 0.0;%mean(Levels(3:4));                       %by default the decission level will be set to 0.65
                                    else
                                        LevelDec2 = LevDec2;
                                    end
                                    %##########################################################################
                                    
                                    if CalcS
                                        [~,~,EyeOpen,EyeOpenHigh,EyeOpenLow] = Olho((Ix),T,NPPB,0);
                                        AberLevS(CurrentTest,ThisCarr)= EyeOpen;
                                        ValsLevS(CurrentTest,ThisCarr)= EyeOpenLow + EyeOpen/2;
                                    else
                                        [~,~,EyeOpen,EyeOpenHigh,EyeOpenLow] = Olho((Ix),T,NPPB,0);
                                    end
                                    %% Ploting the result for qualitative analizes
                                    
                                    if ThisCarr==126
                                        EyeToPlot(CurrentTest,1:length(Ix)) = Ix;
                                        save(['savingforplotingeye' num2str(VetSnr)],'EyeToPlot','VetSnr','IfftOrSum','UsedModula','T','NPPB');
                                    end
                                    if PrintinEye
                                        PrintInfo(Ploting*47,(Ix),T,NPPB);
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
                                        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                        figure;
                                        hold all;
                                        plotpos = zeros(1,length(IxAux));
                                        plot(IxAux,plotpos,'o','color',[1 0.4 0]);
                                        plot(EvmMatRef(ObsCarrUsed(ThisCarr),:),plotpos,'b*');
                                        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                        a=0;
                                    end
                                    AberLev(CurrentTest,ThisCarr)= Inter(SeqFin(LocMax))-Inter(round(SeqFin(LocMax)-MaxVal));
                                    ValsLev(CurrentTest,ThisCarr)= LevDec;
                                    %%         Finding Decission Levels
                                    %The process for decoding the income signal will be based on
                                    %eletronic comparators. Inasmuch as the right decission level
                                    %must be acquired for accurately decide, within a symbol
                                    %periode, what that current leavel means (ones or zeros).
                                    %
                                    %The process hereafter of chosing the  decission levels is not
                                    %deterministic rather it is a statistic process. The main idea
                                    %is to take the decission level from the histogram generated
                                    %from the income signal stream.
                                    %
                                    %This process is realized inside the function Olho.
                                    %
                                    %Basicaly the decission level will be the minimal value of the
                                    %currente eye under evaluation plus the half of the its eye
                                    %opening.The following ilustration better describe this process
                                    %
                                    %Eye Limits:   + (1/2 * Eye Opening:)  =    Comparisson Limit:
                                    %
                                    %UperLevel  ______     ______________     _____
                                    %                 \   /      |       \   /
                                    %                  \ /       |        \ /
                                    %                   \ Half Eye Opening /     Decission Level
                                    %                  / \       |        / \
                                    %LowerLevel ______/   \______|_______/   \_____
                                    %
                                    %
                                    %
                                    %For accurately read the income signal is necessary more than
                                    %just have the decision level, it is also needed to know where
                                    %it occurs. For instance, if we decide to use just to take the
                                    %decision level and measure the mean value across and a portion
                                    %of the time period, which portion should we take? It would be
                                    %logical to take the central portion, but the bit may be
                                    %deformed in such way that the information gets concentrated
                                    %out of the middle part, which means the symbol is not
                                    %symmetric. The symmetry information can be acquired from the
                                    %eye diagram by measuring the longitudinal opening. The
                                    %following sketch better describes this process:
                                    %
                                    %                   Point of Symmetry
                                    %         _____     _______|_______     _____
                                    %              \   /               \   /
                                    %               \ /  Longitudinal   \ /
                                    %                \ ----------------- /
                                    %               / \    Opening      / \
                                    %         _____/   \_______________/   \_____
                                    %
                                    %With those two pieces of information, decision level and point
                                    %of symmetry, we have the X, Y coordinates for the centre of
                                    %the Eye Diagram. Therefore, as long as there is an opening on
                                    %it it will be possible to recover the transmitted information
                                    %without error... theoretically.
                                    %
                                    %As this process is also statistical, first we reshape the
                                    %income vector to analyze all periods at the same time.
                                    
                                    %                 EyeSymMat = reshape(Ix(1+SyncPeriod*NPPB:end-SyncPeriod*...
                                    %                     NPPB),NPPB,(Nb-NumBitDesc)-2*SyncPeriod);
                                    %Then we take the values that compose the decision level
                                    %because they will mark the point of symmetry.
                                    %
                                    %Firstly it was set the interval in which the histogram will be
                                    %build. It is based on the number of samples per bit period.
                                    %                 Interval = linspace(min(Ix),max(Ix),2*NPPB);
                                    %Therefore, the MATLAB hist function returns the number of
                                    %occurrence of each interval.
                                    %                 EyeMax = hist(Ix,Interval);
                                    %                 EyeMaxaux = [0 EyeMax 0];                                      %Zeros are added at the EyeMax to auxiliate the finding peaks process
                                    %                 [~,EyeLoc] = findpeaks(EyeMaxaux,'MinPeakDistance',NPPB/4,...
                                    %                     'SortStr','descend','NPeaks',2);%The peaks on the Eye profile will be the levels at the Eyes limit
                                    
                                    %             figure;findpeaks(EyeMaxaux,'MinPeakDistance',NPPB/4,'SortStr','descend','NPeaks',2);
                                    
                                    %                 ValToSeek = Interval(EyeLoc-1);
                                    %From the location of the max values that occured, which means
                                    %the uper and lower level of the eye diagram it needs to take
                                    %the actual value that those occurences represent that is
                                    %withing the the Interval variable.
                                    %                 ValToSeek = sort(ValToSeek,'ascend');
                                    %The number of ocurrences is a statical measure therefore one
                                    %does not have control which interval will have the highest
                                    %peak, thus it is important to ordenate the values to be seek
                                    %from the lower part of the eye diagram to the uper part of the
                                    %eye diagram.
                                    %                 OccuCount = zeros(1,size(EyeSymMat,1));                        %Auxiliar Variable for accounting.
                                    %                 for kk=1:size(EyeSymMat,1)                                     %For every sample within a symbol period
                                    %                     OccuCount(kk) = OccuCount(kk)+sum((EyeSymMat(kk,:)>=min(...
                                    %                         Ix))&(EyeSymMat(kk,:)<=UpeSymPer*EyeOpenLow));%Account all occurencies of the value 1
                                    %                     OccuCount(kk) = OccuCount(kk)+sum((EyeSymMat(kk,:)>=...
                                    %                         LowSymPer*EyeOpenHigh)&(EyeSymMat(kk,:)<=max(Ix)));%Account all occurencies of the value 2
                                    %                 end
                                    %The point of symmetry of the eye diagram will be where the
                                    %maximum number of occurrences were measured inasmuch as those
                                    %are the points where all the bits go to the center of the
                                    %symbol.  From the maximum number of occurrences, it can happen
                                    %for more than one sample within one symbol period, in the best
                                    %case, all samples would have the same accounting as it is
                                    %shown the ilustration above hence the symmetry will be at the
                                    %middle sample of this group of maximum occurrences. This value
                                    %can be found by the mean of the samples positions within the
                                    %symbol period. The problem with this approach is that the
                                    %signal must be synchronized with the maximum displacement of
                                    %a symbol period minus 25% of the eye Longitudinal opening if
                                    %the displacement is higher than that the point of symmetry
                                    %will be wrongly measured.
                                    %                 [SymLoc] = round(mean(find(ismember(OccuCount,max(OccuCount)...
                                    %                     ))));%The peak on the Eye profile will be the Symmetry level
                                    %%           Ploting for Qualitative Analizes
                                    %             PrintInfo(Ploting*48,NPPB,OccuCount);
                                    %             PrintInfo(Ploting*49,OccuCount,SymLoc);
                                    
                                    %%      Actualy Receiving Data
                                    %Once the signal was processed the next step is through a
                                    %comparator decide the actual information received.
                                    EoutCorr = [];                                                 %Initialization of the vector that will store the income data
                                    EoutCorrD = [];                                                 %Initialization of the vector that will store the income data
                                    EoutCorr2 = [];                                                 %Initialization of the vector that will store the income data
                                    EoutCorrS = [];                                                 %Initialization of the vector that will store the income data
                                    for kk=1:NPPB:length(Ix)                                       %The comparison process will be made for each symbol period
                                        %An small portion of the income signal is take for
                                        %evaluation by measuring the avarage value of the samples
                                        %taken
                                        %                 CalcMean = mean((Ix((kk-1)+SymLoc(1))));
                                        CalcMean = mean((Ix((kk-1)+NPPB/2)));
                                        %Verifying the interval for each symbol received.
                                        if CalcMean >= LevDec%EyeOpenLow+EyeOpen/2                        %If it is the uper level the incoming data
                                            EoutCorr = [EoutCorr 1];                               %is 1
                                        else                                                       %If it is the lowest level the incoming data
                                            EoutCorr = [EoutCorr 0];                               %is 0
                                        end
                                        CalcMean2 = mean((Ix((kk-1)+NPPB/2)));
                                        %Verifying the interval for each symbol received.
                                        if CalcMean2 >= LevDec2%EyeOpenLow+EyeOpen/2                        %If it is the uper level the incoming data
                                            EoutCorr2 = [EoutCorr2 1];                               %is 1
                                        else                                                       %If it is the lowest level the incoming data
                                            EoutCorr2 = [EoutCorr2 0];                               %is 0
                                        end
                                        if CalcMean2 >= EyeOpenLow+EyeOpen/2                        %If it is the uper level the incoming data
                                            EoutCorrD = [EoutCorrD 1];                               %is 1
                                        else                                                       %If it is the lowest level the incoming data
                                            EoutCorrD = [EoutCorrD 0];                               %is 0
                                        end
                                        if CalcS
                                            CalcMeanS = mean((Ix((kk-1)+NPPB/2)));
                                            %Verifying the interval for each symbol received.
                                            if CalcMeanS >= EyeOpenLow + EyeOpen/2%EyeOpenLow+EyeOpen/2                        %If it is the uper level the incoming data
                                                EoutCorrS = [EoutCorrS 1];                               %is 1
                                            else                                                       %If it is the lowest level the incoming data
                                                EoutCorrS = [EoutCorrS 0];                               %is 0
                                            end
                                        end
                                    end
                                    
                                    %PrintInfo(Ploting*50,linspace(0,t(end),length(EoutCorr))...
                                    %                                           ,TxDataMat(ThisCarr,:),EoutCorr);
                                    %%       Calculating the Bit Error Ratio (BER)
                                    %The final process here is to count the number of wrongdoings
                                    %of this whole process upon the transmited data for
                                    %quantitative analizes
                                    if CalcS
                                        BitErrS = sum(xor(TxDataMat(ThisCarr,1+SyncPeriod:end-...
                                            SyncPeriod),EoutCorrS(1+SyncPeriod:end-SyncPeriod)));%Comparison between the Transmited and received and counting the differences
                                        BerOOKS(CurrentTest,ThisCarr) = BitErrS/((Nb-NumBitDesc)-...
                                            2*SyncPeriod);%Calculating the ration of wrong bits and the total number of bits transmited
                                    end
                                    BitErr = sum(xor(TxDataMat(ThisCarr,1+SyncPeriod:end-...
                                        SyncPeriod),EoutCorr(1+SyncPeriod:end-SyncPeriod)));%Comparison between the Transmited and received and counting the differences
                                    BitErr2 = sum(xor(TxDataMat(ThisCarr,1+SyncPeriod:end-...
                                        SyncPeriod),EoutCorr2(1+SyncPeriod:end-SyncPeriod)));%Comparison between the Transmited and received and counting the differences
                                    BitErrD = sum(xor(TxDataMat(ThisCarr,1+SyncPeriod:end-...
                                        SyncPeriod),EoutCorrD(1+SyncPeriod:end-SyncPeriod)));%Comparison between the Transmited and received and counting the differences
                                    if BitErr2<=BitErr
                                        if BitErr2<=BitErrD
                                            BerOOK(CurrentTest,ThisCarr) = BitErr2/((Nb-NumBitDesc)-...
                                                2*SyncPeriod);%Calculating the ration of wrong bits and the total number of bits transmited
                                            AberLev(CurrentTest,ThisCarr)= Inter2(SeqFin2(LocMax2))-Inter2(round(SeqFin2(LocMax2)-MaxVal2));
                                            ValsLev(CurrentTest,ThisCarr)= LevDec2;
                                        else
                                            BerOOK(CurrentTest,ThisCarr) = BitErrD/((Nb-NumBitDesc)-...
                                                2*SyncPeriod);%Calculating the ration of wrong bits and the total number of bits transmited
                                            AberLev(CurrentTest,ThisCarr)= EyeOpen;
                                            ValsLev(CurrentTest,ThisCarr)= EyeOpenLow+EyeOpen/2;
                                        end
                                    else
                                        if BitErr<=BitErrD
                                            BerOOK(CurrentTest,ThisCarr) = BitErr/((Nb-NumBitDesc)-...
                                                2*SyncPeriod);%Calculating the ration of wrong bits and the total number of bits transmited
                                        else
                                            BerOOK(CurrentTest,ThisCarr) = BitErrD/((Nb-NumBitDesc)-...
                                                2*SyncPeriod);%Calculating the ration of wrong bits and the total number of bits transmited
                                            AberLev(CurrentTest,ThisCarr)= EyeOpen;
                                            ValsLev(CurrentTest,ThisCarr)= EyeOpenLow+EyeOpen/2;
                                        end
                                    end
                                    %                             berpos = 1:2:size(BerOOK,2);
                                    %                             BerOOK(size(BerOOK,1),berpos)
                                    %                             a=1;
                                    close all;
                                    %                             BerOOK(CurrentTest,ThisCarr+1) = 1;
                                    
                                end
                                %         end
                            end
                    end
                    
                    if exist('BerDPSK','var')
                        if size(BerDPSK,2)<NumCarr
                            BerDPSK(1,end+1) = 0;
                        end
                    elseif exist('BerDQPSK','var')
                        if size(BerDQPSK,2)<NumCarr
                            BerDQPSK(1,end+1) = 0;
                        end
                    elseif exist('Ber4PAM','var')
                        if size(Ber4PAM,2)<NumCarr
                            Ber4PAM(1,end+1) = 0;
                        end
                    elseif exist('BerOOK','var')
                        if size(BerOOK,2)<NumCarr
                            BerOOK(1,end+1) = 0;
                        end
                    end
                    a=0;
                    
                    % clearvars -except UsedModula TestBundle ThisModula CurrentTest TxDataMat...
                    % CurrentModula ThisTest ThisModula OfcName BerOOK Ber4PAM BerDQPSK OSNRPC...
                    % MaxNumStag EoutUpRec BerDPSK CurrentMedium OcsToTest CurrentOCS t f T...
                    %                                                                   PulseResp
                    
                    
                    %                     end
                else
                    
                    [EoutAux1,~,VetThisCarr]=OpticalFFTN(t,T,MaxNumStag,EoutRec);
                end
            end
            if (UseOfdmElect)&&(UsedModula==5)&&(NumCarr==2)
            else
                if SendingUp
                    %%              Datum Upstream
                    %This section is resposible to generate and to send the
                    %information from the user to the OLT.
                    %% DatumUpstream;
                    
                    
                    
                    [EoutAux1,~,VetThisCarr]=OpticalFFTN(t,T,MaxNumStag,EoutRec);
                    clear EoutRec
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
                        case 'OFDM' %NOT IMPLEMENTED YET
                            for kk=InitCarrUp:CarrPass:NumCarr                             %Generating different data for each carrier
                                TxSig = 0;
                                %                 CpPaFr = 1;
                                %                 for CarrOffSet=-1*(NumFraPar-1)/2:1:(NumFraPar-1)/2
                                for CarDmtM=1:length(DmtMve)
                                    M = DmtMve(CarDmtM);
                                    TxData          = randi([0 M-1],NumFra,1);                %Generation random information
                                    TxDataMat(kk,1+NumFra*(CarDmtM-1):NumFra*CarDmtM) = TxData(:);                               %Saving data for future comparison
                                    switch OfdMod                                              %Sellect which modulation was used
                                        case 'qam'
                                            TxSymbAux(1:NumFra,CarDmtM)  = qammod(TxData,M);                    %Modulating information by QAM
                                        otherwise
                                            TxSymbAux(1:NumFra,CarDmtM)  = dpskmod(TxData,M);                   %Modulating information DPSK as default
                                    end
                                end
                                TxSymb      = TxSymbAux.';
                                TxSigOfdm(CurrentTest,ObsCarrPos==kk,:) = TxSymb(:);
                                EvmMatRef(ObsCarrPos==kk,:) = TxSymb(:);                   %Taking just the middle samples as references
                                %Construction the Hermitian Simetry
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
                                    TxSymbH = [zeros(1,NumFra);TxSymb;zeros(1,NumFra);TxSymbConj];%Hermitian Symetri
                                    TxSymbC = zeros(NFFT,NumFra);
                                    TxSymbC(1+Zt-ZtC:Zt-ZtC+size(TxSymbH,1)/2,:) = TxSymbH(1:size(TxSymbH,1)/2,:);
                                    TxSymbC(end-Zt+ZtC-(size(TxSymbH,1)/2)+1:end-Zt+ZtC,:) = TxSymbH(1+size(TxSymbH,1)/2:end,:);
                                else
                                    TxSymbH = TxSymb;
                                    %     TxSymbC = [zeros((NFFT-size(TxSymbH,1))/2,NumFra);TxSymbH;zeros((NFFT-size(TxSymbH,1))/2,NumFra)];
                                    TxSymbC = [TxSymbH(1:size(TxSymbH,1)/2,:);zeros((NFFT-size(TxSymbH,1)),NumFra);TxSymbH(end+1-size(TxSymbH,1)/2:end,:)];
                                end
                                %                 plot(TxSymbC);
                                %Here was implemented the modulation through FFT process, it basica mix-up
                                %all infomation following Fourier transform.
                                TxSymbModA       = ifft(TxSymbC,NFFT);
                                
                                switch SelModTp
                                    case 'BP'                                              %Base Band transmission
                                        TxSymbModA  = rectpulse(TxSymbModA,OvSam);         %Over sampling
                                        %                         TxSymbMod = TxSymbModA;
                                        tt = linspace(0,1*Te,size(TxSymbModA,1));
                                        for jj=1:NumFra
                                            tta(:,jj) = tt;
                                        end
                                        TxSymbMod   = TxSymbModA;%.*exp(1j*2*pi*CarrOffSet*FramSpac*tta);
                                        ff = time2freq(tt);
                                        if SelecGaus
                                            [ ModFilt ] = FiltroGaussiano(ff,OBw/2,0,1);
                                        else
                                            [ ModFilt ] = Filtro_Retangular(OBw,0,ff);
                                        end
                                        ModFilt = fftshift(ModFilt);
                                        if (kk==2)&&(Ploting)
                                            figure;hold all;
                                            plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                                        end
                                        for jj=1:NumFra
                                            ModFilta(:,jj) = ModFilt;
                                        end
                                        TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
                                        if (kk==2)&&(Ploting)
                                            plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                                            plot(ff,20*log10(abs(fftshift(ModFilt))));
                                            set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                        end
                                    case 'AM'                                              %Out of base band
                                        TxSymbModA  = rectpulse(TxSymbModA,OvSam);         %Over sampling
                                        tt = linspace(0,1*Te,size(TxSymbModA,1));
                                        ff = time2freq(tt);
                                        for jj=1:NumFra
                                            tta(:,jj) = tt;
                                        end
                                        TxSymbMod   = TxSymbModA.*cos(2*pi*Ofc*tta);
                                        if SelecGaus
                                            [ ModFilt ] = FiltroGaussiano(ff,OBw,Ofc,1);
                                        else
                                            [ ModFilt ] = Filtro_Retangular( OBw,Ofc,ff);
                                        end
                                        ModFilt = fftshift(ModFilt);
                                        if (kk==2)&&(Ploting)
                                            figure;hold all;
                                            plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                                        end
                                        for jj=1:NumFra
                                            ModFilta(:,jj) = ModFilt;
                                        end
                                        TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
                                        if (kk==2)&&(Ploting)
                                            plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                                            plot(ff,20*log10(abs(fftshift(ModFilt))));
                                            set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                        end
                                    case 'AMSSB'                                           %Out of base band
                                        TxSymbModA  = rectpulse(TxSymbModA,OvSam);         %Over sampling
                                        tt = linspace(0,1*Te,size(TxSymbModA,1));
                                        ff = time2freq(tt);
                                        for jj=1:NumFra
                                            tta(:,jj) = tt;
                                        end
                                        TxSymbMod   = TxSymbModA.*exp(1j*2*pi*Ofc*tta);
                                        if SelecGaus
                                            [ ModFilt ] = FiltroGaussiano(ff,OBw,Ofc,1);
                                        else
                                            [ ModFilt ] = Filtro_Retangular( OBw,Ofc,ff);
                                        end
                                        ModFilt = fftshift(ModFilt);
                                        if (kk==2)&&(Ploting)
                                            figure;hold all;
                                            plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                                        end
                                        for jj=1:NumFra
                                            ModFilta(:,jj) = ModFilt;
                                        end
                                        TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
                                        if (kk==2)&&(Ploting)
                                            plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                                            plot(ff,20*log10(abs(fftshift(ModFilt))));
                                            set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                        end
                                    otherwise
                                        TxSymbModA  = rectpulse(TxSymbModA,OvSam);         %Over sampling
                                        [TxSymbMod,tt] = modulate(TxSymbModA,Ofc,Ofs,'pm');
                                        ff = time2freq(tt);
                                        if SelecGaus
                                            [ ModFilt ] = FiltroGaussiano(ff,3*Ofc,0,1);
                                        else
                                            [ ModFilt ] = Filtro_Retangular( 6*Ofc,0,ff);
                                        end
                                        ModFilt = fftshift(ModFilt);
                                        if (kk==2)&&(Ploting)
                                            figure;hold all;
                                            plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                                        end
                                        for jj=1:NumFra
                                            ModFilta(:,jj) = ModFilt;
                                        end
                                        TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
                                        if (kk==2)&&(Ploting)
                                            plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                                            plot(ff,20*log10(abs(fftshift(ModFilt))));
                                            set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                        end
                                end
                                
                                
                                AuxTxSymbMod(1,:)  = TxSymbMod(:);
                                TxSymbModA  = [TxSymbMod(1+end-NPOFEX:end,:);TxSymbMod];
                                SigTx  = rectpulse(TxSymbModA,NPPOF);                               %Over sampling
                                
                                TxSigA = SigTx(:).';
                                TxSig = TxSig + [TxSigA(1:NPSTUf) TxSigA];                         %Conforming vector sizes to the same length
                                %                     CpPaFr = CpPaFr + 1;
                                %                 end
                                %Assigning the eletrical signal to one drive of the MZM - The aplitude of
                                %the signal can be controlled by the variable DatGai, which can be
                                %understood as an gain for the eletrical signal or an atenuation. The
                                %second signal will be similar with the only difference a phase shift of pi
                                NormFact = max((TxSig));
                                TxSig = 0.95*(TxSig./NormFact);                     %Normalizing the sinal to mach with the MZM espect to receive.
                                U.U1t = TxSig;
                                U.U2t = exp(-1j*pi).*TxSig;
                                %As both signals will have the mostrly the same characteristics with the
                                %only difference the phase shift of 180 degress. The MZM-I will be working
                                %on the Push-Pull configuration. It is necessary to reduce the Chirp noise
                                %to zero.
                                [EoutModAux,~]=Mach_Zehnder_Modulator_simplificado_2(freqGHz,...
                                    EoutTx(VetThisCarrTx==(RefCarr-1+kk),:),U,MZ_Input_File);%Modulating individual carriers
                                EoutModTem(ObsCarrPos(kk),1:length(EoutModAux)) = ...
                                    EoutModAux;%Adding the current Modulated output to the final OutPut
                                if (kk==2)&&(Ploting)
                                    figure;hold all;grid on;
                                    plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
                                    axis([-25e9 50e9 -200 0]);
                                    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                end
                            end
                            if IfftOrSum
                                if size(EoutModTem,1)>1
                                    [EoutMod,~,~] = OpticalIFFTN(t,T,MaxNumStagT,EoutModTem);
                                else
                                    EoutMod = EoutModTem;
                                end
                            else
                                if size(EoutModTem,1)>1
                                    EoutMod = sum(EoutModTem);
                                else
                                    EoutMod = EoutModTem;
                                end
                            end
                        case 'DPSK'
                            %%            Generate the data DPSK
                            EoutMod    = 0;
                            EoutModTem = 0;
                            for kk=InitCarrUp:CarrPass:NumCarr%Generating different data for each carrier
                                if CarrUsedUp(kk)
                                    % Frist it is chosen to transmite a random information...
                                    TxData = (randi(2,1,NbDPSK)-1);                                %Creating the data stream to be transmited
                                    TxData(1:JusPos) = JusVal;                                     %Adding the Justification simble at the begining of the stream to sincronize received data frame
                                    TxData(end-(JusPos-1):end) = JusValEnd;                        %Adding the Justification simble at the end of the stream to sincronize received data frame
                                    
                                    TxDataMat(kk,:)   = TxData;                                    %Storring the transmited information for latter evaluation
                                    
                                    %Once that we have the data to be transmited it needs to be
                                    %converted to the eletrical signal that hereafter will modulate
                                    %our carrier through the MZM-I
                                    SigTx           = DpskEncodEq( TxData);
                                    SigTx(SigTx==1) =  1;
                                    SigTx(SigTx==0) = -1;
                                    SigTx           = rectpulse(SigTx,NPPB);
                                    
                                    % Adding CP to the data
                                    if AddCP
                                        TxAux = reshape(SigTx,NPPB,NbDPSK);
                                        TxAux = [TxAux(1:NumAmosCP,:);TxAux;TxAux(end-(...
                                            NumAmosCP-1):end,:)];
                                        TxAux = reshape(TxAux,1,(2*NumAmosCP+NPPB)*NbDPSK);
                                        SigTx = [TxAux TxAux(end-(StuffSampels-1):end)];
                                    end
                                    
                                    %Because of reasons, the real world does not have components
                                    %capable instantly change the voltage level. Thus, to emulate
                                    %this behaviour the eletrical PAM signal pass through a
                                    %gaussian filter that will remove the frequencies of higher
                                    %order, which will result in small slope in the level variation
                                    
                                    [BitFilt,~] = FiltroGaussiano(f,BWD,CenFeq,FiltOrd);           %Creating filter for conformation of the input information
                                    BitFilt     = fftshift(BitFilt);                               %Doing a shift on the vector for matching the transmited data
                                    %                 TxSig       = ifft(fft(SigTx).*BitFilt);                       %Conforming the information and Creating the modulation signal
                                    TxSig       = SigTx;                                   %Adding a gain to the eletrical signal
                                    %             PrintInfo(Ploting*2,TxSig,T,NPPB);
                                    %                 TxSigMat(kk,:) = TxSig;                                        %Storing the eletrical signal for further evaluation
                                    %                 a=4;
                                    
                                    EoutT = EoutAux1(VetThisCarr==ObsCarrUsed(kk),:);
                                    %Assigning the eletrical signal to one drive of the MZM -
                                    %The aplitude of the signal can be controlled by the
                                    %variable DatGai, which can be understood as an gain for
                                    %the eletrical signal or an atenuation. The second signal
                                    %will be similar with the only difference a phase shift of
                                    %pi.
                                    U.U1t  = TxSig;
                                    U.U2t  = exp(-1j*pi).*TxSig;
                                    %As both signals will have the mostrly the same
                                    %characteristics with the only difference the phase shift
                                    %of 180 degress. The MZM-I will be working on the Push-Pull
                                    %configuration. It is necessary to reduce the Chirp noise
                                    %to zero.
                                    [EoutModAux,~]=Mach_Zehnder_Modulator_simplificado_2(freqGHz,...
                                        EoutT,U,MZ_Input_File);    %Modulating individual carriers
                                    %                 EoutMod = EoutMod + EoutModAux;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                                    EoutModTem(ObsCarrPos(kk),1:length(EoutModAux)) = EoutModAux;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                                    %% Taking the sampling the EVM meassurement
                                    PhaDel          = 0;
                                    TimDel          = T;
                                    [ESync1,ESync2] = DelayInterf(t,TimDel,PhaDel,EoutModAux); %Coverting phase shift to amplitude variation
                                    ESync1 = ESync1.*conj(ESync1);
                                    ESync2 = ESync2.*conj(ESync2);
                                    Esync  = ESync2-ESync1;
                                    PosAuxEout      = (2*NumAmosCP+NPPB)/2:(2*NumAmosCP+...
                                        NPPB):(length(Esync)-StuffSampels);%Varriable respossible to take just the samples at the middle of the symbol
                                    Esync           = ifft(fft(Esync).*EvmFilt);           %Removing higher signal generated by the receiving process
                                    %                     Esync  = Esync-min(Esync);
                                    Esync  = Esync-mean(Esync);
                                    Esync  = Esync./max(abs(Esync));
                                    IxAux  = Esync(PosAuxEout);
                                    a=0;
                                    TxSymbAmos(CurrentTest,ObsCarrPos==kk,:) = IxAux;%RxSymbAmosUp = [];
                                    EvmMatRef(ObsCarrPos==kk,:) = IxAux;                   %Taking just the middle samples as references
                                    %                 if (kk==2)&&(Ploting)
                                    %                     figure;hold all;grid on;
                                    %                     plot(IxAux,zeros(1,length(IxAux)),'rx');
                                    %                     set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                    %                     a=0;
                                    %                 end
                                else
                                    EoutT = EoutAux1(VetThisCarr==ObsCarrPos(kk),:);
                                    EoutModAux = EoutT;
                                    %                 EoutMod = EoutMod + EoutModAux;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                                    EoutModTem(ObsCarrPos(kk),1:length(EoutModAux)) = EoutModAux;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                                end
                            end
                            if IfftOrSum
                                %             Ekaux=sum(EoutModTem);
                                %             figure;hold all;
                                %             plot(f,20*log10(abs(fftshift(fft(Ekaux)./length(Ekaux)))));
                                %             set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                %             if mod(size(EoutModTem,1),2)
                                %                 EoutModTem(end+1,:) = 0;
                                %             end
                                [EoutMod,~,~] = OpticalIFFTN(t,T,MaxNumStagT,EoutModTem);
                                %             plot(f,20*log10(abs(fftshift(fft(EoutMod)./length(EoutMod)))));
                                %             axis([0.9*fin 1.2*NumCarr*fc+fin -280 -40]);
                                %             a=1;
                            else
                                EoutMod = sum(EoutModTem);
                            end
                        case 'DQPSK'
                            %%            Generate the data DQPSK
                            EoutMod    = 0;
                            EoutModTem = 0;
                            for kk=InitCarrUp:CarrPass:NumCarr%Generating different data for each carrier
                                if CarrUsedUp(kk)
                                    % Frist it is chosen to transmite a random information...
                                    TxData                     = (randi(2,1,NbDQPSK)-1);           %Creating the data stream to be transmited
                                    TxData(1:JusPos)           = JusVal;                           %Adding the Justification simble at the begining of the stream to sincronize received data frame
                                    TxData(end-(JusPos-1):end) = JusValEnd;                        %Adding the Justification simble at the end of the stream to sincronize received data frame
                                    PreBits                    = [0 0];                            %Seting the start point for the DQPSK maping
                                    
                                    TxDataPos       = linspace(1,NbDQPSK,NbDQPSK);                 %Auxiliar variable to split the TxData
                                    [DataI,DataQ]   = DqpskEncodEq(TxData,PreBits);                %Maping the income data to the DQPSK format
                                    %Converting the I and Q components to the polirezed NRZ format
                                    DataI(DataI==1) =  1;
                                    DataI(DataI==0) = -1;
                                    DataQ(DataQ==1) =  1;
                                    DataQ(DataQ==0) = -1;
                                    
                                    TxOdd  = TxData(logical(mod(TxDataPos,2)));                    %Spliting the information of odd positions
                                    TxEven = TxData(~(mod(TxDataPos,2)));                          %Spliting the information of even positions
                                    
                                    TxDataMat(kk,:)         = TxOdd;                               %Storring the transmited information for latter evaluation
                                    TxDataMat(kk+NumCarr,:) = TxEven;                              %Storring the transmited information for latter evaluation
                                    
                                    %Once that we have the data to be transmited it needs to be
                                    %converted to the eletrical signal that hereafter will modulate
                                    %our carrier through the MZM-I
                                    SigI = rectpulse(DataI,NPPB);
                                    SigQ = rectpulse(DataQ,NPPB);
                                    
                                    % Adding CP to the data
                                    if AddCP
                                        TxAux1 = reshape(SigI,NPPB,NbDQPSK/2);
                                        TxAux1 = [TxAux1(1:NumAmosCP,:);TxAux1;TxAux1(end-(...
                                            NumAmosCP-1):end,:)];
                                        TxAux1 = reshape(TxAux1,1,(2*NumAmosCP+NPPB)*NbDQPSK/2);
                                        SigI = [TxAux1 TxAux1(end-(StuffSampels-1):end)];
                                        
                                        TxAux2 = reshape(SigQ,NPPB,NbDQPSK/2);
                                        TxAux2 = [TxAux2(1:NumAmosCP,:);TxAux2;TxAux2(end-(...
                                            NumAmosCP-1):end,:)];
                                        TxAux2 = reshape(TxAux2,1,(2*NumAmosCP+NPPB)*NbDQPSK/2);
                                        SigQ = [TxAux2 TxAux2(end-(StuffSampels-1):end)];
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
                                    TxSigI = SigI;                             %Conforming the information and Creating the modulation signal
                                    TxSigQ = SigQ;                             %Conforming the information and Creating the modulation signal
                                    %             PrintInfo(Ploting*3,TxSigI,T,NPPB,TxSigQ);
                                    %                 TxSigMat(kk,:)         = TxSigI;                               %Storing the eletrical signal for further evaluation
                                    %                 TxSigMat(kk+NumCarr,:) = TxSigQ;                               %Storing the eletrical signal for further evaluation
                                    EoutT = EoutAux1(VetThisCarr==ObsCarrUsed(kk),:);
                                    [EoutModAux] = IqMod(EoutT,TxSigI,TxSigQ,Vpi,V0);
                                    %                 a=4;
                                    EoutModTem(ObsCarrPos(kk),1:length(EoutModAux)) = EoutModAux;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                                    %% Taking the sampling the EVM meassurement
                                    PhaDel = 1*pi/4;
                                    TimDel = T;
                                    [EoutA,EoutB] = DelayInterf(t,TimDel,PhaDel,EoutModAux);%Coverting phase shift to amplitude variation
                                    PhaDel = -1*pi/4;
                                    TimDel = T;
                                    [EoutC,EoutD] = DelayInterf(t,TimDel,PhaDel,EoutModAux);%Coverting phase shift to amplitude variation
                                    EoutA = EoutA.*conj(EoutA);
                                    EoutB = EoutB.*conj(EoutB);
                                    EoutC = EoutC.*conj(EoutC);
                                    EoutD = EoutD.*conj(EoutD);                   %Shifting the filter for matching the received signal
                                    EoutA = ifft(fft(EoutA).*EvmFilt);
                                    EoutB = ifft(fft(EoutB).*EvmFilt);
                                    EoutC = ifft(fft(EoutC).*EvmFilt);
                                    EoutD = ifft(fft(EoutD).*EvmFilt);
                                    EoutI = (EoutB - EoutA);
                                    EoutQ = (EoutD - EoutC);
                                    EoutI = EoutI-mean(EoutI);
                                    EoutQ = EoutQ-mean(EoutQ);
                                    EoutI = EoutI./max(abs(EoutI));                                %Normalizing the signal
                                    EoutQ = EoutQ./max(abs(EoutQ));                                %Normalizing the signal
                                    PosAuxEout      = (2*NumAmosCP+NPPB)/2:(2*NumAmosCP+...
                                        NPPB):(length(EoutModAux)-StuffSampels);%Varriable respossible to take just the samples at the middle of the symbol
                                    IxAux = EoutI(PosAuxEout) + 1j.*EoutQ(PosAuxEout);
                                    TxSymbAmos(CurrentTest,ObsCarrPos==kk,:) = IxAux;%RxSymbAmosUp = [];
                                    EvmMatRef(ObsCarrPos==kk,:) = IxAux;       %Taking just the middle samples as references
                                    %                 if (kk==2)&&(Ploting)
                                    %                     figure;hold all;grid on;
                                    %                     plot(real(IxAux),imag(IxAux),'rx');
                                    %                     set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                    %                     a=0;
                                    %                 end
                                else
                                    EoutT = EoutAux1(VetThisCarr==ObsCarrPos(kk),:);
                                    EoutModAux = EoutT;
                                    %                 EoutMod = EoutMod + EoutModAux;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                                    EoutModTem(ObsCarrPos(kk),1:length(EoutModAux)) = EoutModAux;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                                end
                            end
                            if IfftOrSum
                                %             Ekaux=sum(EoutModTem);
                                %             figure;hold all;
                                %             plot(f,20*log10(abs(fftshift(fft(Ekaux)./length(Ekaux)))));
                                %             set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                %             if mod(size(EoutModTem,1),2)
                                %                 EoutModTem(end+1,:) = 0;
                                %             end
                                [EoutMod,~,~] = OpticalIFFTN(t,T,MaxNumStagT,EoutModTem);
                                %             plot(f,20*log10(abs(fftshift(fft(EoutMod)./length(EoutMod)))));
                                %             axis([0.9*fin 1.2*NumCarr*fc+fin -280 -40]);
                                %             a=1;
                            else
                                EoutMod = sum(EoutModTem);
                            end
                        case '4PAM'
                            %%        Generate the data 4PAM
                            EoutMod    = 0;
                            EoutModTem = 0;
                            for kk=InitCarrUp:CarrPass:NumCarr%Generating different data for each carrier
                                % Frist it is chosen to transmite a random information...
                                if CarrUsedUp(kk)
                                    TxData                     = (randi(2,1,Nb4Pam)-1);        %Creating the data stream to be transmited
                                    TxData(1:JusPos)           = JusVal;                       %Adding the Justification simble at the begining of the stream to sincronize received data frame
                                    TxData(end-(JusPos-1):end) = JusValEnd;                    %Adding the Justification simble at the end of the stream to sincronize received data frame
                                    %             else
                                    %                 TxData                     = zeros(1,Nb4Pam);              %Creating the data stream to be transmited
                                    %             end
                                    TxDataMat(kk,:) = TxData;                                      %Storring the transmited information for latter evaluation
                                    %Once that we have the data to be transmited it needs to be
                                    %converted to the eletrical signal that hereafter will modulate
                                    %our carrier through the MZM-I
                                    
                                    %Selecting if the PAM will be done on Electrical or Optical
                                    %domain
                                    if Which4PAM
                                        %  Good for others kms
                                        [Phi1,Phi2] = Maping4PamIq(TxData,Vmin,Vmax,ModSchem...
                                            ,FiberLength,SetCpSampZer);%Generating the eletrical signal for the optical PAM4
                                        
                                        %  Good for 100 kms
                                        %                 [Phi1,Phi2] = Maping4PamIq2(TxData,0,U_pi2/2);
                                        %                 Vbias=U_pi2/2;
                                    else
                                        [Phi1,Phi2] = Maping4Pam(TxData,VPI,Polirized,MaxAmp4PAM); %Generating the eletrical signal for the electrical PAM4
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
                                    if AddCP
                                        TxAux1 = reshape(TxSig1,NPPB,Nb4Pam/2);
                                        if SetCpSampZer
                                            TxAux1 = [zeros(NumAmosCP,size(TxAux1,2));TxAux1...
                                                ;zeros(NumAmosCP,size(TxAux1,2))];
                                        else
                                            TxAux1 = [TxAux1(1:NumAmosCP,:);TxAux1;TxAux1(...
                                                end-(NumAmosCP-1):end,:)];
                                        end
                                        TxAux1 = reshape(TxAux1,1,(2*NumAmosCP+NPPB)*Nb4Pam/2);
                                        TxSig1 = [TxAux1 TxAux1(end-(StuffSampels-1):end)];
                                        TxAux2 = reshape(TxSig2,NPPB,Nb4Pam/2);
                                        if SetCpSampZer
                                            TxAux2 = [zeros(NumAmosCP,size(TxAux2,2));TxAux2...
                                                ;zeros(NumAmosCP,size(TxAux2,2))];
                                        else
                                            TxAux2 = [TxAux2(1:NumAmosCP,:);TxAux2;TxAux2(...
                                                end-(NumAmosCP-1):end,:)];
                                        end
                                        TxAux2 = reshape(TxAux2,1,(2*NumAmosCP+NPPB)*Nb4Pam/2);
                                        TxSig2 = [TxAux2 TxAux2(end-(StuffSampels-1):end)];
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
                                    
                                    EoutT = EoutAux1(VetThisCarr==ObsCarrUsed(kk),:);
                                    U.U1t  = TxSig1;                                               %Assigning the electrical signal to one drive of the MZM
                                    U.U2t  = TxSig2;                                               %Assigning the electrical signal to another drive of the MZM
                                    if Which4PAM
                                        if ModSchem
                                            %                         [EoutModAux] = IqMod4Pam (EoutT,U.U1t,U.U2t,U_pi2,Vbias);
                                            [EoutModAux] = Mach_Zehnder_Modulator_simplificado_2(freqGHz,EoutT,U,MZ_Input_File);
                                        else
                                            [EoutModAux] = Mach_Zehnder_Modulator_simplificado_2(freqGHz,EoutT,U,MZ_Input_File);
                                        end
                                    else
                                        [EoutModAux] = Mach_Zehnder_Modulator_simplificado_2(freqGHz...
                                            ,EoutT,U,MZ_Input_File);
                                    end
                                    %                 EoutMod = EoutMod + EoutModAux;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                                    EoutModTem(ObsCarrPos(kk),1:length(EoutModAux)) = EoutModAux;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                                    %% Taking the sampling the EVM meassurement
                                    PosAuxEout = (2*NumAmosCP+NPPB)/2:(2*NumAmosCP+NPPB)...
                                        :(length(EoutModAux)-StuffSampels);%Varriable respossible to take just the samples at the middle of the symbol
                                    Ix         = EoutModAux.*conj(EoutModAux);             %Recovering the signal that will be transmited
                                    Ix         = ifft(fft(Ix).*EvmFilt);                   %Removing higher signal generated by the receiving process
                                    Ix         = Ix - min(Ix);
                                    Ix         = Ix./max(abs(Ix));                         %Normalizing the reference
                                    IxAux      = Ix(PosAuxEout);
                                    a=0;
                                    TxSymbAmos(CurrentTest,ObsCarrPos==kk,:) = IxAux;%RxSymbAmosUp = [];
                                    EvmMatRef(ObsCarrPos==kk,:) = IxAux;                   %Taking just the middle samples as references
                                    %                 if (kk==2)&&(Ploting)
                                    %                     figure;hold all;grid on;
                                    %                     plot(t,Ix);plot(t(PosAuxEout),IxAux,'x');
                                    %                     set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                    %                     figure;hold all;grid on;
                                    %                     plot(IxAux,zeros(1,length(IxAux)),'rx');
                                    %                     set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                    %                     a=0;
                                    %                 end
                                else
                                    EoutT = EoutAux1(VetThisCarr==ObsCarrPos(kk),:);
                                    EoutModAux = EoutT;
                                    %                 EoutMod = EoutMod + EoutModAux;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                                    EoutModTem(ObsCarrPos(kk),1:length(EoutModAux)) = EoutModAux;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                                end
                                %             if (kk==2)&&(Ploting)
                                %                 PrintInfo(Ploting*4,abs(EoutModAux(1,1:end-StuffSampels)).^2./max(abs(EoutModAux(1,1:end-StuffSampels)).^2),2*t(2*NumAmosCP+NPPB),2*NumAmosCP+NPPB,abs(EoutModAux(1:end-StuffSampels)).^2./max(abs(EoutModAux(1,1:end-StuffSampels)).^2));
                                %                 a=4;
                                %             end
                            end
                            if IfftOrSum
                                %             Ekaux=sum(EoutModTem);
                                %             figure;hold all;
                                %             plot(f,20*log10(abs(fftshift(fft(Ekaux)./length(Ekaux)))));
                                %             set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                % %             if ~mod(RefCarr,2)
                                % %                 [EoutMod,~,~] = OpticalIFFTNP(t,T,MaxNumStagT,EoutModTem);
                                % %             else
                                if mod(size(EoutModTem,1),2)
                                    EoutModTem(end+1,:) = 0;
                                end
                                [EoutMod,~,~] = OpticalIFFTN(t,T,MaxNumStagT,EoutModTem);
                                % %             end
                                %             plot(f,20*log10(abs(fftshift(fft(EoutMod)./length(EoutMod)))));
                                %             axis([0.9*fin 1.2*NumCarr*fc+fin -280 -40]);
                                %             a=1;
                            else
                                EoutMod = sum(EoutModTem);
                            end
                        otherwise
                            %%        Generate the data OOK
                            EoutMod    = 0;
                            EoutModTem = 0;
                            for kk=InitCarrUp:CarrPass:NumCarr%Generating different data for each carrier
                                % Frist it is chosen to transmite just one high pulse for
                                %testing the channel...
                                if CarrUsedUp(kk)
                                    TxData                     = (randi(2,1,Nb - NumBitDesc)-1);   %Creating Random Information that will be loaded in each individual subcarrier
                                    TxData(1:JusLen)           = JusVal;                           %Making the First 4 bits equal to zero
                                    TxData(end-(JusLen-1):end) = JusVal;                           %Making the Last 4 bits equal to zero
                                    
                                    
                                    TxDataMat(kk,:)            = TxData;                           %Storring the transmited information for latter evaluation
                                    if NRZPolarity
                                        TxData(TxData==0)      = NrzMin;
                                        TxData(TxData==1)      = NrzMax;
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
                                    if AddCP
                                        TxAux = reshape(TxDataRes,NPPB,Nb-NumBitDesc);
                                        TxAux = [TxAux(1:NumAmosCP,:);TxAux;TxAux(end-(NumAmosCP...
                                            -1):end,:)];
                                        TxAux = reshape(TxAux,1,(2*NumAmosCP+NPPB)*(Nb-...
                                            NumBitDesc));
                                        TxDataRes = [TxAux TxAux(end-(StuffSampels-1):end)];
                                    end
                                    
                                    
                                    %Because of reasons, the real world does not have components
                                    %capable instantly change the voltage level. Thus, to emulate
                                    %this behaviour the eletrical PAM signal pass through a
                                    %gaussian filter that will remove the frequencies of higher
                                    %order, which will result in small slope in the level variation
                                    
                                    [BitFilt,~] = FiltroGaussiano(f,BWD,CenFeq,FiltOrd);           %Creating filter for conformation of the input information
                                    BitFilt = fftshift(BitFilt);                                   %Doing a shift on the vector for matching the transmited data
                                    %                 TxSig = ifft(fft(TxDataRes).*BitFilt);                         %Conforming the information and Creating the modulation signal
                                    TxSig = TxDataRes;                         %Conforming the information and Creating the modulation signal
                                    %             TxSigMat(kk,:) = TxSig;                                        %Storring the transmited information for latter evaluation
                                    
                                    EoutT = EoutAux1(VetThisCarr==ObsCarrUsed(kk),:);
                                    %Assigning the eletrical signal to one drive of the MZM -
                                    %The aplitude of the signal can be controlled by the
                                    %variable DatGai, which can be understood as an gain for
                                    %the eletrical signal or an atenuation. The second signal
                                    %will be similar with the only difference a phase shift of
                                    %pi.
                                    U.U1t  = DatGai.*TxSig;
                                    U.U2t  = DatGai.*exp(-1j*pi).*TxSig;
                                    %As both signals will have the mostrly the same
                                    %characteristics with the only difference the phase shift
                                    %of 180 degress. The MZM-I will be working on the Push-Pull
                                    %configuration. It is necessary to reduce the Chirp noise
                                    %to zero.
                                    [EoutModAux,~]=Mach_Zehnder_Modulator_simplificado_2(freqGHz,...
                                        EoutT,U,MZ_Input_File);    %Modulating individual carriers
                                    %                 EoutMod = EoutMod + EoutModAux;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                                    EoutModTem(ObsCarrPos(kk),1:length(EoutModAux)) = EoutModAux;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                                    %% Taking the sampling the EVM meassurement
                                    PosAuxEout = (2*NumAmosCP+NPPB)/2:(2*NumAmosCP+NPPB)...
                                        :(length(EoutModAux)-StuffSampels);%Varriable respossible to take just the samples at the middle of the symbol
                                    Ix         = EoutModAux.*conj(EoutModAux);             %Recovering the signal that will be transmited
                                    Ix         = ifft(fft(Ix).*EvmFilt);                   %Removing higher signal generated by the receiving process
                                    Ix         = Ix - min(Ix);
                                    Ix         = Ix./max(abs(Ix));                         %Normalizing the reference
                                    IxAux      = Ix(PosAuxEout);                           %Normalizing the reference
                                    a=0;
                                    TxSymbAmos(CurrentTest,ObsCarrPos==kk,:) = IxAux;%RxSymbAmosUp = [];
                                    EvmMatRef(ObsCarrPos==kk,:) = IxAux;                   %Taking just the middle samples as references
                                    %                 if (kk==2)&&(Ploting)
                                    %                     figure;hold all;grid on;
                                    %                     plot(IxAux,zeros(1,length(IxAux)),'rx');
                                    %                     set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                    %                     a=0;
                                    %                 end
                                else
                                    EoutT = EoutAux1(VetThisCarr==ObsCarrPos(kk),:);
                                    EoutModAux = EoutT;
                                    %                 EoutMod = EoutMod + EoutModAux;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                                    EoutModTem(ObsCarrPos(kk),1:length(EoutModAux)) = EoutModAux;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                                end
                            end
                            if IfftOrSum
                                %             Ekaux=sum(EoutModTem);
                                %             figure;hold all;
                                %             plot(f,20*log10(abs(fftshift(fft(Ekaux)./length(Ekaux)))));
                                %             set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                if mod(size(EoutModTem,1),2)
                                    EoutModTem(end+1,:) = 0;
                                end
                                [EoutMod,~,~] = OpticalIFFTN(t,T,MaxNumStagT,EoutModTem);
                                %             plot(f,20*log10(abs(fftshift(fft(EoutMod)./length(EoutMod)))));
                                %             axis([0.9*fin 1.2*NumCarr*fc+fin -280 -40]);
                                %             a=1;
                            else
                                EoutMod = sum(EoutModTem);
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
                    
                    
                    if AddingNoiseF
                        EsN0 = 10*log10(snr);
                        [~,sigpower] = MeasPower(EoutRec);
                        EoutRec = awgn(EoutRec,osnrf,sigpower-30);
                        [~,sigpower] = MeasPower(EoutRec);
                    end
                    
                    clear TxSigMat TxSig TxSigI TxSigQ EoutAux1 VetThisCarr EoutMod Ix Esync...
                        Esync1 Esync2 EoutA EoutB EoutC EoutD EoutModTem
                    % clearvars -except fc Rb NPPB fsample T Ta NumberOf_T FinalTime Nb NbDPSK...
                    % TotalSamples t f Eout EoutMod EoutRec TxDataMat NumCarr NbDQPSK ValsLev...
                    % Modulation Polirized MaxAmp4PAM NumbOfTest UsedModula TestBundle Nb4Pam ...
                    % ThisModula CurrentTest CurrentModula ThisTest ThisModula OfcName BerOOK ...
                    % Ber4PAM BerDQPSK BerDPSK FiberDelay Medium CurrentMedium OSNRPC AddCP V0...
                    % OcsToTest CurrentOCS PulseResp NumAmosCP NumBitDesc StuffSampels JusPos ...
                    % contcp Ber4PAM PAM4Type Vpi MZ_Input_File ModAll JusValEnd NRZPolarity ...
                    % VetElecPower AberLev fin FFTSplit PastTime PowRef SnrRef CarSNR ValsLev2...
                    % AberLevI ValsLevI AberLevQ ValsLevQ VetElecPowerQ VetElecPowerI SigAten ...
                    %                                      lambdac Atenuation VetOptiPower JusVal
                    
                    %%                Datum Upstream Reception
                    %This section is responsible to receive the information from
                    %users
                    %% DatumUpstrRec;
                    
                    %%                Main of All optical FFT
                    % The principal importance of the script is to load the All Optical FFT.
                    % Also, it is responsible to receive the transmited data.
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    %%            Recovering the signal
                    % At this point the FFT will be inplemented. The received signal need to be
                    % given as an parameter. This following function recursively implement the
                    % FFT Operation. One may need to investigate how would work a receptor with
                    % filters insted of a OFFT. The second implementation is an array of
                    % filters.
                    if FFTSplit
                        [EoutAux1,~,VetThisCarr]=OpticalFFTN(t,T,MaxNumStag,EoutRec);
                    else
                        VetThisCarr = ObsCarrPos;
                        EoutAux1    = SelectEachCarrier(EoutRec,NumCarr,f,fin,1.0*fc,5,fc);
                    end
                    for kk=1:size(EoutAux1,1)
                        [~,CarrRecPowUp(CurrentTest,kk)] = MeasPower(EoutAux1(VetThisCarr==ObsCarrUsed(kk),:),t);
                    end
                    %%  Ploting some results for qualitative analizes
                    % PrintInfo(Ploting*12,f,db(abs((fftshift(fft(EoutRec))))/length(EoutRec)),EoutAux1,RefCarr*fc);
                    %%      Recovering the Data
                    %This is basicaly the final step, so far, in this process. Here, the
                    %transmited signal will be received and processed to recover the actual
                    %data within it.
                    switch Modulation
                        case 'OFDM'
                            for ThisCarr=InitCarr:CarrPass:NumCarr
                                EoutAux = EoutAux1(VetThisCarr==ObsCarrUsed(ThisCarr),:);
                                switch Medium
                                    case 'Fiber'
                                        EoutAux = ifft(fft(EoutAux).*exp(1j*2*pi*f*(...
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
                                Ix =EoutAux.*conj(EoutAux);
                                [BitFilt,~] = FiltroGaussiano(f,BWD,CenFeq,FiltOrd);           %Creating filter for selection of the received signal
                                BitFilt = fftshift(BitFilt);                                   %Shifting the filter for matching the received signal
                                Ix = ifft(fft(Ix).*BitFilt);                                   %Filtering the signal for better detection
                                Ix = Ix - min(Ix);                                             %Removing any off-set that may exist
                                Ix = Ix./max(abs(Ix));                                         %Normalizing the eletrical signal (amplifying what is needed)
                                if ReceptorNoise                                               %Verify whether the noise is to be added or not
                                    [~,SigPower] = MeasPower(Ix);
                                    SigPower2 = SigPower-30;%20*log10(SigPowerI);
                                    Ix = awgn(Ix,SNR,SigPower2);
                                end
                                if Ploting
                                    figure;hold all;grid on;plot(f,20*log10(abs(fftshift(fft(Ix)./length(Ix)))));axis([-25e9 37.5e9 -200 0]);set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                end
                                Ix = Ix(1+NPSTUf:end);
                                Ix = reshape(Ix,length(Ix)/NumFra,NumFra);
                                Ix = intdump(Ix,NPPOF);                            %Downsampling the income signal
                                SigRecepA = Ix(1+NPOFEX:end,:);
                                
                                %CpPaFr = 1;
                                %for CarrOffSet=-1*(NumFraPar-1)/2:1:(NumFraPar-1)/2
                                switch SelModTp
                                    case 'BP'
                                        %                     SigRecep  = intdump(SigRecepA,OvSam);
                                        if SelecGaus
                                            [ BitFiltEle ] = FiltroGaussiano(ff,OBw/2,0,1);%CarrOffSet*FramSpac,1);
                                        else
                                            [ BitFiltEle ] = Filtro_Retangular(OBw,0,ff);%CarrOffSet*FramSpac,ff);
                                        end
                                        BitFiltEle = fftshift(BitFiltEle);
                                        for jj=1:NumFra
                                            ModFilta(:,jj) = BitFiltEle;
                                        end
                                        SigRecepB   = SigRecepA;%.*exp(-1j*2*pi*CarrOffSet*FramSpac*tta);
                                        if Ploting
                                            figure;
                                            hold all;
                                            plot(ff,20*log10(abs(fftshift(fft(SigRecepB)./length(SigRecepB)))));
                                        end
                                        SigRecepB = ifft(fft(SigRecepB).*ModFilta);
                                        SigRecep  = intdump(SigRecepB,OvSam);
                                        if Ploting
                                            plot(ff,20*log10(abs(fftshift(fft(SigRecepB)./length(SigRecepB)))));
                                            set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                        end
                                    case 'AM'
                                        %                     tt = linspace(0,1*Te,size(SigRecepA,1));
                                        %                     ff = time2freq(tt);
                                        %                     for jj=1:NumFra
                                        %                         tta(:,jj) = tt;
                                        %                     end
                                        [BitFiltEle,~] = FiltroGaussiano(ff,OBw/2,0,1);           %Creating filter for selection of the received signal
                                        BitFiltEle = fftshift(BitFiltEle);
                                        for jj=1:NumFra
                                            ModFilta(:,jj) = BitFiltEle;
                                        end
                                        SigRecepB   = SigRecepA.*cos(2*pi*Ofc*tta);
                                        SigRecepB = ifft(fft(SigRecepB).*ModFilta);
                                        SigRecep  = intdump(SigRecepB,OvSam);
                                        %                         x_ofdm = x_ofdm_ssb;
                                    case 'AMSSB'
                                        %                     tt = linspace(0,1*Te,size(SigRecepA,1));
                                        %                     ff = time2freq(tt);
                                        %                     for jj=1:NumFra
                                        %                         tta(:,jj) = tt;
                                        %                     end
                                        [BitFiltEle,~] = FiltroGaussiano(ff,OBw/2,0,1);           %Creating filter for selection of the received signal
                                        BitFiltEle = fftshift(BitFiltEle);
                                        for jj=1:NumFra
                                            ModFilta(:,jj) = BitFiltEle;
                                        end
                                        SigRecepB   = SigRecepA.*exp(-1j*2*pi*Ofc*tta);
                                        SigRecepB = ifft(fft(SigRecepB).*ModFilta);
                                        SigRecep  = intdump(SigRecepB,OvSam);              %Over sampling
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
                                if Ploting
                                    switch OfdMod                                                  %Modulation Tx signal for future comparison
                                        case 'qam'
                                            TxDataA = TxDataMat(ThisCarr,:);
                                            TxDataA = reshape(TxDataA,NumFra,length(TxDataA)/NumFra);
                                            TxDataA = TxDataA.';
                                            for CarDmtM=1:length(DmtMve)
                                                M = DmtMve(CarDmtM);
                                                TxSigToPlot(CarDmtM,:)  = qammod(TxDataA(CarDmtM,:),M);            %Modulating information
                                            end
                                        otherwise
                                            TxDataA = TxDataMat(ThisCarr,:);
                                            TxDataA = reshape(TxDataA,NumFra,length(TxDataA)/NumFra);
                                            TxDataA = TxDataA.';
                                            %                     SigRecepA = reshape(TxDataMat(ThisCarr,:),NumFra,NFFT);
                                            TxSigToPlot = dpskmod(TxDataA,M);                  %Modulating information
                                            %                     TxSigToPlot = TxSigToPlot(:);
                                    end
                                    %         TxSigToPlot = dpskmod(TxDataMat(ThisCarr,2:end),M);
                                    if Ploting
                                        figure;
                                        txcolor = [0.2 0 1];
                                        rxcolor = [1 0.4 0];
                                        hold all;
                                        plot(TxSigToPlot(:),'*','LineWidth',2,'color',txcolor);
                                        plot(SigRecep3(:),'o','LineWidth',2,'color',rxcolor);
                                        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                    end
                                end
                                
                                switch OfdMod
                                    case 'qam'
                                        RxSigOfdmNoEq(CurrentTest,ThisCarr,:) = SigRecep3(:).';
                                        equa = ChanelEqualizer(ObsCarrUsed(ThisCarr),:);
                                        equ = reshape(equa,length(equa)/NumFra,NumFra);
                                        SigRecep3 = ((SigRecep3)./equ);%Fixing phase rotation and amplitude variation with the equalizer
                                        clear SigRecep4;
                                        for CarDmtM=1:length(DmtMve)
                                            M = DmtMve(CarDmtM);
                                            SigRecep4(CarDmtM,:)  = qamdemod(SigRecep3(CarDmtM,:),M);            %Modulating information
                                        end
                                        %                     SigRecep4  = qamdemod(SigRecep3,M);                    %Modulating information
                                    otherwise
                                        %                     SigRecep3 = reshape(SigRecep3,NumFra,NFFT);            %Reshaping the signal if multiple carriers were transmited
                                        SigRecep4 = dpskdemod(SigRecep3,M);                    %Modulating information
                                        %                     SigRecep4 = SigRecep4(:);                              %Serializing the signal
                                end
                                EvmMatRec(ObsCarrPos==ThisCarr,:) = ...
                                    SigRecep3(:);%Taking just the middle samples as references
                                [EvmDB(CurrentTest,ThisCarr), ...
                                    EvmPer(CurrentTest,ThisCarr), ...
                                    EvmRms(CurrentTest,ThisCarr) ] = ...
                                    EvmCalc( EvmMatRef(ObsCarrPos==ThisCarr,:),SigRecep3(:).' );
                                [EvmDBJ(CurrentTest,ThisCarr),...
                                    EvmPerJ(CurrentTest,ThisCarr),...
                                    EvmRmsJ(CurrentTest,ThisCarr)] = ...
                                    evm1(M,OfdMod,EvmMatRef(ObsCarrPos==ThisCarr,:),SigRecep3(:).');
                                %         SigRecep4 = dpskdemod(SigRecep3,M);
                                SigRecep4 = SigRecep4.';
                                SigRecep4 = SigRecep4(:).';
                                RxSigOfdm(CurrentTest,ThisCarr,:) = SigRecep4;
                                if Ploting
                                    figure(222);
                                    txcolor = [0.2 0 1];
                                    rxcolor = [1 0.4 0];
                                    hold all;
                                    plot(TxSigToPlot(:),'*','LineWidth',2,'color',txcolor);
                                    plot(SigRecep3(:),'o','LineWidth',2,'color',rxcolor);
                                    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                end
                                if Ploting
                                    figure;
                                    hold all;
                                    plot(TxDataMat(ThisCarr,:));
                                    plot(SigRecep4);
                                    SampPos = 1:length(SigRecep4);
                                    plot(SampPos(~(TxDataMat(ThisCarr,:)==SigRecep4)),...
                                        SigRecep4(~(TxDataMat(ThisCarr,:)==SigRecep4)),'o');
                                    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                end
                                SymbErr = ~(TxDataMat(ThisCarr,:)==SigRecep4);%Measuring system bit error ration
                                %             CpPaFr = CpPaFr + 1;
                                %             end
                                DmtMvep = repmat(DmtMve,NumFra,1);
                                DmtKvep = log2(DmtMvep(:).');
                                
                                BerToPlotOfdm(CurrentTest,ThisCarr,:) = SymbErr.*DmtKvep;
                                %             BerOFDM(CurrentTest,ThisCarr) = sum(SymbErr)/length(SymbErr);
                                BerOFDM(CurrentTest,ThisCarr) = sum(SymbErr.*DmtKvep)/sum(DmtKvep);
                                a=0;
                                close all;
                            end
                        case 'DPSK'
                            %%                   Receiver DPSK
                            for ThisCarr=InitCarr:CarrPass:NumCarr                                    %For each carrier the same process of reception will be used.
                                if CarrUsedUp(ThisCarr)
                                    %%  Reception Process: Handling the income optical field
                                    %At the first step the income field will be selected from the
                                    %optical FFT Output with the help of the maping vector
                                    %(previously described).
                                    EoutAux = EoutAux1(VetThisCarr==ObsCarrUsed(ThisCarr),:);
                                    EoutAux = EoutAux1(VetThisCarr==ObsCarrUsed(ThisCarr),:);
                                    if AddingNoiseP
                                        [~,sigpower] = MeasPower(EoutAux);
                                        EoutAux      = awgn(EoutAux,osnrp,sigpower-30);
                                        [~,sigpower] = MeasPower(EoutAux);
                                    end
                                    %%            Fiber Time Delay Compensation
                                    switch Medium
                                        case 'Fiber'
                                            EoutAux = ifft(fft(EoutAux).*exp(1j*2*pi*f*(...
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
                                    [ESync1,ESync2] = DelayInterf(t,TimDel,PhaDel,EoutAux);        %Coverting phase shift to amplitude variation
                                    
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
                                    
                                    if RecFilBanPas
                                        Esync  = ifft(fft(Esync).*BitFilt);                             %Filter is used to remove higher order components
                                        Esync  = Esync-mean(Esync);
                                        Esync  = Esync./max(abs(Esync));
                                    end
                                    %                 ESync1 = ESync1./max(abs(ESync1));                             %Normalizing the signal
                                    %                 ESync2 = ESync2./max(abs(ESync2));                             %Normalizing the signal
                                    %                 Esync  = Esync./max(abs(Esync));                               %Normalizing the signal
                                    
                                    %%         Adding Noise to Received Signal
                                    %Just after the photo diod is where the noise must be added as
                                    %it is the main constrain of the system. Once the signal is
                                    %converted from optical to electrical and it digital values
                                    %were recevered, there is no point to add noise, because the
                                    %information was already received, it just needs decoding.
                                    %The noise here added will be gaussian basically it represents
                                    %the reception sensibility. In another words, how much the
                                    %signal must be about the noise for an acceptable recovering.
                                    if ReceptorNoise                                               %Verify whether the noise is to be added or not
                                        [~,SigPower] = MeasPower(Esync);
                                        SigPower2 = SigPower-30;%20*log10(SigPowerI);
                                        SNR2 = CarSNR + 10*log10(1) - 10*log10(Nsamp) + 10*log10(10^0.36);
                                        if TimeSys
                                            SNR = CarSNR + 10*log10(1) - 10*log10(1*T*BWD2);
                                        else
                                            SNR = CarSNR + 10*log10(1) - 10*log10(1*t(2*NumAmosCP+NPPB)*BWD2);
                                        end
                                        Esync = awgn(Esync,SNR,SigPower2);
                                        %                     Esync = ifft(fft(Esync).*BitFiltN);                             %Filter is used to remove higher order components
                                    end
                                    if ~RecFilBanPas
                                        Esync  = ifft(fft(Esync).*BitFilt);                             %Filter is used to remove higher order components
                                        Esync  = Esync-mean(Esync);
                                        Esync  = Esync./max(abs(Esync));
                                    end
                                    
                                    
                                    VetElecPower(CurrentTest,ThisCarr)= MeasPower(Esync);
                                    EoutAuxF = 20*log10(abs(fftshift(fft(EoutAux)./length(...
                                        EoutAux))));
                                    [VetOptiPower(CurrentTest,ThisCarr),~]= findpeaks(EoutAuxF,...
                                        'SortStr','descend','NPeaks',1);
                                    
                                    SyncAux   = Esync(IniSyncPos:SyncPos);                         %Selecting just the symbol to synchronize
                                    SyncedAux = SyncSymb(IniSyncPos:SyncPos);                      %Selecting just the symbol to synchronize
                                    
                                    %%                   Plot for Qualitative analizes
                                    %                             PrintInfo(Ploting*15,t(IniSyncPos:SyncPos),Esync(IniSyncPos:...
                                    %                             SyncPos),SyncSymb(IniSyncPos:SyncPos),ESync1(IniSyncPos:...
                                    %                                                       SyncPos),ESync2(IniSyncPos:SyncPos));
                                    
                                    %%                   Synchronizing
                                    %This synchronization process is based on the mean expected
                                    %value. That means, the information mean value should be within
                                    %one period of symbol. Thus, the mean value of the received
                                    %signal is acquired and compare of the known sync-word to
                                    %verify if this mean value is at the right possition. Which is
                                    %the midel point (peak) of the highest level at the sync period
                                    SyncAux(SyncAux<0)            = 0;                           %To keep the mean value above zero anything under is neglected
                                    %                 SyncAux(SyncAux>=mean(SyncAux)) = 1;                           %Adding a flag to the first sample of the received mean value
                                    %                 SyncAux(SyncAux<mean(SyncAux))  = -1;                          %All the others samples at set to the lowest level
                                    SyncAux(SyncAux>0) = 1;                           %Adding a flag to the first sample of the received mean value
                                    if PrintinEye
                                        figure; hold all;
                                        plot(t(IniSyncPos:SyncPos),SyncSymb(IniSyncPos:SyncPos),t(IniSyncPos:SyncPos),Esync(IniSyncPos:SyncPos),t(IniSyncPos:SyncPos),SyncAux,t(IniSyncPos:SyncPos),linspace(mean(SyncAux),mean(SyncAux),length(t(IniSyncPos:SyncPos))));
                                        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                    end
                                    SyncAux(SyncAux<mean(SyncAux))  = -1;                          %All the others samples at set to the lowest level
                                    
                                    PosToSyn  = find(ismember(SyncAux,1));                         %Finding where is the location of the samples to synchronize
                                    PosSynced = find(ismember(SyncedAux,1));                       %Finding where is the location of the samples to synchronize
                                    
                                    DiffPos = PosToSyn(round(end/2)) - PosSynced(round(end/2));    %Accounting the peak (midel point) displacement
                                    sn=((1/2)*(((1/2)^MaxNumStag)-1))/((1/2)-1);
                                    DiffPos = round((sn/(3/2.55^IfftOrSum))*T/Ta);
                                    if DiffPos~=0
                                        if DiffPos>0%If the difference is positive, left-shift...
                                            %                     Esync = ifft(fft(Esync).*exp(1j*2*pi*f*(DiffPos*Ta))); %Shift based on time change
                                            Esync = [Esync(DiffPos+1:end) Esync(1:DiffPos)];       %Shift based on sampling sliding
                                        else%... but if the difference is negative, right-shift
                                            %                     Esync = ifft(fft(Esync).*exp(1j*2*pi*f*(DiffPos*Ta))); %Shift based on time change
                                            Esync = [Esync(end+DiffPos+1:end) Esync(1:end + ...
                                                DiffPos)];%Shift based on sampling sliding
                                        end
                                    end
                                    if PrintinEye
                                        figure; hold all;
                                        plot(t(IniSyncPos:SyncPos),SyncSymb(IniSyncPos:SyncPos),t(IniSyncPos:SyncPos),Esync(IniSyncPos:SyncPos));
                                        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                    end
                                    %Because of reasons, sometimes it may be required to make a
                                    %synchronization process with the  end of the data stream as
                                    %well. This following verification check if the user set (or
                                    %not) a second synchronization process to be done.
                                    
                                    %%          Ploting the result for qualitative analizes
                                    %             PrintInfo(Ploting*16,t(end-SyncPos+1:end-IniSyncPos+1),Esync...
                                    %             (end-SyncPos+1:end-IniSyncPos+1),SyncSymb(end-SyncPos+1:end-...
                                    %             IniSyncPos+1),ESync1(end-SyncPos+1:end-IniSyncPos+1),ESync2(...
                                    %                                           end-SyncPos+1:end-IniSyncPos+1));
                                    if SencondAdjust
                                        %%                   Synchronizing
                                        %This synchronization process is based on the mean expected
                                        %value. That means, the information mean value should be
                                        %within one period of symbol. Thus, the mean value of the
                                        %received signal is acquired and compare of the known
                                        %sync-word to verify if this mean value is at the right
                                        %possition.
                                        
                                        SyncAuxEnd   = Esync(end-SyncPos+1:end-IniSyncPos+1);
                                        SyncedAuxEnd = SyncSymbEnd(end-SyncPos+1:end-IniSyncPos+1);
                                        SyncAuxEnd(SyncAuxEnd<0)                 = 0;              %To keep the mean value above zero anything under is neglected
                                        SyncAuxEnd(SyncAuxEnd>=mean(SyncAuxEnd)) = 1;              %Adding a flag to the first sample of the received mean value
                                        SyncAuxEnd(SyncAuxEnd<mean(SyncAuxEnd))  = -1;             %All the others samples at set to the lowest level
                                        
                                        PosToSynEnd  = find(ismember(SyncAuxEnd,1));               %Finding where is the location of the first sample to synchronize
                                        PosSyncedEnd = find(ismember(SyncedAuxEnd,1));             %Finding where is the location of the first sample to synchronize
                                        
                                        DiffPosEnd = PosToSynEnd(round(end/2))-PosSyncedEnd(...
                                            round(end/2));
                                        if DiffPosEnd~=0
                                            if DiffPosEnd>0%If positive difference, left-shift...
                                                %                         Esync = ifft(fft(Esync).*exp(1j*2*pi*f*(...
                                                %                                                           DiffPosEnd*Ta)));%Shift based on time change
                                                Esync = [Esync(DiffPosEnd+1:end) Esync(1:...
                                                    DiffPosEnd)];%Shift based on sampling sliding
                                            else%... but if the difference is negative, right-shift
                                                %                         Esync = ifft(fft(Esync).*exp(1j*2*pi*f*(...
                                                %                                                           DiffPosEnd*Ta)));%Shift based on time change
                                                Esync = [Esync(end+DiffPosEnd+1:end) Esync(...
                                                    1:end+DiffPosEnd)];%Shift based on sampling sliding
                                            end
                                        end
                                    end
                                    
                                    %% Removing CP
                                    %                 if ThisCarr==126
                                    %                     EyeToPlot(CurrentTest,1:length(Esync)) = Esync;
                                    %                     save(['savingforplotingeye' num2str(VetSnr)],'EyeToPlot','VetSnr','IfftOrSum','UsedModula','T','NPPB');
                                    %                 end
                                    if AddCP
                                        IxAux  = Esync(1:end - StuffSampels);
                                        IxAux  = reshape(IxAux,(2*NumAmosCP+NPPB),NbDPSK);
                                        IxAux  = IxAux(1+NumAmosCP:end-NumAmosCP,:);
                                        IxAux  = reshape(IxAux,1,NPPB*NbDPSK);
                                        Esync  = IxAux;
                                    end
                                    
                                    %% Taking the sampling the EVM meassurement
                                    clear IxAux;
                                    PosAuxEout1 = NPPB/2:NPPB:length(Esync);                   %Varriable respossible to take just the samples at the middle of the symbol
                                    PosAuxEout2 = ((NPPB/2)+(NPPB/PerDiv)):NPPB:length(Esync);      %Varriable respossible to take just the samples at the middle of the symbol
                                    PosAuxEout3 = ((NPPB/2)-(NPPB/PerDiv)):NPPB:length(Esync);      %Varriable respossible to take just the samples at the middle of the symbol
                                    IxAux1      = Esync(PosAuxEout1);                          %Normalizing the reference
                                    IxAux2      = Esync(PosAuxEout2);    %Normalizing the reference
                                    IxAux3      = Esync(PosAuxEout3);    %Normalizing the reference
                                    a=0;
                                    EvmMatRecA(1,:) = IxAux1;                       %Taking just the middle samples as references
                                    EvmMatRecA(2,:) = IxAux2;                       %Taking just the middle samples as references
                                    EvmMatRecA(3,:) = IxAux3;                       %Taking just the middle samples as references
                                    EvmMatRecA(4,:) = IxAux1;                       %Taking just the middle samples as references
                                    [EvmDBA(1), EvmPerA(1), EvmRmsA(1) ] = EvmCalc( EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux1 );
                                    [EvmDBA(2), EvmPerA(2), EvmRmsA(2) ] = EvmCalc( EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux2 );
                                    [EvmDBA(3), EvmPerA(3), EvmRmsA(3) ] = EvmCalc( EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux3 );
                                    [EvmDBA(4), EvmPerA(4), EvmRmsA(4) ] = EvmCalc( EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux1 );
                                    [EvmDBJA(1),EvmPerJA(1),EvmRmsJA(1)] = evm1(2,'dpsk',EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux1);
                                    [EvmDBJA(2),EvmPerJA(2),EvmRmsJA(2)] = evm1(2,'dpsk',EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux2);
                                    [EvmDBJA(3),EvmPerJA(3),EvmRmsJA(3)] = evm1(2,'dpsk',EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux3);
                                    [EvmDBJA(4),EvmPerJA(4),EvmRmsJA(4)] = evm1(2,'dpsk',EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux1);
                                    %                 if ThisCarr==2
                                    %                     figure;hold all;grid on;
                                    %                     plot(real(IxAux1),zeros(1,length(IxAux1)),'rx');
                                    %                     plot(real(IxAux2),zeros(1,length(IxAux2)),'yx');
                                    %                     plot(real(IxAux3),zeros(1,length(IxAux3)),'gx');
                                    %                     set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                    %                     a=0;
                                    %                 end
                                    %##########################################################################
                                    PosIx = NPPB/2:NPPB:length(Esync);                                %Possition of the central samples - but the number of samples per symbol is even ????
                                    IxAux = Esync(PosIx);                                             %From the main received signal just a few samples are taken for further evaluation
                                    n     = 100;                                                   %The number of boxes to be filled up on the histogram process
                                    Perce = 0.9;
                                    %At this process each eye will be evaluated separetelly. It is
                                    %important to remember that the PAM4 has 4 levels, which means
                                    %three level of decissions that are exaclty the center of the
                                    %eye diagram.
                                    IxAuxAB = Esync(PosIx);       %Taking just those values relative to the uper eye
                                    InterAB = linspace(Perce*(min(Esync)),Perce*(max(Esync)),n);                     %Building the histogram boxes
                                    EyeAB = hist(IxAuxAB,InterAB);                                 %filling up the boxes with samples that fit on them.
                                    
                                    %The same process described for the uper level will be done at
                                    %the middle and lower eyes levels.
                                    IxAuxCD = Esync(PosIx-NPPB/PerDiv);
                                    InterCD = linspace(Perce*(min(Esync)),Perce*(max(Esync)),n);%NPPB*2^n);
                                    EyeCD = hist(IxAuxCD,InterCD);
                                    
                                    IxAuxEF = Esync(PosIx+NPPB/PerDiv);
                                    InterEF = linspace(Perce*(min(Esync)),Perce*(max(Esync)),n);%NPPB*2^n);
                                    EyeEF = hist(IxAuxEF,InterEF);
                                    
                                    %What we are looking for here are not where exist the
                                    %occurrences of values. The eye is where there is not samples,
                                    %this means, the decission level is a reagion between actual
                                    %levels thus this reagion does not contain any signal.
                                    EyeAB = ~EyeAB;                                                %Changing zeros to one - Zeros compose the eye region
                                    EyeCD = ~EyeCD;
                                    EyeEF = ~EyeEF;
                                    %Starting variables that will be used to identify the eye
                                    %diagram.
                                    CountAB   = 1;
                                    CountCD   = 1;
                                    CountEF   = 1;
                                    
                                    SeqOnesAB = 0;
                                    SeqOnesCD = 0;
                                    SeqOnesEF = 0;
                                    
                                    SeqFinAB  = 0;
                                    SeqFinCD  = 0;
                                    SeqFinEF  = 0;
                                    
                                    SeqIniAB  = 1;
                                    SeqIniCD  = 1;
                                    SeqIniEF  = 1;
                                    %The for loop will take account of every box with ones. It is
                                    %important to take note that the not operator was used in this
                                    %vector, therefore ones means zeros (the eye diagram -
                                    %possibly) and zeros means values abouve zeroa (not the eye).
                                    for kk=1:length(EyeAB)                                         %For every box
                                        if EyeAB(kk)                                               %if it contains "1"
                                            SeqOnesAB(SeqIniAB)=CountAB;                           %count this element as part of a consecutive sequency
                                            CountAB = CountAB + 1;                                 %adds one to the counter of consecutive elements "1"
                                            if kk==length(EyeAB)                                   %if the current box is the last box we got to an end
                                                SeqFinAB(SeqIniAB) = kk;                           %The final sequency element is equal to its possition (kk)
                                            end
                                        else                                                       %else if the current box contains "0"
                                            SeqFinAB(SeqIniAB) = kk-1;                             %the previous element was the last element of a consecutive sequency of ones
                                            SeqIniAB = SeqIniAB + 1;                               %adds one to the counter of consecutive number (take in account how many sequencies there are)
                                            CountAB = 1;                                           %reset the counter of consecutive elements to mark the start of a new consecutive sequence
                                        end
                                        
                                        if EyeCD(kk)
                                            SeqOnesCD(SeqIniCD)=CountCD;
                                            CountCD = CountCD + 1;
                                            if kk==length(EyeCD)
                                                SeqFinCD(SeqIniCD) = kk;
                                            end
                                        else
                                            SeqFinCD(SeqIniCD) = kk-1;
                                            SeqIniCD = SeqIniCD + 1;
                                            CountCD = 1;
                                        end
                                        
                                        if EyeEF(kk)
                                            SeqOnesEF(SeqIniEF)=CountEF;
                                            CountEF = CountEF + 1;
                                            if kk==length(EyeEF)
                                                SeqFinEF(SeqIniEF) = kk;
                                            end
                                        else
                                            SeqFinEF(SeqIniEF) = kk-1;
                                            SeqIniEF = SeqIniEF + 1;
                                            CountEF = 1;
                                        end
                                    end
                                    
                                    
                                    %If the eye is open, which means there is a clear difference
                                    %between adjacent levels, the eye diagram will be the longest
                                    %sequence of ones.
                                    [MaxValAB,LocMaxAB]=max(SeqOnesAB);
                                    if LocMaxAB<2 || MaxValAB<2                                    %if any sequency was found or there is just one sequency it is a error thus
                                        LevDec3 = 0.0;                                            %the decission level will be by default 0.67. Also, the other variables will
                                        LocMaxAB = 1;                                              %will be set with values to not cause errors in the future
                                        SeqFinAB(1)=2;
                                        MaxValAB = 0;
                                        InterAB(1)=LevDec3;
                                    else                                                           %if a sequency was found the middle of the eye will be the middle of the sequency
                                        if (SeqFinAB(LocMaxAB)-MaxValAB/2)<1                       %if for some reason the final element of the sequency minus the half of its
                                            LevDec3 = 0.0;                                         %length results in a negative value, something went very wrong, and by
                                        else                                                       %default it will be set to 0.7
                                            LevDec3 = InterAB(round(SeqFinAB(LocMaxAB)-MaxValAB/...
                                                2));%Otherwise, the decission level is the middle of the sequency
                                        end
                                    end
                                    
                                    [MaxValCD,LocMaxCD]=max(SeqOnesCD);
                                    if LocMaxCD<2 || MaxValCD<2
                                        LevDec2 = 0.0;
                                        LocMaxCD = 1;
                                        SeqFinCD(1)=2;
                                        MaxValCD = 0;
                                        InterCD(1)=LevDec2;
                                    else
                                        if (SeqFinCD(LocMaxCD)-MaxValCD/2)<1
                                            LevDec2 = 0.0;
                                        else
                                            LevDec2 = InterCD(round(SeqFinCD(LocMaxCD)-MaxValCD/2));
                                        end
                                    end
                                    
                                    [MaxValEF,LocMaxEF]=max(SeqOnesEF);
                                    if LocMaxEF<2 || MaxValEF<2
                                        LevDec1 = 0.0;
                                        LocMaxEF = 1;
                                        SeqFinEF(1)=2;
                                        MaxValEF = 0;
                                        InterEF(1)=LevDec1;
                                    else
                                        if (SeqFinEF(LocMaxEF)-MaxValEF/2)<1
                                            LevDec1 = 0.0;
                                        else
                                            LevDec1 = InterEF(round(SeqFinEF(LocMaxEF)-MaxValEF/2));
                                        end
                                        
                                    end
                                    
                                    %another way to measure the eye opening is the get all the
                                    %boxes and find all peaks on it, that will be a plato created
                                    %by the sequences of ones (which was zeros). From thos peaks,
                                    %the eye diagram will be the longer of them hence it will take
                                    %the most part of the vector that store the findpeaks result.
                                    %Thus, the middle of the eye will be basically the middle of
                                    %the peaks vector.
                                    LocAB = find(EyeAB);
                                    if isempty(LocAB)                                              %if for some reason there are no peaks, something went wrong.
                                        LocAB = 1;
                                        LevelDec3 = 0.0;%mean(Levels(3:4));                       %by default the decission level will be set to 0.65
                                    else
                                        LevelDec3 = LevDec3;
                                    end
                                    
                                    LocCD = find(EyeCD);
                                    if isempty(LocCD)
                                        LocCD = 1;
                                        LevelDec2 = 0.0;%mean(Levels(2:3));
                                    else
                                        LevelDec2 = LevDec2;
                                    end
                                    
                                    LocEF = find(EyeEF);
                                    if isempty(LocEF)
                                        LocEF = 1;
                                        LevelDec1 = 0.0;%mean(Levels(1:2));
                                    else
                                        LevelDec1 = LevDec1;
                                    end
                                    %##########################################################################
                                    if CalcS
                                        [~,~,EyeOpenI,EyeOpenHighI,EyeOpenLowI] = Olho(Esync,T,NPPB,0);
                                        AberLevS(CurrentTest,ThisCarr)  = EyeOpenI;
                                        ValsLevS(CurrentTest,ThisCarr)  = EyeOpenLowI + EyeOpenI/2;          %Comparison between the Transmited and received and counting the differences
                                    end
                                    
                                    AberLevAuxI(1)= InterAB(SeqFinAB(LocMaxAB))-InterAB(round(SeqFinAB(LocMaxAB)-MaxValAB));
                                    AberLevAuxI(3)= InterCD(SeqFinCD(LocMaxCD))-InterCD(round(SeqFinCD(LocMaxCD)-MaxValCD));
                                    AberLevAuxI(2)= InterEF(SeqFinEF(LocMaxEF))-InterEF(round(SeqFinEF(LocMaxEF)-MaxValEF));
                                    ValsLevAuxI(1)= LevDec3;
                                    ValsLevAuxI(3)= LevDec2;
                                    ValsLevAuxI(2)= LevDec1;
                                    
                                    %% Ploting the result for qualitative analizes
                                    if PrintinEye
                                        PrintInfo(PlotingThis*17,Esync,T,NPPB);
                                        hold on;
                                        plot(t((NPPB)/2),InterAB(round(SeqFinAB(LocMaxAB)-MaxValAB)),'kx');
                                        plot(t((NPPB)/2),InterAB(SeqFinAB(LocMaxAB)),'ko');
                                        plot(t((NPPB)/2),InterAB(round(SeqFinAB(LocMaxAB)-MaxValAB/2)),'kd');
                                        plot(t(((NPPB)/2)-(NPPB)/PerDiv),InterCD(round(SeqFinCD(LocMaxCD)-MaxValCD)),'gx');
                                        plot(t(((NPPB)/2)-(NPPB)/PerDiv),InterCD(SeqFinCD(LocMaxCD)),'go');
                                        plot(t(((NPPB)/2)-(NPPB)/PerDiv),InterCD(round(SeqFinCD(LocMaxCD)-MaxValCD/2)),'gd');
                                        plot(t(((NPPB)/2)+(NPPB)/PerDiv),InterEF(round(SeqFinEF(LocMaxEF)-MaxValEF)),'bx');
                                        plot(t(((NPPB)/2)+(NPPB)/PerDiv),InterEF(SeqFinEF(LocMaxEF)),'bo');
                                        plot(t(((NPPB)/2)+(NPPB)/PerDiv),InterEF(round(SeqFinEF(LocMaxEF)-MaxValEF/2)),'bd');
                                        if CalcS
                                            plot(t((NPPB)/2),EyeOpenLowI,'mx');
                                            plot(t((NPPB)/2),EyeOpenHighI,'mo');
                                            plot(t((NPPB)/2),EyeOpenLowI + EyeOpenI/2,'md');
                                        end
                                        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                        figure;
                                        hold all;
                                        plotpos = zeros(1,length(IxAux1));
                                        plot(IxAux1,plotpos,'o','color',[1 0.4 0]);
                                        plot(EvmMatRef(ObsCarrUsed(ThisCarr),:),plotpos,'b*');
                                        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                        a=0;
                                    end
                                    %                 [~,~,EyeOpen,EyeOpenHigh,EyeOpenLow] = Olho(Esync,T,NPPB,0);
                                    
                                    %                 AberLev(CurrentTest,ThisCarr)= Inter(SeqFin(LocMax))-Inter(round(SeqFin(LocMax)-MaxVal));
                                    %                 ValsLev(CurrentTest,ThisCarr)= LevDec;
                                    %                 AberLev(CurrentTest,ThisCarr)  = EyeOpen;
                                    %                 ValsLev(CurrentTest,ThisCarr)  = EyeOpenLow + EyeOpen/2;
                                    
                                    %                 Txaux1 = rectpulse(TxDataMat(ThisCarr,:),NPPB);
                                    %                 Txaux1(Txaux1==0) = -1;
                                    %% Ploting the result for qualitative analizes
                                    %                 hold on;
                                    %                 plot(t((NPPB)/2),Inter(round(SeqFin(LocMax)-MaxVal)),'bx');
                                    %                 plot(t((NPPB)/2),Inter(SeqFin(LocMax)),'bo');
                                    %                 plot(t((NPPB)/2),Inter(round(SeqFin(LocMax)-MaxVal/2)),'bd');
                                    %                 plot(t((NPPB)/2),EyeOpenLow,'cx');
                                    %                 plot(t((NPPB)/2),EyeOpenHigh,'co');
                                    %                 plot(t((NPPB)/2),EyeOpenLow + EyeOpen/2,'cd');
                                    %                 a=0;
                                    %                 PrintInfo(Ploting*18,t(1:length(Esync)),Txaux1,Esync);
                                    %%               Recovering the Information
                                    %Once the Income signal was synchronized it is possible to
                                    %recover the signal.
                                    %
                                    %After passing the optical signal to the eletrical domain, for
                                    %actually detect the data withing the signal the following
                                    %steps are needed.
                                    %
                                    %Finding Decission Levels:
                                    %The process for decoding the income signal will be based on
                                    %eletronic comparators. Inasmuch as the right decission level
                                    %must be acquired for accurately decide, within a symbol
                                    %periode, what that current leavel means (ones or zeros).
                                    %
                                    %The process hereafter of chosing the  decission levels is not
                                    %deterministic rather it is a statistic process. The main idea
                                    %is to take the decission level from the histogram generated
                                    %from the income signal stream.
                                    %
                                    %This process is realized inside the function Olho.
                                    %
                                    %Basicaly the decission level will be the minimal value of the
                                    %currente eye under evaluation plus the half of the its eye
                                    %opening.The following ilustration better describe this process
                                    %
                                    %Eye Limits:   + (1/2 * Eye Opening:)  =    Comparisson Limit:
                                    %
                                    %UperLevel  ______     ______________     _____
                                    %                 \   /      |       \   /
                                    %                  \ /       |        \ /
                                    %                   \ Half Eye Opening /     Decission Level
                                    %                  / \       |        / \
                                    %LowerLevel ______/   \______|_______/   \_____
                                    %
                                    %
                                    %Once the signal was processed the next step is through a
                                    %comparator decide the actual information received.
                                    %
                                    %For accurately read the income signal is necessary more than
                                    %just have the decision level, it is also needed to know where
                                    %it occurs. For instance, if we decide to use just to take the
                                    %decision level and measure the mean value across and a portion
                                    %of the time period, which portion should we take? It would be
                                    %logical to take the central portion, but the bit may be
                                    %deformed in such way that the information gets concentrated
                                    %out of the middle part, which means the symbol is not
                                    %symmetric. The symmetry information can be acquired from the
                                    %eye diagram by measuring the longitudinal opening. The
                                    %following sketch better describes this process:
                                    %
                                    %                   Point of Symmetry
                                    %         _____     _______|_______     _____
                                    %              \   /               \   /
                                    %               \ /  Longitudinal   \ /
                                    %                \ ----------------- /
                                    %               / \    Opening      / \
                                    %         _____/   \_______________/   \_____
                                    %
                                    %With those two pieces of information, decision level and point
                                    %of symmetry, we have the X, Y coordinates for the centre of
                                    %the Eye Diagram. Therefore, as long as there is an opening on
                                    %it it will be possible to recover the transmitted information
                                    %without error... theoretically.
                                    %
                                    %As this process is also statistical, first we reshape the
                                    %income vector to analyze all periods at the same time.
                                    %                 EyeSymMat = reshape(Esync(1+SyncPeriod*NPPB:end-SyncPeriod*...
                                    %                     NPPB),NPPB,NbDPSK-2*SyncPeriod);
                                    %Then we take the values that compose the decision level
                                    %because they will mark the point of symmetry.
                                    %
                                    %Firstly it was set the interval in which the histogram will be
                                    %build. It is based on the number of samples per bit period.
                                    %                 Interval = linspace(min(Esync(1+SyncPeriod*NPPB:end-...
                                    %                     SyncPeriod*NPPB)),max(Esync(1+SyncPeriod*NPPB:end-SyncPeriod...
                                    %                     *NPPB)),2*NPPB);
                                    %Therefore, the MATLAB hist function returns the number of
                                    %occurrence of each interval.
                                    %                 EyeMax = hist(Esync,Interval);
                                    %                 EyeMaxaux = [0 EyeMax 0];                                      %Zeros are added at the EyeMax to auxiliate the finding peaks process
                                    %                 [~,EyeLoc] = findpeaks(EyeMaxaux,'MinPeakDistance',NPPB/4,...
                                    %                     'SortStr','descend','NPeaks',2);%The peaks on the Eye profile will be the levels at the Eyes limit
                                    
                                    %From the location of the max values that occured, which means
                                    %the uper and lower level of the eye diagram it needs to take
                                    %the actual value that those occurences represent that is
                                    %withing the the Interval variable.
                                    %                 ValToSeek = Interval(EyeLoc-1);
                                    %The number of ocurrences is a statical measure therefore one
                                    %does not have control which interval will have the highest
                                    %peak, thus it is important to ordenate the values to be seek
                                    %from the lower part of the eye diagram to the uper part of the
                                    %eye diagram.
                                    %                 ValToSeek = sort(ValToSeek,'ascend');
                                    %                 OccuCount = zeros(1,size(EyeSymMat,1));                        %Auxiliar Variable for accounting.
                                    %                 for kk=1:size(EyeSymMat,1)                                     %For every sample within a symbol period
                                    
                                    %                     OccuCount(kk) = OccuCount(kk)+sum((EyeSymMat(kk,:)>=min(...
                                    %                         Esync))&(EyeSymMat(kk,:)<=UpeSymPer*EyeOpenLow)); %Account all occurencies of the value 1
                                    %                     OccuCount(kk) = OccuCount(kk)+sum((EyeSymMat(kk,:)>=...
                                    %                         LowSymPer*EyeOpenHigh)&(EyeSymMat(kk,:)<=max(Esync))); %Account all occurencies of the value 2
                                    %                 end
                                    %The point of symmetry of the eye diagram will be where the
                                    %maximum number of occurrences were measured inasmuch as those
                                    %are the points where all the bits go to the center of the
                                    %symbol.  From the maximum number of occurrences, it can happen
                                    %for more than one sample within one symbol period, in the best
                                    %case, all samples would have the same accounting as it is
                                    %shown the ilustration above hence the symmetry will be at the
                                    %middle sample of this group of maximum occurrences. This value
                                    %can be found by the mean of the samples positions within the
                                    %symbol period. The problem with this approach is that the
                                    %signal must be synchronized with the maximum displacement of
                                    %a symbol period minus 25% of the eye Longitudinal opening if
                                    %the displacement is higher than that the point of symmetry
                                    %will be wrongly measured.
                                    %                 [SymLoc] = round(mean(find(ismember(OccuCount,max(OccuCount)...
                                    %                     ))));
                                    %             figure;findpeaks(OccuCount,'SortStr','descend');
                                    
                                    %% Actualy Receiving Data:
                                    Data = [];                                                  %Initialization of the vector that will store the income data
                                    DataM = [];                                                  %Initialization of the vector that will store the income data
                                    DataL = [];                                                  %Initialization of the vector that will store the income data
                                    DataU = [];                                                  %Initialization of the vector that will store the income data
                                    DataS = [];                                                  %Initialization of the vector that will store the income data
                                    for kk=1:NPPB:length(Esync)                                    %The comparison process will be made for each symbol period
                                        %                 MeanOfData = mean(EoutI((kk-1)+SymLocI(1)));
                                        MeanOfData = mean(Esync((kk-1)+NPPB/2));
                                        if MeanOfData > LevelDec3%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                                            Data = [Data 1];                                 %is 1
                                        else                                                       %If it is the lowest level the incoming data
                                            Data = [Data 0];                                 %is 0
                                        end
                                        MeanOfData = mean(Esync((kk-1)+(NPPB/2)+(+NPPB/PerDiv)));
                                        if MeanOfData > LevelDec1%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                                            DataM = [DataM 1];                                 %is 1
                                        else                                                       %If it is the lowest level the incoming data
                                            DataM = [DataM 0];                                 %is 0
                                        end
                                        MeanOfData = mean(Esync((kk-1)+(NPPB/2)-(+NPPB/PerDiv)));
                                        if MeanOfData > LevelDec2%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                                            DataL = [DataL 1];                                 %is 1
                                        else                                                       %If it is the lowest level the incoming data
                                            DataL = [DataL 0];                                 %is 0
                                        end
                                        MeanOfData = mean(Esync((kk-1)+(NPPB/2)));
                                        if MeanOfData > LevDefDpqsk%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                                            DataU = [DataU 1];                                 %is 1
                                        else                                                       %If it is the lowest level the incoming data
                                            DataU = [DataU 0];                                 %is 0
                                        end
                                        if CalcS
                                            MeanOfData = mean(Esync((kk-1)+(NPPB/2)));
                                            if MeanOfData > EyeOpenLowI + EyeOpenI/2%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                                                DataS = [DataS 1];                                 %is 1
                                            else                                                       %If it is the lowest level the incoming data
                                                DataS = [DataS 0];                                 %is 0
                                            end
                                        end
                                    end
                                    %%       Calculating the Bit Error Ratio (BER)
                                    %The final process here is to count the number of wrongdoings
                                    %of this whole process upon the transmited data for
                                    %quantitative analizes
                                    TxData  = TxDataMat(ThisCarr,1+SyncPeriod:end-SyncPeriod);
                                    
                                    if CalcS
                                        DataS  = DataS(1+SyncPeriod:end-SyncPeriod);
                                        BitErrS    = sum(xor(TxData,DataS));                       %Comparison between the Transmited and received and counting the differences
                                        BerDPSKS(CurrentTest,ThisCarr) = BitErrS/(NbDPSK-(2*SyncPeriod));                  %Calculating the ration of wrong bits and the total number of bits transmited
                                    end
                                    Data  = Data(1+SyncPeriod:end-SyncPeriod);
                                    DataM  = DataM(1+SyncPeriod:end-SyncPeriod);
                                    DataL  = DataL(1+SyncPeriod:end-SyncPeriod);
                                    DataU  = DataU(1+SyncPeriod:end-SyncPeriod);
                                    AberLevAuxI(4) = 0;
                                    ValsLevAuxI(4) = 0.01;
                                    
                                    BitErr(1)  = sum(xor(TxData,Data));                       %Comparison between the Transmited and received and counting the differences
                                    BitErr(2)  = sum(xor(TxData,DataM));                      %Comparison between the Transmited and received and counting the differences
                                    BitErr(3)  = sum(xor(TxData,DataL));                      %Comparison between the Transmited and received and counting the differences
                                    BitErr(4)  = sum(xor(TxData,DataU));                      %Comparison between the Transmited and received and counting the differences
                                    BerDPSK(CurrentTest,ThisCarr) = (min(BitErr))/(NbDPSK-(2*SyncPeriod));%Calculating the ration of wrong bits and the total number of bits transmited
                                    AberIAux = max(AberLevAuxI(BitErr<=min(BitErr)));
                                    AberLev(CurrentTest,ThisCarr)  = AberIAux;
                                    ValsLev(CurrentTest,ThisCarr)  = max(ValsLevAuxI(AberLevAuxI==AberIAux));          %Comparison between the Transmited and received and counting the differences
                                    if BerDPSK(CurrentTest,ThisCarr)==0
                                        LevDefDpqsk = max(ValsLevAuxI(AberLevAuxI==AberIAux));
                                        if LevDefDpqsk==0
                                            LevDefDpqsk = 0.01;
                                        end
                                    end
                                    [~,LocAux] = max(AberLevAuxI==AberIAux);
                                    EvmDB(CurrentTest,ThisCarr)   = EvmDBA(LocAux);
                                    EvmPer(CurrentTest,ThisCarr)  = EvmPerA(LocAux);
                                    EvmRms(CurrentTest,ThisCarr)  = EvmRmsA(LocAux);
                                    EvmDBJ(CurrentTest,ThisCarr)  = EvmDBJA(LocAux);
                                    EvmPerJ(CurrentTest,ThisCarr) = EvmPerJA(LocAux);
                                    EvmRmsJ(CurrentTest,ThisCarr) = EvmRmsJA(LocAux);
                                    
                                    RxSymbAmos(CurrentTest,ObsCarrPos==ThisCarr,:) = EvmMatRecA(LocAux,:);%RxSymbAmos = [];
                                    EvmMatRec(ObsCarrPos==ThisCarr,:) = EvmMatRecA(LocAux,:);
                                    %% Ploting the result for qualitative analizes
                                    %PrintInfo(Ploting*19,TxData,Data);
                                    %%
                                    %             berpos1 = 1:2:size(BerDPSK,2);
                                    %             berpos2 = 2:2:size(BerDPSK,2);
                                    %             b = [BerDPSK(1,berpos1);BerDPSK(1,berpos2)]
                                    a=6;
                                    %                 close all;
                                end
                            end
                        case 'DQPSK'
                            
                            %%                   Receiver DQPSK
                            for ThisCarr=InitCarr:CarrPass:NumCarr                                           %For each carrier the same process of reception will be used.
                                if CarrUsedUp(ThisCarr)
                                    %%  Reception Process: Handling the income optical field
                                    %At the first step the income field will be selected from the
                                    %optical FFT Output with the help of the maping vector
                                    %(previously described).
                                    EoutAux = EoutAux1(VetThisCarr==ObsCarrUsed(ThisCarr),:);
                                    if AddingNoiseP
                                        [~,sigpower] = MeasPower(EoutAux);
                                        EoutAux      = awgn(EoutAux,osnrp,sigpower-30);
                                        [~,sigpower] = MeasPower(EoutAux);
                                    end
                                    %%            Fiber Time Delay Compensation
                                    switch Medium
                                        case 'Fiber'
                                            EoutAux = ifft(fft(EoutAux).*exp(1j*2*pi*f*(...
                                                FiberDelay(ThisCarr)*Ta)));
                                        otherwise
                                    end
                                    %             EoutAux = EoutRec;
                                    %%  Plot for Qualitative analizes
                                    %PrintInfo(Ploting*20,t,abs(EoutAux));
                                    %             PrintInfo(Ploting*21,f,20*log10(abs(fftshift(fft(EoutAux)./...
                                    %                                                        length(EoutAux)))));
                                    
                                    %%                   synchronizing
                                    %Differently from the previuos process at the DQPSK
                                    %demodulation it is needed to make the synchronism before
                                    %interferometric process. Because the recovery of the
                                    %transmited data is extremally dependent of the time delay.
                                    %Thus, for the reception process work properly, it is needed to
                                    %sincronized the recieved signal or the sampling process will
                                    %not work.
                                    
                                    %At this first moment a small part of the income signal needs
                                    %be analized for the sincronism process. As we are not looking
                                    %for the actual data rather for the synchronis symbol, which is
                                    %the bigest phase shift of the income signal. The next delay
                                    %interferometer convert this phase shift to a amplitude
                                    %variation.
                                    %
                                    %             %The phase delay is not needed because the optical field will
                                    %             %be analized as a whole.
                                    %             PhaDel          = 0;
                                    %             %Remember that at DQPSK modulation the information is stored at
                                    %             %the difference of phase between the current and previous
                                    %             %symbol hence this time delay is needed.
                                    %             TimDel          = T;
                                    %             [ESync1,ESync2] = DelayInterf(t,TimDel,PhaDel,EoutAux);        %Coverting phase shift to amplitude variation
                                    %
                                    %             %The second moment is to transfer this information from the
                                    %             %optical domain to the eletrical domain for an eletronic
                                    %             %processing. It is done with the help of an photo diode. The
                                    %             %configuration here used is an balanced receiver as the output
                                    %             %of the Delay Interferometer has two signals resulting from
                                    %             %the constructive and destructive signal interaction.
                                    %             ESync1 = ESync1.*conj(ESync1);
                                    %             ESync2 = ESync2.*conj(ESync2);
                                    %             Esync  = ESync2-ESync1;
                                    %
                                    %             ESync1 = ESync1./max(abs(ESync1));                             %Normalizing the signal
                                    %             ESync2 = ESync2./max(abs(ESync2));                             %Normalizing the signal
                                    %             Esync  = Esync./max(abs(Esync));                               %Normalizing the signal
                                    %
                                    %             SyncAux   = Esync(IniSyncPos:SyncPos);                         %Selecting just the symbol to synchronize
                                    %             SyncedAux = SyncSymb(IniSyncPos:SyncPos);                      %Selecting just the symbol to synchronize
                                    %             %%                   Plot for Qualitative analizes
                                    % %              PrintInfo(Ploting*22,t(IniSyncPos:SyncPos),Esync(IniSyncPos...
                                    % %              :SyncPos),SyncSymb(IniSyncPos:SyncPos),ESync1(IniSyncPos:...
                                    % %                                       SyncPos),ESync2(IniSyncPos:SyncPos));
                                    %             %%                   Synchronizing
                                    %             %This synchronization process is based on the mean expected
                                    %             %value. That means, the information mean value should be within
                                    %             %one period of symbol. Thus, the mean value of the received
                                    %             %signal is acquired and compare of the known sync-word to
                                    %             %verify if this mean value is at the right possition. Which is
                                    %             %the midel point (peak) of the highest level at the sync period
                                    %             SyncAux(SyncAux<0)              = 0;                           %To keep the mean value above zero anything under is neglected
                                    %             SyncAux(SyncAux>=mean(SyncAux)) = 1;                           %Adding a flag to the first sample of the received mean value
                                    %             SyncAux(SyncAux<mean(SyncAux))  = -1;                          %All the others samples at set to the lowest level
                                    %
                                    %             PosToSyn  = find(ismember(SyncAux,1));                         %Finding where is the location of the samples to synchronize
                                    %             PosSynced = find(ismember(SyncedAux,1));                       %Finding where is the location of the samples to synchronize
                                    %
                                    %             DiffPos = PosToSyn(round(end/2)) - PosSynced(round(end/2));    %Accounting the peak (midel point) displacement
                                    %             if DiffPos>=0%If the difference is positive, left-shift...
                                    %                 EoutAux = ifft(fft(EoutAux).*exp(1j*2*pi*f*(DiffPos*Ta))); %Shift based on time change
                                    % %                 EoutAux = [EoutAux(DiffPos+1:end) EoutAux(1:DiffPos)];   %Shift based on sampling sliding
                                    %             else%... but if the difference is negative, right-shift
                                    %                 EoutAux = ifft(fft(EoutAux).*exp(1j*2*pi*f*(DiffPos*Ta))); %Shift based on time change
                                    % %                 EoutAux = [EoutAux(end+DiffPos+1:end) EoutAux(1:end + ...
                                    % %                                                                 DiffPos)];%Shift based on sampling sliding
                                    %             end
                                    %
                                    %                 %                   Synchronizing
                                    %             %Because of reasons, sometimes it may be required to make a
                                    %             %synchronization process with the  end of the data stream as
                                    %             %well. This following verification check if the user set (or
                                    %             %not) a second synchronization process to be done.
                                    %             if SencondAdjust
                                    %                 PhaDel          = 0;
                                    %                 TimDel          = T;
                                    %                 [ESync1,ESync2] = DelayInterf(t,TimDel,PhaDel,EoutAux);
                                    %
                                    %                 ESync1 = ESync1.*conj(ESync1);
                                    %                 ESync2 = ESync2.*conj(ESync2);
                                    %                 Esync  = ESync2-ESync1;
                                    %
                                    %                 ESync1 = ESync1./max(abs(ESync1));
                                    %                 ESync2 = ESync2./max(abs(ESync2));
                                    %                 Esync  = Esync./max(abs(Esync));
                                    %
                                    %
                                    %                 %This synchronization process is based on the mean expected
                                    %                 %value. That means, the information mean value should be
                                    %                 %within one period of symbol. Thus, the mean value of the
                                    %                 %received signal is acquired and compare of the known
                                    %                 %sync-word to verify if this mean value is at the right
                                    %                 %possition.
                                    %
                                    %                 SyncAuxEnd   = Esync(end-SyncPos+1:end-IniSyncPos+1);
                                    %                 SyncedAuxEnd = SyncSymbEnd(end-SyncPos+1:end-IniSyncPos+1);
                                    %                 SyncAuxEnd(SyncAuxEnd<0)                 = 0;              %To keep the mean value above zero anything under is neglected
                                    %                 SyncAuxEnd(SyncAuxEnd>=mean(SyncAuxEnd)) = 1;              %Adding a flag to the first sample of the received mean value
                                    %                 SyncAuxEnd(SyncAuxEnd<mean(SyncAuxEnd))  = -1;             %All the others samples at set to the lowest level
                                    %
                                    %                 PosToSynEnd  = find(ismember(SyncAuxEnd,1));               %Finding where is the location of the first sample to synchronize
                                    %                 PosSyncedEnd = find(ismember(SyncedAuxEnd,1));             %Finding where is the location of the first sample to synchronize
                                    %
                                    %                 DiffPosEnd = PosToSynEnd(round(end/2))-PosSyncedEnd(...
                                    %                                                              round(end/2));
                                    %                 if DiffPosEnd>=0%If positive difference, left-shift...
                                    %                     EoutAux = ifft(fft(EoutAux).*exp(1j*2*pi*f*(...
                                    %                                                           DiffPosEnd*Ta)));%Shift based on time change
                                    % %                     EoutAux = [EoutAux(DiffPosEnd+1:end) EoutAux(1:...
                                    % %                                                             DiffPosEnd)];%Shift based on sampling sliding
                                    %                 else%... but if the difference is negative, right-shift
                                    %                     EoutAux = ifft(fft(EoutAux).*exp(1j*2*pi*f*(...
                                    %                                                           DiffPosEnd*Ta)));%Shift based on time change
                                    % %                     EoutAux = [EoutAux(end+DiffPosEnd+1:end) EoutAux(...
                                    % %                                                       1:end+DiffPosEnd)];%Shift based on sampling sliding
                                    %                 end
                                    %             end
                                    %
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
                                    
                                    E_rec3 = EoutAux./max(abs(EoutAux));                           %Normalizing income signal
                                    %For the interferometric process  take in account just the real
                                    %component it is needed a phase delay of 45° degrees;
                                    PhaDel = 1*pi/4;
                                    %Remember that at DQPSK modulation the information is stored at
                                    %the difference of phase between the current and previous
                                    %symbol hence this time delay of one symbol period is needed.
                                    TimDel = T;
                                    [EoutA,EoutB] = DelayInterf(t,TimDel,PhaDel,EoutAux);          %Coverting phase shift to amplitude variation
                                    %For the interferometric process  take in account just the real
                                    %component it is needed a phase delay of -45° degrees;
                                    PhaDel = -1*pi/4;
                                    %Remember that at DQPSK modulation the information is stored at
                                    %the difference of phase between the current and previous
                                    %symbol hence this time delay of one symbol period is needed.
                                    TimDel = T;
                                    [EoutC,EoutD] = DelayInterf(t,TimDel,PhaDel,EoutAux);          %Coverting phase shift to amplitude variation
                                    
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
                                    
                                    if ReceptorNoise                                               %Verify whether the noise is to be added or not
                                        [~,SigPowerA] = MeasPower(EoutA);                          %it own received signal or
                                        [~,SigPowerB] = MeasPower(EoutB);
                                        [~,SigPowerC] = MeasPower(EoutC);                          %it own received signal or
                                        [~,SigPowerD] = MeasPower(EoutD);
                                        SNR2 = CarSNR + 10*log10(2) - 10*log10(Nsamp) + 10*log10(10^0.3);
                                        
                                        if TimeSys
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
                                    
                                    EoutI = (EoutB - EoutA);
                                    EoutQ = (EoutD - EoutC);
                                    EoutI = EoutI-mean(EoutI);
                                    EoutQ = EoutQ-mean(EoutQ);
                                    EoutI = EoutI./max(abs(EoutI));                                %Normalizing the signal
                                    EoutQ = EoutQ./max(abs(EoutQ));                                %Normalizing the signal
                                    
                                    VetElecPowerI(CurrentTest,ThisCarr)= MeasPower(EoutI);
                                    VetElecPowerQ(CurrentTest,ThisCarr)= MeasPower(EoutQ);
                                    EoutAuxF = 20*log10(abs(fftshift(fft(EoutAux)./length(...
                                        EoutAux))));
                                    [VetOptiPower(CurrentTest,ThisCarr),~]= findpeaks(EoutAuxF,...
                                        'SortStr','descend','NPeaks',1);
                                    
                                    EoutI = EoutI./max(abs(EoutI));                                %Normalizing the signal
                                    EoutQ = EoutQ./max(abs(EoutQ));                                %Normalizing the signal
                                    
                                    
                                    SyncAuxI   = EoutI(IniSyncPos:SyncPos);                        %Selecting just the symbol to synchronize
                                    SyncAuxQ   = EoutQ(IniSyncPos:SyncPos);                        %Selecting just the symbol to synchronize
                                    %                  PrintInfo(Ploting*22,t(IniSyncPos:SyncPos),EoutI(IniSyncPos:SyncPos),SyncSymb(IniSyncPos:SyncPos),EoutA(IniSyncPos:SyncPos),EoutB(IniSyncPos:SyncPos));
                                    SyncedAux  = SyncSymb(IniSyncPos:SyncPos);                     %Selecting just the symbol to synchronize
                                    %%                   Synchronizing
                                    %This synchronization process is based on the mean expected
                                    %value. That means, the information mean value should be within
                                    %one period of symbol. Thus, the mean value of the received
                                    %signal is acquired and compare of the known sync-word to
                                    %verify if this mean value is at the right possition. Which is
                                    %the midel point (peak) of the highest level at the sync period
                                    %                 SyncAuxI(EoutI(IniSyncPos:SyncPos)<0.5*max((EoutI(IniSyncPos:SyncPos))))               = 0;                        %To keep the mean value above zero anything under is neglected
                                    %                 SyncAuxQ(EoutQ(IniSyncPos:SyncPos)<0.5*max((EoutQ(IniSyncPos:SyncPos))))               = 0;                        %To keep the mean value above zero anything under is neglected
                                    SyncAuxI(EoutI(IniSyncPos:SyncPos)<0)               = 0;                        %To keep the mean value above zero anything under is neglected
                                    SyncAuxQ(EoutQ(IniSyncPos:SyncPos)<0)               = 0;                        %To keep the mean value above zero anything under is neglected
                                    SyncAuxI(SyncAuxI>=mean(SyncAuxI)) = 1;                        %Adding a flag to the first sample of the received mean value
                                    SyncAuxQ(SyncAuxQ>=mean(SyncAuxQ)) = 1;                        %Adding a flag to the first sample of the received mean value
                                    if PrintinEye
                                        figure; hold all;
                                        plot(t(IniSyncPos:SyncPos),SyncSymb(IniSyncPos:SyncPos),t(IniSyncPos:SyncPos),EoutI(IniSyncPos:SyncPos),t(IniSyncPos:SyncPos),SyncAuxI,t(IniSyncPos:SyncPos),linspace(mean(SyncAuxI),mean(SyncAuxI),length(t(IniSyncPos:SyncPos))));
                                        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                        figure; hold all;
                                        plot(t(IniSyncPos:SyncPos),SyncSymb(IniSyncPos:SyncPos),t(IniSyncPos:SyncPos),EoutQ(IniSyncPos:SyncPos),t(IniSyncPos:SyncPos),SyncAuxQ,t(IniSyncPos:SyncPos),linspace(mean(SyncAuxQ),mean(SyncAuxQ),length(t(IniSyncPos:SyncPos))));
                                        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                    end
                                    SyncAuxI(SyncAuxI<mean(SyncAuxI))  = -1;                       %All the others samples at set to the lowest level
                                    SyncAuxQ(SyncAuxQ<mean(SyncAuxQ))  = -1;                       %All the others samples at set to the lowest level
                                    
                                    PosToSynI  = find(ismember(SyncAuxI,1));                       %Finding where is the location of the samples to synchronize
                                    PosToSynQ  = find(ismember(SyncAuxQ,1));                       %Finding where is the location of the samples to synchronize
                                    PosSynced = find(ismember(SyncedAux,1));                       %Finding where is the location of the samples to synchronize
                                    
                                    DiffPosI = ExtDel*(PosToSynI(round(end/2)) - PosSynced(round(end/2)));  %Accounting the peak (midel point) displacement
                                    DiffPosQ = ExtDel*(PosToSynQ(round(end/2)) - PosSynced(round(end/2)));  %Accounting the peak (midel point) displacement
                                    sn=((1/2)*(((1/2)^MaxNumStag)-1))/((1/2)-1);
                                    DiffPosI = round((sn/(3/2.55^IfftOrSum))*T/Ta);
                                    DiffPosQ = round((sn/(3/2.55^IfftOrSum))*T/Ta);
                                    if DiffPosI>=0%If the difference is positive, left-shift...
                                        %                 EoutI = ifft(fft(EoutI).*exp(1j*2*pi*f*(DiffPosI*Ta))); %Shift based on time change
                                        EoutI = [EoutI(DiffPosI+1:end) EoutI(1:DiffPosI)];   %Shift based on sampling sliding
                                        EoutQ = [EoutQ(DiffPosI+1:end) EoutQ(1:DiffPosI)];   %Shift based on sampling sliding
                                    else%... but if the difference is negative, right-shift
                                        %                 EoutI = ifft(fft(EoutI).*exp(-1j*2*pi*f*(DiffPosI*Ta))); %Shift based on time change
                                        EoutI = [EoutI(end+DiffPosI+1:end) EoutI(1:end+DiffPosI)]; %Shift based on sampling sliding
                                        EoutQ = [EoutQ(end+DiffPosI+1:end) EoutQ(1:end+DiffPosI)]; %Shift based on sampling sliding
                                    end
                                    if PrintinEye
                                        figure; hold all;
                                        plot(t(IniSyncPos:SyncPos),SyncSymb(IniSyncPos:SyncPos),t(IniSyncPos:SyncPos),EoutI(IniSyncPos:SyncPos));
                                        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                        figure; hold all;
                                        plot(t(IniSyncPos:SyncPos),SyncSymb(IniSyncPos:SyncPos),t(IniSyncPos:SyncPos),EoutQ(IniSyncPos:SyncPos));
                                        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                    end
                                    %             if DiffPosQ>=0%If the difference is positive, left-shift...
                                    % %                 EoutQ = (ifft(fft(EoutQ).*exp(1j*2*pi*f*(DiffPosQ*Ta)))); %Shift based on time change
                                    %                 EoutQ = [EoutQ(DiffPosQ+1:end) EoutQ(1:DiffPosQ)];   %Shift based on sampling sliding
                                    %             else%... but if the difference is negative, right-shift
                                    % %                 EoutQ = (ifft(fft(EoutQ).*exp(-1j*2*pi*f*(DiffPosQ*Ta)))); %Shift based on time change
                                    %                 EoutQ = [EoutQ(end+DiffPosQ+1:end) EoutQ(1:end+DiffPosQ)]; %Shift based on sampling sliding
                                    %             end
                                    %%                   Plot for Qualitative analizes
                                    %                              PrintInfo(Ploting*22,t(IniSyncPos:SyncPos),EoutI(IniSyncPos...
                                    %                              :SyncPos),SyncSymb(IniSyncPos:SyncPos),EoutA(IniSyncPos:...
                                    %                                                       SyncPos),EoutB(IniSyncPos:SyncPos));
                                    %              PrintInfo(Ploting*22,t(IniSyncPos:SyncPos),EoutQ(IniSyncPos...
                                    %              :SyncPos),SyncSymb(IniSyncPos:SyncPos),EoutC(IniSyncPos:...
                                    %                                       SyncPos),EoutD(IniSyncPos:SyncPos));
                                    %% Removing CP
                                    %                 if ThisCarr==126
                                    %                     EyeToPlot(CurrentTest,1:length([EoutI EoutQ])) = [EoutI EoutQ];
                                    %                     save(['savingforplotingeye' num2str(VetSnr)],'EyeToPlot','VetSnr','IfftOrSum','UsedModula','T','NPPB');
                                    %                 end
                                    if AddCP
                                        IxAux1  = EoutI(1:end - StuffSampels);
                                        IxAux1  = reshape(IxAux1,(2*NumAmosCP+NPPB),NbDQPSK/2);
                                        IxAux1  = IxAux1(1+NumAmosCP:end-NumAmosCP,:);
                                        IxAux1  = reshape(IxAux1,1,NPPB*NbDQPSK/2);
                                        EoutI  = IxAux1;
                                        
                                        IxAux2  = EoutQ(1:end - StuffSampels);
                                        IxAux2  = reshape(IxAux2,(2*NumAmosCP+NPPB),NbDQPSK/2);
                                        IxAux2  = IxAux2(1+NumAmosCP:end-NumAmosCP,:);
                                        IxAux2  = reshape(IxAux2,1,NPPB*NbDQPSK/2);
                                        EoutQ  = IxAux2;
                                        
                                        IxAux3  = E_rec3(1:end - StuffSampels);
                                        IxAux3  = reshape(IxAux3,(2*NumAmosCP+NPPB),NbDQPSK/2);
                                        IxAux3  = IxAux3(1+NumAmosCP:end-NumAmosCP,:);
                                        IxAux3  = reshape(IxAux3,1,NPPB*NbDQPSK/2);
                                        E_rec3  = IxAux3;
                                    end
                                    %% Taking the sampling the EVM meassurement
                                    clear IxAux;
                                    PosAuxEout1 = NPPB/2:NPPB:length(EoutI);%Varriable respossible to take just the samples at the middle of the symbol
                                    PosAuxEout2 = ((NPPB/2)+(NPPB/PerDiv)):NPPB:length(EoutI);%Varriable respossible to take just the samples at the middle of the symbol
                                    PosAuxEout3 = ((NPPB/2)-(NPPB/PerDiv)):NPPB:length(EoutI);%Varriable respossible to take just the samples at the middle of the symbol
                                    IxAux1      = EoutI(PosAuxEout1) + 1j.*EoutQ(PosAuxEout1);    %Normalizing the reference
                                    IxAux2      = EoutI(PosAuxEout2) + 1j.*EoutQ(PosAuxEout2);    %Normalizing the reference
                                    IxAux3      = EoutI(PosAuxEout3) + 1j.*EoutQ(PosAuxEout3);    %Normalizing the reference
                                    a=0;
                                    EvmMatRecA(1,:) = IxAux1;                       %Taking just the middle samples as references
                                    EvmMatRecA(2,:) = IxAux2;                       %Taking just the middle samples as references
                                    EvmMatRecA(3,:) = IxAux3;                       %Taking just the middle samples as references
                                    EvmMatRecA(4,:) = IxAux1;                       %Taking just the middle samples as references
                                    [EvmDBA(1), EvmPerA(1), EvmRmsA(1) ] = EvmCalc( EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux1 );
                                    [EvmDBA(2), EvmPerA(2), EvmRmsA(2) ] = EvmCalc( EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux2 );
                                    [EvmDBA(3), EvmPerA(3), EvmRmsA(3) ] = EvmCalc( EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux3 );
                                    [EvmDBA(4), EvmPerA(4), EvmRmsA(4) ] = EvmCalc( EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux1 );
                                    [EvmDBJA(1),EvmPerJA(1),EvmRmsJA(1)] = evm1(4,'dpsk',EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux1);
                                    [EvmDBJA(2),EvmPerJA(2),EvmRmsJA(2)] = evm1(4,'dpsk',EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux2);
                                    [EvmDBJA(3),EvmPerJA(3),EvmRmsJA(3)] = evm1(4,'dpsk',EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux3);
                                    [EvmDBJA(4),EvmPerJA(4),EvmRmsJA(4)] = evm1(4,'dpsk',EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux1);
                                    %                 if ThisCarr==2
                                    %                     figure;hold all;grid on;
                                    %                     plot(real(IxAux1),imag(IxAux1),'rx');
                                    %                     plot(real(IxAux2),imag(IxAux2),'yx');
                                    %                     plot(real(IxAux3),imag(IxAux3),'gx');
                                    %                     set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                    %                     a=0;
                                    %                 end
                                    %One important step on this project was the confirmation of our
                                    %results from a closed and theoretical equation that relates
                                    %the income optical field witht its respective data (as
                                    %descrived above in the receiving process). The result from
                                    %this equation was further compared with the result from the
                                    %MZ-Interferometer as an proof of concept. This equation can be
                                    %found at the book of Optical Fiber Telecommunications V B,
                                    %which one of the authors is Ivan P. Kaminow at the page 144.
                                    taux = t(1:length(E_rec3));
                                    faux = time2freq(taux);
                                    Ui = real(exp(1j*pi/4).*E_rec3.*conj(ifft(fft(E_rec3).*exp(-...
                                        1j*2*pi*faux*T))));%The data at odd position
                                    Ui = Ui./max(abs(Ui));                                         %Normalizing the signal
                                    Uq = imag(exp(1j*pi/4).*E_rec3.*conj(ifft(fft(E_rec3).*exp(-...
                                        1j*2*pi*faux*T))));%The data at the even position
                                    Uq = Uq./max(abs(Uq));                                         %Normalizing the signal
                                    %##########################################################################
                                    PosIx = NPPB/2:NPPB:length(EoutI);                                %Possition of the central samples - but the number of samples per symbol is even ????
                                    IxAux = EoutI(PosIx);                                             %From the main received signal just a few samples are taken for further evaluation
                                    n     = 100;                                                   %The number of boxes to be filled up on the histogram process
                                    Perce = 0.8;
                                    %At this process each eye will be evaluated separetelly. It is
                                    %important to remember that the PAM4 has 4 levels, which means
                                    %three level of decissions that are exaclty the center of the
                                    %eye diagram.
                                    IxAuxAB = EoutI(PosIx);       %Taking just those values relative to the uper eye
                                    InterAB = linspace(Perce*(min(EoutI)),Perce*(max(EoutI)),n);                     %Building the histogram boxes
                                    EyeAB = hist(IxAuxAB,InterAB);                                 %filling up the boxes with samples that fit on them.
                                    
                                    QIxAuxAB = EoutQ(PosIx);       %Taking just those values relative to the uper eye
                                    QInterAB = linspace(Perce*(min(EoutQ)),Perce*(max(EoutQ)),n);                     %Building the histogram boxes
                                    QEyeAB = hist(QIxAuxAB,QInterAB);                                 %filling up the boxes with samples that fit on them.
                                    %The same process described for the uper level will be done at
                                    %the middle and lower eyes levels.
                                    IxAuxCD = EoutI(PosIx-NPPB/PerDiv);
                                    InterCD = linspace(Perce*(min(EoutI)),Perce*(max(EoutI)),n);%NPPB*2^n);
                                    EyeCD = hist(IxAuxCD,InterCD);
                                    
                                    QIxAuxCD = EoutQ(PosIx-NPPB/PerDiv);
                                    QInterCD = linspace(Perce*(min(EoutQ)),Perce*(max(EoutQ)),n);%NPPB*2^n);
                                    QEyeCD = hist(QIxAuxCD,QInterCD);
                                    
                                    IxAuxEF = EoutI(PosIx+NPPB/PerDiv);
                                    InterEF = linspace(Perce*(min(EoutI)),Perce*(max(EoutI)),n);%NPPB*2^n);
                                    EyeEF = hist(IxAuxEF,InterEF);
                                    
                                    QIxAuxEF = EoutQ(PosIx+NPPB/PerDiv);
                                    QInterEF = linspace(Perce*(min(EoutQ)),Perce*(max(EoutQ)),n);%NPPB*2^n);
                                    QEyeEF = hist(QIxAuxEF,QInterEF);
                                    
                                    %What we are looking for here are not where exist the
                                    %occurrences of values. The eye is where there is not samples,
                                    %this means, the decission level is a reagion between actual
                                    %levels thus this reagion does not contain any signal.
                                    EyeAB = ~EyeAB;                                                %Changing zeros to one - Zeros compose the eye region
                                    EyeCD = ~EyeCD;
                                    EyeEF = ~EyeEF;
                                    
                                    QEyeAB = ~QEyeAB;                                                %Changing zeros to one - Zeros compose the eye region
                                    QEyeCD = ~QEyeCD;
                                    QEyeEF = ~QEyeEF;
                                    %Starting variables that will be used to identify the eye
                                    %diagram.
                                    CountAB   = 1;
                                    CountCD   = 1;
                                    CountEF   = 1;
                                    
                                    QCountAB   = 1;
                                    QCountCD   = 1;
                                    QCountEF   = 1;
                                    
                                    SeqOnesAB = 0;
                                    SeqOnesCD = 0;
                                    SeqOnesEF = 0;
                                    
                                    QSeqOnesAB = 0;
                                    QSeqOnesCD = 0;
                                    QSeqOnesEF = 0;
                                    
                                    SeqFinAB  = 0;
                                    SeqFinCD  = 0;
                                    SeqFinEF  = 0;
                                    
                                    QSeqFinAB  = 0;
                                    QSeqFinCD  = 0;
                                    QSeqFinEF  = 0;
                                    
                                    SeqIniAB  = 1;
                                    SeqIniCD  = 1;
                                    SeqIniEF  = 1;
                                    
                                    QSeqIniAB  = 1;
                                    QSeqIniCD  = 1;
                                    QSeqIniEF  = 1;
                                    %The for loop will take account of every box with ones. It is
                                    %important to take note that the not operator was used in this
                                    %vector, therefore ones means zeros (the eye diagram -
                                    %possibly) and zeros means values abouve zeroa (not the eye).
                                    for kk=1:length(EyeAB)                                         %For every box
                                        if EyeAB(kk)                                               %if it contains "1"
                                            SeqOnesAB(SeqIniAB)=CountAB;                           %count this element as part of a consecutive sequency
                                            CountAB = CountAB + 1;                                 %adds one to the counter of consecutive elements "1"
                                            if kk==length(EyeAB)                                   %if the current box is the last box we got to an end
                                                SeqFinAB(SeqIniAB) = kk;                           %The final sequency element is equal to its possition (kk)
                                            end
                                        else                                                       %else if the current box contains "0"
                                            SeqFinAB(SeqIniAB) = kk-1;                             %the previous element was the last element of a consecutive sequency of ones
                                            SeqIniAB = SeqIniAB + 1;                               %adds one to the counter of consecutive number (take in account how many sequencies there are)
                                            CountAB = 1;                                           %reset the counter of consecutive elements to mark the start of a new consecutive sequence
                                        end
                                        
                                        if QEyeAB(kk)                                               %if it contains "1"
                                            QSeqOnesAB(QSeqIniAB)=QCountAB;                           %count this element as part of a consecutive sequency
                                            QCountAB = QCountAB + 1;                                 %adds one to the counter of consecutive elements "1"
                                            if kk==length(QEyeAB)                                   %if the current box is the last box we got to an end
                                                QSeqFinAB(QSeqIniAB) = kk;                           %The final sequency element is equal to its possition (kk)
                                            end
                                        else                                                       %else if the current box contains "0"
                                            QSeqFinAB(QSeqIniAB) = kk-1;                             %the previous element was the last element of a consecutive sequency of ones
                                            QSeqIniAB = QSeqIniAB + 1;                               %adds one to the counter of consecutive number (take in account how many sequencies there are)
                                            QCountAB = 1;                                           %reset the counter of consecutive elements to mark the start of a new consecutive sequence
                                        end
                                        
                                        if EyeCD(kk)
                                            SeqOnesCD(SeqIniCD)=CountCD;
                                            CountCD = CountCD + 1;
                                            if kk==length(EyeCD)
                                                SeqFinCD(SeqIniCD) = kk;
                                            end
                                        else
                                            SeqFinCD(SeqIniCD) = kk-1;
                                            SeqIniCD = SeqIniCD + 1;
                                            CountCD = 1;
                                        end
                                        
                                        if QEyeCD(kk)
                                            QSeqOnesCD(QSeqIniCD)=QCountCD;
                                            QCountCD = QCountCD + 1;
                                            if kk==length(QEyeCD)
                                                QSeqFinCD(QSeqIniCD) = kk;
                                            end
                                        else
                                            QSeqFinCD(QSeqIniCD) = kk-1;
                                            QSeqIniCD = QSeqIniCD + 1;
                                            QCountCD = 1;
                                        end
                                        
                                        if EyeEF(kk)
                                            SeqOnesEF(SeqIniEF)=CountEF;
                                            CountEF = CountEF + 1;
                                            if kk==length(EyeEF)
                                                SeqFinEF(SeqIniEF) = kk;
                                            end
                                        else
                                            SeqFinEF(SeqIniEF) = kk-1;
                                            SeqIniEF = SeqIniEF + 1;
                                            CountEF = 1;
                                        end
                                        
                                        if QEyeEF(kk)
                                            QSeqOnesEF(QSeqIniEF)=QCountEF;
                                            QCountEF = QCountEF + 1;
                                            if kk==length(QEyeEF)
                                                QSeqFinEF(QSeqIniEF) = kk;
                                            end
                                        else
                                            QSeqFinEF(QSeqIniEF) = kk-1;
                                            QSeqIniEF = QSeqIniEF + 1;
                                            QCountEF = 1;
                                        end
                                    end
                                    
                                    
                                    %If the eye is open, which means there is a clear difference
                                    %between adjacent levels, the eye diagram will be the longest
                                    %sequence of ones.
                                    [MaxValAB,LocMaxAB]=max(SeqOnesAB);
                                    if LocMaxAB<2 || MaxValAB<2                                    %if any sequency was found or there is just one sequency it is a error thus
                                        LevDec3 = 0.1;                                            %the decission level will be by default 0.67. Also, the other variables will
                                        LocMaxAB = 1;                                              %will be set with values to not cause errors in the future
                                        SeqFinAB(1)=2;
                                        MaxValAB = 0;
                                        InterAB(1)=LevDec3;
                                    else                                                           %if a sequency was found the middle of the eye will be the middle of the sequency
                                        if (SeqFinAB(LocMaxAB)-MaxValAB/2)<1                       %if for some reason the final element of the sequency minus the half of its
                                            LevDec3 = 0.1;                                         %length results in a negative value, something went very wrong, and by
                                        else                                                       %default it will be set to 0.7
                                            LevDec3 = InterAB(round(SeqFinAB(LocMaxAB)-MaxValAB/...
                                                2));%Otherwise, the decission level is the middle of the sequency
                                        end
                                    end
                                    [QMaxValAB,QLocMaxAB]=max(QSeqOnesAB);
                                    if QLocMaxAB<2 || QMaxValAB<2                                    %if any sequency was found or there is just one sequency it is a error thus
                                        QLevDec3 = 0.1;                                            %the decission level will be by default 0.67. Also, the other variables will
                                        QLocMaxAB = 1;                                              %will be set with values to not cause errors in the future
                                        QSeqFinAB(1)=2;
                                        QMaxValAB = 0;
                                        QInterAB(1)=QLevDec3;
                                    else                                                           %if a sequency was found the middle of the eye will be the middle of the sequency
                                        if (QSeqFinAB(QLocMaxAB)-QMaxValAB/2)<1                       %if for some reason the final element of the sequency minus the half of its
                                            QLevDec3 = 0.1;                                         %length results in a negative value, something went very wrong, and by
                                        else                                                       %default it will be set to 0.7
                                            QLevDec3 = QInterAB(round(QSeqFinAB(QLocMaxAB)-QMaxValAB/...
                                                2));%Otherwise, the decission level is the middle of the sequency
                                        end
                                    end
                                    
                                    [MaxValCD,LocMaxCD]=max(SeqOnesCD);
                                    if LocMaxCD<2 || MaxValCD<2
                                        LevDec2 = 0.1;
                                        LocMaxCD = 1;
                                        SeqFinCD(1)=2;
                                        MaxValCD = 0;
                                        InterCD(1)=LevDec2;
                                    else
                                        if (SeqFinCD(LocMaxCD)-MaxValCD/2)<1
                                            LevDec2 = 0.1;
                                        else
                                            LevDec2 = InterCD(round(SeqFinCD(LocMaxCD)-MaxValCD/2));
                                        end
                                    end
                                    [QMaxValCD,QLocMaxCD]=max(QSeqOnesCD);
                                    if QLocMaxCD<2 || QMaxValCD<2
                                        QLevDec2 = 0.1;
                                        QLocMaxCD = 1;
                                        QSeqFinCD(1)=2;
                                        QMaxValCD = 0;
                                        QInterCD(1)=LevDec2;
                                    else
                                        if (QSeqFinCD(QLocMaxCD)-QMaxValCD/2)<1
                                            QLevDec2 = 0.1;
                                        else
                                            QLevDec2 = QInterCD(round(QSeqFinCD(QLocMaxCD)-QMaxValCD/2));
                                        end
                                    end
                                    
                                    [MaxValEF,LocMaxEF]=max(SeqOnesEF);
                                    if LocMaxEF<2 || MaxValEF<2
                                        LevDec1 = 0.1;
                                        LocMaxEF = 1;
                                        SeqFinEF(1)=2;
                                        MaxValEF = 0;
                                        InterEF(1)=LevDec1;
                                    else
                                        if (SeqFinEF(LocMaxEF)-MaxValEF/2)<1
                                            LevDec1 = 0.1;
                                        else
                                            LevDec1 = InterEF(round(SeqFinEF(LocMaxEF)-MaxValEF/2));
                                        end
                                        
                                    end
                                    [QMaxValEF,QLocMaxEF]=max(QSeqOnesEF);
                                    if QLocMaxEF<2 || QMaxValEF<2
                                        QLevDec1 = 0.1;
                                        QLocMaxEF = 1;
                                        QSeqFinEF(1)=2;
                                        QMaxValEF = 0;
                                        QInterEF(1)=QLevDec1;
                                    else
                                        if (QSeqFinEF(QLocMaxEF)-QMaxValEF/2)<1
                                            QLevDec1 = 0.1;
                                        else
                                            QLevDec1 = QInterEF(round(QSeqFinEF(QLocMaxEF)-QMaxValEF/2));
                                        end
                                        
                                    end
                                    
                                    %another way to measure the eye opening is the get all the
                                    %boxes and find all peaks on it, that will be a plato created
                                    %by the sequences of ones (which was zeros). From thos peaks,
                                    %the eye diagram will be the longer of them hence it will take
                                    %the most part of the vector that store the findpeaks result.
                                    %Thus, the middle of the eye will be basically the middle of
                                    %the peaks vector.
                                    LocAB = find(EyeAB);
                                    if isempty(LocAB)                                              %if for some reason there are no peaks, something went wrong.
                                        LocAB = 1;
                                        LevelDec3 = 0.1;%mean(Levels(3:4));                       %by default the decission level will be set to 0.65
                                    else
                                        LevelDec3 = LevDec3;
                                    end
                                    QLocAB = find(QEyeAB);
                                    if isempty(QLocAB)                                              %if for some reason there are no peaks, something went wrong.
                                        QLocAB = 1;
                                        QLevelDec3 = 0.1;%mean(Levels(3:4));                       %by default the decission level will be set to 0.65
                                    else
                                        QLevelDec3 = QLevDec3;
                                    end
                                    
                                    LocCD = find(EyeCD);
                                    if isempty(LocCD)
                                        LocCD = 1;
                                        LevelDec2 = 0.1;%mean(Levels(2:3));
                                    else
                                        LevelDec2 = LevDec2;
                                    end
                                    QLocCD = find(QEyeCD);
                                    if isempty(QLocCD)
                                        QLocCD = 1;
                                        QLevelDec2 = 0.1;%mean(Levels(2:3));
                                    else
                                        QLevelDec2 = QLevDec2;
                                    end
                                    
                                    LocEF = find(EyeEF);
                                    if isempty(LocEF)
                                        LocEF = 1;
                                        LevelDec1 = 0.1;%mean(Levels(1:2));
                                    else
                                        LevelDec1 = LevDec1;
                                    end
                                    QLocEF = find(QEyeEF);
                                    if isempty(QLocEF)
                                        QLocEF = 1;
                                        QLevelDec1 = 0.1;%mean(Levels(1:2));
                                    else
                                        QLevelDec1 = QLevDec1;
                                    end
                                    %##########################################################################
                                    if CalcS
                                        [~,~,EyeOpenI,EyeOpenHighI,EyeOpenLowI] = Olho(EoutI,T,NPPB,0);
                                        [~,~,EyeOpenQ,EyeOpenHighQ,EyeOpenLowQ] = Olho(EoutQ,T,NPPB,0);
                                        AberLevIS(CurrentTest,ThisCarr)  = EyeOpenI;
                                        ValsLevIS(CurrentTest,ThisCarr)  = EyeOpenLowI + EyeOpenI/2;
                                        AberLevQS(CurrentTest,ThisCarr)  = EyeOpenQ;
                                        ValsLevQS(CurrentTest,ThisCarr)  = EyeOpenLowQ + EyeOpenQ/2;
                                    end
                                    
                                    AberLevAuxI(1)= InterAB(SeqFinAB(LocMaxAB))-InterAB(round(SeqFinAB(LocMaxAB)-MaxValAB));
                                    AberLevAuxI(3)= InterCD(SeqFinCD(LocMaxCD))-InterCD(round(SeqFinCD(LocMaxCD)-MaxValCD));
                                    AberLevAuxI(2)= InterEF(SeqFinEF(LocMaxEF))-InterEF(round(SeqFinEF(LocMaxEF)-MaxValEF));
                                    ValsLevAuxI(1)= LevDec3;
                                    ValsLevAuxI(3)= LevDec2;
                                    ValsLevAuxI(2)= LevDec1;
                                    AberLevAuxQ(1)= QInterAB(QSeqFinAB(QLocMaxAB))-QInterAB(round(QSeqFinAB(QLocMaxAB)-QMaxValAB));
                                    AberLevAuxQ(3)= QInterCD(QSeqFinCD(QLocMaxCD))-QInterCD(round(QSeqFinCD(QLocMaxCD)-QMaxValCD));
                                    AberLevAuxQ(2)= QInterEF(QSeqFinEF(QLocMaxEF))-QInterEF(round(QSeqFinEF(QLocMaxEF)-QMaxValEF));
                                    ValsLevAuxQ(1)= QLevDec3;
                                    ValsLevAuxQ(3)= QLevDec2;
                                    ValsLevAuxQ(2)= QLevDec1;
                                    %% Ploting the result for qualitative analizes
                                    if PrintinEye
                                        PrintInfo(Ploting*24,EoutI,T,NPPB);
                                        hold on;
                                        plot(t((NPPB)/2),InterAB(round(SeqFinAB(LocMaxAB)-MaxValAB)),'kx');
                                        plot(t((NPPB)/2),InterAB(SeqFinAB(LocMaxAB)),'ko');
                                        plot(t((NPPB)/2),InterAB(round(SeqFinAB(LocMaxAB)-MaxValAB/2)),'kd');
                                        plot(t((NPPB)/PerDiv),InterCD(round(SeqFinCD(LocMaxCD)-MaxValCD)),'gx');
                                        plot(t((NPPB)/PerDiv),InterCD(SeqFinCD(LocMaxCD)),'go');
                                        plot(t((NPPB)/PerDiv),InterCD(round(SeqFinCD(LocMaxCD)-MaxValCD/2)),'gd');
                                        plot(t(3*(NPPB)/PerDiv),InterEF(round(SeqFinEF(LocMaxEF)-MaxValEF)),'bx');
                                        plot(t(3*(NPPB)/PerDiv),InterEF(SeqFinEF(LocMaxEF)),'bo');
                                        plot(t(3*(NPPB)/PerDiv),InterEF(round(SeqFinEF(LocMaxEF)-MaxValEF/2)),'bd');
                                        %PrintInfo(Ploting*25,Ui,T,NPPB);
                                        PrintInfo(Ploting*26,EoutQ,T,NPPB);
                                        hold on;
                                        plot(t((NPPB)/2),QInterAB(round(QSeqFinAB(QLocMaxAB)-QMaxValAB)),'kx');
                                        plot(t((NPPB)/2),QInterAB(QSeqFinAB(QLocMaxAB)),'ko');
                                        plot(t((NPPB)/2),QInterAB(round(QSeqFinAB(QLocMaxAB)-QMaxValAB/2)),'kd');
                                        plot(t((NPPB)/PerDiv),QInterCD(round(QSeqFinCD(QLocMaxCD)-QMaxValCD)),'gx');
                                        plot(t((NPPB)/PerDiv),QInterCD(QSeqFinCD(QLocMaxCD)),'go');
                                        plot(t((NPPB)/PerDiv),QInterCD(round(QSeqFinCD(QLocMaxCD)-QMaxValCD/2)),'gd');
                                        plot(t(3*(NPPB)/PerDiv),QInterEF(round(QSeqFinEF(QLocMaxEF)-QMaxValEF)),'bx');
                                        plot(t(3*(NPPB)/PerDiv),QInterEF(QSeqFinEF(QLocMaxEF)),'bo');
                                        plot(t(3*(NPPB)/PerDiv),QInterEF(round(QSeqFinEF(QLocMaxEF)-QMaxValEF/2)),'bd');
                                        if CalcS
                                            plot(t((NPPB)/2),EyeOpenLowI,'mx');
                                            plot(t((NPPB)/2),EyeOpenHighI,'mo');
                                            plot(t((NPPB)/2),EyeOpenLowI + EyeOpenI/2,'md');
                                            plot(t((NPPB)/2),EyeOpenLowQ,'cx');
                                            plot(t((NPPB)/2),EyeOpenHighQ,'co');
                                            plot(t((NPPB)/2),EyeOpenLowQ + EyeOpenQ/2,'cd');
                                        end
                                        figure;
                                        hold all;
                                        plotpos = zeros(1,length(IxAux1));
                                        plot(IxAux1,'o','color',[1 0.4 0]);
                                        plot(EvmMatRef(ObsCarrUsed(ThisCarr),:),'b*');
                                        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                        a=0;
                                    end
                                    %PrintInfo(Ploting*27,Uq,T,NPPB);
                                    %%
                                    
                                    %                 Txaux1 = rectpulse(TxDataMat(ThisCarr,:),NPPB);
                                    %                 Txaux2 = rectpulse(TxDataMat(ThisCarr+NumCarr,:),NPPB);
                                    %                 Txaux1(Txaux1==0) = -1;
                                    %                 Txaux2(Txaux2==0) = -1;
                                    %%  Ploting some results for qualitative analizes
                                    %             PrintInfo(Ploting*28,t(1:length(EoutI)),Txaux1,Txaux2,EoutI,...
                                    %                           EoutA,EoutB,EoutQ,EoutC,EoutD,real(Ui),real(Uq));
                                    %%               Recovering the Information
                                    %After passing the optical signal to the eletrical domain, for
                                    %actually detect the data withing the signal the following
                                    %steps are needed.
                                    %
                                    %Finding Decission Levels:
                                    %The process for decoding the income signal will be based on
                                    %eletronic comparators. Inasmuch as the right decission level
                                    %must be acquired for accurately decide, within a symbol
                                    %periode, what that current leavel means (ones or zeros).
                                    %
                                    %The process hereafter of chosing the  decission levels is not
                                    %deterministic rather it is a statistic process. The main idea
                                    %is to take the decission level from the histogram generated
                                    %from the income signal stream.
                                    %
                                    %This process is realized inside the function Olho.
                                    %
                                    %Basicaly the decission level will be the minimal value of the
                                    %currente eye under evaluation plus the half of the its eye
                                    %opening.The following ilustration better describe this process
                                    %
                                    %Eye Limits:   + (1/2 * Eye Opening:)  =    Comparisson Limit:
                                    %
                                    %UperLevel  ______     ______________     _____
                                    %                 \   /      |       \   /
                                    %                  \ /       |        \ /
                                    %                   \ Half Eye Opening /     Decission Level
                                    %                  / \       |        / \
                                    %LowerLevel ______/   \______|_______/   \_____
                                    %
                                    %
                                    %Actualy Receiving Data:
                                    %Once the signal was processed the next step is through a
                                    %comparator decide the actual information received.
                                    %
                                    %
                                    %For accurately read the income signal is necessary more than
                                    %just have the decision level, it is also needed to know where
                                    %it occurs. For instance, if we decide to use just to take the
                                    %decision level and measure the mean value across and a portion
                                    %of the time period, which portion should we take? It would be
                                    %logical to take the central portion, but the bit may be
                                    %deformed in such way that the information gets concentrated
                                    %out of the middle part, which means the symbol is not
                                    %symmetric. The symmetry information can be acquired from the
                                    %eye diagram by measuring the longitudinal opening. The
                                    %following sketch better describes this process:
                                    %
                                    %                   Point of Symmetry
                                    %         _____     _______|_______     _____
                                    %              \   /               \   /
                                    %               \ /  Longitudinal   \ /
                                    %                \ ----------------- /
                                    %               / \    Opening      / \
                                    %         _____/   \_______________/   \_____
                                    %
                                    %With those two pieces of information, decision level and point
                                    %of symmetry, we have the X, Y coordinates for the centre of
                                    %the Eye Diagram. Therefore, as long as there is an opening on
                                    %it it will be possible to recover the transmitted information
                                    %without error... theoretically.
                                    %
                                    %As this process is also statistical, first we reshape the
                                    %income vector to analyze all periods at the same time.
                                    %                 EyeSymMatI = reshape(EoutI(1+SyncPeriod*NPPB:end-SyncPeriod*...
                                    %                     NPPB),NPPB,NbDQPSK/2-2*SyncPeriod);
                                    %                 %Then we take the values that compose the decision level
                                    %                 %because they will mark the point of symmetry.
                                    %                 %
                                    %                 %Firstly it was set the interval in which the histogram will be
                                    %                 %build. It is based on the number of samples per bit period.
                                    %                 IntervalI = linspace(min(EoutI(1+SyncPeriod*NPPB:end-...
                                    %                     SyncPeriod*NPPB)),max(EoutI(1+SyncPeriod*NPPB:end-...
                                    %                     SyncPeriod*NPPB)),2*NPPB);
                                    %                 %Therefore, the MATLAB hist function returns the number of
                                    %                 %occurrence of each interval.
                                    %                 EyeMaxI = hist(EoutI,IntervalI);
                                    %                 EyeMaxauxI = [0 EyeMaxI 0];                                    %Zeros are added at the EyeMax to auxiliate the finding peaks process
                                    %                 [~,EyeLocI] = findpeaks(EyeMaxauxI,'MinPeakDistance',NPPB/4,...
                                    %                     'SortStr','descend','NPeaks',2);%The peaks on the Eye profile will be the levels at the Eyes limit
                                    %                 %From the location of the max values that occured, which means
                                    %                 %the uper and lower level of the eye diagram it needs to take
                                    %                 %the actual value that those occurences represent that is
                                    %                 %withing the the Interval variable.
                                    %                 ValToSeekI = IntervalI(EyeLocI-1);
                                    %                 ValToSeekI = sort(ValToSeekI,'ascend');
                                    %                 %The number of ocurrences is a statical measure therefore one
                                    %                 %does not have control which interval will have the highest
                                    %                 %peak, thus it is important to ordenate the values to be seek
                                    %                 %from the lower part of the eye diagram to the uper part of the
                                    %                 %eye diagram.
                                    %                 OccuCountI = zeros(1,size(EyeSymMatI,1));                      %Auxiliar Variable for accounting.
                                    %                 for kk=1:size(EyeSymMatI,1)                                    %For every sample within a symbol period
                                    %
                                    %                     OccuCountI(kk) = OccuCountI(kk)+sum((EyeSymMatI(kk,:)>=...
                                    %                         min(EoutI))&(EyeSymMatI(kk,:)<=UpeSymPer*EyeOpenLowI)); %Account all occurencies of the valeu 1
                                    %                     OccuCountI(kk) = OccuCountI(kk)+sum((EyeSymMatI(kk,:)>=...
                                    %                         LowSymPer*EyeOpenHighI)&(EyeSymMatI(kk,:)<=max(EoutI))); %Account all occurencies of the valeu 2
                                    %                 end
                                    %The point of symmetry of the eye diagram will be where the
                                    %maximum number of occurrences were measured inasmuch as those
                                    %are the points where all the bits go to the center of the
                                    %symbol.  From the maximum number of occurrences, it can happen
                                    %for more than one sample within one symbol period, in the best
                                    %case, all samples would have the same accounting as it is
                                    %shown the ilustration above hence the symmetry will be at the
                                    %middle sample of this group of maximum occurrences. This value
                                    %can be found by the mean of the samples positions within the
                                    %symbol period. The problem with this approach is that the
                                    %signal must be synchronized with the maximum displacement of
                                    %a symbol period minus 25% of the eye Longitudinal opening if
                                    %the displacement is higher than that the point of symmetry
                                    %will be wrongly measured.
                                    %             [SymLocI] = round(mean(find(ismember(OccuCountI,max(...
                                    %                                                            OccuCountI)))));
                                    %             PrintInfo(Ploting*24,EoutI,T,NPPB);
                                    %             hold all;
                                    %             plot(t(NPPB/2),EyeOpenLowI+EyeOpenI/2,'bd');plot(t(NPPB/2),EyeOpenLowI,'bo');plot(t(NPPB/2),EyeOpenHighI,'bx');
                                    %             plot([t(SymLocI(1)) t(SymLocI(1))],[-1 1]);
                                    
                                    %             [~,SymLocI] = findpeaks(OccuCountI,'SortStr','descend');       %The peak on the Eye profile will be the Symmetry level
                                    DataOdd = [];                                                  %Initialization of the vector that will store the income data
                                    DataOddM = [];                                                  %Initialization of the vector that will store the income data
                                    DataOddL = [];                                                  %Initialization of the vector that will store the income data
                                    DataOddU = [];                                                  %Initialization of the vector that will store the income data
                                    DataOddU2 = [];                                                  %Initialization of the vector that will store the income data
                                    DataOddU3 = [];                                                  %Initialization of the vector that will store the income data
                                    DataOddS = [];                                                  %Initialization of the vector that will store the income data
                                    for kk=1:NPPB:length(EoutI)                                    %The comparison process will be made for each symbol period
                                        %                 MeanOfData = mean(EoutI((kk-1)+SymLocI(1)));
                                        MeanOfData = mean(EoutI((kk-1)+NPPB/2));
                                        if MeanOfData > LevelDec3%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                                            DataOdd = [DataOdd 1];                                 %is 1
                                        else                                                       %If it is the lowest level the incoming data
                                            DataOdd = [DataOdd 0];                                 %is 0
                                        end
                                        MeanOfData = mean(EoutI((kk-1)+3*NPPB/PerDiv));
                                        if MeanOfData > LevelDec1%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                                            DataOddM = [DataOddM 1];                                 %is 1
                                        else                                                       %If it is the lowest level the incoming data
                                            DataOddM = [DataOddM 0];                                 %is 0
                                        end
                                        MeanOfData = mean(EoutI((kk-1)+NPPB/PerDiv));
                                        if MeanOfData > LevelDec2%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                                            DataOddL = [DataOddL 1];                                 %is 1
                                        else                                                       %If it is the lowest level the incoming data
                                            DataOddL = [DataOddL 0];                                 %is 0
                                        end
                                        MeanOfData = mean(EoutI((kk-1)+(NPPB/2)));
                                        if MeanOfData > LevDefDpqsk%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                                            DataOddU = [DataOddU 1];                                 %is 1
                                        else                                                       %If it is the lowest level the incoming data
                                            DataOddU = [DataOddU 0];                                 %is 0
                                        end
                                        MeanOfData = mean(EoutI((kk-1)+(NPPB/PerDiv)));
                                        if MeanOfData > LevDefDpqsk%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                                            DataOddU2 = [DataOddU2 1];                                 %is 1
                                        else                                                       %If it is the lowest level the incoming data
                                            DataOddU2 = [DataOddU2 0];                                 %is 0
                                        end
                                        MeanOfData = mean(EoutI((kk-1)+(3*NPPB/PerDiv)));
                                        if MeanOfData > LevDefDpqsk%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                                            DataOddU3 = [DataOddU3 1];                                 %is 1
                                        else                                                       %If it is the lowest level the incoming data
                                            DataOddU3 = [DataOddU3 0];                                 %is 0
                                        end
                                        if CalcS
                                            MeanOfData = mean(EoutI((kk-1)+NPPB/2));
                                            if MeanOfData > EyeOpenLowI + EyeOpenI/2%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                                                DataOddS = [DataOddS 1];                                 %is 1
                                            else                                                       %If it is the lowest level the incoming data
                                                DataOddS = [DataOddS 0];                                 %is 0
                                            end
                                        end
                                    end
                                    %The identical process just described above will also be used
                                    %to recover the data at the even positions.
                                    %
                                    %As this process is also statistical, first we reshape the
                                    %income vector to analyze all periods at the same time.
                                    %                 EyeSymMatQ = reshape(EoutQ(1+SyncPeriod*NPPB:end-SyncPeriod*...
                                    %                     NPPB),NPPB,NbDQPSK/2-2*SyncPeriod);
                                    %                 %Then we take the values that compose the decision level
                                    %                 %because they will mark the point of symmetry.
                                    %                 %
                                    %                 %Firstly it was set the interval in which the histogram will be
                                    %                 %build. It is based on the number of samples per bit period.
                                    %                 IntervalQ = linspace(min(EoutQ(1+SyncPeriod*NPPB:end-...
                                    %                     SyncPeriod*NPPB)),max(EoutQ(1+SyncPeriod*NPPB:end-...
                                    %                     SyncPeriod*NPPB)),2*NPPB);
                                    %                 %Therefore, the MATLAB hist function returns the number of
                                    %                 %occurrence of each interval.
                                    %                 EyeMaxQ = hist(EoutQ,IntervalQ);
                                    %                 EyeMaxauxQ = [0 EyeMaxI 0];                                    %Zeros are added at the EyeMax to auxiliate the finding peaks process
                                    %                 [~,EyeLocQ] = findpeaks(EyeMaxauxQ,'MinPeakDistance',NPPB/4,...
                                    %                     'SortStr','descend','NPeaks',2);%The peaks on the Eye profile will be the levels at the Eyes limit
                                    %                 ValToSeekQ = IntervalQ(EyeLocQ-1);
                                    %                 ValToSeekQ = sort(ValToSeekQ,'ascend');
                                    %                 OccuCountQ = zeros(1,size(EyeSymMatQ,1));                      %Auxiliar Variable for accounting.
                                    %                 for kk=1:size(EyeSymMatQ,1)                                    %For every sample within a symbol period
                                    %                     OccuCountQ(kk) = OccuCountQ(kk)+sum((EyeSymMatI(kk,:)>=...
                                    %                         LowSymPer*ValToSeekQ(1))&(EyeSymMatQ(kk,:)<=UpeSymPer*...
                                    %                         ValToSeekQ(1)));%Account all occurencies of the valeu 1
                                    %                     OccuCountQ(kk) = OccuCountQ(kk)+sum((EyeSymMatQ(kk,:)>=...
                                    %                         LowSymPer*ValToSeekQ(2))&(EyeSymMatQ(kk,:)<=UpeSymPer*...
                                    %                         ValToSeekQ(2)));%Account all occurencies of the valeu 2
                                    %                 end
                                    %                 %             [~,SymLocQ] = findpeaks(OccuCountQ,'SortStr','descend');       %The peak on the Eye profile will be the Symmetry level
                                    %                 [SymLocQ] = round(mean(find(ismember(OccuCountQ,max(...
                                    %                     OccuCountQ)))));
                                    %##############################################################
                                    %######################Important###############################
                                    %The ber results for the Data Even is not as good as the
                                    %results of the Data Odd. One possible reason is the decision
                                    %point of symmetry hence, for testing, we change SymLocQ to
                                    %SymLocI for evaluation of improvement. It is expected as bouth
                                    %signal will the same point of symetry will perform in an equal
                                    %way. If it is confirmed the creation of SymLocQ will be
                                    %erased.
                                    %##############################################################
                                    %##############################################################
                                    DataEven = [];                                                 %Initialization of the vector that will store the income data
                                    DataEvenM = [];                                                 %Initialization of the vector that will store the income data
                                    DataEvenL = [];                                                 %Initialization of the vector that will store the income data
                                    DataEvenU = [];                                                 %Initialization of the vector that will store the income data
                                    DataEvenU2 = [];                                                 %Initialization of the vector that will store the income data
                                    DataEvenU3 = [];                                                 %Initialization of the vector that will store the income data
                                    DataEvenS = [];                                                 %Initialization of the vector that will store the income data
                                    for kk=1:NPPB:length(EoutQ)                                    %The comparison process will be made for each symbol period
                                        %                 MeanOfData = mean(EoutQ((kk-1)+SymLocI(1)));
                                        MeanOfData = mean(EoutQ((kk-1)+NPPB/2));
                                        if MeanOfData > QLevelDec3%EyeOpenLowQ+EyeOpenQ/2                     %If it is the uper level the incoming data
                                            DataEven = [DataEven 1];                               %is 1
                                        else                                                       %If it is the lowest level the incoming data
                                            DataEven = [DataEven 0];                               %is 0
                                        end
                                        MeanOfData = mean(EoutQ((kk-1)+3*NPPB/PerDiv));
                                        if MeanOfData > QLevelDec1%EyeOpenLowQ+EyeOpenQ/2                     %If it is the uper level the incoming data
                                            DataEvenM = [DataEvenM 1];                               %is 1
                                        else                                                       %If it is the lowest level the incoming data
                                            DataEvenM = [DataEvenM 0];                               %is 0
                                        end
                                        MeanOfData = mean(EoutQ((kk-1)+NPPB/PerDiv));
                                        if MeanOfData > QLevelDec3%EyeOpenLowQ+EyeOpenQ/2                     %If it is the uper level the incoming data
                                            DataEvenL = [DataEvenL 1];                               %is 1
                                        else                                                       %If it is the lowest level the incoming data
                                            DataEvenL = [DataEvenL 0];                               %is 0
                                        end
                                        MeanOfData = mean(EoutQ((kk-1)+(NPPB/2)));
                                        if MeanOfData > LevDefDpqsk%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                                            DataEvenU = [DataEvenU 1];                                 %is 1
                                        else                                                       %If it is the lowest level the incoming data
                                            DataEvenU = [DataEvenU 0];                                 %is 0
                                        end
                                        MeanOfData = mean(EoutQ((kk-1)+(NPPB/PerDiv)));
                                        if MeanOfData > LevDefDpqsk%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                                            DataEvenU2 = [DataEvenU2 1];                                 %is 1
                                        else                                                       %If it is the lowest level the incoming data
                                            DataEvenU2 = [DataEvenU2 0];                                 %is 0
                                        end
                                        MeanOfData = mean(EoutQ((kk-1)+(3*NPPB/PerDiv)));
                                        if MeanOfData > LevDefDpqsk%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                                            DataEvenU3 = [DataEvenU3 1];                                 %is 1
                                        else                                                       %If it is the lowest level the incoming data
                                            DataEvenU3 = [DataEvenU3 0];                                 %is 0
                                        end
                                        if CalcS
                                            MeanOfData = mean(EoutI((kk-1)+NPPB/2));
                                            if MeanOfData > EyeOpenLowQ + EyeOpenQ/2%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                                                DataEvenS = [DataEvenS 1];                                 %is 1
                                            else                                                       %If it is the lowest level the incoming data
                                                DataEvenS = [DataEvenS 0];                                 %is 0
                                            end
                                        end
                                    end
                                    %%       Calculating the Bit Error Ratio (BER)
                                    %The final process here is to count the number of wrongdoings
                                    %of this whole process upon the transmited data for
                                    %quantitative analizes
                                    TxDataOdd  = TxDataMat(ThisCarr,1+SyncPeriod:end-SyncPeriod);
                                    TxDataEven = TxDataMat(ThisCarr+NumCarr,1+SyncPeriod:end-...
                                        SyncPeriod);
                                    if CalcS
                                        DataOddS  = DataOddS(1+SyncPeriod:end-SyncPeriod);
                                        DataEvenS = DataEvenS(1+SyncPeriod:end-SyncPeriod);
                                        BitErrOddS    = sum(xor(TxDataOdd,DataOddS));                      %Comparison between the Transmited and received and counting the differences
                                        BitErrEvenS   = sum(xor(TxDataEven,DataEvenS));                    %Comparison between the Transmited and received and counting the differences
                                        BerDQPSKS(CurrentTest,ThisCarr) = (BitErrOddS+BitErrEvenS)/...
                                            ((NbDQPSK)-(2*SyncPeriod));%Calculating the ration of wrong bits and the total number of bits transmited
                                    end
                                    DataOdd  = DataOdd(1+SyncPeriod:end-SyncPeriod);
                                    DataEven = DataEven(1+SyncPeriod:end-SyncPeriod);
                                    DataOddM  = DataOddM(1+SyncPeriod:end-SyncPeriod);
                                    DataEvenM = DataEvenM(1+SyncPeriod:end-SyncPeriod);
                                    DataOddL  = DataOddL(1+SyncPeriod:end-SyncPeriod);
                                    DataEvenL = DataEvenL(1+SyncPeriod:end-SyncPeriod);
                                    DataOddU  = DataOddU(1+SyncPeriod:end-SyncPeriod);
                                    DataEvenU = DataEvenU(1+SyncPeriod:end-SyncPeriod);
                                    DataOddU2  = DataOddU2(1+SyncPeriod:end-SyncPeriod);
                                    DataEvenU2 = DataEvenU2(1+SyncPeriod:end-SyncPeriod);
                                    DataOddU3  = DataOddU3(1+SyncPeriod:end-SyncPeriod);
                                    DataEvenU3 = DataEvenU3(1+SyncPeriod:end-SyncPeriod);
                                    AberLevAuxI(4) = 0;
                                    ValsLevAuxI(4) = 0.01;
                                    AberLevAuxQ(4) = 0;
                                    ValsLevAuxQ(4) = 0.01;
                                    AberLevAuxI(5) = 0;
                                    ValsLevAuxI(5) = 0.01;
                                    AberLevAuxQ(5) = 0;
                                    ValsLevAuxQ(5) = 0.01;
                                    AberLevAuxI(6) = 0;
                                    ValsLevAuxI(6) = 0.01;
                                    AberLevAuxQ(6) = 0;
                                    ValsLevAuxQ(6) = 0.01;
                                    
                                    BitErrOdd(1)  = sum(xor(TxDataOdd,DataOdd));                      %Comparison between the Transmited and received and counting the differences
                                    BitErrEven(1) = sum(xor(TxDataEven,DataEven));                    %Comparison between the Transmited and received and counting the differences
                                    BitErrOdd(2)  = sum(xor(TxDataOdd,DataOddM));                      %Comparison between the Transmited and received and counting the differences
                                    BitErrEven(2) = sum(xor(TxDataEven,DataEvenM));                    %Comparison between the Transmited and received and counting the differences
                                    BitErrOdd(3)  = sum(xor(TxDataOdd,DataOddL));                      %Comparison between the Transmited and received and counting the differences
                                    BitErrEven(3) = sum(xor(TxDataEven,DataEvenL));                    %Comparison between the Transmited and received and counting the differences
                                    BitErrOdd(4)  = sum(xor(TxDataOdd,DataOddU));                      %Comparison between the Transmited and received and counting the differences
                                    BitErrEven(4) = sum(xor(TxDataEven,DataEvenU));                    %Comparison between the Transmited and received and counting the differences
                                    BitErrOdd(5)  = sum(xor(TxDataOdd,DataOddU2));                      %Comparison between the Transmited and received and counting the differences
                                    BitErrEven(5) = sum(xor(TxDataEven,DataEvenU2));                    %Comparison between the Transmited and received and counting the differences
                                    BitErrOdd(6)  = sum(xor(TxDataOdd,DataOddU3));                      %Comparison between the Transmited and received and counting the differences
                                    BitErrEven(6) = sum(xor(TxDataEven,DataEvenU3));                    %Comparison between the Transmited and received and counting the differences
                                    BerDQPSK(CurrentTest,ThisCarr) = (min(BitErrOdd)+min(BitErrEven))/...
                                        ((NbDQPSK)-(2*SyncPeriod));%Calculating the ration of wrong bits and the total number of bits transmited
                                    
                                    AberIAux = max(AberLevAuxI(BitErrOdd<=min(BitErrOdd)));
                                    AberQAux = max(AberLevAuxQ(BitErrEven<=min(BitErrEven)));
                                    AberLevI(CurrentTest,ThisCarr)  = AberIAux;
                                    ValsLevI(CurrentTest,ThisCarr)  = max(ValsLevAuxI(AberLevAuxI==AberIAux));
                                    
                                    AberLevQ(CurrentTest,ThisCarr)  = AberQAux;
                                    ValsLevQ(CurrentTest,ThisCarr)  = max(ValsLevAuxQ(AberLevAuxQ==AberQAux));
                                    
                                    [~,LocAux] = max(AberLevAuxI==AberIAux);
                                    EvmDB(CurrentTest,ThisCarr)   = EvmDBA(LocAux);
                                    EvmPer(CurrentTest,ThisCarr)  = EvmPerA(LocAux);
                                    EvmRms(CurrentTest,ThisCarr)  = EvmRmsA(LocAux);
                                    EvmDBJ(CurrentTest,ThisCarr)  = EvmDBJA(LocAux);
                                    EvmPerJ(CurrentTest,ThisCarr) = EvmPerJA(LocAux);
                                    EvmRmsJ(CurrentTest,ThisCarr) = EvmRmsJA(LocAux);
                                    
                                    RxSymbAmos(CurrentTest,ObsCarrPos==ThisCarr,:) = EvmMatRecA(LocAux,:);%RxSymbAmos = [];
                                    EvmMatRec(ObsCarrPos==ThisCarr,:) = EvmMatRecA(LocAux,:);                       %Taking just the middle samples as references
                                    %% Ploting the result for qualitative analizes
                                    %PrintInfo(Ploting*29,TxDataOdd,DataOdd);
                                    %PrintInfo(Ploting*30,TxDataEven,DataEven);
                                    %%
                                    %             berpos1 = 1:2:size(BerDQPSK,2);
                                    %             berpos2 = 2:2:size(BerDQPSK,2);
                                    %             b = [BerDQPSK(1,berpos1);BerDQPSK(1,berpos2)]
                                    %             a=6;
                                    close all;
                                end
                            end
                        case '4PAM'
                            
                            %%              Receiver 4PAM
                            for ThisCarr=InitCarr:CarrPass:NumCarr                                           %For each carrier the same process of reception will be used.
                                if CarrUsedUp(ThisCarr)
                                    %%  Reception Process: Handling the income optical field
                                    %At the first step the income field will be selected from the
                                    %optical FFT Output with the help of the maping vector
                                    %(previously described).
                                    EoutAux = EoutAux1(VetThisCarr==ObsCarrUsed(ThisCarr),:);
                                    if AddingNoiseP
                                        %                 figure;plot(f,20*log10(abs(fftshift(fft(EoutAux)./length(EoutAux)))));
                                        %                 [~,sigpower] = MeasPower(EoutAux);
                                        [p,l] = findpeaks(20*log10(abs(fftshift(fft(EoutAux)./length(EoutAux)))),'SortStr','descend','NPeaks',1);
                                        EoutAuxA = fftshift(fft(EoutAux));
                                        EoutAuxA(f==f(l)) = EoutAuxA(f==f(l-3));
                                        %                 [testfilt,~] = FiltroGaussiano(f,12.5e9,f(l),5);           %Creating filter for selection of the received signal
                                        %                 testfilt = fftshift(testfilt);                                   %Shifting the filter for matching the received signal
                                        %                 EoutAux = ifft(fftshift(EoutAuxA.*testfilt));
                                        EoutAuxA = ifft(fftshift(EoutAuxA));
                                        %                 figure;plot(f,20*log10(abs(fftshift(fft(EoutAuxA)./length(EoutAuxA)))));
                                        %                 [~,sigpower] = MeasPower(EoutAuxA);
                                        [~,sigpower] = MeasPower(EoutAuxA);
                                        EoutAux      = awgn(EoutAux,osnrp,sigpower-30);
                                        [~,sigpower] = MeasPower(EoutAux);
                                    end
                                    %%            Fiber Time Delay Compensation
                                    switch Medium
                                        case 'Fiber'
                                            EoutAux = ifft(fft(EoutAux).*exp(1j*2*pi*f*(...
                                                FiberDelay((ThisCarr))*Ta)));
                                        otherwise
                                    end
                                    %The current incoming signal is them converted from the optical
                                    %domain to the eletrical domain with the help of an photo
                                    %detector.
                                    Ix = EoutAux.*conj(EoutAux);
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
                                    if RecFilBanPas
                                        Ix = ifft(fft(Ix).*BitFilt);
                                        Ix(1:4*(2*NumAmosCP+NPPB))         = Ix(6*(2*NumAmosCP+NPPB));
                                        Ix(end+1-4*(2*NumAmosCP+NPPB):end) = Ix(6*(2*NumAmosCP+NPPB));
                                        Ix = Ix - min(Ix);                                         %Removing the DC component from them Eletrical signal received
                                        Ix = Ix./max(abs(Ix));                                     %Normalizing the eletrical signal (amplifying what is needed)
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
                                    if ReceptorNoise                                               %Verify whether the noise is to be added or not
                                        [~,SigPower] = MeasPower(Ix);                              %it own received signal or
                                        SigPower2 = SigPower-30;%10*log10(SigPower);
                                        %SNR = CarSNR - 10*log10(2) + 10*log10(Nsamp) + 10*log10(1/2);
                                        SNR2 = CarSNR + 10*log10(2) - 10*log10(Nsamp) - 10*log10(1/2) + 10*log10(10^0.6);
                                        
                                        if TimeSys
                                            SNR = CarSNR + 10*log10(2) - 10*log10(0.25*T*BWD2);
                                        else
                                            SNR = CarSNR + 10*log10(2) - 10*log10(0.25*t(2*NumAmosCP+NPPB)*BWD2);
                                        end
                                        
                                        Ix = awgn(Ix,SNR,SigPower2);
                                    end
                                    
                                    if ~RecFilBanPas
                                        Ix = ifft(fft(Ix).*BitFilt);
                                        Ix(1:4*(2*NumAmosCP+NPPB))         = Ix(6*(2*NumAmosCP+NPPB));
                                        Ix(end+1-4*(2*NumAmosCP+NPPB):end) = Ix(6*(2*NumAmosCP+NPPB));
                                        Ix = Ix - min(Ix);                                         %Removing the DC component from them Eletrical signal received
                                        
                                        Ix = Ix./max(abs(Ix));                                     %Normalizing the eletrical signal (amplifying what is needed)
                                    end
                                    
                                    VetElecPower(CurrentTest,ThisCarr)= MeasPower(Ix);
                                    EoutAuxF = 20*log10(abs(fftshift(fft(EoutAux)./length(...
                                        EoutAux))));
                                    [VetOptiPower(CurrentTest,ThisCarr),~]= findpeaks(EoutAuxF,...
                                        'SortStr','descend','NPeaks',1);
                                    
                                    %% Ploting the result for qualitative analizes
                                    %             PrintInfo(Ploting*31,f,20.*log10(abs(fftshift(fft(Ix)./...
                                    %                                                             length(Ix)))));
                                    %                             PrintInfo(Ploting*32,t,Ix);
                                    %                             PrintInfo(Ploting*33,Ix,T,NPPB);
                                    %             PrintInfo(Ploting*34,t(IniSyncPos:SyncPos),Ix(IniSyncPos:...
                                    %                                     SyncPos),SyncSymb(IniSyncPos:SyncPos));
                                    %%        synchronizing
                                    %For the reception process work properly, it is needed to
                                    %sincronized the recieved signal or the sampling process will
                                    %not work.
                                    AuxSync = (Ix(IniSyncPos:SyncPos));                            %Selecting the sync-word within the received signal
                                    AuxSync1 = AuxSync;
                                    %This synchronization process is based on the mean expected
                                    %value. That means, the information mean value should be within
                                    %one period of symbol. Thus, the mean value of the received
                                    %signal is acquired and compare of the known sync-word to
                                    %verify if this mean value is at the right possition
                                    AuxSync2                           = SyncSymb(IniSyncPos:...
                                        SyncPos);%Selecting the sync-word within the known signal
                                    AuxSync2(AuxSync2>=mean(AuxSync2)) = 1;                        %Adding a flag to the first sample of the known mean value
                                    AuxSync2(AuxSync2<mean(AuxSync2))  = -1;                       %All the others samples at set to the lowest level
                                    AuxSync1(AuxSync1<0)               = 0;                        %To keep the mean value above zero anything under is neglected
                                    AuxSync1(AuxSync1>=mean(AuxSync1)) = 1;                        %Adding a flag to the first sample of the received mean value
                                    if PrintinEye
                                        figure; hold all;
                                        plot(t(IniSyncPos:SyncPos),SyncSymb(IniSyncPos:SyncPos),t(IniSyncPos:SyncPos),Ix(IniSyncPos:SyncPos),t(IniSyncPos:SyncPos),AuxSync1,t(IniSyncPos:SyncPos),linspace(mean(AuxSync1),mean(AuxSync1),length(t(IniSyncPos:SyncPos))));
                                        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                    end
                                    AuxSync1(AuxSync1<mean(AuxSync1))  = -1;                       %All the others samples at set to the lowest level
                                    
                                    PosToSyn  = find(ismember(AuxSync1,1));                        %Finding where is the location of the first sample to synchronize
                                    PosSyn = find(ismember(AuxSync2,1));                           %Finding where is the location of the first sample to synchronize
                                    
                                    AuxSyncCorr = PosToSyn(round(end/2))-PosSyn(round(end/2));
                                    
                                    sn=((1/2)*(((1/2)^MaxNumStag)-1))/((1/2)-1);
                                    AuxSyncCorr = round((sn/(2/2^IfftOrSum))*T/Ta);
                                    %The difference between the PossitionTosynchronize and
                                    %Possitionsynchronized will be used to correct the time
                                    %shifting on the transmition and reception process.
                                    if AuxSyncCorr>=0%If the difference is positive, left-shift...
                                        %                 Ix = ifft(fft(Ix).*exp(1j*2*pi*f*AuxSyncCorr));          %Shift based on time change
                                        Ix = [Ix(AuxSyncCorr+1:end) Ix(1:AuxSyncCorr)];            %Shift based on sampling sliding
                                    else%... but if the difference is negative, right-shift
                                        %                 Ix = ifft(fft(Ix).*exp(-1j*2*pi*f*AuxSyncCorr));         %Shift based on time change
                                        Ix = [Ix(end+AuxSyncCorr+1:end) Ix(1:end + AuxSyncCorr)];  %Shift based on sampling sliding
                                    end
                                    
                                    if PrintinEye
                                        figure; hold all;
                                        plot(t(IniSyncPos:SyncPos),SyncSymb(IniSyncPos:SyncPos),t(IniSyncPos:SyncPos),Ix(IniSyncPos:SyncPos));
                                        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                    end
                                    
                                    if SencondAdjust
                                        %For some reason that we could not understand sometimes the
                                        %time (sampling) sliding of the signal is not equal
                                        %throught the data stream. Thus, the second part of the
                                        %synchronism process will be turn ON or OFF according to
                                        %the user's will.
                                        AuxSyncEnd     = (Ix(end-SyncPos+1:end-IniSyncPos-1));     %Selecting the sync-word within the received signal
                                        SyncSymbEndAux = SyncSymbEnd(end-SyncPos+1:end-...
                                            IniSyncPos-1);%Selecting the sync-word within the known signal
                                        AuxSyncEnd1 = AuxSyncEnd;
                                        AuxSyncEnd1(AuxSyncEnd1<0) = 0;                            %To keep the mean value above zero anything under is neglected
                                        AuxSyncEnd1(AuxSyncEnd1>=mean(AuxSyncEnd1)) = 1;           %Adding a flag to the first sample of the received mean value
                                        AuxSyncEnd1(AuxSyncEnd1<mean(AuxSyncEnd1)) = -1;           %All the others samples at set to the lowest level
                                        AuxSyncEnd2 = SyncSymbEndAux;
                                        AuxSyncEnd2(AuxSyncEnd2>=mean(AuxSyncEnd2)) = 1;           %Adding a flag to the first sample of the known mean value
                                        AuxSyncEnd2(AuxSyncEnd2<mean(AuxSyncEnd2)) = -1;           %All the others samples at set to the lowest level
                                        
                                        
                                        PosToSynEnd  = find(ismember(AuxSyncEnd1,1));              %Finding where is the location of the first sample to synchronize
                                        PosSynEnd = find(ismember(AuxSyncEnd2,1));                 %Finding where is the location of the first sample to synchronize
                                        
                                        AuxSyncCorrEnd = PosToSynEnd(round(end/2))-PosSynEnd(...
                                            round(end/2));
                                        
                                        %The difference between the PossitionTosynchronize and
                                        %Possitionsynchronized will be used to correct the time
                                        %shifting on the transmition and reception process.
                                        if AuxSyncCorrEnd>=0%If possitive difference, left-shift...
                                            %                 Ix = ifft(fft(Ix).*exp(1j*2*pi*f*AuxSyncCorrEnd));   %Shift based on time change
                                            Ix = [Ix(AuxSyncCorrEnd+1:end) Ix(1:AuxSyncCorrEnd)];  %Shift based on sampling sliding
                                        else%... but if the difference is negative, right-shift
                                            %                 Ix = ifft(fft(Ix).*exp(-1j*2*pi*f*AuxSyncCorrEnd));  %Shift based on time change
                                            Ix = [Ix(end+AuxSyncCorrEnd+1:end) Ix(1:end + ...
                                                AuxSyncCorrEnd)];%Shift based on sampling sliding
                                        end
                                    end
                                    
                                    %%          Ploting the result for qualitative analizes
                                    %             PrintInfo(Ploting*35,t(end-SyncPos+1:end-IniSyncPos+1),Ix(...
                                    %                 end-SyncPos+1:end-IniSyncPos+1),SyncSymb(end-SyncPos+1:...
                                    %                                                         end-IniSyncPos+1));
                                    %% Removing CP
                                    %                 if ThisCarr==126
                                    %                     EyeToPlot(CurrentTest,1:length(Ix)) = Ix;
                                    %                     save(['savingforplotingeye' num2str(VetSnr)],'EyeToPlot','VetSnr','IfftOrSum','UsedModula','T','NPPB');
                                    %                 end
                                    if AddCP
                                        IxAux = Ix(1:end - StuffSampels);
                                        IxAux = reshape(IxAux,(2*NumAmosCP+NPPB),Nb4Pam/2);
                                        IxAux = IxAux(1+NumAmosCP:end-NumAmosCP,:);
                                        IxAux = reshape(IxAux,1,NPPB*Nb4Pam/2);
                                        Ix    = IxAux;
                                    end
                                    %% Taking the sampling the EVM meassurement
                                    PosAuxEout = NPPB/2:NPPB:length(Ix);%Varriable respossible to take just the samples at the middle of the symbol
                                    IxAux      = Ix(PosAuxEout);                           %Normalizing the reference
                                    RxSymbAmos(CurrentTest,ObsCarrPos==ThisCarr,:) = IxAux;%RxSymbAmos = [];
                                    EvmMatRec(ObsCarrPos==ThisCarr,:) = IxAux;                       %Taking just the middle samples as references
                                    [EvmDB(CurrentTest,ThisCarr), EvmPer(CurrentTest,ThisCarr), EvmRms(CurrentTest,ThisCarr) ] = EvmCalc( EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux );
                                    [EvmDBJ(CurrentTest,ThisCarr),EvmPerJ(CurrentTest,ThisCarr),EvmRmsJ(CurrentTest,ThisCarr)] = evm1(4,'pam',EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux);
                                    a=0;
                                    %                 if ThisCarr==2
                                    %                     figure;hold all;grid on;
                                    %                     plot(IxAux,zeros(1,length(IxAux)),'rx');
                                    %                     set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                    %                     a=0;
                                    %                 end
                                    %                 %%  Measuring the EVM
                                    %%         Finding Decission Levels
                                    %The process for decoding the income signal will be based on
                                    %eletronic comparators. Inasmuch as the right decission level
                                    %must be acquired for accurately decide, within a symbol
                                    %periode, what that current leavel means (ones or zeros).
                                    %
                                    %The process hereafter of chosing the  decission levels is not
                                    %deterministic rather it is a statistic process. The main idea
                                    %is to take the decission level from the histogram generated
                                    %from the income signal stream.
                                    %
                                    %Firstly it was set the interval in which the histogram will be
                                    %build. It is based on the number of samples per bit period.
                                    %             IxToSeek = reshape(Ix(1+SyncPeriod*NPPB:end-SyncPeriod*NPPB...
                                    %                                        ),NPPB,(Nb4Pam/2) - (2*SyncPeriod));
                                    Interval = linspace(min(Ix(1+SyncPeriod:end-SyncPeriod)),max...
                                        (Ix(1+SyncPeriod:end-SyncPeriod)),IntervalStep);
                                    %Therefore, the MATLAB hist function returns the number of
                                    %occurrence of each interval.
                                    EyeMax = hist(Ix,Interval);
                                    EyeMaxaux = [0 EyeMax 0];                                      %Zeros are added at the EyeMax to auxiliate the finding peaks process
                                    [~,EyeLoc] = findpeaks(EyeMaxaux,'MinPeakDistance',MinDist,...
                                        'SortStr','descend','NPeaks',4);%,'MinPeakHeight',mean(EyeMax));%The peaks on the Eye profile will be the levels at the Eyes limit
                                    if length(EyeLoc)<4
                                        [~,EyeLoc] = findpeaks(EyeMaxaux,'MinPeakDistance',...
                                            MinDist*0.8,'SortStr','descend','NPeaks',4);%,'MinPeakHeight',mean(EyeMax)/4);%The peaks on the Eye profile will be the levels at the Eyes limit
                                    end
                                    %This variable will brings the eye profile of the input signal
                                    %that will be used to generate the decision level. Basicaly the
                                    %decission level will be the minimal value of the currente eye
                                    %under evaluation plus the half of the its eye opening. The
                                    %following ilustration better describe this process.
                                    %
                                    %Eye Limits:   + (1/2 * Eye Opening:)  =    Comparisson Limits:
                                    %
                                    %UperLevel  ______     ______________     _____
                                    %                 \   /      |       \   /
                                    %                  \ /       |        \ /
                                    %                   \ Half Eye Opening /     Decission Level 3
                                    %                  / \       |        / \
                                    %LowerLevel 3_____/   \______|_______/   \_____
                                    %                 \   /      |       \   /
                                    %                  \ /       |        \ /
                                    %                   \ Half Eye Opening /     Decission Level 2
                                    %                  / \       |        / \
                                    %LowerLevel 2_____/   \______|_______/   \_____
                                    %                 \   /      |       \   /
                                    %                  \ /       |        \ /
                                    %                   \ Half Eye Opening /     Decission Level 1
                                    %                  / \       |        / \
                                    %LowerLevel 1_____/   \______|_______/   \_____
                                    %
                                    %%           Ploting for Qualitative Analizes
                                    %             PrintInfo(Ploting*37,Interval,EyeMax);
                                    %%         Finding Decission Levels
                                    %It is not always possible to totaly recover the signal.
                                    %Depending on the configuration of the transmition and
                                    %reception system the eye diagram may be nonexistent. Which
                                    %means, there will not be a profile to be found therefore the
                                    %EyeLoc will not return the correct location. Inasmuch as the
                                    %detection process to works limts will be set accordling to the
                                    %amplitude of the received signal.
                                    if length(EyeLoc)<4%If it was not able to find the eye profile.
                                        EyeLoc = [2 3 4 5];
                                        Levels = [0 0 0.35 0.55 0.85];
                                        %                 Levels = sort(Levels);
                                        %                 LevelDec1 = mean(Ix) - 2*mean(Ix)/3 ;
                                        %                 LevelDec2 = mean(Ix) ;
                                        %                 LevelDec3 = mean(Ix) + 2*mean(Ix)/3 ;
                                    else%Whereas, if there is an profile the decission can be found
                                        Levels = [Interval(EyeLoc(1)-1) Interval(EyeLoc(2)-1) ...
                                            Interval(EyeLoc(3)-1) Interval(EyeLoc(4)-1)];
                                        Levels = sort(Levels);
                                        %                 LevelDec1 = Levels(1) + (Levels(2)-Levels(1))/2 ;
                                        %                 LevelDec2 = Levels(2) + (Levels(3)-Levels(2))/2 ;
                                        %                 LevelDec3 = Levels(3) + (Levels(4)-Levels(3))/2 ;
                                    end
                                    %##############################################################
                                    if CalcS
                                        limiar1 = (1/2)*abs(max(Ix)-min(Ix))/3;
                                        limiar2 = ((1/2)*abs(max(Ix)-min(Ix))/3) + abs(max(Ix)-min(Ix))/3;
                                        limiar3 = ((1/2)*abs(max(Ix)-min(Ix))/3) + 2*(abs(max(Ix)-min(Ix))/3);
                                        limiarPos1 = Interval>limiar1;
                                        limiarPos2 = (Interval<=limiar1)&(Interval>limiar2);
                                        limiarPos3 = (Interval<=limiar2)&(Interval>limiar3);
                                        limiarPos4 = Interval<=limiar3;
                                        
                                        EyeSymMa1 = reshape(Ix(1+SyncPeriod*NPPB:end-SyncPeriod*NPPB...
                                            ),NPPB,(Nb4Pam/2) - (2*SyncPeriod));
                                        for kk = 1:size(EyeSymMa1,1)
                                            EyeHi1    = find((EyeSymMa1(kk,:)<Levels(4))&(EyeSymMa1(kk,:)>limiar1));
                                            EyeHi2    = find((EyeSymMa1(kk,:)<Levels(3))&(EyeSymMa1(kk,:)>limiar2));
                                            EyeHi3    = find((EyeSymMa1(kk,:)<Levels(2))&(EyeSymMa1(kk,:)>limiar3));
                                            
                                            EyeLo1    = find((EyeSymMa1(kk,:)<limiar1)&(EyeSymMa1(kk,:)>Levels(3)));
                                            EyeLo2    = find((EyeSymMa1(kk,:)<limiar2)&(EyeSymMa1(kk,:)>Levels(2)));
                                            EyeLo3    = find((EyeSymMa1(kk,:)<limiar3)&(EyeSymMa1(kk,:)>Levels(1)));
                                            
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
                                    end
                                    %             for ii = 1:N1
                                    %               media(ii) = mean(olho(ii,:));
                                    %               nhi       = find(olho(ii,:) > media(ii));
                                    %               nlow      = find(olho(ii,:) < media(ii));
                                    %               hi(ii)    = mean(olho(ii,nhi));
                                    %               low(ii)   =  mean(olho(ii,nlow));
                                    %             %
                                    %               Lhi(ii)   = mean(hi)  - std(olho(ii,nhi));
                                    %               Llow(ii)  = mean(low) + std(olho(ii,nlow));
                                    %               Ab(ii)    = Lhi(ii) - Llow(ii);
                                    %             end
                                    %             Abertura   = mean(Ab);
                                    %             a=1;
                                    %##############################################################
                                    %For accurately read the income signal is necessary more than
                                    %just have the decision level, it is also needed to know where
                                    %it occurs. For instance, if we decide to use just to take the
                                    %decision level and measure the mean value across and a portion
                                    %of the time period, which portion should we take? It would be
                                    %logical to take the central portion, but the bit may be
                                    %deformed in such way that the information gets concentrated
                                    %out of the middle part, which means the symbol is not
                                    %symmetric. The symmetry information can be acquired from the
                                    %eye diagram by measuring the longitudinal opening. The
                                    %following sketch better describes this process:
                                    %
                                    %                   Point of Symmetry
                                    %         _____     _______|_______     _____
                                    %              \   /               \   /
                                    %               \ /  Longitudinal   \ /
                                    %                \ ----------------- /
                                    %               / \    Opening      / \
                                    %         _____/   \_______________/   \_____
                                    %
                                    %With those two pieces of information, decision level and point
                                    %of symmetry, we have the X, Y coordinates for the centre of
                                    %the Eye Diagram. Therefore, as long as there is an opening on
                                    %it it will be possible to recover the transmitted information
                                    %without error, theoretically.
                                    
                                    %As this process is also statistical, first we reshape the
                                    %income vector to analyze all periods at the same time.
                                    %             EyeSymMat = reshape(Ix(1+SyncPeriod*NPPB:end-SyncPeriod*NPPB...
                                    %                 ),NPPB,(Nb4Pam/2) - (2*SyncPeriod));
                                    %             EyeSymMat = EyeSymMat(1+NPPB/4:end-NPPB/4,:);
                                    %Then we take the values that compose the decision level
                                    %because they will mark the point of symmetry.
                                    %From the location of the max values that occured, which means
                                    %the uper and lower level of the eye diagram it needs to take
                                    %the actual value that those occurences represent that is
                                    %withing the the Interval variable.
                                    %             ValToSeek = Interval(EyeLoc-1);
                                    %The number of ocurrences is a statical measure therefore one
                                    %does not have control which interval will have the highest
                                    %peak, thus it is important to ordenate the values to be seek
                                    %from the lower part of the eye diagram to the uper part of the
                                    %eye diagram.
                                    %             ValToSeek = sort(ValToSeek,'ascend');
                                    %             OccuCount = zeros(1,size(EyeSymMat,1));                        %Auxiliar Variable for accounting.
                                    %             for kk=1:size(EyeSymMat,1)                                     %For every sample within a symbol period
                                    %                 if (kk>=NPPB/4)&&(kk<=NPPB-NPPB/4)
                                    %                     OccuCount(kk) = OccuCount(kk) + sum((EyeSymMat(kk,:)...
                                    %                         >=LowSymPer*ValToSeek(1))&(EyeSymMat(kk,:)<=...
                                    %                         UpeSymPer*ValToSeek(1)));%Account all occurencies of the valeu 1
                                    %                     OccuCount(kk) = OccuCount(kk) + sum((EyeSymMat(kk,:)...
                                    %                         >=LowSymPer*ValToSeek(2))&(EyeSymMat(kk,:)<=...
                                    %                         UpeSymPer*ValToSeek(2)));%Account all occurencies of the valeu 2
                                    %                     OccuCount(kk) = OccuCount(kk) + sum((EyeSymMat(kk,:)...
                                    %                         >=LowSymPer*ValToSeek(3))&(EyeSymMat(kk,:)<=...
                                    %                         UpeSymPer*ValToSeek(3)));%Account all occurencies of the valeu 3
                                    %                     OccuCount(kk) = OccuCount(kk) + sum((EyeSymMat(kk,:)...
                                    %                         >=LowSymPer*ValToSeek(4))&(EyeSymMat(kk,:)<=...
                                    %                         UpeSymPer*ValToSeek(4)));%Account all occurencies of the valeu 4
                                    %                 end
                                    %             end
                                    
                                    %The point of symmetry of the eye diagram will be where the
                                    %maximum number of occurrences were measured inasmuch as those
                                    %are the points where all the bits go to the center of the
                                    %symbol.  From the maximum number of occurrences, it can happen
                                    %for more than one sample within one symbol period, in the best
                                    %case, all samples would have the same accounting as it is
                                    %shown the ilustration above hence the symmetry will be at the
                                    %middle sample of this group of maximum occurrences. This value
                                    %can be found at the highest peak founded.
                                    %             [~,SymLoc] = findpeaks(OccuCount,'SortStr','descend');       %The peak on the Eye profile will be the Symmetry level
                                    
                                    %             PosIx = SymLoc(1):NPPB:length(Ix);
                                    
                                    
                                    
                                    
                                    
                                    %To find the Eye opening this script will concider just the
                                    %samples at the center of a symbol. Thus, it is based on the
                                    %fact that the syncronism process was done succefully.
                                    %##############################################################
                                    PosIx = NPPB/2:NPPB:length(Ix);                                %Possition of the central samples - but the number of samples per symbol is even ????
                                    IxAux = Ix(PosIx);                                             %From the main received signal just a few samples are taken for further evaluation
                                    n     = 100;                                                   %The number of boxes to be filled up on the histogram process
                                    %At this process each eye will be evaluated separetelly. It is
                                    %important to remember that the PAM4 has 4 levels, which means
                                    %three level of decissions that are exaclty the center of the
                                    %eye diagram.
                                    IxAuxAB = IxAux((IxAux<=Levels(4))&(IxAux>=Levels(3)));        %Taking just those values relative to the uper eye
                                    InterAB = linspace(Levels(3),Levels(4),n);                     %Building the histogram boxes
                                    EyeAB = hist(IxAuxAB,InterAB);                                 %filling up the boxes with samples that fit on them.
                                    %The same process described for the uper level will be done at
                                    %the middle and lower eyes levels.
                                    IxAuxCD = IxAux((IxAux<=Levels(3))&(IxAux>=Levels(2)));
                                    InterCD = linspace(Levels(2),Levels(3),n);%NPPB*2^n);
                                    EyeCD = hist(IxAuxCD,InterCD);
                                    
                                    IxAuxEF = IxAux((IxAux<=Levels(2))&(IxAux>=Levels(1)));
                                    InterEF = linspace(Levels(1),Levels(2),n);%NPPB*2^n);
                                    EyeEF = hist(IxAuxEF,InterEF);
                                    
                                    %What we are looking for here are not where exist the
                                    %occurrences of values. The eye is where there is not samples,
                                    %this means, the decission level is a reagion between actual
                                    %levels thus this reagion does not contain any signal.
                                    EyeAB = ~EyeAB;                                                %Changing zeros to one - Zeros compose the eye region
                                    EyeCD = ~EyeCD;
                                    EyeEF = ~EyeEF;
                                    %Starting variables that will be used to identify the eye
                                    %diagram.
                                    CountAB   = 1;
                                    CountCD   = 1;
                                    CountEF   = 1;
                                    
                                    SeqOnesAB = 0;
                                    SeqOnesCD = 0;
                                    SeqOnesEF = 0;
                                    
                                    SeqFinAB  = 0;
                                    SeqFinCD  = 0;
                                    SeqFinEF  = 0;
                                    
                                    SeqIniAB  = 1;
                                    SeqIniCD  = 1;
                                    SeqIniEF  = 1;
                                    %The for loop will take account of every box with ones. It is
                                    %important to take note that the not operator was used in this
                                    %vector, therefore ones means zeros (the eye diagram -
                                    %possibly) and zeros means values abouve zeroa (not the eye).
                                    for kk=1:length(EyeAB)                                         %For every box
                                        if EyeAB(kk)                                               %if it contains "1"
                                            SeqOnesAB(SeqIniAB)=CountAB;                           %count this element as part of a consecutive sequency
                                            CountAB = CountAB + 1;                                 %adds one to the counter of consecutive elements "1"
                                            if kk==length(EyeAB)                                   %if the current box is the last box we got to an end
                                                SeqFinAB(SeqIniAB) = kk;                           %The final sequency element is equal to its possition (kk)
                                            end
                                        else                                                       %else if the current box contains "0"
                                            SeqFinAB(SeqIniAB) = kk-1;                             %the previous element was the last element of a consecutive sequency of ones
                                            SeqIniAB = SeqIniAB + 1;                               %adds one to the counter of consecutive number (take in account how many sequencies there are)
                                            CountAB = 1;                                           %reset the counter of consecutive elements to mark the start of a new consecutive sequence
                                        end
                                        
                                        if EyeCD(kk)
                                            SeqOnesCD(SeqIniCD)=CountCD;
                                            CountCD = CountCD + 1;
                                            if kk==length(EyeCD)
                                                SeqFinCD(SeqIniCD) = kk;
                                            end
                                        else
                                            SeqFinCD(SeqIniCD) = kk-1;
                                            SeqIniCD = SeqIniCD + 1;
                                            CountCD = 1;
                                        end
                                        
                                        if EyeEF(kk)
                                            SeqOnesEF(SeqIniEF)=CountEF;
                                            CountEF = CountEF + 1;
                                            if kk==length(EyeEF)
                                                SeqFinEF(SeqIniEF) = kk;
                                            end
                                        else
                                            SeqFinEF(SeqIniEF) = kk-1;
                                            SeqIniEF = SeqIniEF + 1;
                                            CountEF = 1;
                                        end
                                    end
                                    
                                    
                                    %If the eye is open, which means there is a clear difference
                                    %between adjacent levels, the eye diagram will be the longest
                                    %sequence of ones.
                                    [MaxValAB,LocMaxAB]=max(SeqOnesAB);
                                    if LocMaxAB<1 || MaxValAB<2                                    %if any sequency was found or there is just one sequency it is a error thus
                                        LevDec3 = 0.6953;                                            %the decission level will be by default 0.67. Also, the other variables will
                                        LocMaxAB = 1;                                              %will be set with values to not cause errors in the future
                                        SeqFinAB(1)=2;
                                        MaxValAB = 0;
                                        InterAB(1)=LevDec3;
                                    else                                                           %if a sequency was found the middle of the eye will be the middle of the sequency
                                        if (SeqFinAB(LocMaxAB)-MaxValAB/2)<1                       %if for some reason the final element of the sequency minus the half of its
                                            LevDec3 = 0.6953;                                         %length results in a negative value, something went very wrong, and by
                                        else                                                       %default it will be set to 0.7
                                            LevDec3 = InterAB(round(SeqFinAB(LocMaxAB)-MaxValAB/...
                                                2));%Otherwise, the decission level is the middle of the sequency
                                        end
                                    end
                                    [MaxValCD,LocMaxCD]=max(SeqOnesCD);
                                    if LocMaxCD<1 || MaxValCD<2
                                        LevDec2 = 0.4013;
                                        LocMaxCD = 1;
                                        SeqFinCD(1)=2;
                                        MaxValCD = 0;
                                        InterCD(1)=LevDec2;
                                    else
                                        if (SeqFinCD(LocMaxCD)-MaxValCD/2)<1
                                            LevDec2 = 0.4013;
                                        else
                                            LevDec2 = InterCD(round(SeqFinCD(LocMaxCD)-MaxValCD/2));
                                        end
                                    end
                                    
                                    [MaxValEF,LocMaxEF]=max(SeqOnesEF);
                                    if LocMaxEF<1 || MaxValEF<2
                                        LevDec1 = 0.1877;
                                        LocMaxEF = 1;
                                        SeqFinEF(1)=2;
                                        MaxValEF = 0;
                                        InterEF(1)=LevDec1;
                                    else
                                        if (SeqFinEF(LocMaxEF)-MaxValEF/2)<1
                                            LevDec1 = 0.1877;
                                        else
                                            LevDec1 = InterEF(round(SeqFinEF(LocMaxEF)-MaxValEF/2));
                                        end
                                        
                                    end
                                    
                                    %another way to measure the eye opening is the get all the
                                    %boxes and find all peaks on it, that will be a plato created
                                    %by the sequences of ones (which was zeros). From thos peaks,
                                    %the eye diagram will be the longer of them hence it will take
                                    %the most part of the vector that store the findpeaks result.
                                    %Thus, the middle of the eye will be basically the middle of
                                    %the peaks vector.
                                    LocAB = find(EyeAB);
                                    if isempty(LocAB)                                              %if for some reason there are no peaks, something went wrong.
                                        LocAB = 1;
                                        LevelDec3 = 0.6953;%mean(Levels(3:4));                       %by default the decission level will be set to 0.65
                                    else
                                        if DecMod%InterAB(LocAB(round(end/2)))<=LevDec3
                                            LevelDec3 = LevDec3;
                                        else
                                            LevelDec3 = InterAB(LocAB(round(end/2)));
                                        end
                                    end
                                    
                                    LocCD = find(EyeCD);
                                    if isempty(LocCD)
                                        LocCD = 1;
                                        LevelDec2 = 0.4013;%mean(Levels(2:3));
                                    else
                                        if DecMod%InterCD(LocCD(round(end/2)))>=LevDec2
                                            LevelDec2 = LevDec2;
                                        else
                                            LevelDec2 = InterCD(LocCD(round(end/2)));
                                        end
                                    end
                                    
                                    LocEF = find(EyeEF);
                                    if isempty(LocEF)
                                        LocEF = 1;
                                        LevelDec1 = 0.1877;%mean(Levels(1:2));
                                    else
                                        if DecMod%InterEF(LocEF(round(end/2)))<=LevDec1
                                            LevelDec1 = LevDec1;
                                        else
                                            LevelDec1 = InterEF(LocEF(round(end/2)));
                                        end
                                    end
                                    AberLev(1,ThisCarr,CurrentTest)  = abs(InterAB(SeqFinAB(LocMaxAB)-1) - InterAB(SeqFinAB(LocMaxAB)-MaxValAB+1));
                                    AberLev(2,ThisCarr,CurrentTest)  = abs(InterCD(SeqFinCD(LocMaxCD)-1) - InterCD(SeqFinCD(LocMaxCD)-MaxValCD+1));
                                    AberLev(3,ThisCarr,CurrentTest)  = abs(InterEF(SeqFinEF(LocMaxEF)-1) - InterEF(SeqFinEF(LocMaxEF)-MaxValEF+1));
                                    ValsLev(1,ThisCarr,CurrentTest)  = LevDec3;
                                    ValsLev(2,ThisCarr,CurrentTest)  = LevDec2;
                                    ValsLev(3,ThisCarr,CurrentTest)  = LevDec1;
                                    ValsLev2(1,ThisCarr,CurrentTest) = InterAB(LocAB(round(end/2)));
                                    ValsLev2(2,ThisCarr,CurrentTest) = InterCD(LocCD(round(end/2)));
                                    ValsLev2(3,ThisCarr,CurrentTest) = InterEF(LocEF(round(end/2)));
                                    %%           Ploting for Qualitative Analizes
                                    if PrintinEye
                                        PrintInfo(Ploting*36,Ix,T,NPPB);
                                        hold all;
                                        plot(t(NPPB/2),LevDec1,'bd');plot(t(NPPB/2),InterEF(SeqFinEF(LocMaxEF)-1),'bo');plot(t(NPPB/2),InterEF(round(SeqFinEF(LocMaxEF)-MaxValEF+1)),'bx');
                                        plot(t(NPPB/2),LevDec2,'gd');plot(t(NPPB/2),InterCD(SeqFinCD(LocMaxCD)-1),'go');plot(t(NPPB/2),InterCD(round(SeqFinCD(LocMaxCD)-MaxValCD+1)),'gx');
                                        plot(t(NPPB/2),LevDec3,'kd');plot(t(NPPB/2),InterAB(SeqFinAB(LocMaxAB)-1),'ko');plot(t(NPPB/2),InterAB(round(SeqFinAB(LocMaxAB)-MaxValAB+1)),'kx');
                                        if CalcS
                                            plot(t(NPPB/2),ThisDeciLev3,'cd');plot(t(NPPB/2),mean(LevHi3),'co');plot(t(NPPB/2),mean(LevLo3),'cx');
                                            plot(t(NPPB/2),ThisDeciLev2,'rd');plot(t(NPPB/2),mean(LevHi2),'ro');plot(t(NPPB/2),mean(LevLo2),'rx');
                                            plot(t(NPPB/2),ThisDeciLev1,'md');plot(t(NPPB/2),mean(LevHi1),'mo');plot(t(NPPB/2),mean(LevLo1),'mx');
                                        end
                                        figure;
                                        hold all;
                                        plotpos = zeros(1,length(IxAux));
                                        plot(IxAux,plotpos,'o','color',[1 0.4 0]);
                                        plot(EvmMatRef(ObsCarrUsed(ThisCarr),:),plotpos,'b*');
                                        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                        a=1;
                                    end
                                    %             plot(t(NPPB/2),InterAB(LocAB(round(end/2))),'kx');
                                    %             plot(t(NPPB/2),InterCD(LocCD(round(end/2))),'gx');
                                    %             plot(t(NPPB/2),InterEF(LocEF(round(end/2))),'bx');
                                    %             plot([t(SymLoc(1)) t(SymLoc(1))],[0 1]);
                                    %
                                    %                             drawnow;
                                    %             PrintInfo(Ploting*38,EyeMax,EyeLoc-1,Interval);
                                    % %                         figure;plot(SeqOnesEF);figure;plot(SeqOnesCD);figure;plot(SeqOnesAB)
                                    %             LevDec3
                                    %             a1=InterAB(LocAB(round(end/2)))
                                    %             b1=ThisDeciLev1
                                    %             LevDec2
                                    %             a2=InterCD(LocCD(round(end/2)))
                                    %             b2=ThisDeciLev2
                                    %             LevDec1
                                    %             c1=InterEF(LocEF(round(end/2)))
                                    %             c2=ThisDeciLev3
                                    %             PrintInfo(Ploting*39,NPPB,OccuCount);
                                    %             PrintInfo(Ploting*40,OccuCount,SymLoc);
                                    %%      Actualy Receiving Data
                                    %Once the signal was processed the next step is through a
                                    %comparator decide the actual information received.
                                    IxRec = [];                                                    %Initialization of the vector that will store the income data
                                    IxRecDef = [];                                                 %Initialization of the vector that will store the income data
                                    IxRecDeS = [];                                                 %Initialization of the vector that will store the income data
                                    for kk=1:NPPB:length(Ix)                                       %The comparison process will be made for each symbol period
                                        %                 midaux = round(mean(SymLoc(1:round(end/2))));
                                        midaux = NPPB/2;%SymLoc(1);
                                        aux1 = Ix((kk-1)+midaux);     %An small portion of the income signal is take for evaluation
                                        MeanRec = mean(aux1);                                      %Measuring the avarage value of the samples taken
                                        %Verifying the interval for each symbol received.
                                        if MeanRec <= LevelDec1                                    %If it is the lowest level the incoming data
                                            IxRec = [IxRec 0 0];                                   %is 01 (1)
                                        elseif (MeanRec <= LevelDec2)&&(MeanRec > LevelDec1)       %If it is the second level the incoming data
                                            IxRec = [IxRec 0 1];                                   %is 00 (0)
                                        elseif (MeanRec <= LevelDec3)&&(MeanRec > LevelDec2)       %If it is the tird level the incoming data
                                            IxRec = [IxRec 1 1];                                   %is 10 (2)
                                        elseif MeanRec > LevelDec3                                 %If it is the uper level the incoming data
                                            IxRec = [IxRec 1 0];                                   %is 11 (3)
                                        else                                                       %If for some misteriose reason neither of previous verification were sucedded
                                            IxRec = [IxRec 0 0];                                   %by default the current data is set to be 00 (0)
                                        end
                                        
                                        %Verifying the interval for each symbol received.
                                        if MeanRec <= DecLevDef1                                   %If it is the lowest level the incoming data
                                            IxRecDef = [IxRecDef 0 0];                             %is 01 (1)
                                        elseif (MeanRec <= DecLevDef2)&&(MeanRec > DecLevDef1)     %If it is the second level the incoming data
                                            IxRecDef = [IxRecDef 0 1];                             %is 00 (0)
                                        elseif (MeanRec <= DecLevDef3)&&(MeanRec > DecLevDef2)     %If it is the tird level the incoming data
                                            IxRecDef = [IxRecDef 1 1];                             %is 10 (2)
                                        elseif MeanRec > DecLevDef3                                %If it is the uper level the incoming data
                                            IxRecDef = [IxRecDef 1 0];                             %is 11 (3)
                                        else                                                       %If for some misteriose reason neither of previous verification were sucedded
                                            IxRecDef = [IxRecDef 0 0];                             %by default the current data is set to be 00 (0)
                                        end
                                        
                                        if CalcS
                                            if MeanRec <= ThisDeciLev1                                   %If it is the lowest level the incoming data
                                                IxRecDeS = [IxRecDeS 0 0];                             %is 01 (1)
                                            elseif (MeanRec <= ThisDeciLev2)&&(MeanRec > ThisDeciLev1)     %If it is the second level the incoming data
                                                IxRecDeS = [IxRecDeS 0 1];                             %is 00 (0)
                                            elseif (MeanRec <= ThisDeciLev3)&&(MeanRec > ThisDeciLev2)     %If it is the tird level the incoming data
                                                IxRecDeS = [IxRecDeS 1 1];                             %is 10 (2)
                                            elseif MeanRec > ThisDeciLev3                                %If it is the uper level the incoming data
                                                IxRecDeS = [IxRecDeS 1 0];                             %is 11 (3)
                                            else                                                       %If for some misteriose reason neither of previous verification were sucedded
                                                IxRecDeS = [IxRecDeS 0 0];                             %by default the current data is set to be 00 (0)
                                            end
                                        end
                                    end
                                    %%           Ploting for Qualitative Analizes
                                    %             PrintInfo(Ploting*41,t(length(TxDataMat(ThisCarr,:))),Nb4Pam...
                                    %                                            /2,TxDataMat(ThisCarr,:),IxRec);
                                    %%       Calculating the Bit Error Ratio (BER)
                                    %The final process here is to count the number of wrongdoings
                                    %of this whole process upon the transmited data for
                                    %quantitative analizes
                                    BitErr = sum(xor(TxDataMat(ThisCarr,1+2*SyncPeriod:end-2*...
                                        SyncPeriod),IxRec(1+2*SyncPeriod:end-2*SyncPeriod)));%Comparison between the Transmited and received and counting the differences
                                    BitErrAux1 = BitErr;
                                    BitErrAux2 = sum(xor(TxDataMat(ThisCarr,1+2*SyncPeriod:end-...
                                        2*SyncPeriod),IxRecDef(1+2*SyncPeriod:end-2*SyncPeriod)));%Comparison between the Transmited and received and counting the differences
                                    if BitErr ~= 0
                                        if BitErrAux2<BitErrAux1
                                            BitErr = BitErrAux2;
                                        end
                                    else
                                        DecLevDef1 = LevelDec1;
                                        DecLevDef2 = LevelDec2;
                                        DecLevDef3 = LevelDec3;
                                    end
                                    if CalcS
                                        BitErrS = sum(xor(TxDataMat(ThisCarr,1+2*SyncPeriod:end-2*...
                                            SyncPeriod),IxRecDeS(1+2*SyncPeriod:end-2*SyncPeriod)));%Comparison between the Transmited and received and counting the differences
                                        Ber4PAMS(CurrentTest,ThisCarr) = BitErrS/(Nb4Pam-(2*SyncPeriod));%Calculating the ration of wrong bits and the total number of bits transmited
                                    end
                                    if BitErrS<BitErr
                                        Ber4PAM(CurrentTest,ThisCarr) = BitErrS/(Nb4Pam-(2*SyncPeriod));%Calculating the ration of wrong bits and the total number of bits transmited
                                    else
                                        Ber4PAM(CurrentTest,ThisCarr) = BitErr/(Nb4Pam-(2*SyncPeriod));%Calculating the ration of wrong bits and the total number of bits transmited
                                    end
                                    %             if BitErr>0
                                    %                 berpos = 1:2:size(Ber4PAM,2);
                                    %                 Ber4PAM(size(Ber4PAM,1),berpos)
                                    %                 drawnow;%pause(1);
                                    %                 a=2;
                                    %             end
                                    %             berpos1 = 1:2:size(Ber4PAM,2);
                                    %             berpos2 = 2:2:size(Ber4PAM,2);
                                    %             b = [Ber4PAM(1,berpos1);Ber4PAM(1,berpos2)]
                                    %             drawnow;%pause(1);
                                    %             a=3;
                                    %             drawnow;
                                    %             pause(0.5);
                                    close all;
                                else
                                    Ber4PAM(CurrentTest,ThisCarr) = 1;
                                end
                            end
                        otherwise
                            %%              Receiver OOK
                            for ThisCarr=InitCarr:CarrPass:NumCarr                                           %For each carrier the same process of reception will be used.
                                if CarrUsedUp(ThisCarr)
                                    %%  Reception Process: Handling the income optical field
                                    %At the first step the income field will be selected from the
                                    %optical FFT Output with the help of the maping vector
                                    %(previously described).
                                    EoutAux = EoutAux1(VetThisCarr==ObsCarrUsed(ThisCarr),:);
                                    if AddingNoiseP
                                        [p,l] = findpeaks(20*log10(abs(fftshift(fft(EoutAux)./length(EoutAux)))),'SortStr','descend','NPeaks',1);
                                        EoutAuxA = fftshift(fft(EoutAux));
                                        EoutAuxA(f==f(l)) = EoutAuxA(f==f(l-3));
                                        EoutAuxA = ifft(fftshift(EoutAuxA));
                                        [~,sigpower] = MeasPower(EoutAuxA);
                                        EoutAux      = awgn(EoutAux,osnrp,sigpower-30);
                                        [~,sigpower] = MeasPower(EoutAux);
                                    end
                                    %%            Fiber Time Delay Compensation
                                    switch Medium
                                        case 'Fiber'
                                            EoutAux = ifft(fft(EoutAux).*exp(1j*2*pi*f*(...
                                                FiberDelay(ThisCarr)*Ta)));
                                        otherwise
                                    end
                                    
                                    %The current incoming signal is them converted from the optical
                                    %domain to the eletrical domain with the help of an photo
                                    %detector.
                                    Ix =EoutAux.*conj(EoutAux);
                                    %%           Creating the Reception Filter
                                    [BitFilt,~] = FiltroGaussiano(f,BWD,CenFeq,FiltOrd);           %Creating filter for selection of the received signal
                                    BitFilt = fftshift(BitFilt);                                   %Shifting the filter for matching the received signal
                                    %%
                                    if RecFilBanPas
                                        Ix = ifft(fft(Ix).*BitFilt);
                                        Ix = Ix - min(Ix);                                         %Removing the DC component from them Eletrical signal received
                                        Ix = Ix./max(abs(Ix));                                     %Normalizing the eletrical signal (amplifying what is needed)
                                    end
                                    
                                    if ReceptorNoise
                                        [~,SigPower] = MeasPower(Ix);
                                        SigPower2 = SigPower-30;%10*log10(SigPower);
                                        SNR2 = CarSNR + 10*log10(1) - 10*log10(Nsamp) - 10*log10(1/2) + 10*log10(10^0.3);
                                        
                                        if TimeSys
                                            SNR = CarSNR + 10*log10(1) - 10*log10(0.25*T*BWD);
                                        else
                                            SNR = CarSNR + 10*log10(1) - 10*log10(0.25*t(2*NumAmosCP+NPPB)*BWD);
                                        end
                                        Ix = awgn(Ix,SNR,SigPower2);
                                    end
                                    if ~RecFilBanPas
                                        Ix = ifft(fft(Ix).*BitFilt);
                                        Ix = Ix - min(Ix);                                         %Removing the DC component from them Eletrical signal received
                                        Ix = Ix./max(abs(Ix));                                     %Normalizing the eletrical signal (amplifying what is needed)
                                    end
                                    
                                    if CurrentModula==2
                                        VetElecPowerUp(CurrentTest,ThisCarr)= MeasPower(Ix);
                                        EoutAuxF = 20*log10(abs(fftshift(fft(EoutAux)./length(...
                                            EoutAux))));
                                        [VetOptiPowerUp(CurrentTest,ThisCarr),~]= findpeaks(EoutAuxF,...
                                            'SortStr','descend','NPeaks',1);
                                    else
                                        VetElecPower(CurrentTest,ThisCarr)= MeasPower(Ix);
                                        EoutAuxF = 20*log10(abs(fftshift(fft(EoutAux)./length(...
                                            EoutAux))));
                                        [VetOptiPower(CurrentTest,ThisCarr),~]= findpeaks(EoutAuxF,...
                                            'SortStr','descend','NPeaks',1);
                                    end
                                    
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
                                    
                                    AuxSync = (Ix(IniSyncPos:SyncPos));                            %Selecting the sync-word within the received signal
                                    AuxSync1 = AuxSync;
                                    %This synchronization process is based on the mean expected
                                    %value. That means, the information mean value should be within
                                    %one period of symbol. Thus, the mean value of the received
                                    %signal is acquired and compare of the known sync-word to
                                    %verify if this mean value is at the right possition
                                    AuxSync1(AuxSync1<0) = 0;                                      %To keep the mean value above zero anything under is neglected
                                    AuxSync1(AuxSync1>=mean(AuxSync1)) = 1;                        %Adding a flag to the first sample of the received mean value
                                    AuxSync1(AuxSync1<mean(AuxSync1)) = -1;                        %All the others samples at set to the lowest level
                                    AuxSync2 = SyncSymb(IniSyncPos:SyncPos);                       %Selecting the sync-word within the known signal
                                    AuxSync2(AuxSync2>=mean(AuxSync2)) = 1;                        %Adding a flag to the first sample of the known mean value
                                    AuxSync2(AuxSync2<mean(AuxSync2)) = -1;                        %All the others samples at set to the lowest level
                                    
                                    PosToSyn  = find(ismember(AuxSync1,1));
                                    PosSyn    = find(ismember(AuxSync2,1));
                                    
                                    AuxSyncCorr = PosToSyn(round(end/2))-PosSyn(round(end/2));
                                    sn=((1/2)*(((1/2)^MaxNumStag)-1))/((1/2)-1);
                                    AuxSyncCorr = round((sn/(2/2^IfftOrSum))*T/Ta);
                                    %The difference between the PossitionTosynchronize and
                                    %Possitionsynchronized will be used to correct the time
                                    %shifting on the transmition and reception process.
                                    
                                    if AuxSyncCorr>=0%If the difference is positive, left-shift...
                                        %                 Ix = ifft(fft(Ix).*exp(1j*2*pi*f*(AuxSyncCorr*Ta)));     %Shift based on time change
                                        Ix = [Ix(AuxSyncCorr+1:end) Ix(1:AuxSyncCorr)];            %Shift based on sampling sliding
                                    else%... but if the difference is negative, right-shift
                                        %                 Ix = ifft(fft(Ix).*exp(-1j*2*pi*f*(AuxSyncCorr*Ta)));    %Shift based on time change
                                        Ix = [Ix(end+AuxSyncCorr+1:end) Ix(1:end + AuxSyncCorr)];  %Shift based on sampling sliding
                                    end
                                    
                                    
                                    %%        synchronizing
                                    if SencondAdjust
                                        %For some reason that we could not understand sometimes the
                                        %time (sampling) sliding of the signal is not equal
                                        %throught the data stream. Thus, the second part of the
                                        %synchronism process will be turn ON or OFF according to the
                                        %user's will.
                                        AuxSyncEnd = (Ix(end-(SyncPos+1*abs(AuxSyncCorr))+1:end-...
                                            (2*NPPB+1*abs(AuxSyncCorr))));%Selecting the sync-word within the received signal
                                        SyncSymbEndAux = SyncSymbEnd(end-SyncPos+1:end-2*NPPB);    %Selecting the sync-word within the known signal
                                        AuxSyncEnd1 = AuxSyncEnd;
                                        AuxSyncEnd1(AuxSyncEnd1<0) = 0;                            %To keep the mean value above zero anything under is neglected
                                        AuxSyncEnd1(AuxSyncEnd1>=mean(AuxSyncEnd1)) = 1;           %Adding a flag to the first sample of the received mean value
                                        AuxSyncEnd1(AuxSyncEnd1<mean(AuxSyncEnd1)) = -1;           %All the others samples at set to the lowest level
                                        AuxSyncEnd2 = SyncSymbEndAux;
                                        AuxSyncEnd2(AuxSyncEnd2>=mean(AuxSyncEnd2)) = 1;           %Adding a flag to the first sample of the known mean value
                                        AuxSyncEnd2(AuxSyncEnd2<mean(AuxSyncEnd2)) = -1;           %All the others samples at set to the lowest level
                                        
                                        PosToSynEnd  = find(ismember(AuxSyncEnd1,1));
                                        PosSynEnd = find(ismember(AuxSyncEnd2,1));
                                        
                                        AuxSyncCorrEnd = PosToSynEnd(round(end/2))-PosSynEnd(...
                                            round(end/2));
                                        
                                        %The difference between the PossitionTosynchronize and
                                        %Possitionsynchronized will be used to correct the time
                                        %shifting on the transmition and reception process.
                                        
                                        if AuxSyncCorrEnd>=0%If possitive difference, left-shift...
                                            %                     Ix = ifft(fft(Ix).*exp(1j*2*pi*f*(AuxSyncCorrEnd*...
                                            %                                                                    Ta)));%Shift based on time change
                                            Ix = [Ix(AuxSyncCorrEnd+1:end) Ix(1:AuxSyncCorrEnd)];
                                        else%... but if the difference is negative, right-shift
                                            %                     Ix = ifft(fft(Ix).*exp(-1j*2*pi*f*(AuxSyncCorrEnd*...
                                            %                                                                    Ta)));%Shift based on time change
                                            Ix = [Ix(end+AuxSyncCorrEnd+1:end) Ix(1:end + ...
                                                AuxSyncCorrEnd)];%Shift based on sampling sliding
                                        end
                                    end
                                    
                                    %% Ploting the result for qualitative analizes
                                    
                                    %                             PrintInfo(Ploting*45,t(IniSyncPos:SyncPos),Ix(IniSyncPos:...
                                    %                                                     SyncPos),SyncSymb(IniSyncPos:SyncPos));
                                    %                                             PrintInfo(Ploting*46,t(1:length(Ix)),Ix);
                                    %%                  Removing CP
                                    %                 if ThisCarr==126
                                    %                     EyeToPlot(CurrentTest,1:length(Ix)) = Ix;
                                    %                     save(['savingforplotingeye' num2str(VetSnr)],'EyeToPlot','VetSnr','IfftOrSum','UsedModula','T','NPPB');
                                    %                 end
                                    if AddCP
                                        IxAux = Ix(1:end - StuffSampels);
                                        IxAux = reshape(IxAux,(2*NumAmosCP+NPPB),Nb-NumBitDesc);
                                        IxAux = IxAux(1+NumAmosCP:end-NumAmosCP,:);
                                        IxAux = reshape(IxAux,1,NPPB*(Nb-NumBitDesc));
                                        Ix    = IxAux;
                                    end
                                    %% Taking the sampling the EVM meassurement
                                    PosAuxEout = NPPB/2:NPPB:length(Ix);                       %Varriable respossible to take just the samples at the middle of the symbol
                                    IxAux      = Ix(PosAuxEout);                               %Normalizing the reference
                                    a=0;
                                    RxSymbAmos(CurrentTest,ObsCarrPos==ThisCarr,:) = IxAux;%RxSymbAmos = [];
                                    EvmMatRec(ObsCarrPos==ThisCarr,:) = IxAux;                       %Taking just the middle samples as references
                                    [EvmDB(CurrentTest,ThisCarr), EvmPer(CurrentTest,ThisCarr), EvmRms(CurrentTest,ThisCarr) ] = EvmCalc( EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux );
                                    [EvmDBJ(CurrentTest,ThisCarr),EvmPerJ(CurrentTest,ThisCarr),EvmRmsJ(CurrentTest,ThisCarr)] = evm1(2,'pam',EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux);
                                    %                 if ThisCarr==2
                                    %                     figure;hold all;grid on;
                                    %                     plot(IxAux,zeros(1,length(IxAux)),'rx');
                                    %                     set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                    %                     a=0;
                                    %                 end
                                    %% Ploting the result for qualitative analizes
                                    
                                    %##########################################################################
                                    PosIx = NPPB/2:NPPB:length(Ix);                                %Possition of the central samples - but the number of samples per symbol is even ????
                                    IxAux = Ix(PosIx);                                             %From the main received signal just a few samples are taken for further evaluation
                                    n     = 100;                                                   %The number of boxes to be filled up on the histogram process
                                    MinDist = n/3;
                                    Perce = 0.5;
                                    %At this process each eye will be evaluated separetelly. It is
                                    %important to remember that the PAM4 has 4 levels, which means
                                    %three level of decissions that are exaclty the center of the
                                    %eye diagram.
                                    Interval = linspace(Perce*(min(Ix)),Perce*(max(Ix)),n);                     %Building the histogram boxes
                                    
                                    %Therefore, the MATLAB hist function returns the number of
                                    %occurrence of each interval.
                                    EyeMax = hist(Ix,Interval);
                                    EyeMaxaux = [0 EyeMax 0];                                      %Zeros are added at the EyeMax to auxiliate the finding peaks process
                                    [~,EyeLoc] = findpeaks(EyeMaxaux,'MinPeakDistance',MinDist,'SortStr','descend','NPeaks',2);%The peaks on the Eye profile will be the levels at the Eyes limit
                                    
                                    if length(EyeLoc)<2%If it was not able to find the eye profile.
                                        [~,EyeLoc] = findpeaks(EyeMaxaux,'MinPeakDistance',MinDist/2,'SortStr','descend','NPeaks',2);%The peaks on the Eye profile will be the levels at the Eyes limit
                                    end
                                    
                                    if length(EyeLoc)<2%If it was not able to find the eye profile.
                                        EyeLoc = [2 3];
                                        Levels = [0 -0.2 0.35];
                                    else%Whereas, if there is an profile the decission can be found
                                        Levels = [Interval(EyeLoc(1)-1) Interval(EyeLoc(2)-1)];
                                        Levels = sort(Levels);
                                        %                 LevelDec1 = Levels(1) + (Levels(2)-Levels(1))/2 ;
                                        %                 LevelDec2 = Levels(2) + (Levels(3)-Levels(2))/2 ;
                                        %                 LevelDec3 = Levels(3) + (Levels(4)-Levels(3))/2 ;
                                    end
                                    IxAuxAB = IxAux((IxAux<=Levels(2))&(IxAux>=Levels(1)));
                                    Inter = linspace(Levels(1),Levels(2),n);                     %Building the histogram boxes
                                    Inter2 = linspace(min(IxAux),max(IxAux),n);                     %Building the histogram boxes
                                    EyeAB = hist(IxAuxAB,Inter);                                 %filling up the boxes with samples that fit on them.
                                    EyeCD = hist(IxAux,Inter2);                                 %filling up the boxes with samples that fit on them.
                                    
                                    %What we are looking for here are not where exist the
                                    %occurrences of values. The eye is where there is not samples,
                                    %this means, the decission level is a reagion between actual
                                    %levels thus this reagion does not contain any signal.
                                    EyeAB = ~EyeAB;                                                %Changing zeros to one - Zeros compose the eye region
                                    EyeCD = ~EyeCD;                                                %Changing zeros to one - Zeros compose the eye region
                                    %Starting variables that will be used to identify the eye
                                    %diagram.
                                    Count   = 1;
                                    Count2   = 1;
                                    
                                    SeqOnes = 0;
                                    SeqOnes2 = 0;
                                    
                                    SeqIni  = 1;
                                    SeqIni2  = 1;
                                    
                                    SeqFin  = 0;
                                    SeqFin2  = 0;
                                    %The for loop will take account of every box with ones. It is
                                    %important to take note that the not operator was used in this
                                    %vector, therefore ones means zeros (the eye diagram -
                                    %possibly) and zeros means values abouve zeroa (not the eye).
                                    for kk=1:length(EyeAB)                                         %For every box
                                        if EyeAB(kk)                                               %if it contains "1"
                                            SeqOnes(SeqIni)=Count;                           %count this element as part of a consecutive sequency
                                            Count = Count + 1;                                 %adds one to the counter of consecutive elements "1"
                                            if kk==length(EyeAB)                                   %if the current box is the last box we got to an end
                                                SeqFin(SeqIni) = kk;                           %The final sequency element is equal to its possition (kk)
                                            end
                                        else                                                       %else if the current box contains "0"
                                            SeqFin(SeqIni) = kk-1;                             %the previous element was the last element of a consecutive sequency of ones
                                            SeqIni = SeqIni + 1;                               %adds one to the counter of consecutive number (take in account how many sequencies there are)
                                            Count = 1;                                           %reset the counter of consecutive elements to mark the start of a new consecutive sequence
                                        end
                                        
                                        if EyeCD(kk)                                               %if it contains "1"
                                            SeqOnes2(SeqIni2)=Count2;                           %count this element as part of a consecutive sequency
                                            Count2 = Count2 + 1;                                 %adds one to the counter of consecutive elements "1"
                                            if kk==length(EyeCD)                                   %if the current box is the last box we got to an end
                                                SeqFin2(SeqIni2) = kk;                           %The final sequency element is equal to its possition (kk)
                                            end
                                        else                                                       %else if the current box contains "0"
                                            SeqFin2(SeqIni2) = kk-1;                             %the previous element was the last element of a consecutive sequency of ones
                                            SeqIni2 = SeqIni2 + 1;                               %adds one to the counter of consecutive number (take in account how many sequencies there are)
                                            Count2 = 1;                                           %reset the counter of consecutive elements to mark the start of a new consecutive sequence
                                        end
                                    end
                                    
                                    
                                    %If the eye is open, which means there is a clear difference
                                    %between adjacent levels, the eye diagram will be the longest
                                    %sequence of ones.
                                    [MaxVal,LocMax]=max(SeqOnes);
                                    if LocMax<2 || MaxVal<2                                    %if any sequency was found or there is just one sequency it is a error thus
                                        LevDec = 0.0;                                            %the decission level will be by default 0.67. Also, the other variables will
                                        LocMax = 1;                                              %will be set with values to not cause errors in the future
                                        SeqFin(1)=2;
                                        MaxVal = 0;
                                        Inter(1)=LevDec;
                                    else                                                           %if a sequency was found the middle of the eye will be the middle of the sequency
                                        if (SeqFin(LocMax)-MaxVal/2)<1                       %if for some reason the final element of the sequency minus the half of its
                                            LevDec = 0.0;                                         %length results in a negative value, something went very wrong, and by
                                        else                                                       %default it will be set to 0.7
                                            LevDec = Inter(round(SeqFin(LocMax)-MaxVal/2));%Otherwise, the decission level is the middle of the sequency
                                        end
                                    end
                                    [MaxVal2,LocMax2]=max(SeqOnes2);
                                    if LocMax2<2 || MaxVal2<2                                    %if any sequency was found or there is just one sequency it is a error thus
                                        LevDec2 = 0.0;                                            %the decission level will be by default 0.67. Also, the other variables will
                                        LocMax2 = 1;                                              %will be set with values to not cause errors in the future
                                        SeqFin2(1)=2;
                                        MaxVal2 = 0;
                                        Inter2(1)=LevDec2;
                                    else                                                           %if a sequency was found the middle of the eye will be the middle of the sequency
                                        if (SeqFin2(LocMax2)-MaxVal2/2)<1                       %if for some reason the final element of the sequency minus the half of its
                                            LevDec2 = 0.0;                                         %length results in a negative value, something went very wrong, and by
                                        else                                                       %default it will be set to 0.7
                                            LevDec2 = Inter2(round(SeqFin2(LocMax2)-MaxVal2/2));%Otherwise, the decission level is the middle of the sequency
                                        end
                                    end
                                    
                                    %another way to measure the eye opening is the get all the
                                    %boxes and find all peaks on it, that will be a plato created
                                    %by the sequences of ones (which was zeros). From thos peaks,
                                    %the eye diagram will be the longer of them hence it will take
                                    %the most part of the vector that store the findpeaks result.
                                    %Thus, the middle of the eye will be basically the middle of
                                    %the peaks vector.
                                    Loc = find(EyeAB);
                                    if isempty(Loc)                                              %if for some reason there are no peaks, something went wrong.
                                        Loc = 1;
                                        LevelDec = 0.0;%mean(Levels(3:4));                       %by default the decission level will be set to 0.65
                                    else
                                        LevelDec = LevDec;
                                    end
                                    Loc2 = find(EyeCD);
                                    if isempty(Loc2)                                              %if for some reason there are no peaks, something went wrong.
                                        Loc2 = 1;
                                        LevelDec2 = 0.0;%mean(Levels(3:4));                       %by default the decission level will be set to 0.65
                                    else
                                        LevelDec2 = LevDec2;
                                    end
                                    %##########################################################################
                                    if CalcS
                                        [~,~,EyeOpen,EyeOpenHigh,EyeOpenLow] = Olho((Ix),T,NPPB,0);
                                        AberLevS(CurrentTest,ThisCarr)= EyeOpen;
                                        ValsLevS(CurrentTest,ThisCarr)= EyeOpenLow + EyeOpen/2;
                                    else
                                        [~,~,EyeOpen,EyeOpenHigh,EyeOpenLow] = Olho((Ix),T,NPPB,0);
                                    end
                                    %% Ploting the result for qualitative analizes
                                    if PrintinEye
                                        PrintInfo(PlotingThis*47,(Ix),T,NPPB);
                                        if ThisCarr==126
                                            a=1;
                                        end
                                        hold on;
                                        plot(t((NPPB)/2),Inter(round(SeqFin(LocMax)-MaxVal)),'kx');
                                        plot(t((NPPB)/2),Inter(SeqFin(LocMax)),'ko');
                                        plot(t((NPPB)/2),Inter(round(SeqFin(LocMax)-MaxVal/2)),'kd');
                                        plot(t((NPPB)/2),Inter2(round(SeqFin2(LocMax2)-MaxVal2)),'bx');
                                        plot(t((NPPB)/2),Inter2(SeqFin2(LocMax2)),'bo');
                                        plot(t((NPPB)/2),Inter2(round(SeqFin2(LocMax2)-MaxVal2/2)),'bd');
                                        if CalcS
                                            plot(t((NPPB)/2),EyeOpenLow,'mx');
                                            plot(t((NPPB)/2),EyeOpenHigh,'mo');
                                            plot(t((NPPB)/2),EyeOpenLow + EyeOpen/2,'md');
                                        end
                                        figure;
                                        hold all;
                                        plotpos = zeros(1,length(IxAux));
                                        plot(IxAux,plotpos,'o','color',[1 0.4 0]);
                                        plot(EvmMatRef(ObsCarrUsed(ThisCarr),:),plotpos,'b*');
                                        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                        a=0;
                                    end
                                    %                 AberLev(CurrentTest,ThisCarr) = EyeOpen;
                                    %                 ValsLev(CurrentTest,ThisCarr) = EyeOpenLow+EyeOpen/2;
                                    %                 [~,~,EyeOpen,EyeOpenHigh,EyeOpenLow] = Olho((Ix),T,NPPB,0);
                                    %                 hold on;
                                    %                 plot(t((NPPB)/2),EyeOpenLow,'kx');
                                    %                 plot(t((NPPB)/2),EyeOpenHigh,'ko');
                                    %                 plot(t((NPPB)/2),EyeOpenLow + EyeOpen/2,'kd');
                                    if CurrentModula==2
                                        AberLevUp(CurrentTest,ThisCarr) = Inter(SeqFin(LocMax))-Inter(round(SeqFin(LocMax)-MaxVal));
                                        ValsLevUp(CurrentTest,ThisCarr) = LevDec;
                                    else
                                        AberLev(CurrentTest,ThisCarr) = Inter(SeqFin(LocMax))-Inter(round(SeqFin(LocMax)-MaxVal));
                                        ValsLev(CurrentTest,ThisCarr) = LevDec;
                                    end
                                    %%         Finding Decission Levels
                                    %The process for decoding the income signal will be based on
                                    %eletronic comparators. Inasmuch as the right decission level
                                    %must be acquired for accurately decide, within a symbol
                                    %periode, what that current leavel means (ones or zeros).
                                    %
                                    %The process hereafter of chosing the  decission levels is not
                                    %deterministic rather it is a statistic process. The main idea
                                    %is to take the decission level from the histogram generated
                                    %from the income signal stream.
                                    %
                                    %This process is realized inside the function Olho.
                                    %
                                    %Basicaly the decission level will be the minimal value of the
                                    %currente eye under evaluation plus the half of the its eye
                                    %opening.The following ilustration better describe this process
                                    %
                                    %Eye Limits:   + (1/2 * Eye Opening:)  =    Comparisson Limit:
                                    %
                                    %UperLevel  ______     ______________     _____
                                    %                 \   /      |       \   /
                                    %                  \ /       |        \ /
                                    %                   \ Half Eye Opening /     Decission Level
                                    %                  / \       |        / \
                                    %LowerLevel ______/   \______|_______/   \_____
                                    %
                                    %
                                    %
                                    %For accurately read the income signal is necessary more than
                                    %just have the decision level, it is also needed to know where
                                    %it occurs. For instance, if we decide to use just to take the
                                    %decision level and measure the mean value across and a portion
                                    %of the time period, which portion should we take? It would be
                                    %logical to take the central portion, but the bit may be
                                    %deformed in such way that the information gets concentrated
                                    %out of the middle part, which means the symbol is not
                                    %symmetric. The symmetry information can be acquired from the
                                    %eye diagram by measuring the longitudinal opening. The
                                    %following sketch better describes this process:
                                    %
                                    %                   Point of Symmetry
                                    %         _____     _______|_______     _____
                                    %              \   /               \   /
                                    %               \ /  Longitudinal   \ /
                                    %                \ ----------------- /
                                    %               / \    Opening      / \
                                    %         _____/   \_______________/   \_____
                                    %
                                    %With those two pieces of information, decision level and point
                                    %of symmetry, we have the X, Y coordinates for the centre of
                                    %the Eye Diagram. Therefore, as long as there is an opening on
                                    %it it will be possible to recover the transmitted information
                                    %without error... theoretically.
                                    %
                                    %As this process is also statistical, first we reshape the
                                    %income vector to analyze all periods at the same time.
                                    
                                    %                 EyeSymMat = reshape(Ix(1+SyncPeriod*NPPB:end-SyncPeriod*...
                                    %                     NPPB),NPPB,(Nb-NumBitDesc)-2*SyncPeriod);
                                    %Then we take the values that compose the decision level
                                    %because they will mark the point of symmetry.
                                    %
                                    %Firstly it was set the interval in which the histogram will be
                                    %build. It is based on the number of samples per bit period.
                                    %                 Interval = linspace(min(Ix),max(Ix),2*NPPB);
                                    %Therefore, the MATLAB hist function returns the number of
                                    %occurrence of each interval.
                                    %                 EyeMax = hist(Ix,Interval);
                                    %                 EyeMaxaux = [0 EyeMax 0];                                      %Zeros are added at the EyeMax to auxiliate the finding peaks process
                                    %                 [~,EyeLoc] = findpeaks(EyeMaxaux,'MinPeakDistance',NPPB/4,...
                                    %                     'SortStr','descend','NPeaks',2);%The peaks on the Eye profile will be the levels at the Eyes limit
                                    
                                    %             figure;findpeaks(EyeMaxaux,'MinPeakDistance',NPPB/4,'SortStr','descend','NPeaks',2);
                                    
                                    %                 ValToSeek = Interval(EyeLoc-1);
                                    %From the location of the max values that occured, which means
                                    %the uper and lower level of the eye diagram it needs to take
                                    %the actual value that those occurences represent that is
                                    %withing the the Interval variable.
                                    %                 ValToSeek = sort(ValToSeek,'ascend');
                                    %The number of ocurrences is a statical measure therefore one
                                    %does not have control which interval will have the highest
                                    %peak, thus it is important to ordenate the values to be seek
                                    %from the lower part of the eye diagram to the uper part of the
                                    %eye diagram.
                                    %                 OccuCount = zeros(1,size(EyeSymMat,1));                        %Auxiliar Variable for accounting.
                                    %                 for kk=1:size(EyeSymMat,1)                                     %For every sample within a symbol period
                                    %                     OccuCount(kk) = OccuCount(kk)+sum((EyeSymMat(kk,:)>=min(...
                                    %                         Ix))&(EyeSymMat(kk,:)<=UpeSymPer*EyeOpenLow));%Account all occurencies of the value 1
                                    %                     OccuCount(kk) = OccuCount(kk)+sum((EyeSymMat(kk,:)>=...
                                    %                         LowSymPer*EyeOpenHigh)&(EyeSymMat(kk,:)<=max(Ix)));%Account all occurencies of the value 2
                                    %                 end
                                    %The point of symmetry of the eye diagram will be where the
                                    %maximum number of occurrences were measured inasmuch as those
                                    %are the points where all the bits go to the center of the
                                    %symbol.  From the maximum number of occurrences, it can happen
                                    %for more than one sample within one symbol period, in the best
                                    %case, all samples would have the same accounting as it is
                                    %shown the ilustration above hence the symmetry will be at the
                                    %middle sample of this group of maximum occurrences. This value
                                    %can be found by the mean of the samples positions within the
                                    %symbol period. The problem with this approach is that the
                                    %signal must be synchronized with the maximum displacement of
                                    %a symbol period minus 25% of the eye Longitudinal opening if
                                    %the displacement is higher than that the point of symmetry
                                    %will be wrongly measured.
                                    %                 [SymLoc] = round(mean(find(ismember(OccuCount,max(OccuCount)...
                                    %                     ))));%The peak on the Eye profile will be the Symmetry level
                                    %%           Ploting for Qualitative Analizes
                                    %             PrintInfo(Ploting*48,NPPB,OccuCount);
                                    %             PrintInfo(Ploting*49,OccuCount,SymLoc);
                                    
                                    %%      Actualy Receiving Data
                                    %Once the signal was processed the next step is through a
                                    %comparator decide the actual information received.
                                    EoutCorr = [];                                                 %Initialization of the vector that will store the income data
                                    EoutCorrD = [];                                                 %Initialization of the vector that will store the income data
                                    EoutCorr2 = [];                                                 %Initialization of the vector that will store the income data
                                    EoutCorrS = [];                                                 %Initialization of the vector that will store the income data
                                    for kk=1:NPPB:length(Ix)                                       %The comparison process will be made for each symbol period
                                        %An small portion of the income signal is take for
                                        %evaluation by measuring the avarage value of the samples
                                        %taken
                                        %                 CalcMean = mean((Ix((kk-1)+SymLoc(1))));
                                        CalcMean = mean((Ix((kk-1)+NPPB/2)));
                                        %Verifying the interval for each symbol received.
                                        if CalcMean >= LevDec%EyeOpenLow+EyeOpen/2                        %If it is the uper level the incoming data
                                            EoutCorr = [EoutCorr 1];                               %is 1
                                        else                                                       %If it is the lowest level the incoming data
                                            EoutCorr = [EoutCorr 0];                               %is 0
                                        end
                                        CalcMean2 = mean((Ix((kk-1)+NPPB/2)));
                                        %Verifying the interval for each symbol received.
                                        if CalcMean2 >= LevDec2%EyeOpenLow+EyeOpen/2                        %If it is the uper level the incoming data
                                            EoutCorr2 = [EoutCorr2 1];                               %is 1
                                        else                                                       %If it is the lowest level the incoming data
                                            EoutCorr2 = [EoutCorr2 0];                               %is 0
                                        end
                                        if CalcMean2 >= EyeOpenLow+EyeOpen/2                        %If it is the uper level the incoming data
                                            EoutCorrD = [EoutCorrD 1];                               %is 1
                                        else                                                       %If it is the lowest level the incoming data
                                            EoutCorrD = [EoutCorrD 0];                               %is 0
                                        end
                                        if CalcS
                                            CalcMeanS = mean((Ix((kk-1)+NPPB/2)));
                                            %Verifying the interval for each symbol received.
                                            if CalcMeanS >= EyeOpenLow + EyeOpen/2%EyeOpenLow+EyeOpen/2                        %If it is the uper level the incoming data
                                                EoutCorrS = [EoutCorrS 1];                               %is 1
                                            else                                                       %If it is the lowest level the incoming data
                                                EoutCorrS = [EoutCorrS 0];                               %is 0
                                            end
                                        end
                                    end
                                    
                                    %PrintInfo(Ploting*50,linspace(0,t(end),length(EoutCorr))...
                                    %                                           ,TxDataMat(ThisCarr,:),EoutCorr);
                                    %%       Calculating the Bit Error Ratio (BER)
                                    %The final process here is to count the number of wrongdoings
                                    %of this whole process upon the transmited data for
                                    %quantitative analizes
                                    if CalcS
                                        BitErrS = sum(xor(TxDataMat(ThisCarr,1+SyncPeriod:end-...
                                            SyncPeriod),EoutCorrS(1+SyncPeriod:end-SyncPeriod)));%Comparison between the Transmited and received and counting the differences
                                        BerOOKS(CurrentTest,ThisCarr) = BitErrS/((Nb-NumBitDesc)-...
                                            2*SyncPeriod);%Calculating the ration of wrong bits and the total number of bits transmited
                                    end
                                    BitErr = sum(xor(TxDataMat(ThisCarr,1+SyncPeriod:end-...
                                        SyncPeriod),EoutCorr(1+SyncPeriod:end-SyncPeriod)));%Comparison between the Transmited and received and counting the differences
                                    BitErr2 = sum(xor(TxDataMat(ThisCarr,1+SyncPeriod:end-...
                                        SyncPeriod),EoutCorr2(1+SyncPeriod:end-SyncPeriod)));%Comparison between the Transmited and received and counting the differences
                                    BitErrD = sum(xor(TxDataMat(ThisCarr,1+SyncPeriod:end-...
                                        SyncPeriod),EoutCorrD(1+SyncPeriod:end-SyncPeriod)));%Comparison between the Transmited and received and counting the differences
                                    if BitErr2<=BitErr
                                        if BitErr2<=BitErrD
                                            BerOOK(CurrentTest,ThisCarr) = BitErr2/((Nb-NumBitDesc)-...
                                                2*SyncPeriod);%Calculating the ration of wrong bits and the total number of bits transmited
                                            AberLev(CurrentTest,ThisCarr)= Inter2(SeqFin2(LocMax2))-Inter2(round(SeqFin2(LocMax2)-MaxVal2));
                                            ValsLev(CurrentTest,ThisCarr)= LevDec2;
                                        else
                                            BerOOK(CurrentTest,ThisCarr) = BitErrD/((Nb-NumBitDesc)-...
                                                2*SyncPeriod);%Calculating the ration of wrong bits and the total number of bits transmited
                                            AberLev(CurrentTest,ThisCarr)= EyeOpen;
                                            ValsLev(CurrentTest,ThisCarr)= EyeOpenLow+EyeOpen/2;
                                        end
                                    else
                                        if BitErr<=BitErrD
                                            BerOOK(CurrentTest,ThisCarr) = BitErr/((Nb-NumBitDesc)-...
                                                2*SyncPeriod);%Calculating the ration of wrong bits and the total number of bits transmited
                                        else
                                            BerOOK(CurrentTest,ThisCarr) = BitErrD/((Nb-NumBitDesc)-...
                                                2*SyncPeriod);%Calculating the ration of wrong bits and the total number of bits transmited
                                            AberLev(CurrentTest,ThisCarr)= EyeOpen;
                                            ValsLev(CurrentTest,ThisCarr)= EyeOpenLow+EyeOpen/2;
                                        end
                                    end
                                    %PrintInfo(Ploting*50,linspace(0,t(end),length(EoutCorr))...
                                    %                                           ,TxDataMat(ThisCarr,:),EoutCorr);
                                    
                                    %             berpos1 = 1:2:size(BerOOK,2);
                                    %             berpos2 = 2:2:size(BerOOK,2);
                                    %             b = [BerOOK(1,berpos1);BerOOK(1,berpos2)]
                                    %             a=1;
                                    close all;
                                else
                                    BerOOK(CurrentTest,ThisCarr) = 1;
                                end
                                %         end
                            end
                    end
                    a=0;
                    % if CurrentModula == 1
                    %     if size(BerOOK,2)<NumCarr
                    %         BerOOK(:,end:end+NumCarr-size(BerOOK,2)) = 1;
                    %     end
                    % end
                    clear EoutAux1 VetThisCarr Esync Esync1 Esync2 EoutA EoutB EoutC EoutD
                    % clearvars -except UsedModula TestBundle ThisModula CurrentTest TxDataMat...
                    % CurrentModula ThisTest ThisModula OfcName BerOOK Ber4PAM BerDQPSK OSNRPC...
                    % MaxNumStag EoutUpRec BerDPSK CurrentMedium OcsToTest CurrentOCS t f T...
                    %                                                                   PulseResp
                    
                    
                end
            end
            %%              Saving
            
            clearvars -except ThisPath UofcsLocal OfcName NumbOfTest PowRef ...
                PAM4Type ThisTest TestBundle CurrentMedium UsedModula UsedModula...
                Ber4PAM SnrRef CarSNR FiberDelay ValsLev2 ValsLev contcp ExRa nn...
                PastTime SigAten CurrentModula Atenuation VetOptiPower AberLev T...
                VetElecPower FFTSplit ThisModula SendingDowStr VetElecPowerQ Tcp...
                VetElecPowerI BerDPSK BerDQPSK CurrentTest VetOptiPowerUp mm Nb ...
                AberLevUp VetElecPowerUp ValsLevUp EvmRms EvmDB EvmPer osnr snr ...
                DecLevDef1 DecLevDef2 DecLevDef3 UpStModula EvmRmsUp EvmDBUp Ns ...
                EvmPerUp UsingAllCarriers ObsCarrUsed ll CarrToSendInt3 Bref Tu ...
                CarrToSend6 CarrToSend5 CarrToSend4 CarrToSend3 CarrToSendInt2  ...
                CarrToSend2 ObsCarr CarrToSendInt CarrToSend CarrEvenPos planck ...
                CarrOddPos CarrPosInt ObsCarrPos NumCarrUsed VetNumCar VetOpen k...
                FiberLength lambdac NumCarr CarrUsed SavingAT RefCarr VetOpenO1 ...
                VetOpenO2 VetOpenO3 VetOpenA1 VetOpenA2 VetOpenA3 SavingAT2 NFFT...
                VetNumCar2 VetOpen2 IfftOrSum MaxNumStagT OSNRPC EoutTx BerOOK N...
                BERt4PAM BERtOOK InitCarrUp InitCarrDo CarrPass CarrUsedDo snrp ...
                CarrUsedUp SelecSetUP VetThisCarrTx CarrRecPowUp CarrRecPowDo Tg...
                AberLevQ ValsLevQ AberLevI ValsLevI LevDefDpqsk LightSpeed Dtf g...
                AddingNoiseO AddingNoiseE AsePower osnrf AddingNoiseF BerOOKS Ts...
                AddingNoiseO AddingNoiseE AddingNoiseP osnrp AberLevS ValsLevS M...
                BerDPSKS BerDQPSKS AberLevIS ValsLevIS AberLevQS ValsLevQS CalcS...
                Ber4PAMS PrintinEye SetCpSampZer ReceptorNoise EvmDBJ EvmRmsJ BW...
                EvmPerJ EvmMatRec EvmMatRef Nsamp RecFilBanPas VetSnr SendingUp ...
                TimeSys TxSymbAmos RxSymbAmos nRb NumFra ChanelEqualizer DifNsN ...
                CpLe OfdMod BerOFDM DatSiz ChanelEqualizer2 ChanelEqualizer3 Ofs...
                Tap3 ExtSam NuSaTs NuSaTu NuSaTg SaveToScatter Zt EleCar Ofc OBw...
                SelModTp OvSam NPOFEX NPPOF NPSTUf NPPB SelecGaus ControlVars ...
                GatingAtReception ZtC Te tta ff NumFraPar DmtMve RxSigOfdmNoEq ...
                RxSigOfdm TxSigOfdm nRbAgg DmtMinM DmtPas NuDmSuDi TapN NpEl ...
                ControlZt OveT MZ_Input_File NuPart NuSuSe freqGHz EyeToPlot ...
                UsingHermitian UseOfdmElect ChanAwgn SaveRxNotEq SaveRxEq SNR ...
                VetSnrIni VetSnrPass VerSnrEnd BerTheo CompEbn0 LegendThis ...
                ThisSimuBer ThisPlotCont VerThisMod ThisM Modulation DmtKvep...
                SaveBerEbn0Path MinTesNumb MinSavTim MzFilesGenerated...
                MzFilesDpGenerated BerToPlotOfdm
            
            if ~exist('CarrRecPowUp','var')
                CarrRecPowUp = 0;
            end
            if ~exist('CarrRecPowDo','var')
                CarrRecPowDo = 0;
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
            if CurrentTest > MinTesNumb
                if (((CurrentTest)*ceil(0.8*Nb))>(5/min(mean(b,1))))
                    break;
                end
            end
            baux = mean(b,1);
            bdis = [baux(berpos1);baux(berpos2)];
            display(bdis);
            display(sprintf(['Elapsed Time(s) = %f | Current Test = %d '...
                '| \nTested = %0.1f | Threshold = %0.2f | Decision = %0.2f '...
                '\n'],ElapsedTime,CurrentTest,VetSnr,((CurrentTest)*ceil(...
                0.8*Nb)),(5/min(mean(b,1)))));
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
        a=1;
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
a=a+1;

