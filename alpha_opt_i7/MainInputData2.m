%c
%c                                                       ..'Ž`'..'Ž`..'Ž`..                                                   
%c       File: MainInputDAta File With the Main Script Parameters 
%c(Main file to call each idividual parameters)
%c
%c     This main code is resposible to  load all necessary parameters of 
%c the MAIN script, just in this file changes should be done. Do not change
%c the MAIN file. 
%c
%c      
%c
%c                                           by P.Marciano LG
%c                                           28/10/2017
%c                                           last updated
%c                                           09/05/2018
%c                                           pablorafael.mcx@gmail.com
%c 
%c     References:
%c@article{essiambre2010capacity,
%c  title={Capacity limits of optical fiber networks},
%c  author={Essiambre, Ren{\'e}-Jean and Kramer, Gerhard and Winzer, Peter J and Foschini, Gerard J and Goebel, Bernhard},
%c  journal={Journal of Lightwave Technology},
%c  volume={28},
%c  number={4},
%c  pages={662--701},
%c  year={2010},
%c  publisher={IEEE}
%c}
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

%%                    OFC Generator Parameters
ThisPath = pwd;                                                            %Path for the folter where the OFC is saved
UofcsLocal = '\save_files\';                                               %Files name for the location of the OFC
if exist('CurrentOCS','var')                                               %This was not been used, this code works when different OFCS need to be tested and it was called by another script
    OfcName = [ThisPath UofcsLocal 'OCS_' num2str(OcsToTest(CurrentOCS))];
else
%     OfcName = [ThisPath UofcsLocal 'OCS1_08-Mar-2018'];                  %OFCS with 128 carriers NPPB 512 Nb 1024 fc 12.5e9   EDFA Power 40 db
    OfcName = [ThisPath UofcsLocal 'OCS_fc_125_16-May-2018'];              %OFCS with 128 carriers NPPB 512 Nb 1024 fc 12.5e9   EDFA Power 20 db
%     OfcName = [ThisPath UofcsLocal 'OCS_fc_250_27-Apr-2018'];            %OFCS with 128 carriers NPPB 512 Nb 1024 fc 25e9     EDFA Power 20 db
%     OfcName = [ThisPath UofcsLocal 'OCS_fc_15625_16-May-2018'];          %OFCS with 128 carriers NPPB 512 Nb 1024 fc 15.625e9 EDFA Power 40 db
%     OfcName = [ThisPath UofcsLocal 'OCS_fc_125_Long_28-May-2018'];              %OFCS with 128 carriers NPPB 512 Nb 1024 fc 12.5e9   EDFA Power 20 db
end

%%                 Control Variable
SelecSetUP    = 1;                                                         %Select the simulated setup 1 for original 0 for IFFT/OFFT at receptor
SendingDowStr = 1;                                                         %Set 1 to send data downstream 0 to not send.
SigAten       = 0;                                                         %Set 1 to atenuate the OFCS signal. Set to zero to turn off attenuation
SnrRef        = 0;                                                         %Set 1 to select the own signal for the noise reference at receptor. 
                                                                           %Set 0 for the pilote carrier as noise reference at receptor.
FFTSplit      = 1;                                                         %Set 1 to use the OFFT to split the income signal on receptor. Set 0 to 
                                                                           %use filter to split the income signal at receptor.
CarSNR        = -15;                                                         %Set the snr to be applied at receptor.
PastTime      = 0;                                                         %Storage for the time flow on simulation
% Atenuation    = 10^(4.0);                                                  %Set how much the OFCS will be attenuated
NumbOfTest    = 750-120;                                                      %The number of repetition that one test will be done
TestBundle    = 1;                                                         %Set if the scrip will be run once or many times
CurrentMedium = 0;                                                         %Sellection of channel 1 for Fiber other of B2B
PAM4Type      = [1 0];                                                     %Sellection the type of 4PAM that will be tested
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
ReceptorNoise = 1;
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
LevDefDpqsk = 0.13;                                                        %Default level for the DPSK modulation

%This will be the variable that will select the type of modulation used by
%this script. So far, this the options goes from 1 to 4, where:
%       1   -->    OOK
%       2   -->    4PAM
%       3   -->    DQPSK
%       4   -->    DPSK
% ThisModula = 1;                                            
UsedModula = [4];% 3 4];
UpStModula = [4];
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
    NumCarr      = 8;                                                    %Total number of carriers (DownStream and UpStream) - This variable must be even
    lambdac      = 1550e-9;                                                %Central wave length of the signal throughout the fiber
    FiberLength  = 270;                                                    %Setting the length of the fiber in this simulation
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
    else                                                                   %whereas, when the first carrier was not multiple of the number of carrier used
        ObsCarrPos  = [mod(RefCarr,NumCarr):NumCarr 1:mod(RefCarr,...
                                                               NumCarr)-1];%The first carrier is what left from the divission Reference carrier by 
                                                                           %the number of carriers used and this sequency goes untill the last carrier 
                                                                           %selected. Then, the following carrier will be from 1 untill the divission result minus 1
        ObsCarrUsed = ObsCarrPos;                                          %As in this simulation all carrier will be used the observed carriers is equal 
                                                                           %to the possition of the observed carriers
    end
%     SavingAT     = 'DpskCp122DemuxFiltOifftB2BOfftSnr-15_32';                                               
    SavingAT     = 'test_0';
end

%The relation between OSNR and SNR was taken from the first article 
%mentioned in this script.

AsePower      = 10^((0-30)/10);                                            %Ase Power(dB)
LightSpeed    = 3e8;                                                       %Ligth Speed (ms)
planck        = 6.626069e-34;                                              %Plank constant
Bref          = 12.5e9;%0.1e-9;                                            %Reference Bandwidth for measuring the OSNR (0.1 nm at 1550 nm)
BrefAll       = NumCarr*12.5e9;%0.1e-9;                                    %Reference Bandwidth for measuring the OSNR (0.1 nm at 1550 nm)
snr           = -12;                                                       %Signal to Noite Ratio of the eletrical signal that will modulate the optical carrier in dB
osnr          = snr + 10*log10(10e9) - 10*log10(2) - 10*log10(Bref);       %Optical Signal to Noite Ratio of the optical signal after modulation in dB
snradded      = -15;
osnrf         = snradded + 10*log10(10e9) - 10*log10(2) -10*log10(BrefAll);%Optical Signal to Noite Ratio of the optical signal after fiber in dB 
snrp          = 0;
osnrp         = snrp + 10*log10(10e9) - 10*log10(2) - 10*log10(Bref);      %Optical Signal to Noite Ratio of the optical signal before photodetector in dB
a=0;

