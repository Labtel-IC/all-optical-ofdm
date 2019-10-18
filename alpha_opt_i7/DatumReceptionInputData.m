%c
%c                                                       ..'Ž`'..'Ž`..'Ž`..                                                   
%c       File: DatumReceptionInputDAta File With the Main Script Parameters 
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
%c                                           Last UpDate
%c                                           02/04/2018
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       Variables for Syncronism

%%    Control Variables of this program
ON            =  1; %DO NOT CHANGE THIS PARAMETER
OFF           =  0; %DO NOT CHANGE THIS PARAMETER
Ploting       = ON; %It is used to select an ploting.
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




