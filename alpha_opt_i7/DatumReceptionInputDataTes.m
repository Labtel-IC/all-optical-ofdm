%c
%c                                                       ..'�`'..'�`..'�`..                                                   
%c       File: DatumReceptionInputDAta File With the Main Script Parameters 
%c(Main fail to call each idividual parameters)
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
%c                                           23/12/2017
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
Ploting       = OFF; %It is used to select an ploting.
PlotingThis   = ON; %It is used to select an especific ploting among others
Elucidation   = OFF;%It is used to enable the OFFT example
UsePhotDiod   = OFF;%To select or not the usage of the photodiod
BitFilBanWid  = 0.9*Rb;%Badnwidth for the bit conformation
BitFiltCenFre = 0;  %Center frequency for the bit conformation
BitFilOrd     = 3;  %The order for the filter which will conform data
%%    Constants to controu the Optical FFT
%The maximum number of carriers that the actual Optical FFT can implement
%is 32. Which means, a FFT of order 5. Further worked is needed to change
%this optical FFT to order 'n'.
% if NumCarr>33                                                              %Cheking the maximum number of carriers given
%     NumCarr = 32;                                                          %If the number of carrier were biger than the maximum value the 
%     MaxNumCarr = 32;                                                       %variables here presented are set to its maximum corespondent value
    %Maximum Number of sub carriers. 32 is the maximum allowed, I did not 
    %adjust the the FFT to work with more than 32
%     MaxNumStag = ceil(log2(MaxNumCarr)); 
% else
    %If the number of carriers is less than 32, this scrip will set the
    %OFFT accordingly.
    %The order of this FFT is defined by the number of carriers. But the 
    %number of carriers must be one value of the power of 2
    MaxNumStag = nextpow2(NumCarr);                                        
    MaxNumCarr = 2^MaxNumStag; 
% end

%As result of the OFFT process the carriers get mingled at the output ports
%of this sistem. Thus, those vector were created as a map for the
%simulation. Where the possition of the vector is the output of the OFFT
%and the value within tha possition is the respective carrier filtered. The
%folloing ilustration better shows how each carrier get mingled in this
%process:
%               Ein : Is the input fild
%               MZI : Is the Mach-Zehnder Interferometre
%               N�  : The number are each carrier from 1 to 32 within Ein
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
%         VetThisCarr = [1 17 9 25 5 21 13 29 3 19 11 27 7 23 15 31 2 18 ...
%                                    10 26 6 22 14 30 4 20 12 28 8 24 16 32];
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
%                                 ' cange the number Maximum of carriers.']);
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
    case 'DPSK'
        %%             Syncing DPSK signal
        
        BWD2          = 1*12.5e9;                                          %Badnwidth for the bit conformation
        CenFeq2       = 0;                                                 %Center frequency for the bit conformation
        FiltOrd2      = 1;                                                 %The order for the filter which will conform data
        SyncSymb      = [zeros(1,10) ones(1,4) zeros(1,10)];              %Generating the sync-word at the beginning of the data stream
        IniSyncPos    = 1+5*NPPB;
        SyncPos       = 19*NPPB;                                           %Seting the possition where the sync-symble is
        SyncSymbEnd   = SyncSymb;                                          %Generating the sync-word at the end of the data stream
        SyncPeriod    = length(SyncSymb);                                  %Seting how many symble periods will be within the sync-word
        AuxSyncLength = SyncPeriod*NPPB;                                   %auxiliar variable to set where is the sync-symble
        tsync         = linspace(0,SyncPeriod*T,NPPB*SyncPeriod);          %create a time vetor for the sync-word to help on the ploting of this signal
        fsync         = time2freq(tsync);                                  %creating the the frequency vector of the sync-word for hereafter filtering
        BWD           = 2*12.5e9;                                          %Badnwidth for the bit conformation
        CenFeq        = 0;                                                 %Center frequency for the bit conformation
        FiltOrd       = 1;                                                 %The order for the filter which will conform data
        BitFilBanWid  = 10e9;                                              %Badnwidth for the bit conformation
        BitFiltCenFre = 0;                                                 %Center frequency for the bit conformation
        BitFilOrd     = 1;                                                 %The order for the filter which will conform data
        [BitFilt,~]   = FiltroGaussiano(fsync,BWD,CenFeq,FiltOrd);         %Creating filter for conformation of the input information
        BitFilt       = fftshift(BitFilt);
        
        SyncSymb      = rectpulse(SyncSymb,NPPB);
        SyncSymbEnd   = rectpulse(SyncSymbEnd,NPPB);
        SencondAdjust = OFF;                                               %Setting whether is needed to sync the information by the Syncronism Symble at the End
        LowSymPer     = 1;
        UpeSymPer     = 1;
    case 'DQPSK'
        %%             Syncing DQPSK signal
        
        BWD2          = 1.2*12.5e9;                                          %Badnwidth for the bit conformation
        CenFeq2       = 0;                                                 %Center frequency for the bit conformation
        FiltOrd2      = 5;                                                 %The order for the filter which will conform data
        SyncSymb      = [zeros(1,10) ones(1,4) zeros(1,10)];               %Generating the sync-word at the beginning of the data stream
        IniSyncPos    = 1+5*NPPB;
        SyncPos       = 19*NPPB;                                           %Seting the possition where the sync-symble is
        SyncSymbEnd   = SyncSymb;                                          %Generating the sync-word at the end of the data stream
        SyncPeriod    = length(SyncSymb);                                  %Seting how many symble periods will be within the sync-word
        AuxSyncLength = SyncPeriod*NPPB;                                   %auxiliar variable to set where is the sync-symble
        tsync         = linspace(0,SyncPeriod*T,NPPB*SyncPeriod);          %create a time vetor for the sync-word to help on the ploting of this signal
        fsync         = time2freq(tsync);                                  %creating the the frequency vector of the sync-word for hereafter filtering
        BWD           = 2*12.5e9;                                          %Badnwidth for the bit conformation
        CenFeq        = 0;                                                 %Center frequency for the bit conformation
        FiltOrd       = 1;                                                 %The order for the filter which will conform data
        BitFilBanWid  = 10e9;                                              %Badnwidth for the bit conformation
        BitFiltCenFre = 0;                                                 %Center frequency for the bit conformation
        BitFilOrd     = 1;                                                 %The order for the filter which will conform data
        [BitFilt,~]   = FiltroGaussiano(fsync,BWD,CenFeq,FiltOrd);         %Creating filter for conformation of the input information
        BitFilt       = fftshift(BitFilt);
        
        SyncSymb      = rectpulse(SyncSymb,NPPB);
        SyncSymbEnd   = rectpulse(SyncSymbEnd,NPPB);
        SencondAdjust = ON;                                                %Setting whether is needed to sync the information by the Syncronism Symble at the End
        LowSymPer     = 1;
        UpeSymPer     = 1;
    case '4PAM'
        %%             Syncing 4PAM signal
        SyncSymb      = [zeros(1,10) ones(1,4) zeros(1,10)];               %Generating the sync-word at the beginning of the data stream
        IniSyncPos    = 1+5*NPPB;
        SyncPos       = 19*NPPB;                                           %Seting the possition where the sync-symble is
        SyncSymbEnd   = SyncSymb;                                          %Generating the sync-word at the end of the data stream
        SyncPeriod    = length(SyncSymb);                                  %Seting how many symble periods will be within the sync-word
        AuxSyncLength = SyncPeriod*NPPB;                                   %auxiliar variable to set where is the sync-symble
        tsync         = linspace(0,SyncPeriod*T,NPPB*SyncPeriod);          %create a time vetor for the sync-word to help on the ploting of this signal
        fsync         = time2freq(tsync);                                  %creating the the frequency vector of the sync-word for hereafter filtering
        BWD           = 2*12.5e9;                                          %Badnwidth for the bit conformation
        CenFeq        = 0;                                                 %Center frequency for the bit conformation
        FiltOrd       = 1;                                                 %The order for the filter which will conform data
        BWD2          = 2.0*12.5e9;                                        %Badnwidth for the bit conformation
        CenFeq2       = 0;                                                 %Center frequency for the bit conformation
        FiltOrd2      = 3;                                                 %The order for the filter which will conform data
        [BitFilt,~]   = FiltroGaussiano(fsync,BWD,CenFeq,FiltOrd);         %Creating filter for conformation of the input information
        BitFilt       = fftshift(BitFilt);
        SyncSymb      = rectpulse(SyncSymb,NPPB);                          %Resampling the sync-word accordingly with the current system
        SyncSymbEnd   = rectpulse(SyncSymbEnd,NPPB);                       %Resampling the sync-word accordingly with the current system
        IntervalStep  = 40;
        MinDist       = IntervalStep/6;
        %Because of reasons, the real world does not have components
        %capable instantly change the voltage level. Thus, to emulate this 
        %behaviour the eletrical PAM signal pass through a gaussian filter 
        %that will remove the frequencies of higher order, which will 
        %result in small slope in the level variation
        SyncSymb      = ifft(fft(SyncSymb).*BitFilt);
        SyncSymbEnd   = ifft(fft(SyncSymbEnd).*BitFilt);
        SencondAdjust = ON;                                                %Setting whether is needed to sync the information by the Syncronism Symble at the End
        LowSymPer     = 0.9;
        UpeSymPer     = 1.1;
    otherwise
        %%             Syncing OOK signal
        
        SyncSymb      = [zeros(1,10) ones(1,4) zeros(1,10)];               %Generating the sync-word at the beginning of the data stream
        SyncPos       = 19*NPPB;                                           %Seting the possition where the sync-symble is
        IniSyncPos    = 1+5*NPPB;
        SyncSymbEnd   = [zeros(1,10) ones(1,4) zeros(1,10)];               %Generating the sync-word at the end of the data stream
        SyncPeriod    = length(SyncSymb);                                  %Seting how many symble periods will be within the sync-word
        AuxSyncLength = SyncPeriod*NPPB;                                   %auxiliar variable to set where is the sync-symble
        tsync         = linspace(0,SyncPeriod*T,NPPB*SyncPeriod);          %create a time vetor for the sync-word to help on the ploting of this signal
        fsync         = time2freq(tsync);                                  %creating the the frequency vector of the sync-word for hereafter filtering
        BWD           = 1.2*12.5e9;                                        %Badnwidth for the bit conformation
        CenFeq        = 0;                                                 %Center frequency for the bit conformation
        FiltOrd       = 3;                                                 %The order for the filter which will conform data
        BitFilBanWid  = 10e9;                                              %Badnwidth for the bit conformation
        BitFiltCenFre = 0;                                                 %Center frequency for the bit conformation
        BitFilOrd     = 1;                                                 %The order for the filter which will conform data
        [BitFilt,~]   = FiltroGaussiano(fsync,BitFilBanWid,BitFiltCenFre...
                                                               ,BitFilOrd);%Creating filter for conformation of the input information
        BitFilt       = fftshift(BitFilt);
        SyncSymb      = rectpulse(SyncSymb,NPPB);                          %Resampling the sync-word accordingly with the current system
        SyncSymbEnd   = rectpulse(SyncSymbEnd,NPPB);                       %Resampling the sync-word accordingly with the current system
        
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




%%   Control Variables for the All Optical FFT Example
% Time dellay used on the example
%Stage 1
Del0   = T/2;%Time delay at the frist stage
%Stage 2
Del1a  = T/4;%Time delay at the second stage
Del1b  = T/4;%Time delay at the second stage
%Stage 3
Del2aa = T/8;%Time delay at the final stage
Del2ab = T/8;%Time delay at the final stage
Del2ba = T/8;%Time delay at the final stage
Del2bb = T/8;%Time delay at the final stage

% Phase dellay used on the example given
%Stage 1
Pha0   = (0)*(pi/180); %Phase dalay at the first stage
%Stage 2
Pha1a  = (90)*(pi/180);%Phase dely at the second stage
Pha1b  = 0;            %Phase delay at the second stage
%Stage 3
Pha2aa = 3*pi/4;       %Phase delay at final stage
Pha2ab = pi/4;         %Phase delay at final stage
Pha2ba = pi/2;         %Phase delay at final stage
Pha2bb = 0;            %Phase delay at final stage

% Variables for the optical coupler
% it was fundamental to use a teoretical model for the couple. Just add or
% subtract results will not work. The model used at this script was taken
% from Optical Networks by Ramaswami
kl = (135)*pi/180;    %Coupling Constant. It is the argument of the function
MutCoef = 1;          %Multiplier, it is used to sete the polarization 
ExpCoef = exp(-1j*1); %Exponential that multiplies all frequency componente

%%      Filters For Channel Selection
% Finaly, it was created filter responsible to better receive the signal
% and eliminate any non  linearities.O
SelecBand = 2*fc; %Bandwidth of the filter
CentFreq  = fc;   %Center Frequency
SelFilOrd = 5;    %rden for the selected Filter

[ SelecFilt3aaa ] = FiltroGaussiano(f,SelecBand,CentFreq,SelFilOrd);
% [ SelecFilt1 ] = Filtro_Retangular( SelecBand,CentFreq,fb);
SelecFilt3aaa = fftshift(SelecFilt3aaa);

[ SelecFilt3baa ] = FiltroGaussiano(f,SelecBand,2*CentFreq,SelFilOrd);
% [ SelecFilt2 ] = Filtro_Retangular( SelecBand,2*CentFreq,fb);
SelecFilt3baa = fftshift(SelecFilt3baa);

[ SelecFilt3aba ] = FiltroGaussiano(f,SelecBand,3*CentFreq,SelFilOrd);
% [ SelecFilt3 ] = Filtro_Retangular( SelecBand,3*CentFreq,fb);
SelecFilt3aba = fftshift(SelecFilt3aba);

[ SelecFilt3bba ] = FiltroGaussiano(f,SelecBand,4*CentFreq,SelFilOrd);
% [ SelecFilt4 ] = Filtro_Retangular( SelecBand,4*CentFreq,fb);
SelecFilt3bba = fftshift(SelecFilt3bba);

[ SelecFilt3aab ] = FiltroGaussiano(f,SelecBand,5*CentFreq,SelFilOrd);
% [ SelecFilt5 ] = Filtro_Retangular( SelecBand,5*CentFreq,fb);
SelecFilt3aab = fftshift(SelecFilt3aab);

[ SelecFilt3bab ] = FiltroGaussiano(f,SelecBand,6*CentFreq,SelFilOrd);
% [ SelecFilt6 ] = Filtro_Retangular( SelecBand,6*CentFreq,fb);
SelecFilt3bab = fftshift(SelecFilt3bab);

[ SelecFilt3abb ] = FiltroGaussiano(f,SelecBand,7*CentFreq,SelFilOrd);
% [ SelecFilt7 ] = Filtro_Retangular( SelecBand,7*CentFreq,fb);
SelecFilt3abb = fftshift(SelecFilt3abb);

[ SelecFilt3bbb ] = FiltroGaussiano(f,SelecBand,8*CentFreq,SelFilOrd);
% [ SelecFilt8 ] = Filtro_Retangular( SelecBand,8*CentFreq,fb);
SelecFilt3bbb = fftshift(SelecFilt3bbb);



