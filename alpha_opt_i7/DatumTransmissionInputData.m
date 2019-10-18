%c
%c                                                       ..'Ž`'..'Ž`..'Ž`..                                                   
%c File: DatumTransmissionInputDAta File With the Main Script Parameters 
%c (Main file to call each idividual parameters)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%       General Parameter for the simulation
ON          =  1;                                                          %DO NOT CHANGE THIS PARAMETER
OFF         =  0;                                                          %DO NOT CHANGE THIS PARAMETER
Ploting     = ON;                                                         %Use On or Off to say if the figures should be plot or not
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
MZ_Input_Data;                                                             %Loading the basic data to be saved
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
    Make_MZ_Input_Files_Simp;                                                  %Create the files for the MZM input data
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
MZ_Input_Data_DP;                                                          %Loading the basic data to be saved
if ~exist('MzFilesDpGenerated','var')
    Make_MZ_Input_Files_DP;                                                    %Create the files for the MZM input data
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
        OfdmInputData;
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
    ControlVars.FeqC          = FeqC;
    ControlVars.L             = L;
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











