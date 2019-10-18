%c
%c                                                       ..'Ž`'..'Ž`..'Ž`..                                                   
%c File: DatumUpstreamInputDAta File With the Main Script Parameters 
%c (Main file to call each idividual parameters)
%c
%c     This main code is resposible to  load all necessary parameters of 
%c the MAIN script, just in this file changes should be done. Do not change
%c the MAIN file. 
%c
%c      
%c
%c                                           by P.Marciano LG
%c                                           05/04/2018
%c                                           Last UpDate
%c                                           05/04/2018
%c                                           pablorafael.mcx@gmail.com
%c 
%c     References:
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
ON          =  1; %DO NOT CHANGE THIS PARAMETER
OFF         =  0; %DO NOT CHANGE THIS PARAMETER
Ploting     = OFF;%Use On or Off to say if the figures should be plot or not
PlotingThis = OFF;%Use On or Off to say if the figures should be plot or not
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
















