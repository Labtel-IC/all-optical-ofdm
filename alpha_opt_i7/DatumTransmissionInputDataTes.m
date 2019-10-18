%c
%c                                                       ..'Ž`'..'Ž`..'Ž`..                                                   
%c       File: DatumTransmissionInputDAta File With the Main Script Parameters 
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

%%       General Parameter for the simulation
ON          =  1; %DO NOT CHANGE THIS PARAMETER
OFF         =  0; %DO NOT CHANGE THIS PARAMETER
Ploting     = OFF;%Use On or Off to say if the figures should be plot or not
PlotingThis = ON;%Use On or Off to say if the figures should be plot or not
Selecting   = OFF;%Use to filter or not the carriers to be transmited
ModAll      = OFF;%Set if all carriers will be modulate at once or not
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
switch CurrentModula
    case 1
        Modulation = 'NoN';
    case 2
        Modulation = '4PAM';
    case 3
        Modulation = 'DQPSK';
    case 4
        Modulation = 'DPSK';
    otherwise
        Modulation = 'NoN';
end
                                     
%%          Variables for the Medium Fiber
NumCarr      = 8;                                                          %Total number of carriers (DownStream and UpStream) - This variable must be even
lambdac      = 1550e-9;                                                    %Central wave length of the signal throughout the fiber
FiberLength  = 50;                                                         %Length of the optical fiber [km]
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
    [FiberDelay] = ProgationDelay(NumCarr,t,lambdac,T,FiberLength,0,f,...
                                                                  Ploting);
end
%%         Collecting the impulse response
% PulseSig          = Eout;%zeros(1,length(t));
% CW                = ones(1,length(t));
% PulseSig(1)       = 1;
% U.U1t             = PulseSig;
% U.U2t             = exp(-1j*pi).*PulseSig;
% [PulseEout,~]     = Mach_Zehnder_Modulator_simplificado(t,CW,U,2);
% [PulseRec,~,~]    = Fibra_Monomodo1(t,PulseSig,lambdac,T,FiberLength,0,f);
% [PulseAux,~,~]    = OpticalFFTN(t,T,nextpow2(NumCarr),PulseRec);
% PulseResp         = zeros(size(PulseAux,1),size(PulseAux,2));
% 
% for kk=1:size(PulseAux,1)
%     PulseResp(kk,:) = abs(PulseAux(kk,:)).^2;
% end
%%          Variables for different types of Modulations
switch Modulation
    case 'DPSK' %Tested to fiberlength of 80 -> worked 
        %%      Variables for the DPSK Datum Modification
        %This variable allow the user to chose to send an random 
        %information or just a high pulse within the midle of the datastram
        SendingData   = 1;                                                 %Selecting the nature of the information to be send                                                
        NbDPSK        = Nb;                                                %The total number of bits to be transmited
        DatGai        = 1;                                                 %Seting the maximum range for the eletrical signal
        JusVal        = [zeros(1,10) ones(1,4) zeros(1,10)];               %Value that will be load into the justificaiton space
        JusValEnd     = JusVal;                                            %Value that will be load into the justificaiton space
        JusPos        = length(JusVal);
        BWD           = 2*12.5e9;                                          %Badnwidth for the bit conformation
        CenFeq        = 0;                                                 %Center frequency for the bit conformation
        FiltOrd       = 1;                                                 %The order for the filter which will conform data
        MZ_Input_File = 3;                                                 %Telling the script witch file will be used for modulation
        %Variables for the Selection of carriers - Wave Divided.
        fin           = 12.5e9;                                            %The frequency of the first carrier of the optical COMB
        FBWD          = 1.2*Rb;                                            %The filter BandWidth to split each individual carrier
        Order         = 5;                                                 %The order of the filter used to split each carrier
    case 'DQPSK'%Tested to fiberlength of 40 -> work in progress
        %%      Variables for the DQPSK Datum Modification
        %This variable allow the user to chose to send an random 
        %information or just a high pulse within the midle of the datastram
        SendingData = 1;                                                   %Selecting the nature of the information to be send
        NbDQPSK     = 2*Nb;                                                %The total number of bits to be transmited
        JusVal      = [zeros(1,20) ones(1,8) zeros(1,20)];                 %Value that will be load into the justificaiton space
        JusValEnd   = JusVal;                                              %Value that will be load into the justificaiton space
        JusPos      = length(JusVal);
        BWD         = 4*12.5e9;                                            %Badnwidth for the bit conformation
        CenFeq      = 0;                                                   %Center frequency for the bit conformation
        FiltOrd     = 1;                                                   %The order for the filter which will conform data
        DatGai      = 1;                                                   %Gain to increase the modulatin signal amplitude
        V0          = 3.8;
        Vpi         = 3.8;
        %Variables for the Selection of carriers - Wave Divided.
        fin         = 12.5e9;                                              %The frequency of the first carrier of the optical COMB
        FBWD        = 1.2*Rb;                                              %The filter BandWidth to split each individual carrier
        Order       = 5;                                                   %The order of the filter used to split each carrier
    case '4PAM'%Tested to fiberlength of 1 -> work in progress 
        %%      Variables for the 4PAM Datum Modification
        %This variable allow the user to chose to send an random 
        %information or just a high pulse within the midle of the datastram
        SendingData = 1;                                                   %Selecting the nature of the information to be send                                                
        Nb4Pam      = 2*Nb;                                                %The total number of bits to be transmited
        %This variable controls if the electrical signal will vary from 
        %0[V] to max[V] (signal unipolar) or from -max[V] to max[V] 
        %(signal bipolar).
        Polirized     = 0; 
        MaxAmp4PAM    = U_pi1/8;                                           %Seting the maximum range for the eletrical signal
        JusVal        = [zeros(1,20) 1 0 1 0 1 0 1 0 zeros(1,20)];         %
        JusValEnd     = JusVal;                                            %Value that will be load into the justificaiton space
        JusPos        = length(JusVal);                                    %
        BWD           = 2.0*12.5e9;                                        %Badnwidth for the bit conformation
        CenFeq        = 0;                                                 %Center frequency for the bit conformation
        FiltOrd       = 3;                                                 %The order for the filter which will conform data
        MZ_Input_File = 4;                                                 %Telling the script witch file will be used for modulation
        U_pi2         = 3.8;
        Vbias         = -1.9;
        %Variables for the Selection of carriers - Wave Divided.
        fin           = 12.5e9;                                            %The frequency of the first carrier of the optical COMB
        FBWD          = 1.2*Rb;                                            %The filter BandWidth to split each individual carrier
        Order         = 5;                                                 %The order of the filter used to split each carrier
    otherwise %Tested to fiberlength of 80 -> worked 
        %%      Variables for the OOK Datum Modification
        AdjusData     = 0;                                                 %Selecting the nature of the information to be send
        JusVal        = [zeros(1,10) ones(1,4) zeros(1,10)];               %Value that will be load into the justificaiton space
        JusValEnd     = [zeros(1,10) ones(1,4) zeros(1,10)];               %Value that will be load into the justificaiton space
        JusLen        = length(JusVal);                                    %Length for the information justification
        DatMin        = 0;                                                 %Minimum value that the Data can assume
        DatMax        = 1;                                                 %Maximum value that the DAta can assume
        NrzMin        = -1;                                                %Minimum value for the NRZ bit
        NrzMax        = 1;                                                 %Maximum value for the NRZ bit
        IQAmp         = 1;                                                 %Amplitude of the signal IQ
        DatGai        = 0.6;                                               %Gain to increase the modulatin signal amplitude
        BWD           = 11e9;                                              %Badnwidth for the bit conformation
        CenFeq        = 0;                                                 %Center frequency for the bit conformation
        FiltOrd       = 1;                                                 %The order for the filter which will conform data
        MZ_Input_File = 2;                                                 %Telling the script witch file will be used for modulation
        %Variables for the Selection of carriers - Wave Divided.
        fin           = 12.5e9;                                            %The frequency of the first carrier of the optical COMB
        FBWD          = 1.2*Rb;                                            %The filter BandWidth to split each individual carrier
        Order         = 5;                                                 %The order of the filter used to split each carrier
end
%%      Filter to Select the Carriers
[ SelecFilt ] = FiltroGaussiano(f,18*fc,10*fc,9);
% [ SelecFilt ] = Filtro_Retangular( 34*fc,16*fc,f);
SelecFilt = fftshift(SelecFilt);
















