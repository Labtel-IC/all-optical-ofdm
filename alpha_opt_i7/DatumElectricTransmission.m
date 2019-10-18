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
MZ_Input_File_Dp = 1;                                                      %Telling the script witch file will be used for modulation
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
% if CurrentMedium
%     if ~exist('FiberDelay','var')
%         [FiberDelay,PowRef] = ProgationDelay(RefCarr,NumCarr,t,lambdac,T,...
%                                                   FiberLength,0,f,PlotingThis);
%     end
% end
%% ########################################################################
% ########################################################################
% ########################################################################
% ########################################################################
% ########################################################################
% ########################################################################
%% ########################################################################
OfdmInputData;

%% ########################################################################
% ########################################################################
% ########################################################################
% ########################################################################
% ########################################################################
% ########################################################################
% ########################################################################
%% ########################################################################
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