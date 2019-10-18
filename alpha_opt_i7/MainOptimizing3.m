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
%Checking if this main file was not called by another script
%%
close all;clc;clear;tic;a=1;
MainOptmInputData;
%clear MzFilesGenerated MzFilesDpGenerated
tic;                                                                       %Starting timer to check the simulation time
for VetSnr=VetSnrIni:VetSnrPass:VerSnrEnd
    clear BerOOK BerPAM4 BerDQPSK BerDPSK BerOOKS BerPAM4S BerDQPSKS ...
        BerDPSKS BerOFDM ChanelEqualizer
    
    FiberLength  = VetSnr;                                                     %Setting the length of the fiber in this simulation
    CarSNR       = VetSnr;                                                    %Set the snr to be applied at receptor.
    switch Modulation
        case 'OFDM'
%             SNR = CarSNR + 10*log10(log2(max(DmtMve))) - 10*log10(1*Tu*OBw);
            SNR = CarSNR + 10*log10(log2(max(DmtMve))) - 10*log10(OvSam);
            TxSymbAmos = zeros(NumbOfTest,NumCarr,NumFra*length(DmtMve));
        case 'DPSK' %Tested to fiberlength of 80 -> worked
            if TimeSys==1
                SNR = CarSNR + 10*log10(1) - 10*log10(1*T*BWD2);
            else
                SNR = CarSNR + 10*log10(1) - 10*log10(1*t(2*NumAmosCP+NPPB)*BWD2);
            end
            TxSymbAmos = zeros(NumbOfTest,NumCarr,Nb - NumBitDesc);
        case 'DQPSK'%Tested to fiberlength of 40 -> work in progress
            if TimeSys==1
                SNR = CarSNR + 10*log10(2) - 10*log10(1*T*BWD2);
            else
                SNR = CarSNR + 10*log10(2) - 10*log10(1*t(2*NumAmosCP+NPPB)*BWD2);
            end
            TxSymbAmos = zeros(NumbOfTest,NumCarr,Nb - NumBitDesc);
        case '4PAM'%Tested to fiberlength of 1 -> work in progress
            if TimeSys==1
                SNR = CarSNR + 10*log10(2) - 10*log10(1*T*BWD2);
            else
                SNR = CarSNR + 10*log10(2) - 10*log10(1*t(2*NumAmosCP+NPPB)*BWD2);
            end
            TxSymbAmos = zeros(NumbOfTest,NumCarr,Nb - NumBitDesc);
        otherwise %Tested to fiberlength of 80 -> worked
            if TimeSys==1
                SNR = CarSNR + 10*log10(1) - 10*log10(1*T*BWD2);
            else
                SNR = CarSNR + 10*log10(1) - 10*log10(1*t(2*NumAmosCP+NPPB)*BWD2);
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
            [FiberDelay,PowRef] = ProgationDelay(RefCarr,NumCarr,t,lambdac,T,...
                FiberLength,0,f,Ploting);
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
            RxSigOfdmNoEq = zeros(Ns*NumFra,CurTesExtSiz);
            RxSigOfdm = zeros(Ns*NumFra,CurTesExtSiz);
            EvmMatRef = zeros(1,CurTesExtSiz*length(DmtMve)*NumFra);
        else
            RxSigOfdmNoEq = zeros(NumbOfTest,NumCarr,Ns*NumFra);
            RxSigOfdm = zeros(NumbOfTest,NumCarr,Ns*NumFra);
            EvmMatRef = zeros(NumCarr,CurTesExtSiz*length(DmtMve)*NumFra);
        end
        SigRecep4 = zeros(length(DmtMve),NumFra);
        
        TxSigOfdm = zeros(NumbOfTest,NumCarr,NumFra*length(DmtMve));
        BerToPlotOfdm = zeros(NumbOfTest,NumCarr,length(DmtMve)*NumFra);
    else
        RxSigOfdmNoEq = 0;
        SigRecep4 = 0;
        RxSigOfdm = 0;
        TxSigOfdm = 0;
        BerToPlotOfdm = 0;
        EvmMatRef = zeros(NumCarr,Nb - NumBitDesc);
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
    VetElecPower = zeros(NumbOfTest,NumCarr);
    VetOptiPower = zeros(NumbOfTest,NumCarr);
    AberLevS = zeros(NumbOfTest,NumCarr);
    ValsLevS = zeros(NumbOfTest,NumCarr);
    EyeToPlot = zeros(NumbOfTest,Nb*NPPB);
    BerDPSK = zeros(NumbOfTest,NumCarr);
    BerDPSKS = zeros(NumbOfTest,NumCarr);
    AberLev = zeros(NumbOfTest,NumCarr);
    ValsLev = zeros(NumbOfTest,NumCarr);
    RxSymbAmos = zeros(NumbOfTest,NumCarr,Nb - NumBitDesc);
    VetElecPowerI = zeros(NumbOfTest,NumCarr);
    VetElecPowerQ = zeros(NumbOfTest,NumCarr);
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
    VetElecPowerUp = zeros(NumbOfTest,NumCarr);
    VetOptiPowerUp = zeros(NumbOfTest,NumCarr);
    AberLevUp = zeros(NumbOfTest,NumCarr);
    ValsLevUp = zeros(NumbOfTest,NumCarr);
    SaveRxNotEq=[];
    SaveRxEq=[];
    CurrentModula=ThisModula;
    for CurrentTest=1:CurTesPas:ThisTest
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
            EletricOfdmTransm;
        else
            %% DatumTransmission;
            DownStreamOptm;
            
        end
        if ~exist('ChanelEqualizer','var')
            if exist('OfdMod','var')
				
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
            ElectricReceiverOptm;
           
        else
            DownStreamREceiverOptm;
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
                    EoutGpu = gpuArray(EoutRec);
                    [EoutAux1Gpu,~,VetThisCarrGpu]=OpticalFFTG(fgpu,gpuArray(T),gpuArray(MaxNumStag),EoutGpu);
                    EoutAux1 = gather(EoutAux1Gpu);
                    VetThisCarr = gather(VetThisCarrGpu);
                    clear EoutGpu EoutAux1Gpu;
                else
                    [EoutAux1,~,VetThisCarr]=OpticalFFTN(f,T,MaxNumStag,EoutRec);
                end
                clear EoutRec
                % PrintInfo(Ploting*12,f,db(abs((fftshift(fft(EoutRec))))/length(EoutRec)),EoutAux1,RefCarr*fc);
                UpStreamTransmiterOptm;
                
                UpStreamReceptionOptm;
            end
        end
        %%              Saving
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
        a=a+1;
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
















