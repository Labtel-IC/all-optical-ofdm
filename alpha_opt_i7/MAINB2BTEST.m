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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Checking if this main file was not called by another script
if ~exist('UsingAllCarriers','var')
    clear ; clc; close all;
end
tic;                                                                       %Starting timer to check the simulation time
for VetSnr=31:40
    clear BerOOK BerPAM4 BerDQPSK BerDPSK BerOOKS BerPAM4S BerDQPSKS BerDPSKS
MainInputData;
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
for CurrentModula=ThisModula
    for CurrentTest=1:ThisTest
        if TestBundle                                                      %When many tests will be run, the OFCS will loaded from a file
            load([OfcName '.mat']);                                        %Loading a OFC from memory.
        else                                                               %When it will be just a general test it will ask if user wants a new OFCS
            Check = input(['Do you want to load a new OFC?\n If not pre'...%Fristly the script check if the user want to 
            'ss Enter to continue. The OFC will be loaded from file.\n']...
                                                                     ,'s');%load an OFC from file or to generate a new one
            if ~isempty(Check)
                GenerateOFC;                                               %Generate a new OFC if one wishes so. 
            else
                load([OfcName '.mat']);                                    %Otherwise load a OFC from memory.
            end
        end
        % Ploting the OFC that will be used in this simulation
        %opticalnoiseadded;
        Ploting = 0;
        PrintInfo(Ploting*1,f,Eout);
        %xmax = max(20*log10(abs(fftshift(fft(Eout)./length(Eout)))));
        %axis([-1e11 2e12 xmax-60 xmax+0.1]);
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
        
        DatumTransmissionB2bTest;

        %%              Datum Reception
        % The following step is to recover the information transmited 
        %throughout an medium. It is the most important part of this work, 
        %because it will deal with an All-Optical FFT implementation. It 
        %aims to preciselly separate each sub carrier that compose the OFDM 
        %Symble. It will apply methods previouly presented in the 
        %literature. Looking to do all the FFT processing at the optical 
        %domain.
        
        if SendingDowStr
            DatumReceptionB2bTest;
        else
            DatumReceptionInputData;
            [EoutAux1,~,VetThisCarr]=OpticalFFTN(t,T,MaxNumStag,EoutRec);
        end
        
        %%              Datum Upstream
        %This section is resposible to generate and to send the information
        %from the user to the OLT.
        DatumUpstreamB2bTest;
        
        %%                Datum Upstream Reception
        %This section is responsible to receive the information from users
        DatumUpstrRecB2bTest;
        %%              Saving
        
        clearvars -except ThisPath UofcsLocal OfcName NumbOfTest PowRef ...
        PAM4Type ThisTest TestBundle CurrentMedium UsedModula UsedModula...
        Ber4PAM SnrRef CarSNR FiberDelay ValsLev2 ValsLev contcp ExRa nn...
        PastTime SigAten CurrentModula Atenuation VetOptiPower AberLev ...
        VetElecPower FFTSplit ThisModula SendingDowStr VetElecPowerQ Tcp...
        VetElecPowerI BerDPSK BerDQPSK CurrentTest VetOptiPowerUp mm Nb...
        AberLevUp VetElecPowerUp ValsLevUp EvmRms EvmDB EvmPer osnr snr...
        DecLevDef1 DecLevDef2 DecLevDef3 UpStModula EvmRmsUp EvmDBUp ...
        EvmPerUp UsingAllCarriers ObsCarrUsed ll CarrToSendInt3 Bref...
        CarrToSend6 CarrToSend5 CarrToSend4 CarrToSend3 CarrToSendInt2 ...
        CarrToSend2 ObsCarr CarrToSendInt CarrToSend CarrEvenPos planck...
        CarrOddPos CarrPosInt ObsCarrPos NumCarrUsed VetNumCar VetOpen ...
        FiberLength lambdac NumCarr CarrUsed SavingAT RefCarr VetOpenO1 ...
        VetOpenO2 VetOpenO3 VetOpenA1 VetOpenA2 VetOpenA3 SavingAT2 ...
        VetNumCar2 VetOpen2 IfftOrSum MaxNumStagT OSNRPC EoutTx BerOOK ...
        BERt4PAM BERtOOK InitCarrUp InitCarrDo CarrPass CarrUsedDo snrp ...
        CarrUsedUp SelecSetUP VetThisCarrTx CarrRecPowUp CarrRecPowDo ...
        AberLevQ ValsLevQ AberLevI ValsLevI LevDefDpqsk LightSpeed ...
        AddingNoiseO AddingNoiseE AsePower osnrf AddingNoiseF BerOOKS ...
        AddingNoiseO AddingNoiseE AddingNoiseP osnrp AberLevS ValsLevS ...
        BerDPSKS BerDQPSKS AberLevIS ValsLevIS AberLevQS ValsLevQS CalcS...
        Ber4PAMS PrintinEye SetCpSampZer ReceptorNoise EvmDBJ EvmRmsJ ...
        EvmPerJ EvmMatRec EvmMatRef Nsamp RecFilBanPas VetSnr SendingUp
                  
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
        if CurrentModula == 2
            berpos1 = 1:2:size(Ber4PAM,2);
            berpos2 = 2:2:size(Ber4PAM,2);
            ber = [Ber4PAM(CurrentTest,berpos1);Ber4PAM(CurrentTest,berpos2)]
            b = Ber4PAM;
            
            toc - PastTime
            PastTime = toc;
            if PastTime/60 > 15
                if CurrentModula == UpStModula
                    save(SavingAT,'NumbOfTest',...
'Ber4PAM','ValsLev','VetOptiPower','VetElecPower','AberLev',...
'EvmRms','EvmDB','EvmPer','DecLevDef1','DecLevDef2','DecLevDef3',...
'OSNRPC','OfcName','RefCarr', 'CarrUsed','SelecSetUP','FFTSplit',...
'CurrentMedium','UsedModula','NumCarr','PAM4Type','UpStModula',...
'IfftOrSum','NumCarr','lambdac','FiberLength','InitCarrUp', 'InitCarrDo'...
,'CarrPass', 'CarrUsedDo','CarrUsedUp','ObsCarrUsed',...
'osnrf','osnr','snr','AddingNoiseF','AddingNoiseE',...
'AddingNoiseO','Ber4PAMS','AberLevS','ValsLevS','EvmMatRec','EvmMatRef',...
'EvmDB','EvmPer','EvmRms','EvmDBJ','EvmPerJ','EvmRmsJ');
                else
                    save(SavingAT,'NumbOfTest',...
'Ber4PAM','ValsLev2','ValsLev','VetOptiPower','VetElecPower','AberLev',...
'EvmRms','EvmDB','EvmPer','DecLevDef1','DecLevDef2','DecLevDef3',...
'OSNRPC','OfcName','VetOptiPowerUp','VetElecPowerUp','ValsLevUp',...
'AberLevUp','EvmRmsUp','EvmDBUp','EvmPerUp','RefCarr', 'CarrUsed',...
'SelecSetUP','FFTSplit','CurrentMedium','UsedModula','NumCarr',...
'PAM4Type','UpStModula','IfftOrSum','NumCarr','lambdac','FiberLength',...
'InitCarrUp', 'InitCarrDo', 'CarrPass', 'CarrUsedDo','CarrUsedUp',...
'ObsCarrUsed','CarrRecPowUp', 'CarrRecPowDo','osnrf','osnr','snr',...
'AddingNoiseF','AddingNoiseE','AddingNoiseO','Ber4PAMS','AberLevS',...
'ValsLevS','EvmMatRec','EvmMatRef','EvmDB','EvmPer','EvmRms','EvmDBJ',...
'EvmPerJ','EvmRmsJ');
                end
                PastTime = 0;
                tic;
            end
        elseif CurrentModula == 1
            berpos1 = 1:2:size(BerOOK,2);
            berpos2 = 2:2:size(BerOOK,2);
            ber = [BerOOK(CurrentTest,berpos1);BerOOK(CurrentTest,berpos2)]
            b = BerOOK;
            
            toc - PastTime
            PastTime = toc;
            if PastTime/60 > 15
                save(SavingAT,'NumbOfTest',...
'BerOOK','ValsLev','VetOptiPower','VetElecPower','AberLev',...
'EvmRms','EvmDB','EvmPer','DecLevDef1','DecLevDef2','DecLevDef3',...
'OSNRPC','OfcName','RefCarr', 'CarrUsed','SelecSetUP','FFTSplit',...
'CurrentMedium','UsedModula','NumCarr','PAM4Type','UpStModula',...
'IfftOrSum','NumCarr','lambdac','FiberLength','InitCarrUp', 'InitCarrDo'...
,'CarrPass', 'CarrUsedDo','CarrUsedUp','ObsCarrUsed','CarrRecPowUp',...
'CarrRecPowDo','osnrf','osnr','snr','AddingNoiseF','AddingNoiseE',...
'AddingNoiseO','BerOOKS','AberLevS','ValsLevS','EvmMatRec','EvmMatRef',...
'EvmDB','EvmPer','EvmRms','EvmDBJ','EvmPerJ','EvmRmsJ');
                PastTime = 0;
                tic;
            end
        elseif CurrentModula == 3
            berpos1 = 1:2:size(BerDQPSK,2);
            berpos2 = 2:2:size(BerDQPSK,2);
            ber=[BerDQPSK(CurrentTest,berpos1);BerDQPSK(CurrentTest,berpos2)]
            b = BerDQPSK;
            toc - PastTime
            PastTime = toc;
            if PastTime/60 > 15
                save(SavingAT,'NumbOfTest','BerDQPSK','VetOptiPower',...
'AberLevQ','ValsLevQ','AberLevI','ValsLevI','DecLevDef1','DecLevDef2',...
'DecLevDef3','OSNRPC','OfcName','RefCarr','CarrUsed','SelecSetUP',...
'FFTSplit','CurrentMedium','UsedModula','NumCarr','PAM4Type',...
'UpStModula','IfftOrSum','NumCarr','lambdac','FiberLength',...
'InitCarrUp','InitCarrDo','CarrPass','CarrUsedDo','CarrUsedUp',...
'ObsCarrUsed','osnrf','osnr','snr',...
'AddingNoiseF','AddingNoiseE','AddingNoiseO','BerDQPSKS','AberLevIS',...
'ValsLevIS','AberLevQS','ValsLevQS','EvmMatRec','EvmMatRef','EvmDB',...
'EvmPer','EvmRms','EvmDBJ','EvmPerJ','EvmRmsJ');
                PastTime = 0;
                tic;
            end
        elseif CurrentModula == 4
            berpos1 = 1:2:size(BerDPSK,2);
            berpos2 = 2:2:size(BerDPSK,2);
            ber = [BerDPSK(CurrentTest,berpos1);BerDPSK(CurrentTest,berpos2)]
            b = BerDPSK;
            
            toc - PastTime
            PastTime = toc;
            if PastTime/60 > 15
                save(SavingAT,'NumbOfTest','BerDPSK','ValsLev',...
'VetOptiPower','VetElecPower','AberLev','DecLevDef1','DecLevDef2',...
'DecLevDef3','OSNRPC','OfcName','RefCarr', 'CarrUsed','SelecSetUP',...
'FFTSplit','CurrentMedium','UsedModula','NumCarr','PAM4Type',...
'UpStModula','IfftOrSum','NumCarr','lambdac','FiberLength','InitCarrUp',...
'InitCarrDo','CarrPass', 'CarrUsedDo','CarrUsedUp','ObsCarrUsed',...
'CarrRecPowUp','CarrRecPowDo','osnrf','osnr','snr','AddingNoiseF',...
'AddingNoiseE','AddingNoiseO','BerDPSKS','AberLevS','ValsLevS',...
'EvmMatRec','EvmMatRef','EvmDB','EvmPer','EvmRms','EvmDBJ','EvmPerJ',...
'EvmRmsJ');
                PastTime = 0;
                tic;
            end
        end
        CurrentTest
        VetSnr
        if CurrentTest > 2
            if (((CurrentTest)*ceil(0.8*Nb))>(5/min(mean(b))))||(CurrentTest==3000)
                break;
            end
        end
    end
    if CurrentModula == 2
        if CurrentModula == UpStModula
                    save(SavingAT,'NumbOfTest',...
'Ber4PAM','ValsLev','VetOptiPower','VetElecPower','AberLev',...
'EvmRms','EvmDB','EvmPer','DecLevDef1','DecLevDef2','DecLevDef3',...
'OSNRPC','OfcName','RefCarr', 'CarrUsed','FFTSplit',...
'CurrentMedium','UsedModula','NumCarr','PAM4Type','UpStModula',...
'IfftOrSum','NumCarr','lambdac','FiberLength'...
,'CarrPass','ObsCarrUsed','osnrf','osnr','snr','AddingNoiseF','AddingNoiseE',...
'AddingNoiseO','Ber4PAMS','AberLevS','ValsLevS','EvmMatRec','EvmMatRef',...
'EvmDB','EvmPer','EvmRms','EvmDBJ','EvmPerJ','EvmRmsJ');
                else
                    save(SavingAT,'NumbOfTest',...
'Ber4PAM','ValsLev2','ValsLev','VetOptiPower','VetElecPower','AberLev',...
'EvmRms','EvmDB','EvmPer','DecLevDef1','DecLevDef2','DecLevDef3',...
'OSNRPC','OfcName','VetOptiPowerUp','VetElecPowerUp','ValsLevUp',...
'AberLevUp','EvmRmsUp','EvmDBUp','EvmPerUp','RefCarr', 'CarrUsed',...
'SelecSetUP','FFTSplit','CurrentMedium','UsedModula','NumCarr',...
'PAM4Type','UpStModula','IfftOrSum','NumCarr','lambdac','FiberLength',...
'InitCarrUp', 'InitCarrDo', 'CarrPass', 'CarrUsedDo','CarrUsedUp',...
'ObsCarrUsed','CarrRecPowUp', 'CarrRecPowDo','osnrf','osnr','snr',...
'AddingNoiseF','AddingNoiseE','AddingNoiseO','Ber4PAMS','AberLevS',...
'ValsLevS','EvmMatRec','EvmMatRef','EvmDB','EvmPer','EvmRms','EvmDBJ',...
'EvmPerJ','EvmRmsJ');
        end
    elseif CurrentModula == 1
                save(SavingAT,'NumbOfTest',...
'BerOOK','ValsLev','VetOptiPower','VetElecPower','AberLev',...
'EvmRms','EvmDB','EvmPer','DecLevDef1','DecLevDef2','DecLevDef3',...
'OSNRPC','OfcName','RefCarr','FFTSplit',...
'CurrentMedium','UsedModula','NumCarr','PAM4Type','UpStModula',...
'IfftOrSum','NumCarr','lambdac','FiberLength'...
,'CarrPass','ObsCarrUsed','osnrf','osnr','snr','AddingNoiseF','AddingNoiseE',...
'AddingNoiseO','BerOOKS','AberLevS','ValsLevS','EvmMatRec','EvmMatRef',...
'EvmDB','EvmPer','EvmRms','EvmDBJ','EvmPerJ','EvmRmsJ');
    elseif CurrentModula == 3
                save(SavingAT,'NumbOfTest','BerDQPSK','VetOptiPower',...
'AberLevQ','ValsLevQ','AberLevI','ValsLevI','DecLevDef1','DecLevDef2',...
'DecLevDef3','OSNRPC','OfcName','RefCarr','CarrUsed',...
'FFTSplit','CurrentMedium','UsedModula','NumCarr','PAM4Type',...
'UpStModula','IfftOrSum','NumCarr','lambdac','FiberLength',...
'InitCarrUp','InitCarrDo','CarrPass','CarrUsedDo','CarrUsedUp',...
'ObsCarrUsed','osnrf','osnr','snr',...
'AddingNoiseF','AddingNoiseE','AddingNoiseO','BerDQPSKS','AberLevIS',...
'ValsLevIS','AberLevQS','ValsLevQS','EvmMatRec','EvmMatRef','EvmDB',...
'EvmPer','EvmRms','EvmDBJ','EvmPerJ','EvmRmsJ');
    elseif CurrentModula == 4
                save(SavingAT,'NumbOfTest','BerDPSK','ValsLev',...
'VetOptiPower','VetElecPower','AberLev','DecLevDef1','DecLevDef2',...
'DecLevDef3','OSNRPC','OfcName','RefCarr', 'CarrUsed','SelecSetUP',...
'FFTSplit','CurrentMedium','UsedModula','NumCarr','PAM4Type',...
'UpStModula','IfftOrSum','NumCarr','lambdac','FiberLength','InitCarrUp',...
'InitCarrDo','CarrPass', 'CarrUsedDo','CarrUsedUp','ObsCarrUsed',...
'CarrRecPowUp','CarrRecPowDo','osnrf','osnr','snr','AddingNoiseF',...
'AddingNoiseE','AddingNoiseO','BerDPSKS','AberLevS','ValsLevS',...
'EvmMatRec','EvmMatRef','EvmDB','EvmPer','EvmRms','EvmDBJ','EvmPerJ',...
'EvmRmsJ');
    end
    a=1;
end
end
a=1;



















