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
%c                                           02/04/2018
%c                                           pablorafael.mcx@gmail.com
%c 
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
if ~exist('CurrentOCS','var')
    clear ; clc; close all;tic;
end
PARMainInputData;
%%              Loading the Ultra OFC Source
% Initially this scrip will load the Ultra Optical Flat Comb Source, which
% will be the carriers equaly spaced of our sistem. It is important to
% remember that each carrier was originated from an optical ring in which a
% single laser CW (Continuos Wave) was used. The main point of this scrip
% is to demonstrate how it is possible to implement an All Optical FFT.
% Therefore, the UOFCS needs to be precisely spaced according to the bit
% ratio Rs.That is directly related to the Symbol ration of the sistem.

% clc;
%The two following loops are to control the test bundle. The first loop is 
%responsible to select the type of modulation that will be under test and 
%the second loop will determinate how many tests will be running.
for CurrentModula=ThisModula
    for CurrentTest=1:ThisTest
        if TestBundle
            load([OfcName '.mat']);                                        %Loading a OFC from memory.
        else
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
        if SigAten
            Eout = Eout/Atenuation;
        end
        Ploting = 0;
        PrintInfo(Ploting*1,f,Eout);
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
        
        PARDatumTransmission;

        %%              Datum Reception
        % The following step is to recover the information transmited 
        %throughout an medium. It is the most important part of this work, 
        %because it will deal with an All-Optical FFT implementation. It 
        %aims to preciselly separate each sub carrier that compose the OFDM 
        %Symble. It will apply methods previouly presented in the 
        %literature. Looking to do all the FFT processing at the optical 
        %domain.

        if SendingDowStr
            PARDatumReception;
        else
            PARDatumReceptionInputData;
            [EoutAux1,~,VetThisCarr]=OpticalFFTN(t,T,MaxNumStag,EoutRec);
        end
        
        %%              Datum Upstream
        %This section is resposible to generate and to send the information
        %from the user to the OLT.
        PARDatumUpstream;
        
        %%                Datum Upstream Reception
        %This section is responsible to receive the information from users
        PARDatumUpstrRec;
        %%              Saving
        clearvars -except ThisPath UofcsLocal OfcName NumbOfTest PowRef...
PAM4Type ThisTest TestBundle CurrentMedium UsedModula UsedModula SnrRef ...
BerDQPSKF BerDPSKF BerOOKF Ber4PAMF FiberDelay ValsLev2 ValsLev PastTime...
VetOptiPower VetElecPower AberLev FFTSplit CarSNR CurrentModula BerOOK ...
VetOptiPowerF VetElecPowerF AberLev1 AberLev2 AberLev3 ValsLev1 SigAten...
ValsLev02 ValsLev3 ValsLev21 ValsLev22 ValsLev23 ValsLev21 Ber4PAM ...
VetElecPowerIF VetElecPowerQF BerDQPSK BerDPSK Atenuation  ThisModula...
SendingDowStr EvmRms EvmDB EvmPer CurrentTest VetElecPowerI EvmRmsF ...
ValsLevUp AberLevUp VetOptiPowerUp VetElecPowerUp EvmRmsUp EvmDBUp ...
DecLevDef1 DecLevDef2 DecLevDef3 ValsLev01 AberLev01 EvmPerUp UpStModula...
                                               EvmDBF EvmPerF VetElecPowerQ
                                                        
        if exist('BerDPSK','var')
            VetOptiPower(CurrentTest,:)            = VetOptiPowerF;
            VetElecPower(CurrentTest,:)            = VetElecPowerF;
            AberLev(CurrentTest,:)                 = AberLev1;
            ValsLev(CurrentTest,:)                 = ValsLev1;
            berpos1 = 1:2:size(BerDPSK,2);
            berpos2 = 2:2:size(BerDPSK,2);
            b = [BerDPSK(1,berpos1);BerDPSK(1,berpos2)]
            BerDPSKF(CurrentTest,1:length(BerDPSK))= BerDPSK;
        end
        
        if exist('BerDQPSK','var')
            VetOptiPower(CurrentTest,:)      = VetOptiPowerF;
            VetElecPowerI(CurrentTest,:)       = VetElecPowerIF;
            VetElecPowerQ(CurrentTest,:)       = VetElecPowerQF;
            AberLevI(CurrentTest,:)           = AberLev1;
            AberLevQ(CurrentTest,:)           = AberLev2;
            ValsLevI(CurrentTest,:)           = ValsLev1;
            ValsLevQ(CurrentTest,:)           = ValsLev02;
            berpos1 = 1:2:size(BerDQPSK,2);
            berpos2 = 2:2:size(BerDQPSK,2);
            b = [BerDQPSK(1,berpos1);BerDQPSK(1,berpos2)]
            BerDQPSKF(CurrentTest,:)          = BerDQPSK;
        end
        if exist('Ber4PAM','var')
            if CurrentModula == UpStModula
                EvmRms(CurrentTest,:) = EvmRmsF;
                EvmDB(CurrentTest,:)  = EvmDBF;
                EvmPer(CurrentTest,:) = EvmPerF;
                VetOptiPower(CurrentTest,:)=VetOptiPowerF;
                VetElecPower(CurrentTest,:)=VetElecPowerF;
                AberLev(1,1:length(AberLev1),CurrentTest)=AberLev1;
                AberLev(2,1:length(AberLev2),CurrentTest)=AberLev2;
                AberLev(3,1:length(AberLev3),CurrentTest)=AberLev3;
                ValsLev(1,1:length(ValsLev1),CurrentTest)=ValsLev1;
                ValsLev(2,1:length(ValsLev02),CurrentTest)=ValsLev02;
                ValsLev(3,1:length(ValsLev3),CurrentTest)=ValsLev3;
                ValsLev2(1,1:length(ValsLev21),CurrentTest)=ValsLev21;
                ValsLev2(2,1:length(ValsLev22),CurrentTest)=ValsLev22;
                ValsLev2(3,1:length(ValsLev23),CurrentTest)=ValsLev23;
                berpos1 = 1:2:size(Ber4PAM,2);
                berpos2 = 2:2:size(Ber4PAM,2);
                b = [Ber4PAM(1,berpos1);Ber4PAM(1,berpos2)]
                Ber4PAMF(CurrentTest,1:length(Ber4PAM))= Ber4PAM;
            else
                AberLevUp(CurrentTest,:)= AberLev01;
                ValsLevUp(CurrentTest,:)= ValsLev01;
                EvmRms(CurrentTest,:) = EvmRmsF;
                EvmDB(CurrentTest,:)  = EvmDBF;
                EvmPer(CurrentTest,:) = EvmPerF;
                VetOptiPower(CurrentTest,:)=VetOptiPowerF;
                VetElecPower(CurrentTest,:)=VetElecPowerF;
                AberLev(1,1:length(AberLev1),CurrentTest)=AberLev1;
                AberLev(2,1:length(AberLev2),CurrentTest)=AberLev2;
                AberLev(3,1:length(AberLev3),CurrentTest)=AberLev3;
                ValsLev(1,1:length(ValsLev1),CurrentTest)=ValsLev1;
                ValsLev(2,1:length(ValsLev02),CurrentTest)=ValsLev02;
                ValsLev(3,1:length(ValsLev3),CurrentTest)=ValsLev3;
                ValsLev2(1,1:length(ValsLev21),CurrentTest)=ValsLev21;
                ValsLev2(2,1:length(ValsLev22),CurrentTest)=ValsLev22;
                ValsLev2(3,1:length(ValsLev23),CurrentTest)=ValsLev23;
                berpos1 = 1:2:size(Ber4PAM,2);
                berpos2 = 2:2:size(Ber4PAM,2);
                b = [Ber4PAM(1,berpos1);Ber4PAM(1,berpos2)]
                Ber4PAMF(CurrentTest,1:length(Ber4PAM))= Ber4PAM;end
        end
        if exist('BerOOK','var')
            EvmRms(CurrentTest,:) = EvmRmsF;
            EvmDB(CurrentTest,:)  = EvmDBF;
            EvmPer(CurrentTest,:) = EvmPerF;
            VetOptiPower(CurrentTest,:)            = VetOptiPowerF;
            VetElecPower(CurrentTest,:)            = VetElecPowerF;
            AberLev(CurrentTest,:)                 = AberLev1;
            ValsLev(CurrentTest,:)                 = ValsLev1;
            berpos1 = 1:2:size(BerOOK,2);
            berpos2 = 2:2:size(BerOOK,2);
            b = [BerOOK(1,berpos1);BerOOK(1,berpos2)]
            BerOOKF(CurrentTest,1:length(BerOOK))= BerOOK;
        end
%         DatumUpstreamReception;
        clear BerOOK Ber4PAM BerDQPSK BerDPSK VetOptiPowerF AberLev1 ...
VetElecPowerF AberLev3 AberLev2 ValsLev1 ValsLev02 ValsLev3 ValsLev21 ...
EvmRmsF EvmDBF EvmPerF VetElecPowerIF VetElecPowerQF ValsLev22 ValsLev23... 
                                                        ValsLev01 AberLev01 
        CurrentTest
%         contcp = contcp + 1;
%         save('Teste0_4PAM_B2B');
        toc - PastTime
        PastTime = toc;
        if PastTime/60 > 15
            if CurrentModula == 2
                save('Teste1_Do4PAM_UpOOK_40km');
            elseif CurrentModula == 1
                save('Teste0_OOK');
            elseif CurrentModula == 3
                save('Teste0_DQPSK');
            elseif CurrentModula == 4
                save('Teste0_DPSK');
            end
            PastTime = 0;
            tic;
        end
    end
    if CurrentModula == 2
        save('Teste1_Do4PAM_UpOOK_40km');
    elseif CurrentModula == 1
        save('Teste0_OOK');
    elseif CurrentModula == 3
        save('Teste0_DQPSK');
    elseif CurrentModula == 4
        save('Teste0_DPSK');
    end
end

a=1;



















