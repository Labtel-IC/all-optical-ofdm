%c
%c                                                       ..'Ž`'..'Ž`..'Ž`..
%c       File: DatumTransmission File For the All Optical OFDM
%c(Main file to call each idividual parameters)
%c
%c     This main code is resposible to call and run the all the fuction to
%c related to this simulation, which is responsible to generate a random
%c strem of data, modulated it and transmite it accross a SSMF.
%c
%c
%c
%c                                           by P.Marciano LG
%c                                           29/10/2017
%c                                           Last UpDate
%c                                           02/01/2018
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
% The first step here is to configure all the parameters that will used
% hereafter by the transmiter. Hence, implementing functionalities such as
% generation of random data to fill the pay load, create the syncronism
% symble, modulating the incoming data acordingly, transmiting this
% information throughout a medium (fiber or B2B). All those functionalities
% will need variables to be first initialyzed in this following Input Data.
DatumTransmissionInputData;
%%               Creating Data for transmission
% Creating the information accordingly with the type of transmition

%It is important to mention that an stream of data will be formed for each
%individual carrier. For the last but not the least, the DownStream and the
%UpStream carriers will be interleaved within the same transmition pass
%band.For example, if the user choses carriers 1,3,5 and 7 as DownStream
%carriers, the UpStream will be formed by carriers 2,4,6, and 8. This
%manuver was suggested by Segatto to address the problem of ICI (Inter
%Carrier Interference). The transmited signal still an OFDM once the
%carriers from 1 to 8, for instance, still orthogonal to each other.
%Therefore, the for loop to generate the information to be transmited will
%have the passe of 2.

if SendingDowStr
    
    %%  Selecting part of the UOFCS
    %The generation fo the Ultra Optical Flat Comb Source is not perfect, which
    %may drive the user to select a specific part of it for its actual use.
    %As a result, this condition verify if the user wants to select a part of
    %the UOFCS or use it as an whole.
    if Selecting
        EoutT = ifft(fft(Eout).*SelecFilt);
    else
        EoutT = Eout;
    end
    
    if ~exist('EoutTx','var')
        EoutTx=SelectEachCarrier(EoutT,NumCarr,f,fin,FBWD,Order,fc);          %This function is resposible to split each carrier
        %         figure;
        %         hold all;
        %         for kk=1:size(EoutTx,1)
        %             plot(f,20*log10(abs(fftshift(fft(EoutTx(kk,:)./length(EoutTx(kk,:)))))));
        %         end
    end
    
    switch Modulation
        case 'DPSK'
            %%            Generate the data DPSK
            for kk=InitCarrDo:CarrPass:NumCarr%Generating different data for each carrier
                % Frist it is chosen to transmite a random information...
                if SendingData
                    TxData = (randi(2,1,NbDPSK)-1);                            %Creating the data stream to be transmited
                    TxData(1:JusPos) = JusVal;                                 %Adding the Justification simble at the begining of the stream to sincronize received data frame
                    TxData(end-(JusPos-1):end) = JusValEnd;                    %Adding the Justification simble at the end of the stream to sincronize received data frame
                else %... or just one high pulse for testing the channel
                    TxData        = zeros(1,NbDPSK);                           %Creating the base fo the data stream
                    TxData(end/2) = 1;                                         %Adding the code in the data stream to produce the highest level possible
                end
                
                TxDataMat(kk,:)   = TxData;                                    %Storring the transmited information for latter evaluation
                
                %Once that we have the data to be transmited it needs to be
                %converted to the eletrical signal that hereafter will modulate
                %our carrier through the MZM-I
                SigTx           = DpskEncodEq( TxData);
                SigTx(SigTx==1) =  1;
                SigTx(SigTx==0) = -1;
                SigTx           = rectpulse(SigTx,NPPB);
                
                % Adding CP to the data
                if AddCP
                    TxAux = reshape(SigTx,NPPB,NbDPSK);
                    TxAux = [TxAux(1:NumAmosCP,:);TxAux;TxAux(end-(...
                        NumAmosCP-1):end,:)];
                    TxAux = reshape(TxAux,1,(2*NumAmosCP+NPPB)*NbDPSK);
                    SigTx = [TxAux TxAux(end-(StuffSampels-1):end)];
                end
                
                %Because of reasons, the real world does not have components
                %capable instantly change the voltage level. Thus, to emulate
                %this behaviour the eletrical PAM signal pass through a
                %gaussian filter that will remove the frequencies of higher
                %order, which will result in small slope in the level variation
                
                [BitFilt,~] = FiltroGaussiano(f,BWD,CenFeq,FiltOrd);           %Creating filter for conformation of the input information
                BitFilt     = fftshift(BitFilt);                               %Doing a shift on the vector for matching the transmited data
                TxSig       = ifft(fft(SigTx).*BitFilt);                       %Conforming the information and Creating the modulation signal
                TxSig       = DatGai.*TxSig;                                   %Adding a gain to the eletrical signal
                %             PrintInfo(Ploting*2,TxSig,T,NPPB);
                TxSigMat(kk,:) = TxSig;                                        %Storing the eletrical signal for further evaluation
                a=4;
                %If the user choses to modulate all carriers at the same time
                %there is no need of this script to generate data for each
                %individual carrier. Therefore, this For loop can be halted
                if ModAll                                                      %Check if ModAll is set OFF
                    break;                                                     %If not the loop is halted
                end
            end
        case 'DQPSK'
            %%            Generate the data DQPSK
            for kk=InitCarrDo:CarrPass:NumCarr%Generating different data for each carrier
                % Frist it is chosen to transmite a random information...
                if SendingData == 1
                    TxData = (randi(2,1,NbDQPSK)-1);                           %Creating the data stream to be transmited
                    TxData(1:JusPos) = JusVal;                                 %Adding the Justification simble at the begining of the stream to sincronize received data frame
                    TxData(end-(JusPos-1):end) = JusValEnd;                    %Adding the Justification simble at the end of the stream to sincronize received data frame
                    PreBits  = [0 0];                                          %Seting the start point for the DQPSK maping
                elseif SendingData == 2
                    %For previous test it was needed a known bit sequence. This
                    %section creates it. The [I Q] output will be:
                    %[1 1 0 0 0 1 0 0 0 0 1 0 0 1 0 1 1 1 1 1 0 1 1 0 1 1 1 0 1
                    % 0 0 0]
                    %This sequency covers all possible combinations of I and Q
                    %which is:
                    %I = 1 0 0 0 0 1 0 0 1 1 0 1 1 1 1 0
                    %Q = 1 0 1 0 0 0 1 1 1 1 1 0 1 0 0 0
                    TxDataSeg = [0 0; %01
                        0 0; %02
                        0 1; %03
                        1 0; %04
                        1 1; %05
                        1 0; %06
                        0 0; %07
                        1 1; %08
                        0 1; %09
                        1 1; %10
                        1 0; %11
                        0 0; %12
                        1 0; %13
                        0 1; %14
                        1 1; %15
                        0 1];%16
                    TxData   = [];
                    PreBits  = [0 0];
                    countseg = 1;
                    for jj=1:Nb
                        TxData   = [TxData TxDataSeg(countseg,:)];
                        countseg = countseg + 1;
                        if countseg > 16
                            countseg = 1;
                        end
                    end
                else %... or just one high pulse for testing the channel
                    TxData = zeros(1,NbDQPSK);                                 %Creating the base fo the data stream
                    TxData((end/2)-1:end/2) = [1 1];                           %Adding the code in the data stream to produce the highest level possible
                    PreBits  = [0 0];
                end
                TxDataPos       = linspace(1,NbDQPSK,NbDQPSK);                 %Auxiliar variable to split the TxData
                [DataI,DataQ]   = DqpskEncodEq(TxData,PreBits);                %Maping the income data to the DQPSK format
                %Converting the I and Q components to the polirezed NRZ format
                DataI(DataI==1) =  1;
                DataI(DataI==0) = -1;
                DataQ(DataQ==1) =  1;
                DataQ(DataQ==0) = -1;
                
                TxOdd  = TxData(logical(mod(TxDataPos,2)));                    %Spliting the information of odd positions
                TxEven = TxData(~(mod(TxDataPos,2)));                          %Spliting the information of even positions
                
                TxDataMat(kk,:)         = TxOdd;                               %Storring the transmited information for latter evaluation
                TxDataMat(kk+NumCarr,:) = TxEven;                              %Storring the transmited information for latter evaluation
                
                %Once that we have the data to be transmited it needs to be
                %converted to the eletrical signal that hereafter will modulate
                %our carrier through the MZM-I
                SigI = rectpulse(DataI,NPPB);
                SigQ = rectpulse(DataQ,NPPB);
                
                % Adding CP to the data
                if AddCP
                    TxAux1 = reshape(SigI,NPPB,NbDQPSK/2);
                    TxAux1 = [TxAux1(1:NumAmosCP,:);TxAux1;TxAux1(end-(...
                        NumAmosCP-1):end,:)];
                    TxAux1 = reshape(TxAux1,1,(2*NumAmosCP+NPPB)*NbDQPSK/2);
                    SigI = [TxAux1 TxAux1(end-(StuffSampels-1):end)];
                    
                    TxAux2 = reshape(SigQ,NPPB,NbDQPSK/2);
                    TxAux2 = [TxAux2(1:NumAmosCP,:);TxAux2;TxAux2(end-(...
                        NumAmosCP-1):end,:)];
                    TxAux2 = reshape(TxAux2,1,(2*NumAmosCP+NPPB)*NbDQPSK/2);
                    SigQ = [TxAux2 TxAux2(end-(StuffSampels-1):end)];
                end
                
                %Because of reasons, the real world does not have components
                %capable instantly change the voltage level. Thus, to emulate
                %this behaviour the eletrical PAM signal pass through a
                %gaussian filter that will remove the frequencies of higher
                %order, which will result in small slope in the level variation
                
                [BitFilt,~] = FiltroGaussiano(f,BWD,CenFeq,FiltOrd);           %Creating filter for conformation of the input information
                BitFilt = fftshift(BitFilt);                                   %Doing a shift on the vector for matching the transmited data
                TxSigI = ifft(fft(SigI).*BitFilt);                             %Conforming the information and Creating the modulation signal
                TxSigQ = ifft(fft(SigQ).*BitFilt);                             %Conforming the information and Creating the modulation signal
                %             PrintInfo(Ploting*3,TxSigI,T,NPPB,TxSigQ);
                TxSigMat(kk,:)         = TxSigI;                               %Storing the eletrical signal for further evaluation
                TxSigMat(kk+NumCarr,:) = TxSigQ;                               %Storing the eletrical signal for further evaluation
                a=4;
                %If the user choses to modulate all carriers at the same time
                %there is no need of this script to generate data for each
                %individual carrier. Therefore, this For loop can be halted
                if ModAll                                                      %Check if ModAll is set OFF
                    break;                                                     %If not the loop is halted
                end
            end
        case '4PAM'
            %%        Generate the data 4PAM
            EoutMod    = 0;
            EoutModTem = 0;
            for kk=InitCarrDo:CarrPass:NumCarr%Generating different data for each carrier
                % Frist it is chosen to transmite a random information...
                if CarrUsedDo(kk)
                    if SendingData
                        TxData = (randi(2,1,Nb4Pam)-1);                            %Creating the data stream to be transmited
                        TxData(1:JusPos) = JusVal;                                 %Adding the Justification simble at the begining of the stream to sincronize received data frame
                        TxData(end-(JusPos-1):end) = JusValEnd;                    %Adding the Justification simble at the end of the stream to sincronize received data frame
                    else%... or just one high pulse for testing the channel
                        TxData = zeros(1,Nb4Pam);                                  %Creating the base fo the data stream
                        TxDataPos = linspace(1,Nb4Pam,Nb4Pam);                     %Creating a vector to auxilliate the addresing process
                        TxData(~mod(TxDataPos,2))=1;                               %Adding the code in the data stream to produce the lowest level possible
                        TxData(end/2 - 1) = 1;                                     %Adding the code in the data stream to produce the highest level possible
                    end
                    
                    TxDataMat(kk,:) = TxData;                                      %Storring the transmited information for latter evaluation
                    %Once that we have the data to be transmited it needs to be
                    %converted to the eletrical signal that hereafter will modulate
                    %our carrier through the MZM-I
                    
                    %Selecting if the PAM will be done on Electrical or Optical
                    %domain
                    if Which4PAM
                        %  Good for others kms
                        [Phi1,Phi2] = Maping4PamIq(TxData,0,U_pi2/2,ModSchem...
                            ,FiberLength);%Generating the eletrical signal for the optical PAM4
                        Vbias=U_pi2/2;
                        %  Good for 100 kms
                        %                 [Phi1,Phi2] = Maping4PamIq2(TxData,0,U_pi2/2);
                        %                 Vbias=U_pi2/2;
                    else
                        [Phi1,Phi2] = Maping4Pam(TxData,VPI,Polirized,MaxAmp4PAM); %Generating the eletrical signal for the electrical PAM4
                    end
                    %The signal generated are not yet with the same number of
                    %samples as the OFCS loaded. These nexte lines do oversampling
                    %of the electrical signal.
                    TxSig1 = rectpulse(Phi1,NPPB);
                    TxSig2 = rectpulse(Phi2,NPPB);
                    
                    %Thus, if it would be required to add cycle prefix the number
                    %of samples per symbol needs to change as well as some
                    %adjustments needs to be done for the new signal match in size
                    %with the size of the vector time. This problem just exist on
                    %simulation, at practice the main point is the syncronism of
                    %the signals.
                    % Adding CP to the data
                    if AddCP
                        TxAux1 = reshape(TxSig1,NPPB,Nb4Pam/2);
                        TxAux1 = [TxAux1(1:NumAmosCP,:);TxAux1;TxAux1(end-(...
                            NumAmosCP-1):end,:)];
                        TxAux1 = reshape(TxAux1,1,(2*NumAmosCP+NPPB)*Nb4Pam/2);
                        TxSig1 = [TxAux1 TxAux1(end-(StuffSampels-1):end)];
                        
                        TxAux2 = reshape(TxSig2,NPPB,Nb4Pam/2);
                        TxAux2 = [TxAux2(1:NumAmosCP,:);TxAux2;TxAux2(end-(...
                            NumAmosCP-1):end,:)];
                        TxAux2 = reshape(TxAux2,1,(2*NumAmosCP+NPPB)*Nb4Pam/2);
                        TxSig2 = [TxAux2 TxAux2(end-(StuffSampels-1):end)];
                    end
                    
                    %Because of reasons, the real world does not have components
                    %capable instantly change the voltage level. Thus, to emulate
                    %this behaviour the eletrical PAM signal pass through a
                    %gaussian filter that will remove the frequencies of higher
                    %order, which will result in small slope in the level variation
                    
                    [BitFilt,~]            = FiltroGaussiano(f,BWD,CenFeq,FiltOrd);%Creating filter for conformation of the input information
                    BitFilt                = fftshift(BitFilt);                    %Doing a shift on the vector for matching the transmited data
                    TxSig1                 = ifft(fft(TxSig1).*BitFilt);           %Conforming the information and Creating the modulation signal
                    TxSig2                 = ifft(fft(TxSig2).*BitFilt);           %Conforming the information and Creating the modulation signal
                    %                 TxSigMat(kk,:)         = TxSig1;                               %Storing the eletrical signal for further evaluation
                    %                 TxSigMat(kk+NumCarr,:) = TxSig2;                               %Storing the eletrical signal for further evaluation
                                if kk==1
                                    PrintInfo(PlotingThis*4,TxSig1(kk,1:end-StuffSampels),t(2*NumAmosCP+NPPB),2*NumAmosCP+NPPB,TxSig2(1:end-StuffSampels));
                                    a=4;
                                end
                    %If the user choses to modulate all carriers at the same time
                    %there is no need of this script to generate data for each
                    %individual carrier. Therefore, this For loop can be halted
                    
                    U.U1t = TxSig1;                                        %Assigning the electrical signal to one drive of the MZM
                    U.U2t = TxSig2;                                        %Assigning the electrical signal to another drive of the MZM
                    if Which4PAM
                        if ModSchem
                            [EoutModAux] = IqMod4Pam (EoutTx(kk,:),U.U1t,...
                                U.U2t,U_pi2,Vbias);
                        else
                            [EoutModAux] = ...
                                Mach_Zehnder_Modulator_simplificado(t,EoutTx(kk,...
                                :),U,MZ_Input_File);
                        end
                    else
                        [EoutModAux] = Mach_Zehnder_Modulator_simplificado(t...
                            ,EoutTx(kk,:),U,MZ_Input_File);
                    end
                    %EoutMod = EoutMod + EoutModAux + EoutTx(kk+1,:);       %Adding the current Modulated output to the final OutPut
                    EoutModTem(ObsCarrPos(kk),1:length(EoutModAux))   = EoutModAux;    %Adding the current Modulated output to the final OutPut
                    if ~mod(CarrPass,2)
                        EoutModTem(ObsCarrPos(kk+1),1:length(EoutModAux)) = EoutTx(kk+1,:);%Adding the current Modulated output to the final OutPut
                    end
                else
                    EoutModAux = EoutTx(kk,:);
                    %EoutMod = EoutMod + EoutModAux + EoutTx(kk+1,:);       %Adding the current Modulated output to the final OutPut
                    EoutModTem(ObsCarrPos(kk),1:length(EoutModAux))   = EoutModAux;    %Adding the current Modulated output to the final OutPut
                    if ~mod(CarrPass,2)
                        EoutModTem(ObsCarrPos(kk+1),1:length(EoutModAux)) = EoutTx(kk+1,:);%Adding the current Modulated output to the final OutPut
                    end
                end
                %% Taking the sampling the EVM meassurement
                PosAuxEout      = (2*NumAmosCP+NPPB)/2:(2*NumAmosCP+...
                    NPPB):(length(EoutModAux)-StuffSampels);%Varriable respossible to take just the samples at the middle of the symbol
                IxAux           = EoutModAux.*conj(EoutModAux);        %Recovering the signal that will be transmited
                IxAux           = ifft(fft(IxAux).*EvmFilt);           %Removing higher signal generated by the receiving process
                IxAux           = IxAux./max(abs(IxAux));              %Normalizing the reference
                EvmMatRef(ObsCarrPos==kk,:) = IxAux(PosAuxEout);                   %Taking just the middle samples as references
            end
            if IfftOrSum
                Ekaux=sum(EoutModTem);
                figure;hold all;
                plot(f,20*log10(abs(fftshift(fft(Ekaux)./length(Ekaux)))));
                set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%                 if ~mod(RefCarr,2)
%                     [EoutMod,~,~] = OpticalIFFTNP(t,T,MaxNumStagT,EoutModTem);
%                 else
                    [EoutMod,~,~] = OpticalIFFTN(t,T,MaxNumStagT,EoutModTem);
%                 end
                plot(f,20*log10(abs(fftshift(fft(EoutMod)./length(EoutMod)))));
                axis([0.5*fin 1.2*NumCarr*fc+fin -140 -20]);
                a=1;
            else
                EoutMod = sum(EoutModTem);
            end
        otherwise
            %%        Generate the data OOK
            EoutMod    = 0;
            EoutModTem = 0;
            for kk=InitCarrDo:CarrPass:NumCarr%Generating different data for each carrier
                if CarrUsedDo(kk)
                    % Frist it is chosen to transmite just one high pulse for
                    %testing the channel...
                    if AdjusData
                        TxData        = zeros(1,Nb - NumBitDesc);                  %creating vector of zerros for correcting the delay caused by the transmission
                        TxData(end/2) = 1;                                         %Adding an strategical pulse for measuring the result.
                        TxDataMat(kk,:) = TxData;                                  %Storring the transmited information for latter evaluation
                        if NRZPolarity
                            TxData(TxData==0) = NrzMin;
                            TxData(TxData==1) = NrzMax;
                        end
                    else%... or just a random information
                        
                        TxData                     = (randi(2,1,Nb - NumBitDesc)...
                            -1);%Creating Random Information that will be loaded in each individual subcarrier
                        TxData(1:JusLen)           = JusVal;                       %Making the First 4 bits equal to zero
                        TxData(end-(JusLen-1):end) = JusVal;                       %Making the Last 4 bits equal to zero
                        TxDataMat(kk,:) = TxData;                                  %Storring the transmited information for latter evaluation
                        if NRZPolarity
                            TxData(TxData==0) = NrzMin;
                            TxData(TxData==1) = NrzMax;
                        end
                    end
                    
                    %The signal generated are not yet with the same number of
                    %samples as the OFCS loaded. These nexte lines do oversampling
                    TxDataRes = rectpulse(TxData,NPPB);                            %Changing the length of the Data acordingly with the time vector
                    
                    %Thus, if it would be required to add cycle prefix the number
                    %of samples per symbol needs to change as well as some
                    %adjustments needs to be done for the new signal match in size
                    %with the size of the vector time. This problem just exist on
                    %simulation, at practice the main point is the syncronism of
                    %the signals.
                    % Adding CP to the data
                    if AddCP
                        TxAux = reshape(TxDataRes,NPPB,Nb-NumBitDesc);
                        TxAux = [TxAux(1:NumAmosCP,:);TxAux;TxAux(end-(NumAmosCP...
                            -1):end,:)];
                        TxAux = reshape(TxAux,1,(2*NumAmosCP+NPPB)*(Nb-...
                            NumBitDesc));
                        TxDataRes = [TxAux TxAux(end-(StuffSampels-1):end)];
                    end
                    
                    
                    %Because of reasons, the real world does not have components
                    %capable instantly change the voltage level. Thus, to emulate
                    %this behaviour the eletrical PAM signal pass through a
                    %gaussian filter that will remove the frequencies of higher
                    %order, which will result in small slope in the level variation
                    
                    [BitFilt,~] = FiltroGaussiano(f,BWD,CenFeq,FiltOrd);           %Creating filter for conformation of the input information
                    BitFilt = fftshift(BitFilt);                                   %Doing a shift on the vector for matching the transmited data
                    TxSig = ifft(fft(TxDataRes).*BitFilt);                         %Conforming the information and Creating the modulation signal
                    %                 TxSigMat(kk,:) = TxSig;                                        %Storring the transmited information for latter evaluation
                    %             if kk==1
                    %                 PrintInfo(Ploting*5,t,TxSigMat(kk,:),TxDataRes);
                    %             end
                    %If the user choses to modulate all carriers at the same time
                    %there is no need of this script to generate data for each
                    %individual carrier. Therefore, this For loop can be halted
                    
                    
                    %Assigning the eletrical signal to one drive of the MZM -
                    %The aplitude of the signal can be controlled by the
                    %variable DatGai, which can be understood as an gain for
                    %the eletrical signal or an atenuation. The second signal
                    %will be similar with the only difference a phase shift of
                    %pi.
                    U.U1t = DatGai.*TxSig;
                    U.U2t = DatGai.*exp(-1j*pi).*TxSig;
                    %As both signals will have the mostrly the same
                    %characteristics with the only difference the phase shift
                    %of 180 degress. The MZM-I will be working on the Push-Pull
                    %configuration. It is necessary to reduce the Chirp noise
                    %to zero.
                    [EoutModAux,~]=Mach_Zehnder_Modulator_simplificado(t,...
                        EoutTx(kk,:),U,MZ_Input_File);%Modulating individual carriers
                    %EoutMod = EoutMod + EoutModAux + EoutTx(kk+1,:);       %Adding the current Modulated output to the final OutPut
                    EoutModTem(ObsCarrPos(kk),1:length(EoutModAux))   = EoutModAux;    %Adding the current Modulated output to the final OutPut
                    if ~mod(CarrPass,2)
                        EoutModTem(ObsCarrPos(kk+1),1:length(EoutModAux)) = EoutTx(kk+1,:);%Adding the current Modulated output to the final OutPut
                    end
                else
                    EoutModAux = EoutTx(kk,:);
                    %EoutMod = EoutMod + EoutModAux + EoutTx(kk+1,:);       %Adding the current Modulated output to the final OutPut
                    EoutModTem(ObsCarrPos(kk),1:length(EoutModAux))   = EoutModAux;    %Adding the current Modulated output to the final OutPut
                    if ~mod(CarrPass,2)
                        EoutModTem(ObsCarrPos(kk+1),1:length(EoutModAux)) = EoutTx(kk+1,:);%Adding the current Modulated output to the final OutPut
                    end
                end
                %% Taking the sampling the EVM meassurement
                PosAuxEout      = (2*NumAmosCP+NPPB)/2:(2*NumAmosCP+...
                    NPPB):(length(EoutModAux)-StuffSampels);%Varriable respossible to take just the samples at the middle of the symbol
                IxAux           = EoutModAux.*conj(EoutModAux);        %Recovering the signal that will be transmited
                IxAux           = ifft(fft(IxAux).*EvmFilt);           %Removing higher signal generated by the receiving process
                IxAux           = IxAux./max(abs(IxAux));              %Normalizing the reference
                EvmMatRef(ObsCarrPos==kk,:) = IxAux(PosAuxEout);                   %Taking just the middle samples as references
            end
            if IfftOrSum
                Ekaux=sum(EoutModTem);
                figure;hold all;
                plot(f,20*log10(abs(fftshift(fft(Ekaux)./length(Ekaux)))));
                set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%                 if ~mod(RefCarr,2)
%                     [EoutMod,~,~] = OpticalIFFTNP(t,T,MaxNumStagT,EoutModTem);
%                 else
                    [EoutMod,~,~] = OpticalIFFTN(t,T,MaxNumStagT,EoutModTem);
%                 end
                plot(f,20*log10(abs(fftshift(fft(EoutMod)./length(EoutMod)))));
                axis([0.5*fin 1.2*NumCarr*fc+fin -140 -20]);
                a=1;
            else
                EoutMod = sum(EoutModTem);
            end
    end
    
    
    %%   Modulation of Given UOFCS
    % After the OFC acquizition and the generation of the modulation signal.
    % The folloiwing step is to insert the Data inside each sub carriers. It
    % process can be individual or all carriers at once.
    
    %%  Actualy modulating the Data into each carrier
    switch Modulation
        case 'DPSK'
            %%                  DPSK Modulation
            %Just before sent the carriers and the modulating signal this
            %script verify if the user set to modulate individual carriers or
            %all at once.
            EoutMod = 0;
            if ModAll%Modulating all carrier at once...
                U.U1t = TxSigMat(1,:);                                         %Assigning the eletrical signal to one drive of the MZM
                U.U2t = exp(-1j*pi).*TxSigMat(1,:);                            %Assigning the same eletrica signal but with a phase shift of pi
                %As both signals will have the mostrly the same characteristics
                %with the only difference the phase shift of 180 degress. The
                %MZM-I will be working on the Push-Pull configuration. It is
                %necessary to reduce the Chirp noise to zero.
                [EoutModAux] = Mach_Zehnder_Modulator_simplificado(t,EoutT,U...
                    ,MZ_Input_File);%Modulating all carriers at the same time for simplicity
                %             PrintInfo(Ploting*6,real(EoutModAux),T,NPPB);
                EoutMod = EoutMod + EoutModAux;                                %Adding the current Modulated output to the final OutPut
            else%... or each one at the time
                %                 EoutTx=SelectEachCarrier(EoutT,NumCarr,f,fin,FBWD,Order);      %This function is resposible to split each carrier
                for kk=1:2:NumCarr
                    U.U1t = TxSigMat(kk,:);                                    %Assigning the eletrical signal to one drive of the MZM
                    U.U2t = exp(-1j*pi).*TxSigMat(kk,:);                       %Assigning the same eletrica signal but with a phase shift of pi
                    %As both signals will have the mostrly the same
                    %characteristics with the only difference the phase shift
                    %of 180 degress. The MZM-I will be working on the Push-Pull
                    %configuration. It is necessary to reduce the Chirp noise
                    %to zero.
                    [EoutModAux,~]=Mach_Zehnder_Modulator_simplificado(t,...
                        EoutTx(kk,:),U,MZ_Input_File);%Modulating individual carriers
                    EoutMod = EoutMod + EoutModAux + EoutTx(kk+1,:);           %Adding the current Modulated output to the final OutPut
                    %                 if kk == 1
                    %                     PrintInfo(Ploting*6,abs(EoutModAux),T,NPPB);
                    %                 end
                end
            end
        case 'DQPSK'
            %%            DQPSK Modulation
            %Just before sent the carriers and the modulating signal this
            %script verify if the user set to modulate individual carriers or
            %all at once.
            EoutMod = 0;
            if ModAll%Modulating all carrier at once...
                [EoutModAux] = IqMod(EoutT,DatGai.*TxSigMat(1,:),DatGai.*...
                    TxSigMat(1+NumCarr,:),Vpi,V0);
                %             PrintInfo(Ploting*7,abs(EoutModAux),T,NPPB);
                EoutMod = EoutMod + EoutModAux;                                %Adding the current Modulated output to the final OutPut
            else%... or each one at the time
                %                 EoutTx=SelectEachCarrier(EoutT,NumCarr,f,fin,FBWD,Order);      %This function is resposible to split each carrier
                for kk=1:2:NumCarr
                    [EoutModAux] = IqMod(EoutTx(kk,:),DatGai.*TxSigMat(kk,:)...
                        ,DatGai.*TxSigMat(kk+NumCarr,:),Vpi,V0);
                    %                 if kk == 1
                    %                     PrintInfo(Ploting*7,abs(EoutModAux),T,NPPB);
                    %                 end
                    EoutMod = EoutMod + EoutModAux + EoutTx(kk+1,:);           %Adding the current Modulated output to the final OutPut
                end
            end
            %         case '4PAM'
            %             %%            4PAM Modulation
            %             %Just before sent the carriers and the modulating signal this
            %             %script verify if the user set to modulate individual carriers or
            %             %all at once.
            %             EoutMod = 0;
            %             if ModAll%Modulating all carrier at once...
            %                 U.U1t = TxSigMat(1,:);                                         %Assigning the electrical signal to one drive of the MZM
            %                 U.U2t = TxSigMat(1+NumCarr,:);                                 %Assigning the electrical signal to another drive of the MZM
            %     %             EoutModAux = PamMod(EoutTx(kk,:),U.U1t,U.U2t,U_pi2,Vbias,t);
            %                 [EoutModAux,~]=Mach_Zehnder_Modulator_simplificado(t,EoutT,U...
            %                                                                ,MZ_Input_File);%Modulating all carriers at the same time for simplicity
            %     %             PrintInfo(Ploting*8,EoutModAux.*conj(EoutModAux)./...
            %     %                             max(abs(EoutModAux.*conj(EoutModAux))),T,NPPB);
            %                 EoutMod = EoutMod + EoutModAux;                                %Adding the current Modulated output to the final OutPut
            %             else%... or each one at the time
            %             end
        otherwise
            %             %%            OOK Modulation
            %             %Just before sent the carriers and the modulating signal this
            %             %script verify if the user set to modulate individual carriers or
            %             %all at once.
            %             EoutMod = 0;
            %             if ModAll%Modulating all carrier at once...
            %                 %Assigning the eletrical signal to one drive of the MZM - The
            %                 %aplitude of the signal can be controlled by the variable
            %                 %DatGai, which can be understood as an gain for the eletrical
            %                 %signal or an atenuation. The second signal will be similar
            %                 %with the only difference a phase shift of pi.
            %                 U.U1t = DatGai.*TxSigMat(1,:);
            %                 U.U2t = DatGai.*exp(-1j*pi).*TxSigMat(1,:);
            %                     %As both signals will have the mostrly the same
            %                     %characteristics with the only difference the phase shift
            %                     %of 180 degress. The MZM-I will be working on the Push-Pull
            %                     %configuration. It is necessary to reduce the Chirp noise
            %                     %to zero.
            %                 [EoutModAux,~]=Mach_Zehnder_Modulator_simplificado(t,EoutT...
            %                                                              ,U,MZ_Input_File);%Modulating all carriers at the same time for simplicity
            %     %             PrintInfo(Ploting*9,abs(EoutModAux),T,NPPB);
            %                 EoutMod = EoutMod + EoutModAux;                                %Adding the current Modulated output to the final OutPut
            %             else%... or each one at the time
            % %                 EoutTx=SelectEachCarrier(EoutT,NumCarr,f,fin,FBWD,Order);      %This function is resposible to split each carrier
            %                 for kk=1:2:NumCarr
            %                     %Assigning the eletrical signal to one drive of the MZM -
            %                     %The aplitude of the signal can be controlled by the
            %                     %variable DatGai, which can be understood as an gain for
            %                     %the eletrical signal or an atenuation. The second signal
            %                     %will be similar with the only difference a phase shift of
            %                     %pi.
            %                     U.U1t = DatGai.*TxSigMat(kk,:);
            %                     U.U2t = DatGai.*exp(-1j*pi).*TxSigMat(kk,:);
            %                     %As both signals will have the mostrly the same
            %                     %characteristics with the only difference the phase shift
            %                     %of 180 degress. The MZM-I will be working on the Push-Pull
            %                     %configuration. It is necessary to reduce the Chirp noise
            %                     %to zero.
            %                     [EoutModAux,~]=Mach_Zehnder_Modulator_simplificado(t,...
            %                                                  EoutTx(kk,:),U,MZ_Input_File);%Modulating individual carriers
            %                     EoutMod = EoutMod + EoutModAux + EoutTx(kk+1,:);           %Adding the current Modulated output to the final OutPut
            %                     %% Taking the sampling the EVM meassurement
            %                     PosAuxEout      = (2*NumAmosCP+NPPB)/2:(2*NumAmosCP+...
            %                                    NPPB):(length(EoutModAux)-StuffSampels);%Varriable respossible to take just the samples at the middle of the symbol
            %                     IxAux           = EoutModAux.*conj(EoutModAux);        %Recovering the signal that will be transmited
            %                     IxAux           = ifft(fft(IxAux).*EvmFilt);           %Removing higher signal generated by the receiving process
            %                     IxAux           = IxAux./max(abs(IxAux));              %Normalizing the reference
            %                     EvmMatRef(kk,:) = IxAux(PosAuxEout);                   %Taking just the middle samples as references
            %     %                 if kk == 1
            %     %                     figure;plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));axis([kk*fc-1.1*fc kk*fc+1.1*fc max(20*log10(abs(fft(EoutModAux)./length(EoutModAux))))-100 max(20*log10(abs(fft(EoutModAux)./length(EoutModAux))))+10]);set(gcf,'units','normalized','outerposition',[0 0 1 1]);%,f,20*log10(abs(fftshift(fft(AA)./length(AA)))))
            %     %                     PrintInfo(PlotingThis*9,abs(EoutModAux(1:end-StuffSampels)),t(2*NumAmosCP+NPPB),(2*NumAmosCP+NPPB));
            %     %                     a=0;
            %     %                 end
            %                 end
            %             end
    end
else
    EoutTx = SelectEachCarrier(Eout,NumCarr,f,fin,FBWD,Order);
    EoutMod = 0;
    for kk=1:NumCarr
        EoutMod = EoutMod + EoutTx(kk,:);
    end
end
% figure;
% plot(f,20*log10(abs(fftshift(fft(EoutMod)./length(EoutMod)))));
% set(gcf,'units','normalized','outerposition',[0 0 1 1]);
% Ploting result for an qualitative analizes
% PrintInfo(Ploting*10,f,db(abs((fftshift(fft(Eout))))/length(Eout)),db(...
% abs((fftshift(fft(EoutT))))/length(EoutT)),db(abs((fftshift(fft(EoutMod)...
%                              )))/length(EoutMod)),db(fftshift(SelecFilt)));

%%   Transmission of the OFDM Symble through a channel
% Having data stored and ready to be sent to end user. At the stage this
% script is responsible to chose the medium where this signal will travel.
% It may be withing an optical fiber or Back-toBack transmission.
% [~,PdBm] = MeasPower(EoutMod)
EoutMod = EoutMod/SplitRatio;
% [~,PdBm] = MeasPower(EoutMod)
switch Medium
    case 'B2B'
        EoutRec = EoutMod;
    case 'Fiber'
        [EoutRec,~,~]=Fibra_Monomodo1(t,EoutMod,lambdac,T,FiberLength,0,f);
        %         PrintInfo(Ploting*11,EoutRec.*conj(EoutRec)./max(abs(EoutRec.*...
        %                                                    conj(EoutRec))),T,NPPB);
    otherwise
        EoutRec = EoutMod;
end


clear TxSigMat TxSig TxSigI TxSigQ EoutAux1 VetThisCarr EoutMod Ix...
    EoutModAux Esync Esync1 Esync2 EoutA EoutB EoutC EoutD EoutModTem
% clearvars -except fc Rb NPPB fsample T Ta NumberOf_T FinalTime Nb NbDPSK...
% TotalSamples t f Eout EoutMod EoutRec TxDataMat NumCarr NbDQPSK ValsLev...
% Modulation Polirized MaxAmp4PAM NumbOfTest UsedModula TestBundle Nb4Pam ...
% ThisModula CurrentTest CurrentModula ThisTest ThisModula OfcName BerOOK ...
% Ber4PAM BerDQPSK BerDPSK FiberDelay Medium CurrentMedium OSNRPC AddCP V0...
% OcsToTest CurrentOCS PulseResp NumAmosCP NumBitDesc StuffSampels JusPos ...
% contcp Ber4PAM PAM4Type Vpi MZ_Input_File ModAll JusValEnd NRZPolarity ...
% VetElecPower AberLev fin FFTSplit PastTime PowRef SnrRef CarSNR ValsLev2...
% AberLevI ValsLevI AberLevQ ValsLevQ VetElecPowerQ VetElecPowerI SigAten ...
%                                      lambdac Atenuation VetOptiPower JusVal