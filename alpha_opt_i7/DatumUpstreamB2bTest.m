%c
%c                                                       ..'Ž`'..'Ž`..'Ž`..
%c       File: DatumUpstream File For the All Optical OFDM
%c(Main file to call each idividual parameters)
%c
%c     This main code is resposible to call and run the all the fuction to
%c related to this simulation. Once the signal was received by the user the
%c next step is to modulate this information and send it back to the OLT.
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
DatumUpstreamInputData;
clear TxData TxDataMat EoutAux1 VetThisCarr

[EoutAux1,~,VetThisCarr]=OpticalFFTN(t,T,MaxNumStag,EoutRec);
clear EoutRec
% PrintInfo(PlotingThis*12,f,db(abs((fftshift(fft(EoutRec))))/length(EoutRec)),EoutAux1,RefCarr*fc);
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

switch Modulation
    case 'DPSK'
        %%            Generate the data DPSK
        EoutMod    = 0;
        EoutModTem = 0;
        for kk=InitCarrUp:CarrPass:NumCarr%Generating different data for each carrier
            if CarrUsedUp(kk)
                % Frist it is chosen to transmite a random information...
                TxData = (randi(2,1,NbDPSK)-1);                                %Creating the data stream to be transmited
                TxData(1:JusPos) = JusVal;                                     %Adding the Justification simble at the begining of the stream to sincronize received data frame
                TxData(end-(JusPos-1):end) = JusValEnd;                        %Adding the Justification simble at the end of the stream to sincronize received data frame

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
%                 TxSig       = ifft(fft(SigTx).*BitFilt);                       %Conforming the information and Creating the modulation signal
                TxSig       = SigTx;                                   %Adding a gain to the eletrical signal
                %             PrintInfo(Ploting*2,TxSig,T,NPPB);
%                 TxSigMat(kk,:) = TxSig;                                        %Storing the eletrical signal for further evaluation
%                 a=4;
                
                EoutT = EoutAux1(VetThisCarr==ObsCarrUsed(kk),:);
                %Assigning the eletrical signal to one drive of the MZM -
                %The aplitude of the signal can be controlled by the
                %variable DatGai, which can be understood as an gain for
                %the eletrical signal or an atenuation. The second signal
                %will be similar with the only difference a phase shift of
                %pi.
                U.U1t  = TxSig;
                U.U2t  = exp(-1j*pi).*TxSig;
                %As both signals will have the mostrly the same
                %characteristics with the only difference the phase shift
                %of 180 degress. The MZM-I will be working on the Push-Pull
                %configuration. It is necessary to reduce the Chirp noise
                %to zero.
                [EoutModAux,~]=Mach_Zehnder_Modulator_simplificado(t,...
                    EoutT,U,MZ_Input_File);    %Modulating individual carriers
%                 EoutMod = EoutMod + EoutModAux;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                EoutModTem(ObsCarrPos(kk),1:length(EoutModAux)) = EoutModAux;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                %% Taking the sampling the EVM meassurement
                PhaDel          = 0;
                TimDel          = T;
                [ESync1,ESync2] = DelayInterf(t,TimDel,PhaDel,EoutModAux); %Coverting phase shift to amplitude variation
                ESync1 = ESync1.*conj(ESync1);
                ESync2 = ESync2.*conj(ESync2);
                Esync  = ESync2-ESync1;
                PosAuxEout      = (2*NumAmosCP+NPPB)/2:(2*NumAmosCP+...
                               NPPB):(length(Esync)-StuffSampels);%Varriable respossible to take just the samples at the middle of the symbol
                Esync           = ifft(fft(Esync).*EvmFilt);           %Removing higher signal generated by the receiving process
%                     Esync  = Esync-min(Esync);
                Esync  = Esync-mean(Esync);
                Esync  = Esync./max(abs(Esync));
                IxAux  = Esync(PosAuxEout);  
                a=0;
                EvmMatRef(ObsCarrPos==kk,:) = IxAux;                   %Taking just the middle samples as references
%                 if kk==2
%                     figure;hold all;grid on;
%                     plot(IxAux,zeros(1,length(IxAux)),'rx');
%                     set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%                     a=0;
%                 end
            else
                EoutT = EoutAux1(VetThisCarr==ObsCarrPos(kk),:);
                EoutModAux = EoutT;
%                 EoutMod = EoutMod + EoutModAux;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                EoutModTem(ObsCarrPos(kk),1:length(EoutModAux)) = EoutModAux;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
            end
        end
        if IfftOrSum
%             Ekaux=sum(EoutModTem);
%             figure;hold all;
%             plot(f,20*log10(abs(fftshift(fft(Ekaux)./length(Ekaux)))));
%             set(gcf,'units','normalized','outerposition',[0 0 1 1]);
            if mod(size(EoutModTem,1),2)
                EoutModTem(end+1,:) = 0;
            end
            if size(EoutModTem,1)>1
                [EoutMod,~,~] = OpticalIFFTN(t,T,MaxNumStagT,EoutModTem);
            else
                EoutMod = EoutModTem;
            end
%             [EoutMod,~,~] = OpticalIFFTN(t,T,MaxNumStagT,EoutModTem);
%             plot(f,20*log10(abs(fftshift(fft(EoutMod)./length(EoutMod)))));
%             axis([0.9*fin 1.2*NumCarr*fc+fin -280 -40]);
%             a=1;
        else
            EoutMod = sum(EoutModTem);
        end
    case 'DQPSK'
        %%            Generate the data DQPSK
        EoutMod    = 0;
        EoutModTem = 0;
        for kk=InitCarrUp:CarrPass:NumCarr%Generating different data for each carrier
            if CarrUsedUp(kk)
                % Frist it is chosen to transmite a random information...
                TxData                     = (randi(2,1,NbDQPSK)-1);           %Creating the data stream to be transmited
                TxData(1:JusPos)           = JusVal;                           %Adding the Justification simble at the begining of the stream to sincronize received data frame
                TxData(end-(JusPos-1):end) = JusValEnd;                        %Adding the Justification simble at the end of the stream to sincronize received data frame
                PreBits                    = [0 0];                            %Seting the start point for the DQPSK maping

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
%                 TxSigI = ifft(fft(SigI).*BitFilt);                             %Conforming the information and Creating the modulation signal
%                 TxSigQ = ifft(fft(SigQ).*BitFilt);                             %Conforming the information and Creating the modulation signal
                TxSigI = SigI;                             %Conforming the information and Creating the modulation signal
                TxSigQ = SigQ;                             %Conforming the information and Creating the modulation signal
                %             PrintInfo(Ploting*3,TxSigI,T,NPPB,TxSigQ);
%                 TxSigMat(kk,:)         = TxSigI;                               %Storing the eletrical signal for further evaluation
%                 TxSigMat(kk+NumCarr,:) = TxSigQ;                               %Storing the eletrical signal for further evaluation
                EoutT = EoutAux1(VetThisCarr==ObsCarrUsed(kk),:);
                [EoutModAux] = IqMod(EoutT,TxSigI,TxSigQ,Vpi,V0);
%                 a=4;
                EoutModTem(ObsCarrPos(kk),1:length(EoutModAux)) = EoutModAux;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                %% Taking the sampling the EVM meassurement
                PhaDel = 1*pi/4;
                TimDel = T;
                [EoutA,EoutB] = DelayInterf(t,TimDel,PhaDel,EoutModAux);%Coverting phase shift to amplitude variation
                PhaDel = -1*pi/4;
                TimDel = T;
                [EoutC,EoutD] = DelayInterf(t,TimDel,PhaDel,EoutModAux);%Coverting phase shift to amplitude variation
                EoutA = EoutA.*conj(EoutA);
                EoutB = EoutB.*conj(EoutB);
                EoutC = EoutC.*conj(EoutC);
                EoutD = EoutD.*conj(EoutD);                   %Shifting the filter for matching the received signal
                EoutA = ifft(fft(EoutA).*EvmFilt);
                EoutB = ifft(fft(EoutB).*EvmFilt);
                EoutC = ifft(fft(EoutC).*EvmFilt);
                EoutD = ifft(fft(EoutD).*EvmFilt);
                EoutI = (EoutB - EoutA);
                EoutQ = (EoutD - EoutC);
                EoutI = EoutI-mean(EoutI);
                EoutQ = EoutQ-mean(EoutQ);
                EoutI = EoutI./max(abs(EoutI));                                %Normalizing the signal
                EoutQ = EoutQ./max(abs(EoutQ));                                %Normalizing the signal
                PosAuxEout      = (2*NumAmosCP+NPPB)/2:(2*NumAmosCP+...
                               NPPB):(length(EoutModAux)-StuffSampels);%Varriable respossible to take just the samples at the middle of the symbol
                IxAux = EoutI(PosAuxEout) + 1j.*EoutQ(PosAuxEout);
                EvmMatRef(ObsCarrPos==kk,:) = IxAux;       %Taking just the middle samples as references
%                 if kk==2
%                     figure;hold all;grid on;
%                     plot(real(IxAux),imag(IxAux),'rx');
%                     set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%                     a=0;
%                 end
            else
                EoutT = EoutAux1(VetThisCarr==ObsCarrPos(kk),:);
                EoutModAux = EoutT;
%                 EoutMod = EoutMod + EoutModAux;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                EoutModTem(ObsCarrPos(kk),1:length(EoutModAux)) = EoutModAux;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
            end
        end
        if IfftOrSum
%             Ekaux=sum(EoutModTem);
%             figure;hold all;
%             plot(f,20*log10(abs(fftshift(fft(Ekaux)./length(Ekaux)))));
%             set(gcf,'units','normalized','outerposition',[0 0 1 1]);
            if mod(size(EoutModTem,1),2)
                EoutModTem(end+1,:) = 0;
            end
%             [EoutMod,~,~] = OpticalIFFTN(t,T,MaxNumStagT,EoutModTem);
            if size(EoutModTem,1)>1
                [EoutMod,~,~] = OpticalIFFTN(t,T,MaxNumStagT,EoutModTem);
            else
                EoutMod = EoutModTem;
            end
%             plot(f,20*log10(abs(fftshift(fft(EoutMod)./length(EoutMod)))));
%             axis([0.9*fin 1.2*NumCarr*fc+fin -280 -40]);
%             a=1;
        else
            EoutMod = sum(EoutModTem);
        end
    case '4PAM'
        %%        Generate the data 4PAM
        EoutMod    = 0;
        EoutModTem = 0;
        for kk=InitCarrUp:CarrPass:NumCarr%Generating different data for each carrier
            % Frist it is chosen to transmite a random information...
            if CarrUsedUp(kk)
                TxData                     = (randi(2,1,Nb4Pam)-1);        %Creating the data stream to be transmited
                TxData(1:JusPos)           = JusVal;                       %Adding the Justification simble at the begining of the stream to sincronize received data frame
                TxData(end-(JusPos-1):end) = JusValEnd;                    %Adding the Justification simble at the end of the stream to sincronize received data frame
                %             else
                %                 TxData                     = zeros(1,Nb4Pam);              %Creating the data stream to be transmited
                %             end
                TxDataMat(kk,:) = TxData;                                      %Storring the transmited information for latter evaluation
                %Once that we have the data to be transmited it needs to be
                %converted to the eletrical signal that hereafter will modulate
                %our carrier through the MZM-I
                
                %Selecting if the PAM will be done on Electrical or Optical
                %domain
                if Which4PAM
                    %  Good for others kms
                    [Phi1,Phi2] = Maping4PamIq(TxData,Vmin,Vmax,ModSchem...
                                                ,FiberLength,SetCpSampZer);%Generating the eletrical signal for the optical PAM4
                    
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
                    if SetCpSampZer
                        TxAux1 = [zeros(NumAmosCP,size(TxAux1,2));TxAux1...
                                         ;zeros(NumAmosCP,size(TxAux1,2))];
                    else
                        TxAux1 = [TxAux1(1:NumAmosCP,:);TxAux1;TxAux1(...
                                                 end-(NumAmosCP-1):end,:)];
                    end
                    TxAux1 = reshape(TxAux1,1,(2*NumAmosCP+NPPB)*Nb4Pam/2);
                    TxSig1 = [TxAux1 TxAux1(end-(StuffSampels-1):end)];
                    TxAux2 = reshape(TxSig2,NPPB,Nb4Pam/2);
                    if SetCpSampZer
                        TxAux2 = [zeros(NumAmosCP,size(TxAux2,2));TxAux2...
                                         ;zeros(NumAmosCP,size(TxAux2,2))];
                    else
                        TxAux2 = [TxAux2(1:NumAmosCP,:);TxAux2;TxAux2(...
                                                 end-(NumAmosCP-1):end,:)];
                    end
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
%                 TxSig1                 = ifft(fft(TxSig1).*BitFilt);           %Conforming the information and Creating the modulation signal
%                 TxSig2                 = ifft(fft(TxSig2).*BitFilt);           %Conforming the information and Creating the modulation signal
                %             TxSigMat(kk,:)         = TxSig1;                               %Storing the eletrical signal for further evaluation
                %             TxSigMat(kk+NumCarr,:) = TxSig2;                               %Storing the eletrical signal for further evaluation
                
                EoutT = EoutAux1(VetThisCarr==ObsCarrUsed(kk),:);
                U.U1t  = TxSig1;                                               %Assigning the electrical signal to one drive of the MZM
                U.U2t  = TxSig2;                                               %Assigning the electrical signal to another drive of the MZM
                if Which4PAM
                    if ModSchem
%                         [EoutModAux] = IqMod4Pam (EoutT,U.U1t,U.U2t,U_pi2,Vbias);
                        [EoutModAux] = Mach_Zehnder_Modulator_simplificado(t,EoutT,U,MZ_Input_File);
                    else
                        [EoutModAux] = Mach_Zehnder_Modulator_simplificado(t,EoutT,U,MZ_Input_File);
                    end
                else
                    [EoutModAux] = Mach_Zehnder_Modulator_simplificado(t...
                        ,EoutT,U,MZ_Input_File);
                end
%                 EoutMod = EoutMod + EoutModAux;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                EoutModTem(ObsCarrPos(kk),1:length(EoutModAux)) = EoutModAux;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                %% Taking the sampling the EVM meassurement
                PosAuxEout = (2*NumAmosCP+NPPB)/2:(2*NumAmosCP+NPPB)...
                                    :(length(EoutModAux)-StuffSampels);%Varriable respossible to take just the samples at the middle of the symbol
                Ix         = EoutModAux.*conj(EoutModAux);             %Recovering the signal that will be transmited
                Ix         = ifft(fft(Ix).*EvmFilt);                   %Removing higher signal generated by the receiving process
                Ix         = Ix - min(Ix);
                Ix         = Ix./max(abs(Ix));                         %Normalizing the reference
                IxAux      = Ix(PosAuxEout);
                a=0;
                EvmMatRef(ObsCarrPos==kk,:) = IxAux;                   %Taking just the middle samples as references
%                 if kk==2
%                     figure;hold all;grid on;
%                     plot(t,Ix);plot(t(PosAuxEout),IxAux,'x');
%                     set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%                     figure;hold all;grid on;
%                     plot(IxAux,zeros(1,length(IxAux)),'rx');
%                     set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%                     a=0;
%                 end
            else
                EoutT = EoutAux1(VetThisCarr==ObsCarrPos(kk),:);
                EoutModAux = EoutT;
%                 EoutMod = EoutMod + EoutModAux;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                EoutModTem(ObsCarrPos(kk),1:length(EoutModAux)) = EoutModAux;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
            end
%             if kk==2
%                 PrintInfo(PlotingThis*4,abs(EoutModAux(1,1:end-StuffSampels)).^2./max(abs(EoutModAux(1,1:end-StuffSampels)).^2),2*t(2*NumAmosCP+NPPB),2*NumAmosCP+NPPB,abs(EoutModAux(1:end-StuffSampels)).^2./max(abs(EoutModAux(1,1:end-StuffSampels)).^2));
%                 a=4;
%             end
        end
        if IfftOrSum
%             Ekaux=sum(EoutModTem);
%             figure;hold all;
%             plot(f,20*log10(abs(fftshift(fft(Ekaux)./length(Ekaux)))));
%             set(gcf,'units','normalized','outerposition',[0 0 1 1]);
% %             if ~mod(RefCarr,2)
% %                 [EoutMod,~,~] = OpticalIFFTNP(t,T,MaxNumStagT,EoutModTem);
% %             else
            if mod(size(EoutModTem,1),2)
                EoutModTem(end+1,:) = 0;
            end
            if size(EoutModTem,1)>1
                [EoutMod,~,~] = OpticalIFFTN(t,T,MaxNumStagT,EoutModTem);
            else
                EoutMod = EoutModTem;
            end
% %             end
%             plot(f,20*log10(abs(fftshift(fft(EoutMod)./length(EoutMod)))));
%             axis([0.9*fin 1.2*NumCarr*fc+fin -280 -40]);
%             a=1;
        else
            EoutMod = sum(EoutModTem);
        end
    otherwise
        %%        Generate the data OOK
        EoutMod    = 0;
        EoutModTem = 0;
        for kk=InitCarrUp:CarrPass:NumCarr%Generating different data for each carrier
            % Frist it is chosen to transmite just one high pulse for
            %testing the channel...
            if CarrUsedUp(kk)
                TxData                     = (randi(2,1,Nb - NumBitDesc)-1);   %Creating Random Information that will be loaded in each individual subcarrier
                TxData(1:JusLen)           = JusVal;                           %Making the First 4 bits equal to zero
                TxData(end-(JusLen-1):end) = JusVal;                           %Making the Last 4 bits equal to zero
                
                
                TxDataMat(kk,:)            = TxData;                           %Storring the transmited information for latter evaluation
                if NRZPolarity
                    TxData(TxData==0)      = NrzMin;
                    TxData(TxData==1)      = NrzMax;
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
%                 TxSig = ifft(fft(TxDataRes).*BitFilt);                         %Conforming the information and Creating the modulation signal
                TxSig = TxDataRes;                         %Conforming the information and Creating the modulation signal
                %             TxSigMat(kk,:) = TxSig;                                        %Storring the transmited information for latter evaluation
                
                EoutT = EoutAux1(VetThisCarr==ObsCarrUsed(kk),:);
                %Assigning the eletrical signal to one drive of the MZM -
                %The aplitude of the signal can be controlled by the
                %variable DatGai, which can be understood as an gain for
                %the eletrical signal or an atenuation. The second signal
                %will be similar with the only difference a phase shift of
                %pi.
                U.U1t  = DatGai.*TxSig;
                U.U2t  = DatGai.*exp(-1j*pi).*TxSig;
                %As both signals will have the mostrly the same
                %characteristics with the only difference the phase shift
                %of 180 degress. The MZM-I will be working on the Push-Pull
                %configuration. It is necessary to reduce the Chirp noise
                %to zero.
                [EoutModAux,~]=Mach_Zehnder_Modulator_simplificado(t,...
                    EoutT,U,MZ_Input_File);    %Modulating individual carriers
%                 EoutMod = EoutMod + EoutModAux;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                EoutModTem(ObsCarrPos(kk),1:length(EoutModAux)) = EoutModAux;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                %% Taking the sampling the EVM meassurement
                PosAuxEout = (2*NumAmosCP+NPPB)/2:(2*NumAmosCP+NPPB)...
                                    :(length(EoutModAux)-StuffSampels);%Varriable respossible to take just the samples at the middle of the symbol
                Ix         = EoutModAux.*conj(EoutModAux);             %Recovering the signal that will be transmited
                Ix         = ifft(fft(Ix).*EvmFilt);                   %Removing higher signal generated by the receiving process
                Ix         = Ix - min(Ix);
                Ix         = Ix./max(abs(Ix));                         %Normalizing the reference
                IxAux      = Ix(PosAuxEout);                           %Normalizing the reference
                a=0;
                EvmMatRef(ObsCarrPos==kk,:) = IxAux;                   %Taking just the middle samples as references
%                 if kk==2
%                     figure;hold all;grid on;
%                     plot(IxAux,zeros(1,length(IxAux)),'rx');
%                     set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%                     a=0;
%                 end
            else
                EoutT = EoutAux1(VetThisCarr==ObsCarrPos(kk),:);
                EoutModAux = EoutT;
%                 EoutMod = EoutMod + EoutModAux;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                EoutModTem(ObsCarrPos(kk),1:length(EoutModAux)) = EoutModAux;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
            end
        end
        if IfftOrSum
%             Ekaux=sum(EoutModTem);
%             figure;hold all;
%             plot(f,20*log10(abs(fftshift(fft(Ekaux)./length(Ekaux)))));
%             set(gcf,'units','normalized','outerposition',[0 0 1 1]);
            if mod(size(EoutModTem,1),2)
                EoutModTem(end+1,:) = 0;
            end
            if size(EoutModTem,1)>1
                [EoutMod,~,~] = OpticalIFFTN(t,T,MaxNumStagT,EoutModTem);
            else
                EoutMod = EoutModTem;
            end
%             plot(f,20*log10(abs(fftshift(fft(EoutMod)./length(EoutMod)))));
%             axis([0.9*fin 1.2*NumCarr*fc+fin -280 -40]);
%             a=1;
        else
            EoutMod = sum(EoutModTem);
        end
end

%%   Transmission of the OFDM Symble through a channel
% Having data stored and ready to be sent to end user. At the stage this
% script is responsible to chose the medium where this signal will travel.
% It may be withing an optical fiber or Back-toBack transmission.
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


if AddingNoiseF
    EsN0 = 10*log10(snr);
    [~,sigpower] = MeasPower(EoutRec);
    EoutRec = awgn(EoutRec,osnrf,sigpower-30);
    [~,sigpower] = MeasPower(EoutRec);
end

clear TxSigMat TxSig TxSigI TxSigQ EoutAux1 VetThisCarr EoutMod Ix Esync...
    Esync1 Esync2 EoutA EoutB EoutC EoutD EoutModTem
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