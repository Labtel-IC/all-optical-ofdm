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
                    case 'OFDM' %NOT IMPLEMENTED YET
                        EoutModTem = zeros(NumCarr,NPPB*Nb);
                        for kk=InitCarrUp:CarrPass:NumCarr                             %Generating different data for each carrier
                            %%U1t = 0;
                            %                 CpPaFr = 1;
                            %                 for CarrOffSet=-1*(NumFraPar-1)/2:1:(NumFraPar-1)/2
                            UniDmtMve = unique(DmtMve);
                            UniDmtMve = fliplr(UniDmtMve);
                            for CarDmtM=1:length(UniDmtMve)
                                M = UniDmtMve(CarDmtM);
                                DaSiAu = sum(ismember(DmtMve,M));
                                TxData          = randi([0 M-1],NumFra,DaSiAu);                %Generation random information
                                TxDataMat(kk,1+DaSiAu*NumFra*(CarDmtM-1):DaSiAu*NumFra*CarDmtM) = TxData(:);                               %Saving data for future comparison
%                             for CarDmtM=1:length(DmtMve)
%                                 M = DmtMve(CarDmtM);
%                                 TxData          = randi([0 M-1],NumFra,1);                %Generation random information
%                                 TxDataMat(kk,1+NumFra*(CarDmtM-1):NumFra*CarDmtM) = TxData(:);                               %Saving data for future comparison
                                switch OfdMod                                              %Sellect which modulation was used
                                    case 'qam'
                                        TxSymbAux(1:NumFra,1+DaSiAu*(CarDmtM-1):DaSiAu*CarDmtM)  = qammod(TxData,M);                    %Modulating information by QAM
                                    otherwise
                                        TxSymbAux(1:NumFra,1+DaSiAu*(CarDmtM-1):DaSiAu*CarDmtM)  = dpskmod(TxData,M);                   %Modulating information DPSK as default
                                end
                            end
                            TxSymb      = TxSymbAux.';
                            TxSigOfdm(CurrentTest,ObsCarrPos==kk,:) = TxSymb(:);
                            EvmMatRef(ObsCarrPos==kk,:) = TxSymb(:);                   %Taking just the middle samples as references
                            %Construction the Hermitian Simetry
                            TxSymbConj      = conj(flipud(TxSymb));                    %Singnal conjugate format at reverse order
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
                            %                 plot(TxSymbC);
                            %Here was implemented the modulation through FFT process, it basica mix-up
                            %all infomation following Fourier transform.
                            TxSymbMod       = ifft(TxSymbC,NFFT);
                            TxSymbMod       = rectpulse(TxSymbMod,OvSam);         %Over sampling
                            tt = linspace(0,1*Te,size(TxSymbMod,1));
                            ff = time2freq(tt);
                            tta = repmat(tt.',1,NumFra);
                            switch SelModTp
                                case 'BP'                                              %Base Band transmission
                                    
%                                     TxSymbMod   = TxSymbMod;%.*exp(1j*2*pi*CarrOffSet*FramSpac*tta);
                                    ff = time2freq(tt);
                                    if SelecGaus==1
                                        [ ModFilt ] = FiltroGaussiano(ff,OBw/2,0,1);
                                    else
                                        [ ModFilt ] = Filtro_Retangular(OBw,0,ff);
                                    end
                                    ModFilt = fftshift(ModFilt);
                                    if (kk==2)&&(Ploting)
                                        figure;hold all;
                                        plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                                    end
                                    ModFilta = repmat(ModFilt.',1,NumFra);
                                    TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
                                    if (kk==2)&&(Ploting)
                                        plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                                        plot(ff,20*log10(abs(fftshift(ModFilt))));
                                        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                    end
                                case 'AM'                                              %Out of base band
%                                     TxSymbModA  = rectpulse(TxSymbModA,OvSam);         %Over sampling
                                    TxSymbMod   = TxSymbMod.*cos(2*pi*Ofc*tta);
                                    if SelecGaus==1
                                        [ ModFilt ] = FiltroGaussiano(ff,OBw,Ofc,1);
                                    else
                                        [ ModFilt ] = Filtro_Retangular( OBw,Ofc,ff);
                                    end
                                    ModFilt = fftshift(ModFilt);
                                    if (kk==2)&&(Ploting)
                                        figure;hold all;
                                        plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                                    end
                                    ModFilta = repmat(ModFilt.',1,NumFra);
                                    TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
                                    if (kk==2)&&(Ploting)
                                        plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                                        plot(ff,20*log10(abs(fftshift(ModFilt))));
                                        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                    end
                                case 'AMSSB'                                           %Out of base band
%                                     TxSymbModA  = rectpulse(TxSymbModA,OvSam);         %Over sampling
                                    TxSymbMod   = TxSymbMod.*exp(1j*2*pi*Ofc*tta);
                                    if SelecGaus==1
                                        [ ModFilt ] = FiltroGaussiano(ff,OBw,Ofc,1);
                                    else
                                        [ ModFilt ] = Filtro_Retangular( OBw,Ofc,ff);
                                    end
                                    ModFilt = fftshift(ModFilt);
                                    if (kk==2)&&(Ploting)
                                        figure;hold all;
                                        plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                                    end
                                    ModFilta = repmat(ModFilt.',1,NumFra);
                                    TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
                                    if (kk==2)&&(Ploting)
                                        plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                                        plot(ff,20*log10(abs(fftshift(ModFilt))));
                                        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                    end
                                otherwise
%                                     TxSymbModA  = rectpulse(TxSymbModA,OvSam);         %Over sampling
                                    [TxSymbMod,tt] = modulate(TxSymbMod,Ofc,Ofs,'pm');
                                    if SelecGaus==1
                                        [ ModFilt ] = FiltroGaussiano(ff,3*Ofc,0,1);
                                    else
                                        [ ModFilt ] = Filtro_Retangular( 6*Ofc,0,ff);
                                    end
                                    ModFilt = fftshift(ModFilt);
                                    if (kk==2)&&(Ploting)
                                        figure;hold all;
                                        plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                                    end
                                    ModFilta = repmat(ModFilt.',1,NumFra);
                                    TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
                                    if (kk==2)&&(Ploting)
                                        plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                                        plot(ff,20*log10(abs(fftshift(ModFilt))));
                                        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                    end
                            end
                            
                            
                            %AuxTxSymbMod(1,:)  = TxSymbMod(:);
                            TxSymbModA  = [TxSymbMod(1+end-NPOFEX:end,:);TxSymbMod];
%                             SigTx  = rectpulse(TxSymbModA,NPPOF);                               %Over sampling
                            
                            TxSigA = TxSymbModA(:).';
%                             TxSig = TxSig + [TxSigA(1:NPSTUf) TxSigA];                         %Conforming vector sizes to the same length
                            U1t = [TxSigA(1:NPSTUf) TxSigA];                         %Conforming vector sizes to the same length
                            %                     CpPaFr = CpPaFr + 1;
                            %                 end
                            %Assigning the eletrical signal to one drive of the MZM - The aplitude of
                            %the signal can be controlled by the variable DatGai, which can be
                            %understood as an gain for the eletrical signal or an atenuation. The
                            %second signal will be similar with the only difference a phase shift of pi
                            NormFact = max((U1t));
                            U1t = 0.95*(U1t./NormFact);                     %Normalizing the sinal to mach with the MZM espect to receive.
                            U.U1t = U1t;
                            U.U2t = exp(-1j*pi).*U1t;
                            %As both signals will have the mostrly the same characteristics with the
                            %only difference the phase shift of 180 degress. The MZM-I will be working
                            %on the Push-Pull configuration. It is necessary to reduce the Chirp noise
                            %to zero.
                            [EoutMod,~]=MZM(freqGHz,EoutTx(VetThisCarrTx==(RefCarr-1+kk),:),U,L,U0,U_pi1,U_pi2,nopt,nel,C,alfa0);%Modulating individual carriers
                            EoutModTem(ObsCarrPos(kk),1:length(EoutMod)) = ...
                                EoutMod;%Adding the current Modulated output to the final OutPut
                            if (kk==2)&&(Ploting)
                                figure;hold all;grid on;
                                plot(f,20*log10(abs(fftshift(fft(EoutMod)./length(EoutMod)))));
                                axis([-25e9 50e9 -200 0]);
                                set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                            end
                        end
                        if IfftOrSum==1
                            if size(EoutModTem,1)>1
                                if UsingGpu==1
                                    EoutModTemGpu = gpuArray(EoutModTem);
                                    [EoutAux1Gpu,~,~]=OpticalIFFTG(fgpu,gpuArray(T),gpuArray(MaxNumStagT),EoutModTemGpu);
                                    EoutMod = gather(EoutAux1Gpu);
                                    clear EoutModTemGpu EoutAux1Gpu;
                                else
                                    [EoutMod,~,~] = OpticalIFFTN(f,T,MaxNumStagT,EoutModTem);
                                end
                            else
                                EoutMod = EoutModTem;
                            end
                        else
                            if size(EoutModTem,1)>1
                                EoutMod = sum(EoutModTem);
                            else
                                EoutMod = EoutModTem;
                            end
                        end
                    case 'DPSK'
                        %%            Generate the data DPSK
                        EoutModTem = zeros(NumCarr,NPPB*Nb);
                        % 						EoutModTem = 0;
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
                                U1t           = DpskEncodEq( TxData);
                                U1t(U1t==1) =  1;
                                U1t(U1t==0) = -1;
                                U1t           = rectpulse(U1t,NPPB);
                                
                                % Adding CP to the data
                                if AddCP==1
                                    TxAux = reshape(U1t,NPPB,NbDPSK);
                                    TxAux1 = [TxAux(1:NumAmosCP,:);TxAux;TxAux(end-(...
                                        NumAmosCP-1):end,:)];
                                    TxAux2 = reshape(TxAux1,1,(2*NumAmosCP+NPPB)*NbDPSK);
                                    U1t = [TxAux2 TxAux2(end-(StuffSampels-1):end)];
                                end
                                
                                %Because of reasons, the real world does not have components
                                %capable instantly change the voltage level. Thus, to emulate
                                %this behaviour the eletrical PAM signal pass through a
                                %gaussian filter that will remove the frequencies of higher
                                %order, which will result in small slope in the level variation
                                
                                [BitFilt,~] = FiltroGaussiano(f,BWD,CenFeq,FiltOrd);           %Creating filter for conformation of the input information
                                BitFilt     = fftshift(BitFilt);                               %Doing a shift on the vector for matching the transmited data
                                %                 TxSig       = ifft(fft(SigTx).*BitFilt);                       %Conforming the information and Creating the modulation signal
%                                 U1t       = U1t;                                   %Adding a gain to the eletrical signal
                                %             PrintInfo(Ploting*2,TxSig,T,NPPB);
                                %                 TxSigMat(kk,:) = TxSig;                                        %Storing the eletrical signal for further evaluation
                                %                 a=a+4;
                                
                                EoutMod = EoutAux1(VetThisCarr==ObsCarrUsed(kk),:);
                                %Assigning the eletrical signal to one drive of the MZM -
                                %The aplitude of the signal can be controlled by the
                                %variable DatGai, which can be understood as an gain for
                                %the eletrical signal or an atenuation. The second signal
                                %will be similar with the only difference a phase shift of
                                %pi.
                                U.U1t  = U1t;
                                U.U2t  = exp(-1j*pi).*U1t;
                                %As both signals will have the mostrly the same
                                %characteristics with the only difference the phase shift
                                %of 180 degress. The MZM-I will be working on the Push-Pull
                                %configuration. It is necessary to reduce the Chirp noise
                                %to zero.
                                [EoutMod,~]=MZM(freqGHz,EoutMod,U,L,U0,U_pi1,U_pi2,nopt,nel,C,alfa0);    %Modulating individual carriers
                                %                 EoutMod = EoutMod + EoutModAux;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                                EoutModTem(ObsCarrPos(kk),1:length(EoutMod)) = EoutMod;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                                %% Taking the sampling the EVM meassurement
                                PhaDel          = 0;
                                TimDel          = T;
                                [ESync1,ESync2] = DelayInterf(t,TimDel,PhaDel,EoutMod); %Coverting phase shift to amplitude variation
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
                                a=a+0;
                                TxSymbAmos(CurrentTest,ObsCarrPos==kk,:) = IxAux;%RxSymbAmosUp = [];
                                EvmMatRef(ObsCarrPos==kk,:) = IxAux;                   %Taking just the middle samples as references
                           else
                                EoutMod = EoutAux1(VetThisCarr==ObsCarrPos(kk),:);
%                                 EoutMod = EoutMod;
                                %                 EoutMod = EoutMod + EoutModAux;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                                EoutModTem(ObsCarrPos(kk),1:length(EoutMod)) = EoutMod;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                            end
                        end
                        if IfftOrSum==1
                            if UsingGpu==1
                                EoutModTemGpu = gpuArray(EoutModTem);
                                [EoutAux1Gpu,~,~]=OpticalIFFTG(fgpu,gpuArray(T),gpuArray(MaxNumStagT),EoutModTemGpu);
                                EoutMod = gather(EoutAux1Gpu);
                                clear EoutModTemGpu EoutAux1Gpu;
                            else
                                [EoutMod,~,~] = OpticalIFFTN(f,T,MaxNumStagT,EoutModTem);
                            end
                        else
                            EoutMod = sum(EoutModTem);
                        end
                    case 'DQPSK'
                        %%            Generate the data DQPSK
                        EoutModTem = zeros(NumCarr,NPPB*Nb);
                        TxDataMat = zeros(2*NumCarr,NbDQPSK/2);
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
                                if AddCP==1
                                    TxAux1 = reshape(SigI,NPPB,NbDQPSK/2);
                                    TxAux2 = [TxAux1(1:NumAmosCP,:);TxAux1;TxAux1(end-(...
                                        NumAmosCP-1):end,:)];
                                    TxAux3 = reshape(TxAux2,1,(2*NumAmosCP+NPPB)*NbDQPSK/2);
                                    SigI = [TxAux3 TxAux3(end-(StuffSampels-1):end)];
                                    
                                    TxAux4 = reshape(SigQ,NPPB,NbDQPSK/2);
                                    TxAux5 = [TxAux4(1:NumAmosCP,:);TxAux4;TxAux4(end-(...
                                        NumAmosCP-1):end,:)];
                                    TxAux6 = reshape(TxAux5,1,(2*NumAmosCP+NPPB)*NbDQPSK/2);
                                    SigQ = [TxAux6 TxAux6(end-(StuffSampels-1):end)];
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
                                EoutMod = EoutAux1(VetThisCarr==ObsCarrUsed(kk),:);
                                [EoutMod] = IqMod(EoutMod,TxSigI,TxSigQ,Vpi,V0);
                                %                 a=a+4;
                                EoutModTem(ObsCarrPos(kk),1:length(EoutMod)) = EoutMod;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                                %% Taking the sampling the EVM meassurement
                                PhaDel = 1*pi/4;
                                TimDel = T;
                                [EoutA,EoutB] = DelayInterf(t,TimDel,PhaDel,EoutMod);%Coverting phase shift to amplitude variation
                                PhaDel = -1*pi/4;
                                TimDel = T;
                                [EoutC,EoutD] = DelayInterf(t,TimDel,PhaDel,EoutMod);%Coverting phase shift to amplitude variation
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
                                    NPPB):(length(EoutMod)-StuffSampels);%Varriable respossible to take just the samples at the middle of the symbol
                                IxAux = EoutI(PosAuxEout) + 1j.*EoutQ(PosAuxEout);
                                TxSymbAmos(CurrentTest,ObsCarrPos==kk,:) = IxAux;%RxSymbAmosUp = [];
                                EvmMatRef(ObsCarrPos==kk,:) = IxAux;       %Taking just the middle samples as references
                            else
                                EoutMod = EoutAux1(VetThisCarr==ObsCarrPos(kk),:);
%                                 EoutMod = EoutMod;
                                %                 EoutMod = EoutMod + EoutModAux;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                                EoutModTem(ObsCarrPos(kk),1:length(EoutMod)) = EoutMod;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                            end
                        end
                        if IfftOrSum==1
                            if UsingGpu==1
                                EoutModTemGpu = gpuArray(EoutModTem);
                                [EoutAux1Gpu,~,~]=OpticalIFFTG(fgpu,gpuArray(T),gpuArray(MaxNumStagT),EoutModTemGpu);
                                EoutMod = gather(EoutAux1Gpu);
                                clear EoutModTemGpu EoutAux1Gpu;
                            else
                                [EoutMod,~,~] = OpticalIFFTN(f,T,MaxNumStagT,EoutModTem);
                            end
                        else
                            EoutMod = sum(EoutModTem);
                        end
                    case '4PAM'
                        %%        Generate the data 4PAM
                        EoutModTem = zeros(NumCarr,NPPB*Nb);
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
                                if AddCP==1
                                    TxAux1 = reshape(TxSig1,NPPB,Nb4Pam/2);
                                    if SetCpSampZer==1
                                        TxAux2 = [zeros(NumAmosCP,size(TxAux1,2));TxAux1...
                                            ;zeros(NumAmosCP,size(TxAux1,2))];
                                    else
                                        TxAux2 = [TxAux1(1:NumAmosCP,:);TxAux1;TxAux1(...
                                            end-(NumAmosCP-1):end,:)];
                                    end
                                    TxAux3 = reshape(TxAux2,1,(2*NumAmosCP+NPPB)*Nb4Pam/2);
                                    TxSig1 = [TxAux3 TxAux3(end-(StuffSampels-1):end)];
                                    TxAux4 = reshape(TxSig2,NPPB,Nb4Pam/2);
                                    if SetCpSampZer==1
                                        TxAux5 = [zeros(NumAmosCP,size(TxAux4,2));TxAux4...
                                            ;zeros(NumAmosCP,size(TxAux4,2))];
                                    else
                                        TxAux5 = [TxAux4(1:NumAmosCP,:);TxAux4;TxAux4(...
                                            end-(NumAmosCP-1):end,:)];
                                    end
                                    TxAux6 = reshape(TxAux5,1,(2*NumAmosCP+NPPB)*Nb4Pam/2);
                                    TxSig2 = [TxAux6 TxAux6(end-(StuffSampels-1):end)];
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
                                
                                EoutMod = EoutAux1(VetThisCarr==ObsCarrUsed(kk),:);
                                U.U1t  = TxSig1;                                               %Assigning the electrical signal to one drive of the MZM
                                U.U2t  = TxSig2;                                               %Assigning the electrical signal to another drive of the MZM
                                if Which4PAM
                                    if ModSchem
                                        %                         [EoutModAux] = IqMod4Pam (EoutT,U.U1t,U.U2t,U_pi2,Vbias);
                                        [EoutMod] = MZM(freqGHz,EoutMod,U,L,U0Up,U_pi1,U_pi2,nopt,nel,C,alfa0);
                                    else
                                        [EoutMod] = MZM(freqGHz,EoutMod,U,L,U0Up,U_pi1,U_pi2,nopt,nel,C,alfa0);
                                    end
                                else
                                    [EoutMod] = MZM(freqGHz,EoutMod,U,L,U0Up,U_pi1,U_pi2,nopt,nel,C,alfa0);
                                end
                                %                 EoutMod = EoutMod + EoutModAux;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                                EoutModTem(ObsCarrPos(kk),1:length(EoutMod)) = EoutMod;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                                %% Taking the sampling the EVM meassurement
                                PosAuxEout = (2*NumAmosCP+NPPB)/2:(2*NumAmosCP+NPPB)...
                                    :(length(EoutMod)-StuffSampels);%Varriable respossible to take just the samples at the middle of the symbol
                                Ix         = EoutMod.*conj(EoutMod);             %Recovering the signal that will be transmited
                                Ix         = ifft(fft(Ix).*EvmFilt);                   %Removing higher signal generated by the receiving process
                                Ix         = Ix - min(Ix);
                                Ix         = Ix./max(abs(Ix));                         %Normalizing the reference
                                IxAux      = Ix(PosAuxEout);
                                a=a+0;
                                TxSymbAmos(CurrentTest,ObsCarrPos==kk,:) = IxAux;%RxSymbAmosUp = [];
                                EvmMatRef(ObsCarrPos==kk,:) = IxAux;                   %Taking just the middle samples as references
                            else
                                EoutMod = EoutAux1(VetThisCarr==ObsCarrPos(kk),:);
%                                 EoutMod = EoutMod;
                                %                 EoutMod = EoutMod + EoutModAux;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                                EoutModTem(ObsCarrPos(kk),1:length(EoutMod)) = EoutMod;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                            end
                        end
                        if IfftOrSum==1
                            if UsingGpu==1
                                EoutModTemGpu = gpuArray(EoutModTem);
                                [EoutAux1Gpu,~,~]=OpticalIFFTG(fgpu,gpuArray(T),gpuArray(MaxNumStagT),EoutModTemGpu);
                                EoutMod = gather(EoutAux1Gpu);
                                clear EoutModTemGpu EoutAux1Gpu;
                            else
                                [EoutMod,~,~] = OpticalIFFTN(f,T,MaxNumStagT,EoutModTem);
                            end
                        else
                            EoutMod = sum(EoutModTem);
                        end
                    otherwise
                        %%        Generate the data OOK
                        EoutModTem = zeros(NumCarr,NPPB*Nb);
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
                                if AddCP==1
                                    TxAux = reshape(TxDataRes,NPPB,Nb-NumBitDesc);
                                    TxAux1 = [TxAux(1:NumAmosCP,:);TxAux;TxAux(end-(NumAmosCP...
                                        -1):end,:)];
                                    TxAux2 = reshape(TxAux1,1,(2*NumAmosCP+NPPB)*(Nb-...
                                        NumBitDesc));
                                    TxDataRes = [TxAux2 TxAux2(end-(StuffSampels-1):end)];
                                end
                                
                                
                                %Because of reasons, the real world does not have components
                                %capable instantly change the voltage level. Thus, to emulate
                                %this behaviour the eletrical PAM signal pass through a
                                %gaussian filter that will remove the frequencies of higher
                                %order, which will result in small slope in the level variation
                                
                                [BitFilt,~] = FiltroGaussiano(f,BWD,CenFeq,FiltOrd);           %Creating filter for conformation of the input information
                                BitFilt = fftshift(BitFilt);                                   %Doing a shift on the vector for matching the transmited data
                                %                 TxSig = ifft(fft(TxDataRes).*BitFilt);                         %Conforming the information and Creating the modulation signal
                                U1t = TxDataRes;                         %Conforming the information and Creating the modulation signal
                                %             TxSigMat(kk,:) = TxSig;                                        %Storring the transmited information for latter evaluation
                                
                                EoutMod = EoutAux1(VetThisCarr==ObsCarrUsed(kk),:);
                                %Assigning the eletrical signal to one drive of the MZM -
                                %The aplitude of the signal can be controlled by the
                                %variable DatGai, which can be understood as an gain for
                                %the eletrical signal or an atenuation. The second signal
                                %will be similar with the only difference a phase shift of
                                %pi.
                                U.U1t  = DatGai.*U1t;
                                U.U2t  = DatGai.*exp(-1j*pi).*U1t;
                                %As both signals will have the mostrly the same
                                %characteristics with the only difference the phase shift
                                %of 180 degress. The MZM-I will be working on the Push-Pull
                                %configuration. It is necessary to reduce the Chirp noise
                                %to zero.
                                [EoutMod,~]=MZM(freqGHz,EoutMod,U,L,U0,U_pi1,U_pi2,nopt,nel,C,alfa0);    %Modulating individual carriers
                                %                 EoutMod = EoutMod + EoutModAux;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                                EoutModTem(ObsCarrPos(kk),1:length(EoutMod)) = EoutMod;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                                %% Taking the sampling the EVM meassurement
                                PosAuxEout = (2*NumAmosCP+NPPB)/2:(2*NumAmosCP+NPPB)...
                                    :(length(EoutMod)-StuffSampels);%Varriable respossible to take just the samples at the middle of the symbol
                                Ix         = EoutMod.*conj(EoutMod);             %Recovering the signal that will be transmited
                                Ix         = ifft(fft(Ix).*EvmFilt);                   %Removing higher signal generated by the receiving process
                                Ix         = Ix - min(Ix);
                                Ix         = Ix./max(abs(Ix));                         %Normalizing the reference
                                IxAux      = Ix(PosAuxEout);                           %Normalizing the reference
                                a=a+0;
                                TxSymbAmos(CurrentTest,ObsCarrPos==kk,:) = IxAux;%RxSymbAmosUp = [];
                                EvmMatRef(ObsCarrPos==kk,:) = IxAux;                   %Taking just the middle samples as references
                            else
                                EoutMod = EoutAux1(VetThisCarr==ObsCarrPos(kk),:);
%                                 EoutMod = EoutMod;
                                %                 EoutMod = EoutMod + EoutModAux;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                                EoutModTem(ObsCarrPos(kk),1:length(EoutMod)) = EoutMod;% + EoutAux1(VetThisCarr==ObsCarrPos(kk-1),:);
                            end
                        end
                        if IfftOrSum==1
                            if UsingGpu==1
                                EoutModTemGpu = gpuArray(EoutModTem);
                                [EoutAux1Gpu,~,~]=OpticalIFFTG(fgpu,gpuArray(T),gpuArray(MaxNumStagT),EoutModTemGpu);
                                EoutMod = gather(EoutAux1Gpu);
                                clear EoutModTemGpu EoutAux1Gpu;
                            else
                                [EoutMod,~,~] = OpticalIFFTN(f,T,MaxNumStagT,EoutModTem);
                            end
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