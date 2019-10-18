%c
%c                                                       ..'Ž`'..'Ž`..'Ž`..
%c       File: DatumReception File For the All Optical OFDM
%c(Main file to call each idividual parameters)
%c
%c     This main code is resposible to call and run the all the fuction to
%c related to this simulation. Here it is possible to change any
%c configuration that was previouly stated on the Input data file of this
%c simulation.
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
%%                Main of All optical FFT
% The principal importance of the script is to load the All Optical FFT.
% Also, it is responsible to receive the transmited data.

DatumReceptionInputData;                                                   % Configuration File of DatumReception.

%%            Recovering the signal
% At this point the FFT will be inplemented. The received signal need to be
% given as an parameter. This following function recursively implement the
% FFT Operation. One may need to investigate how would work a receptor with
% filters insted of a OFFT. The second implementation is an array of
% filters.
if FFTSplit
    [EoutAux1,~,VetThisCarr]=OpticalFFTN(t,T,MaxNumStag,EoutRec);
%     [EoutAux1,~,VetThisCarr]=OpticalFFTN(t,T,1,EoutRec);
%     for kk=1:size(EoutAux1,1)
%         [RecFilt,~]   = FiltroGaussiano(f,2*fc,(kk-1)*fc+fin,5);%Creating filter for conformation of the input information
%         RecFilt       = fftshift(RecFilt);
%         EoutAux1(kk,:) = ifft(fft(EoutAux1(kk,:)).*RecFilt);
%     end
else
    VetThisCarr = ObsCarrPos;
    EoutAux1    = EoutRec;
end
for kk=1:size(EoutAux1,1)
    [~,CarrRecPowDo(CurrentTest,kk)] = MeasPower(EoutAux1(VetThisCarr==...
                                                     ObsCarrUsed(kk),:),t);
end
%%  Ploting some results for qualitative analizes
% PrintInfo(PlotingThis*12,f,db(abs((fftshift(fft(EoutRec))))/length(EoutRec)),EoutAux1,RefCarr*fc);
%%      Recovering the Data
%This is basicaly the final step, so far, in this process. Here, the
%transmited signal will be received and processed to recover the actual
%data within it.
switch Modulation
    case 'OFDM'
        SigRecep = zeros(NumCarr,NbOFDM);
        CountOfdm = 1;
        for ThisCarr=InitCarr:CarrPass:NumCarr
            EoutAux = EoutAux1(VetThisCarr==ObsCarrUsed(ThisCarr),:);
            switch Medium
                case 'Fiber'
                    EoutAux = ifft(fft(EoutAux).*exp(1j*2*pi*f*(...
                        FiberDelay(ThisCarr)*Ta)));
                otherwise
            end
            %The current incoming signal is them converted from the optical
            %domain to the eletrical domain with the help of an photo
            %detector.
            Ix =EoutAux.*conj(EoutAux);
            [BitFilt,~] = FiltroGaussiano(f,BWD,CenFeq,FiltOrd);           %Creating filter for selection of the received signal
            BitFilt = fftshift(BitFilt);                                   %Shifting the filter for matching the received signal

            Ix = ifft(fft(Ix).*BitFilt);                                   %Filtering the signal for better detection

            if ReceptorNoise
                [~,SigPower] = MeasPower(Ix);
                SigPower2 = SigPower-30;%10*log10(SigPower);
                Ix = awgn(Ix,CarSNR,SigPower2);
            end
            Ix = Ix - min(Ix);                                             %Removing any off-set that may exist
            Ix = Ix./max(abs(Ix));                                         %Normalizing the eletrical signal (amplifying what is needed)
            
            %%        synchronizing
            %For the reception process work properly, it is needed to
            %sincronized the recieved signal or the sampling process will
            %not work.

            AuxSync = (Ix(IniSyncPos:SyncPos));                            %Selecting the sync-word within the received signal
            AuxSync1 = AuxSync;
            %This synchronization process is based on the mean expected
            %value. That means, the information mean value should be within
            %one period of symbol. Thus, the mean value of the received
            %signal is acquired and compare of the known sync-word to
            %verify if this mean value is at the right possition
            AuxSync1(AuxSync1<0) = 0;                                      %To keep the mean value above zero anything under is neglected
            AuxSync1(AuxSync1>=mean(AuxSync1)) = 1;                        %Adding a flag to the first sample of the received mean value
            AuxSync1(AuxSync1<mean(AuxSync1)) = -1;                        %All the others samples at set to the lowest level
            AuxSync2 = SyncSymb(IniSyncPos:SyncPos);                       %Selecting the sync-word within the known signal
            AuxSync2(AuxSync2>=mean(AuxSync2)) = 1;                        %Adding a flag to the first sample of the known mean value
            AuxSync2(AuxSync2<mean(AuxSync2)) = -1;                        %All the others samples at set to the lowest level

            PosToSyn  = find(ismember(AuxSync1,1));
            PosSyn    = find(ismember(AuxSync2,1));

            AuxSyncCorr = PosToSyn(round(end/2))-PosSyn(round(end/2));
            %The difference between the PossitionTosynchronize and
            %Possitionsynchronized will be used to correct the time
            %shifting on the transmition and reception process.

            if AuxSyncCorr>=0%If the difference is positive, left-shift...
                %                 Ix = ifft(fft(Ix).*exp(1j*2*pi*f*(AuxSyncCorr*Ta)));     %Shift based on time change
                Ix = [Ix(AuxSyncCorr+1:end) Ix(1:AuxSyncCorr)];            %Shift based on sampling sliding
            else%... but if the difference is negative, right-shift
                %                 Ix = ifft(fft(Ix).*exp(-1j*2*pi*f*(AuxSyncCorr*Ta)));    %Shift based on time change
                Ix = [Ix(end+AuxSyncCorr+1:end) Ix(1:end + AuxSyncCorr)];  %Shift based on sampling sliding
            end

            %% Ploting the result for qualitative analizes

%             PrintInfo(PlotingThis*45,t(IniSyncPos:SyncPos),Ix(IniSyncPos:SyncPos),SyncSymb(IniSyncPos:SyncPos));
%             PrintInfo(PlotingThis*46,t(1:length(Ix)),Ix);
            %%                  Removing CP
            if AddCP
                IxAux = Ix(1:end - StuffSampels);
                IxAux = reshape(IxAux,(2*NumAmosCP+NPPB),Nb-NumBitDesc);
                IxAux = IxAux(1+NumAmosCP:end-NumAmosCP,:);
                IxAux = reshape(IxAux,1,NPPB*(Nb-NumBitDesc));
                Ix    = IxAux;
            end
            contaux1 = 1;
            for kk=1:NPPB:length(Ix)                                       %The comparison process will be made for each symbol period
                IxDownSamp(contaux1) = Ix(kk+NPPB/2);
                contaux1 = contaux1 + 1;
            end
            SigRecep(CountOfdm,:) = IxDownSamp(JustSize/2+1:JustSize/2+NbOFDM);
            CountOfdm             = CountOfdm + 1;
            
        end
        SigRecep = fft(SigRecep,EleNfft);
        SigRecep = qamdemod(SigRecep,M);
%         figure;
%         hold on;
%         for kk=1:size(TxDataMat,1)
%             plot(TxDataMat(kk,:));
%             plot(SigRecep(kk+1,:));
%             set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%             a=0;
%         end
%         a=0;
    case 'DPSK'
        %%                   Receiver DPSK
        for ThisCarr=InitCarr:CarrPass:NumCarr                                    %For each carrier the same process of reception will be used.
            if CarrUsedDo(ThisCarr)
                %%  Reception Process: Handling the income optical field
                %At the first step the income field will be selected from the
                %optical FFT Output with the help of the maping vector
                %(previously described).
                EoutAux = EoutAux1(VetThisCarr==ObsCarrUsed(ThisCarr),:);
                if AddingNoiseP
                    [~,sigpower] = MeasPower(EoutAux);
                    EoutAux      = awgn(EoutAux,osnrp,sigpower-30);
                    [~,sigpower] = MeasPower(EoutAux);
                end
                %%            Fiber Time Delay Compensation
                switch Medium
                    case 'Fiber'
                        EoutAux = ifft(fft(EoutAux).*exp(1j*2*pi*f*(...
                            FiberDelay(ThisCarr)*Ta)));
                    otherwise
                end
                %             EoutAux = EoutRec;
                %%  Plot for Qualitative analizes
                %PrintInfo(Ploting*13,t,abs(EoutAux));
                %             PrintInfo(Ploting*14,f,20*log10(abs(fftshift(fft(EoutAux)./...
                %                                                        length(EoutAux)))));

                %%                   synchronizing
                %For the reception process work properly, it is needed to
                %sincronized the recieved signal or the sampling process will
                %not work.

                %At this first moment a small part of the income signal needs
                %be analized for the sincronism process. As we are looking for
                %the actual synchronis symbol, which is the bigest phase shift
                %of the income signal. The next delay interferometer convert
                %this phase shift to a amplitude variation.

                %The phase delay is not needed because the optical field will
                %be analized as a whole as it is composed only by the real part
                PhaDel          = 0;
                %Remember that at DPSK modulation the information is stored at
                %the phase difference of the current symbol with the previous
                %one hence this time delay is needed to analyse each symbol by
                %analizes of its interaction of the current symbel with the
                %previous one.
                TimDel          = T;
                [ESync1,ESync2] = DelayInterf(t,TimDel,PhaDel,EoutAux);        %Coverting phase shift to amplitude variation

                %The second moment is to transfer this information from the
                %optical domain to the eletrical domain for an eletronic
                %processing. It is done with the help of an photo diode. The
                %configuration here used is an balanced receiver as the output
                %of the Delay Interferometer has two signals resulting from
                %the constructive and destructive signal interaction.
                ESync1 = ESync1.*conj(ESync1);
                %             ESync1 = ifft(fft(ESync1)./fft(PulseResp(VetThisCarr==ThisCarr,:)));
                ESync2 = ESync2.*conj(ESync2);
                %             ESync2 = ifft(fft(ESync2)./fft(PulseResp(VetThisCarr==ThisCarr,:)));

                %%           Creating the Reception Filter
                
                
                [BitFilt,~] = FiltroGaussiano(f,BWD2,CenFeq2,FiltOrd2);        %Creating filter for selection of the received signal
                BitFilt = fftshift(BitFilt);                                   %Shifting the filter for matching the received signal

                %%
                Esync  = ESync2-ESync1;
                
                if RecFilBanPas
                    Esync  = ifft(fft(Esync).*BitFilt);                             %Filter is used to remove higher order components
                    Esync  = Esync-mean(Esync);
                    Esync  = Esync./max(abs(Esync));
                end
%                 ESync1 = ESync1./max(abs(ESync1));                             %Normalizing the signal
%                 ESync2 = ESync2./max(abs(ESync2));                             %Normalizing the signal
%                 Esync  = Esync./max(abs(Esync));                               %Normalizing the signal

                %%         Adding Noise to Received Signal
                %Just after the photo diod is where the noise must be added as
                %it is the main constrain of the system. Once the signal is
                %converted from optical to electrical and it digital values
                %were recevered, there is no point to add noise, because the
                %information was already received, it just needs decoding.
                %The noise here added will be gaussian basically it represents
                %the reception sensibility. In another words, how much the
                %signal must be about the noise for an acceptable recovering.
%                 if ReceptorNoise                                               %Verify whether the noise is to be added or not
%                     [BitFiltN,~] = FiltroGaussiano(f,2*BWD2,CenFeq2,FiltOrd2);        %Creating filter for selection of the received signal
%                     BitFiltN = fftshift(BitFiltN);                                   %Shifting the filter for matching the received signal
%                     [~,SigPower] = MeasPower(Esync);                           %it own received signal or
%                     SigPower2    = SigPower-30;%10*log10(SigPower);
%                     SigAux       = linspace(SigPower2,SigPower2,length(t));
%                     IrecNoise = awgn(SigAux,CarSNR,'measured');
% %                     IrecNoise = ifft(fft(IrecNoise).*BitFiltN);                             %Filter is used to remove higher order components
%                     Esync     = Esync + IrecNoise;
%                 end
                if ReceptorNoise                                               %Verify whether the noise is to be added or not
                    [BitFiltN,~] = FiltroGaussiano(f,2*BWD2,CenFeq2,FiltOrd2);        %Creating filter for selection of the received signal
                    BitFiltN = fftshift(BitFiltN);                                   %Shifting the filter for matching the received signal
                    SNR2 = CarSNR + 10*log10(1) - 10*log10(Nsamp) + 10*log10(10^0.36);
                    SNR = CarSNR + 10*log10(1) - 10*log10((t(NPPB))/(2*(t(2*NumAmosCP+NPPB))/(NPPB)));
                    Esync = awgn(Esync,SNR,'measured');
%                     Esync = ifft(fft(Esync).*BitFiltN);                             %Filter is used to remove higher order components
                end
                if ~RecFilBanPas
                    Esync  = ifft(fft(Esync).*BitFilt);                             %Filter is used to remove higher order components
                    Esync  = Esync-mean(Esync);
                    Esync  = Esync./max(abs(Esync));
                end
                
                VetElecPower(CurrentTest,ThisCarr)= MeasPower(Esync);
                EoutAuxF = 20*log10(abs(fftshift(fft(EoutAux)./length(...
                    EoutAux))));
                [VetOptiPower(CurrentTest,ThisCarr),~]= findpeaks(EoutAuxF,...
                    'SortStr','descend','NPeaks',1);
                
                SyncAux   = Esync(IniSyncPos:SyncPos);                         %Selecting just the symbol to synchronize
                SyncedAux = SyncSymb(IniSyncPos:SyncPos);                      %Selecting just the symbol to synchronize

                %%                   Plot for Qualitative analizes
                %             PrintInfo(Ploting*15,t(IniSyncPos:SyncPos),Esync(IniSyncPos:...
                %             SyncPos),SyncSymb(IniSyncPos:SyncPos),ESync1(IniSyncPos:...
                %                                       SyncPos),ESync2(IniSyncPos:SyncPos));

                %%                   Synchronizing
                %This synchronization process is based on the mean expected
                %value. That means, the information mean value should be within
                %one period of symbol. Thus, the mean value of the received
                %signal is acquired and compare of the known sync-word to
                %verify if this mean value is at the right possition. Which is
                %the midel point (peak) of the highest level at the sync period
                SyncAux(SyncAux<0)              = 0;                           %To keep the mean value above zero anything under is neglected
                SyncAux(SyncAux>=mean(SyncAux)) = 1;                           %Adding a flag to the first sample of the received mean value
                SyncAux(SyncAux<mean(SyncAux))  = -1;                          %All the others samples at set to the lowest level

                PosToSyn  = find(ismember(SyncAux,1));                         %Finding where is the location of the samples to synchronize
                PosSynced = find(ismember(SyncedAux,1));                       %Finding where is the location of the samples to synchronize

                DiffPos = PosToSyn(round(end/2)) - PosSynced(round(end/2));    %Accounting the peak (midel point) displacement
                if DiffPos~=0
                    if DiffPos>0%If the difference is positive, left-shift...
                        %                     Esync = ifft(fft(Esync).*exp(1j*2*pi*f*(DiffPos*Ta))); %Shift based on time change
                        Esync = [Esync(DiffPos+1:end) Esync(1:DiffPos)];       %Shift based on sampling sliding
                    else%... but if the difference is negative, right-shift
                        %                     Esync = ifft(fft(Esync).*exp(1j*2*pi*f*(DiffPos*Ta))); %Shift based on time change
                        Esync = [Esync(end+DiffPos+1:end) Esync(1:end + ...
                            DiffPos)];%Shift based on sampling sliding
                    end
                end

                %Because of reasons, sometimes it may be required to make a
                %synchronization process with the  end of the data stream as
                %well. This following verification check if the user set (or
                %not) a second synchronization process to be done.

                %%          Ploting the result for qualitative analizes
                %             PrintInfo(Ploting*16,t(end-SyncPos+1:end-IniSyncPos+1),Esync...
                %             (end-SyncPos+1:end-IniSyncPos+1),SyncSymb(end-SyncPos+1:end-...
                %             IniSyncPos+1),ESync1(end-SyncPos+1:end-IniSyncPos+1),ESync2(...
                %                                           end-SyncPos+1:end-IniSyncPos+1));
                if SencondAdjust
                    %%                   Synchronizing
                    %This synchronization process is based on the mean expected
                    %value. That means, the information mean value should be
                    %within one period of symbol. Thus, the mean value of the
                    %received signal is acquired and compare of the known
                    %sync-word to verify if this mean value is at the right
                    %possition.

                    SyncAuxEnd   = Esync(end-SyncPos+1:end-IniSyncPos+1);
                    SyncedAuxEnd = SyncSymbEnd(end-SyncPos+1:end-IniSyncPos+1);
                    SyncAuxEnd(SyncAuxEnd<0)                 = 0;              %To keep the mean value above zero anything under is neglected
                    SyncAuxEnd(SyncAuxEnd>=mean(SyncAuxEnd)) = 1;              %Adding a flag to the first sample of the received mean value
                    SyncAuxEnd(SyncAuxEnd<mean(SyncAuxEnd))  = -1;             %All the others samples at set to the lowest level

                    PosToSynEnd  = find(ismember(SyncAuxEnd,1));               %Finding where is the location of the first sample to synchronize
                    PosSyncedEnd = find(ismember(SyncedAuxEnd,1));             %Finding where is the location of the first sample to synchronize

                    DiffPosEnd = PosToSynEnd(round(end/2))-PosSyncedEnd(...
                        round(end/2));
                    if DiffPosEnd~=0
                        if DiffPosEnd>0%If positive difference, left-shift...
                            %                         Esync = ifft(fft(Esync).*exp(1j*2*pi*f*(...
                            %                                                           DiffPosEnd*Ta)));%Shift based on time change
                            Esync = [Esync(DiffPosEnd+1:end) Esync(1:...
                                DiffPosEnd)];%Shift based on sampling sliding
                        else%... but if the difference is negative, right-shift
                            %                         Esync = ifft(fft(Esync).*exp(1j*2*pi*f*(...
                            %                                                           DiffPosEnd*Ta)));%Shift based on time change
                            Esync = [Esync(end-DiffPosEnd+1:end) Esync(...
                                1:end-DiffPosEnd)];%Shift based on sampling sliding
                        end
                    end
                end

                %% Removing CP
                if AddCP
                    IxAux  = Esync(1:end - StuffSampels);
                    IxAux  = reshape(IxAux,(2*NumAmosCP+NPPB),NbDPSK);
                    IxAux  = IxAux(1+NumAmosCP:end-NumAmosCP,:);
                    IxAux  = reshape(IxAux,1,NPPB*NbDPSK);
                    Esync  = IxAux;
                end
                
                %% Taking the sampling the EVM meassurement
                clear IxAux;
                PosAuxEout1 = NPPB/2:NPPB:length(Esync);                   %Varriable respossible to take just the samples at the middle of the symbol
                PosAuxEout2 = ((NPPB/2)+(NPPB/16)):NPPB:length(Esync);      %Varriable respossible to take just the samples at the middle of the symbol
                PosAuxEout3 = ((NPPB/2)-(NPPB/16)):NPPB:length(Esync);      %Varriable respossible to take just the samples at the middle of the symbol
                IxAux1      = Esync(PosAuxEout1);                          %Normalizing the reference
                IxAux2      = Esync(PosAuxEout2);    %Normalizing the reference
                IxAux3      = Esync(PosAuxEout3);    %Normalizing the reference
                a=0;
                EvmMatRecA(1,:) = IxAux1;                       %Taking just the middle samples as references
                EvmMatRecA(2,:) = IxAux2;                       %Taking just the middle samples as references
                EvmMatRecA(3,:) = IxAux3;                       %Taking just the middle samples as references
                EvmMatRecA(4,:) = IxAux1;                       %Taking just the middle samples as references
                [EvmDBA(1), EvmPerA(1), EvmRmsA(1) ] = EvmCalc( EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux1 );
                [EvmDBA(2), EvmPerA(2), EvmRmsA(2) ] = EvmCalc( EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux2 );
                [EvmDBA(3), EvmPerA(3), EvmRmsA(3) ] = EvmCalc( EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux3 );
                [EvmDBA(4), EvmPerA(4), EvmRmsA(4) ] = EvmCalc( EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux1 );
                [EvmDBJA(1),EvmPerJA(1),EvmRmsJA(1)] = evm1(2,'dpsk',EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux1);
                [EvmDBJA(2),EvmPerJA(2),EvmRmsJA(2)] = evm1(2,'dpsk',EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux2);
                [EvmDBJA(3),EvmPerJA(3),EvmRmsJA(3)] = evm1(2,'dpsk',EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux3);
                [EvmDBJA(4),EvmPerJA(4),EvmRmsJA(4)] = evm1(2,'dpsk',EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux1);
%                 if ThisCarr==1
%                     figure;hold all;grid on;
%                     plot(real(IxAux1),zeros(1,length(IxAux1)),'rx');
%                     plot(real(IxAux2),zeros(1,length(IxAux2)),'yx');
%                     plot(real(IxAux3),zeros(1,length(IxAux3)),'gx');
%                     set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%                     a=0;
%                 end
%##########################################################################
%                 Esync = Esync(1+SyncPeriod:end-SyncPeriod);
                PosIx = NPPB/2:NPPB:length(Esync);                                %Possition of the central samples - but the number of samples per symbol is even ????
                IxAux = Esync(PosIx);                                             %From the main received signal just a few samples are taken for further evaluation
                n     = 100;                                                   %The number of boxes to be filled up on the histogram process
                Perce = 0.7;
                %At this process each eye will be evaluated separetelly. It is
                %important to remember that the PAM4 has 4 levels, which means
                %three level of decissions that are exaclty the center of the
                %eye diagram.
                IxAuxAB = Esync(PosIx);       %Taking just those values relative to the uper eye
                InterAB = linspace(Perce*(min(Esync)),Perce*(max(Esync)),n);                     %Building the histogram boxes
                EyeAB = hist(IxAuxAB,InterAB);                                 %filling up the boxes with samples that fit on them.
                
                %The same process described for the uper level will be done at
                %the middle and lower eyes levels.
                IxAuxCD = Esync(PosIx-NPPB/16);
                InterCD = linspace(Perce*(min(Esync)),Perce*(max(Esync)),n);%NPPB*2^n);
                EyeCD = hist(IxAuxCD,InterCD);
                
                IxAuxEF = Esync(PosIx+NPPB/16);
                InterEF = linspace(Perce*(min(Esync)),Perce*(max(Esync)),n);%NPPB*2^n);
                EyeEF = hist(IxAuxEF,InterEF);
                
                %What we are looking for here are not where exist the
                %occurrences of values. The eye is where there is not samples,
                %this means, the decission level is a reagion between actual
                %levels thus this reagion does not contain any signal.
                EyeAB = ~EyeAB;                                                %Changing zeros to one - Zeros compose the eye region
                EyeCD = ~EyeCD;
                EyeEF = ~EyeEF;
                %Starting variables that will be used to identify the eye
                %diagram.
                CountAB   = 1;
                CountCD   = 1;
                CountEF   = 1;
                
                SeqOnesAB = 0;
                SeqOnesCD = 0;
                SeqOnesEF = 0;
                
                SeqFinAB  = 0;
                SeqFinCD  = 0;
                SeqFinEF  = 0;
                
                SeqIniAB  = 1;
                SeqIniCD  = 1;
                SeqIniEF  = 1;
                %The for loop will take account of every box with ones. It is
                %important to take note that the not operator was used in this
                %vector, therefore ones means zeros (the eye diagram -
                %possibly) and zeros means values abouve zeroa (not the eye).
                for kk=1:length(EyeAB)                                         %For every box
                    if EyeAB(kk)                                               %if it contains "1"
                        SeqOnesAB(SeqIniAB)=CountAB;                           %count this element as part of a consecutive sequency
                        CountAB = CountAB + 1;                                 %adds one to the counter of consecutive elements "1"
                        if kk==length(EyeAB)                                   %if the current box is the last box we got to an end
                            SeqFinAB(SeqIniAB) = kk;                           %The final sequency element is equal to its possition (kk)
                        end
                    else                                                       %else if the current box contains "0"
                        SeqFinAB(SeqIniAB) = kk-1;                             %the previous element was the last element of a consecutive sequency of ones
                        SeqIniAB = SeqIniAB + 1;                               %adds one to the counter of consecutive number (take in account how many sequencies there are)
                        CountAB = 1;                                           %reset the counter of consecutive elements to mark the start of a new consecutive sequence
                    end
                    
                    if EyeCD(kk)
                        SeqOnesCD(SeqIniCD)=CountCD;
                        CountCD = CountCD + 1;
                        if kk==length(EyeCD)
                            SeqFinCD(SeqIniCD) = kk;
                        end
                    else
                        SeqFinCD(SeqIniCD) = kk-1;
                        SeqIniCD = SeqIniCD + 1;
                        CountCD = 1;
                    end
                    
                    if EyeEF(kk)
                        SeqOnesEF(SeqIniEF)=CountEF;
                        CountEF = CountEF + 1;
                        if kk==length(EyeEF)
                            SeqFinEF(SeqIniEF) = kk;
                        end
                    else
                        SeqFinEF(SeqIniEF) = kk-1;
                        SeqIniEF = SeqIniEF + 1;
                        CountEF = 1;
                    end
                end
                
                
                %If the eye is open, which means there is a clear difference
                %between adjacent levels, the eye diagram will be the longest
                %sequence of ones.
                [MaxValAB,LocMaxAB]=max(SeqOnesAB);
                if LocMaxAB<2 || MaxValAB<2                                    %if any sequency was found or there is just one sequency it is a error thus
                    LevDec3 = 0.0;                                            %the decission level will be by default 0.67. Also, the other variables will
                    LocMaxAB = 1;                                              %will be set with values to not cause errors in the future
                    SeqFinAB(1)=2;
                    MaxValAB = 0;
                    InterAB(1)=LevDec3;
                else                                                           %if a sequency was found the middle of the eye will be the middle of the sequency
                    if (SeqFinAB(LocMaxAB)-MaxValAB/2)<1                       %if for some reason the final element of the sequency minus the half of its
                        LevDec3 = 0.0;                                         %length results in a negative value, something went very wrong, and by
                    else                                                       %default it will be set to 0.7
                        LevDec3 = InterAB(round(SeqFinAB(LocMaxAB)-MaxValAB/...
                            2));%Otherwise, the decission level is the middle of the sequency
                    end
                end
                
                [MaxValCD,LocMaxCD]=max(SeqOnesCD);
                if LocMaxCD<2 || MaxValCD<2
                    LevDec2 = 0.0;
                    LocMaxCD = 1;
                    SeqFinCD(1)=2;
                    MaxValCD = 0;
                    InterCD(1)=LevDec2;
                else
                    if (SeqFinCD(LocMaxCD)-MaxValCD/2)<1
                        LevDec2 = 0.0;
                    else
                        LevDec2 = InterCD(round(SeqFinCD(LocMaxCD)-MaxValCD/2));
                    end
                end
                
                [MaxValEF,LocMaxEF]=max(SeqOnesEF);
                if LocMaxEF<2 || MaxValEF<2
                    LevDec1 = 0.0;
                    LocMaxEF = 1;
                    SeqFinEF(1)=2;
                    MaxValEF = 0;
                    InterEF(1)=LevDec1;
                else
                    if (SeqFinEF(LocMaxEF)-MaxValEF/2)<1
                        LevDec1 = 0.0;
                    else
                        LevDec1 = InterEF(round(SeqFinEF(LocMaxEF)-MaxValEF/2));
                    end
                    
                end
                
                %another way to measure the eye opening is the get all the
                %boxes and find all peaks on it, that will be a plato created
                %by the sequences of ones (which was zeros). From thos peaks,
                %the eye diagram will be the longer of them hence it will take
                %the most part of the vector that store the findpeaks result.
                %Thus, the middle of the eye will be basically the middle of
                %the peaks vector.
                LocAB = find(EyeAB);
                if isempty(LocAB)                                              %if for some reason there are no peaks, something went wrong.
                    LocAB = 1;
                    LevelDec3 = 0.0;%mean(Levels(3:4));                       %by default the decission level will be set to 0.65
                else
                    LevelDec3 = LevDec3;
                end
                
                LocCD = find(EyeCD);
                if isempty(LocCD)
                    LocCD = 1;
                    LevelDec2 = 0.0;%mean(Levels(2:3));
                else
                    LevelDec2 = LevDec2;
                end
                
                LocEF = find(EyeEF);
                if isempty(LocEF)
                    LocEF = 1;
                    LevelDec1 = 0.0;%mean(Levels(1:2));
                else
                    LevelDec1 = LevDec1;
                end
%##########################################################################
                if CalcS
                    [~,~,EyeOpenI,EyeOpenHighI,EyeOpenLowI] = Olho(Esync,T,NPPB,0);
                    AberLevS(CurrentTest,ThisCarr)  = EyeOpenI;
                    ValsLevS(CurrentTest,ThisCarr)  = EyeOpenLowI + EyeOpenI/2;          %Comparison between the Transmited and received and counting the differences
                end
%                 [~,~,EyeOpenQ,EyeOpenHighQ,EyeOpenLowQ] = Olho(EoutQ,T,NPPB,0);
                AberLevAuxI(1)= InterAB(SeqFinAB(LocMaxAB))-InterAB(round(SeqFinAB(LocMaxAB)-MaxValAB));
                AberLevAuxI(3)= InterCD(SeqFinCD(LocMaxCD))-InterCD(round(SeqFinCD(LocMaxCD)-MaxValCD));
                AberLevAuxI(2)= InterEF(SeqFinEF(LocMaxEF))-InterEF(round(SeqFinEF(LocMaxEF)-MaxValEF));
                ValsLevAuxI(1)= LevDec3;
                ValsLevAuxI(3)= LevDec2;
                ValsLevAuxI(2)= LevDec1;

                %% Ploting the result for qualitative analizes
                if PrintinEye
                    PrintInfo(PlotingThis*17,Esync,T,NPPB);
                    hold on;
                    plot(t((NPPB)/2),InterAB(round(SeqFinAB(LocMaxAB)-MaxValAB)),'kx');
                    plot(t((NPPB)/2),InterAB(SeqFinAB(LocMaxAB)),'ko');
                    plot(t((NPPB)/2),InterAB(round(SeqFinAB(LocMaxAB)-MaxValAB/2)),'kd');
                    plot(t(((NPPB)/2)-(NPPB)/16),InterCD(round(SeqFinCD(LocMaxCD)-MaxValCD)),'gx');
                    plot(t(((NPPB)/2)-(NPPB)/16),InterCD(SeqFinCD(LocMaxCD)),'go');
                    plot(t(((NPPB)/2)-(NPPB)/16),InterCD(round(SeqFinCD(LocMaxCD)-MaxValCD/2)),'gd');
                    plot(t(((NPPB)/2)+(NPPB)/16),InterEF(round(SeqFinEF(LocMaxEF)-MaxValEF)),'bx');
                    plot(t(((NPPB)/2)+(NPPB)/16),InterEF(SeqFinEF(LocMaxEF)),'bo');
                    plot(t(((NPPB)/2)+(NPPB)/16),InterEF(round(SeqFinEF(LocMaxEF)-MaxValEF/2)),'bd');
                    if CalcS
                        plot(t((NPPB)/2),EyeOpenLowI,'mx');
                        plot(t((NPPB)/2),EyeOpenHighI,'mo');
                        plot(t((NPPB)/2),EyeOpenLowI + EyeOpenI/2,'md');
                    end
                    a=0;
                end
%                 [~,~,EyeOpen,EyeOpenHigh,EyeOpenLow] = Olho(Esync,T,NPPB,0);

%                 AberLev(CurrentTest,ThisCarr)= Inter(SeqFin(LocMax))-Inter(round(SeqFin(LocMax)-MaxVal));
%                 ValsLev(CurrentTest,ThisCarr)= LevDec;
%                 AberLev(CurrentTest,ThisCarr)  = EyeOpen;
%                 ValsLev(CurrentTest,ThisCarr)  = EyeOpenLow + EyeOpen/2;

%                 Txaux1 = rectpulse(TxDataMat(ThisCarr,:),NPPB);
%                 Txaux1(Txaux1==0) = -1;
                %% Ploting the result for qualitative analizes
%                 hold on;
%                 plot(t((NPPB)/2),Inter(round(SeqFin(LocMax)-MaxVal)),'bx');
%                 plot(t((NPPB)/2),Inter(SeqFin(LocMax)),'bo');
%                 plot(t((NPPB)/2),Inter(round(SeqFin(LocMax)-MaxVal/2)),'bd');
%                 plot(t((NPPB)/2),EyeOpenLow,'cx');
%                 plot(t((NPPB)/2),EyeOpenHigh,'co');
% %                 plot(t((NPPB)/2),EyeOpenLow + EyeOpen/2,'cd');
%                 a=0;
%                 PrintInfo(PlotingThis*18,t(1:length(Esync)),Txaux1,Esync);
                %%               Recovering the Information
                %Once the Income signal was synchronized it is possible to
                %recover the signal.
                %
                %After passing the optical signal to the eletrical domain, for
                %actually detect the data withing the signal the following
                %steps are needed.
                %
                %Finding Decission Levels:
                %The process for decoding the income signal will be based on
                %eletronic comparators. Inasmuch as the right decission level
                %must be acquired for accurately decide, within a symbol
                %periode, what that current leavel means (ones or zeros).
                %
                %The process hereafter of chosing the  decission levels is not
                %deterministic rather it is a statistic process. The main idea
                %is to take the decission level from the histogram generated
                %from the income signal stream.
                %
                %This process is realized inside the function Olho.
                %
                %Basicaly the decission level will be the minimal value of the
                %currente eye under evaluation plus the half of the its eye
                %opening.The following ilustration better describe this process
                %
                %Eye Limits:   + (1/2 * Eye Opening:)  =    Comparisson Limit:
                %
                %UperLevel  ______     ______________     _____
                %                 \   /      |       \   /
                %                  \ /       |        \ /
                %                   \ Half Eye Opening /     Decission Level
                %                  / \       |        / \
                %LowerLevel ______/   \______|_______/   \_____
                %
                %
                %Once the signal was processed the next step is through a
                %comparator decide the actual information received.
                %
                %For accurately read the income signal is necessary more than
                %just have the decision level, it is also needed to know where
                %it occurs. For instance, if we decide to use just to take the
                %decision level and measure the mean value across and a portion
                %of the time period, which portion should we take? It would be
                %logical to take the central portion, but the bit may be
                %deformed in such way that the information gets concentrated
                %out of the middle part, which means the symbol is not
                %symmetric. The symmetry information can be acquired from the
                %eye diagram by measuring the longitudinal opening. The
                %following sketch better describes this process:
                %
                %                   Point of Symmetry
                %         _____     _______|_______     _____
                %              \   /               \   /
                %               \ /  Longitudinal   \ /
                %                \ ----------------- /
                %               / \    Opening      / \
                %         _____/   \_______________/   \_____
                %
                %With those two pieces of information, decision level and point
                %of symmetry, we have the X, Y coordinates for the centre of
                %the Eye Diagram. Therefore, as long as there is an opening on
                %it it will be possible to recover the transmitted information
                %without error... theoretically.
                %
                %As this process is also statistical, first we reshape the
                %income vector to analyze all periods at the same time.
%                 EyeSymMat = reshape(Esync(1+SyncPeriod*NPPB:end-SyncPeriod*...
%                     NPPB),NPPB,NbDPSK-2*SyncPeriod);
                %Then we take the values that compose the decision level
                %because they will mark the point of symmetry.
                %
                %Firstly it was set the interval in which the histogram will be
                %build. It is based on the number of samples per bit period.
%                 Interval = linspace(min(Esync(1+SyncPeriod*NPPB:end-...
%                     SyncPeriod*NPPB)),max(Esync(1+SyncPeriod*NPPB:end-SyncPeriod...
%                     *NPPB)),2*NPPB);
                %Therefore, the MATLAB hist function returns the number of
                %occurrence of each interval.
%                 EyeMax = hist(Esync,Interval);
%                 EyeMaxaux = [0 EyeMax 0];                                      %Zeros are added at the EyeMax to auxiliate the finding peaks process
%                 [~,EyeLoc] = findpeaks(EyeMaxaux,'MinPeakDistance',NPPB/4,...
%                     'SortStr','descend','NPeaks',2);%The peaks on the Eye profile will be the levels at the Eyes limit

                %From the location of the max values that occured, which means
                %the uper and lower level of the eye diagram it needs to take
                %the actual value that those occurences represent that is
                %withing the the Interval variable.
%                 ValToSeek = Interval(EyeLoc-1);
                %The number of ocurrences is a statical measure therefore one
                %does not have control which interval will have the highest
                %peak, thus it is important to ordenate the values to be seek
                %from the lower part of the eye diagram to the uper part of the
                %eye diagram.
%                 ValToSeek = sort(ValToSeek,'ascend');
%                 OccuCount = zeros(1,size(EyeSymMat,1));                        %Auxiliar Variable for accounting.
%                 for kk=1:size(EyeSymMat,1)                                     %For every sample within a symbol period

%                     OccuCount(kk) = OccuCount(kk)+sum((EyeSymMat(kk,:)>=min(...
%                         Esync))&(EyeSymMat(kk,:)<=UpeSymPer*EyeOpenLow)); %Account all occurencies of the value 1
%                     OccuCount(kk) = OccuCount(kk)+sum((EyeSymMat(kk,:)>=...
%                         LowSymPer*EyeOpenHigh)&(EyeSymMat(kk,:)<=max(Esync))); %Account all occurencies of the value 2
%                 end
                %The point of symmetry of the eye diagram will be where the
                %maximum number of occurrences were measured inasmuch as those
                %are the points where all the bits go to the center of the
                %symbol.  From the maximum number of occurrences, it can happen
                %for more than one sample within one symbol period, in the best
                %case, all samples would have the same accounting as it is
                %shown the ilustration above hence the symmetry will be at the
                %middle sample of this group of maximum occurrences. This value
                %can be found by the mean of the samples positions within the
                %symbol period. The problem with this approach is that the
                %signal must be synchronized with the maximum displacement of
                %a symbol period minus 25% of the eye Longitudinal opening if
                %the displacement is higher than that the point of symmetry
                %will be wrongly measured.
%                 [SymLoc] = round(mean(find(ismember(OccuCount,max(OccuCount)...
%                     ))));
                %             figure;findpeaks(OccuCount,'SortStr','descend');

                %% Actualy Receiving Data:
                Data = [];                                                  %Initialization of the vector that will store the income data
                DataM = [];                                                  %Initialization of the vector that will store the income data
                DataL = [];                                                  %Initialization of the vector that will store the income data
                DataU = [];                                                  %Initialization of the vector that will store the income data
                DataS = [];                                                  %Initialization of the vector that will store the income data
                for kk=1:NPPB:length(Esync)                                    %The comparison process will be made for each symbol period
                    %                 MeanOfData = mean(EoutI((kk-1)+SymLocI(1)));
                    MeanOfData = mean(Esync((kk-1)+NPPB/2));
                    if MeanOfData > LevelDec3%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                        Data = [Data 1];                                 %is 1
                    else                                                       %If it is the lowest level the incoming data
                        Data = [Data 0];                                 %is 0
                    end
                    MeanOfData = mean(Esync((kk-1)+(NPPB/2)+(+NPPB/16)));
                    if MeanOfData > LevelDec1%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                        DataM = [DataM 1];                                 %is 1
                    else                                                       %If it is the lowest level the incoming data
                        DataM = [DataM 0];                                 %is 0
                    end
                    MeanOfData = mean(Esync((kk-1)+(NPPB/2)-(+NPPB/16)));
                    if MeanOfData > LevelDec2%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                        DataL = [DataL 1];                                 %is 1
                    else                                                       %If it is the lowest level the incoming data
                        DataL = [DataL 0];                                 %is 0
                    end
                    MeanOfData = mean(Esync((kk-1)+(NPPB/2)));
                    if MeanOfData > LevDefDpqsk%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                        DataU = [DataU 1];                                 %is 1
                    else                                                       %If it is the lowest level the incoming data
                        DataU = [DataU 0];                                 %is 0
                    end
                    if CalcS
                        MeanOfData = mean(Esync((kk-1)+(NPPB/2)));
                        if MeanOfData > EyeOpenLowI + EyeOpenI/2%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                            DataS = [DataS 1];                                 %is 1
                        else                                                       %If it is the lowest level the incoming data
                            DataS = [DataS 0];                                 %is 0
                        end
                    end
                end
                %%       Calculating the Bit Error Ratio (BER)
                %The final process here is to count the number of wrongdoings
                %of this whole process upon the transmited data for
                %quantitative analizes
                TxData  = TxDataMat(ThisCarr,1+SyncPeriod:end-SyncPeriod);
               
                if CalcS
                    DataS  = DataS(1+SyncPeriod:end-SyncPeriod);
                    BitErrS    = sum(xor(TxData,DataS));                       %Comparison between the Transmited and received and counting the differences
                    BerDPSKS(CurrentTest,ThisCarr) = BitErrS/(NbDPSK-(2*SyncPeriod));                  %Calculating the ration of wrong bits and the total number of bits transmited
                end
                Data  = Data(1+SyncPeriod:end-SyncPeriod);
                DataM  = DataM(1+SyncPeriod:end-SyncPeriod);
                DataL  = DataL(1+SyncPeriod:end-SyncPeriod);
                DataU  = DataU(1+SyncPeriod:end-SyncPeriod);
                AberLevAuxI(4) = 0;
                ValsLevAuxI(4) = 0.01;

                BitErr(1)  = sum(xor(TxData,Data));                        %Comparison between the Transmited and received and counting the differences
                BitErr(2)  = sum(xor(TxData,DataM));                       %Comparison between the Transmited and received and counting the differences
                BitErr(3)  = sum(xor(TxData,DataL));                       %Comparison between the Transmited and received and counting the differences
                BitErr(4)  = sum(xor(TxData,DataU));                       %Comparison between the Transmited and received and counting the differences
                BerDPSK(CurrentTest,ThisCarr) = (min(BitErr))/(NbDPSK-(2*SyncPeriod));%Calculating the ration of wrong bits and the total number of bits transmited
                AberIAux = max(AberLevAuxI(BitErr<=min(BitErr)));
                AberLev(CurrentTest,ThisCarr)  = AberIAux;
                ValsLev(CurrentTest,ThisCarr)  = max(ValsLevAuxI(AberLevAuxI==AberIAux));          %Comparison between the Transmited and received and counting the differences
                if BerDPSK(CurrentTest,ThisCarr)==0
                    LevDefDpqsk = max(ValsLevAuxI(AberLevAuxI==AberIAux));
                    if LevDefDpqsk==0
                        LevDefDpqsk = 0.01;
                    end
                end
                [~,LocAux] = max(AberLevAuxI==AberIAux);
                EvmDB(CurrentTest,ThisCarr)   = EvmDBA(LocAux);
                EvmPer(CurrentTest,ThisCarr)  = EvmPerA(LocAux);
                EvmRms(CurrentTest,ThisCarr)  = EvmRmsA(LocAux);
                EvmDBJ(CurrentTest,ThisCarr)  = EvmDBJA(LocAux);
                EvmPerJ(CurrentTest,ThisCarr) = EvmPerJA(LocAux);
                EvmRmsJ(CurrentTest,ThisCarr) = EvmRmsJA(LocAux);
                
                EvmMatRec(ObsCarrPos==ThisCarr,:) = EvmMatRecA(LocAux,:); 
                %% Ploting the result for qualitative analizes
                %PrintInfo(Ploting*19,TxData,Data);
                %%
                %             berpos = 1:2:size(BerDPSK,2);
                %             BerDPSK(size(BerDPSK,1),berpos)
            %             a=6;
                        close all;
            end
        end
    case 'DQPSK'
        
        %%                   Receiver DQPSK
        for ThisCarr=InitCarr:CarrPass:NumCarr                                           %For each carrier the same process of reception will be used.
            if CarrUsedDo(ThisCarr)
                %%  Reception Process: Handling the income optical field
                %At the first step the income field will be selected from the
                %optical FFT Output with the help of the maping vector
                %(previously described).
                EoutAux = EoutAux1(VetThisCarr==ObsCarrUsed(ThisCarr),:);
                if AddingNoiseP
                    [~,sigpower] = MeasPower(EoutAux);
                    EoutAux      = awgn(EoutAux,osnrp,sigpower-30);
                    [~,sigpower] = MeasPower(EoutAux);
                end
                %%            Fiber Time Delay Compensation
                switch Medium
                    case 'Fiber'
                        EoutAux = ifft(fft(EoutAux).*exp(1j*2*pi*f*(...
                            FiberDelay(ThisCarr)*Ta)));
                    otherwise
                end
                %             EoutAux = EoutRec;
                %%  Plot for Qualitative analizes
                %PrintInfo(Ploting*20,t,abs(EoutAux));
                %             PrintInfo(Ploting*21,f,20*log10(abs(fftshift(fft(EoutAux)./...
                %                                                        length(EoutAux)))));

                %%                   synchronizing
                %Differently from the previuos process at the DQPSK
                %demodulation it is needed to make the synchronism before
                %interferometric process. Because the recovery of the
                %transmited data is extremally dependent of the time delay.
                %Thus, for the reception process work properly, it is needed to
                %sincronized the recieved signal or the sampling process will
                %not work.

                %At this first moment a small part of the income signal needs
                %be analized for the sincronism process. As we are not looking
                %for the actual data rather for the synchronis symbol, which is
                %the bigest phase shift of the income signal. The next delay
                %interferometer convert this phase shift to a amplitude
                %variation.

                %             %The phase delay is not needed because the optical field will
                %             %be analized as a whole.
                %             PhaDel          = 0;
                %             %Remember that at DQPSK modulation the information is stored at
                %             %the difference of phase between the current and previous
                %             %symbol hence this time delay is needed.
                %             TimDel          = T;
                %             [ESync1,ESync2] = DelayInterf(t,TimDel,PhaDel,EoutAux);        %Coverting phase shift to amplitude variation
                %
                %             %The second moment is to transfer this information from the
                %             %optical domain to the eletrical domain for an eletronic
                %             %processing. It is done with the help of an photo diode. The
                %             %configuration here used is an balanced receiver as the output
                %             %of the Delay Interferometer has two signals resulting from
                %             %the constructive and destructive signal interaction.
                %             ESync1 = ESync1.*conj(ESync1);
                %             ESync2 = ESync2.*conj(ESync2);
                %             Esync  = ESync2-ESync1;
                %
                %             ESync1 = ESync1./max(abs(ESync1));                             %Normalizing the signal
                %             ESync2 = ESync2./max(abs(ESync2));                             %Normalizing the signal
                %             Esync  = Esync./max(abs(Esync));                               %Normalizing the signal
                %
                %             SyncAux   = Esync(IniSyncPos:SyncPos);                         %Selecting just the symbol to synchronize
                %             SyncedAux = SyncSymb(IniSyncPos:SyncPos);                      %Selecting just the symbol to synchronize
                %             %%                   Plot for Qualitative analizes
                % %              PrintInfo(Ploting*22,t(IniSyncPos:SyncPos),Esync(IniSyncPos...
                % %              :SyncPos),SyncSymb(IniSyncPos:SyncPos),ESync1(IniSyncPos:...
                % %                                       SyncPos),ESync2(IniSyncPos:SyncPos));
                %             %%                   Synchronizing
                %             %This synchronization process is based on the mean expected
                %             %value. That means, the information mean value should be within
                %             %one period of symbol. Thus, the mean value of the received
                %             %signal is acquired and compare of the known sync-word to
                %             %verify if this mean value is at the right possition. Which is
                %             %the midel point (peak) of the highest level at the sync period
                %             SyncAux(SyncAux<0)              = 0;                           %To keep the mean value above zero anything under is neglected
                %             SyncAux(SyncAux>=mean(SyncAux)) = 1;                           %Adding a flag to the first sample of the received mean value
                %             SyncAux(SyncAux<mean(SyncAux))  = -1;                          %All the others samples at set to the lowest level
                %
                %             PosToSyn  = find(ismember(SyncAux,1));                         %Finding where is the location of the samples to synchronize
                %             PosSynced = find(ismember(SyncedAux,1));                       %Finding where is the location of the samples to synchronize
                %
                %             DiffPos = PosToSyn(round(end/2)) - PosSynced(round(end/2));    %Accounting the peak (midel point) displacement
                %             if DiffPos>=0%If the difference is positive, left-shift...
                %                 EoutAux = ifft(fft(EoutAux).*exp(1j*2*pi*f*(DiffPos*Ta))); %Shift based on time change
                % %                 EoutAux = [EoutAux(DiffPos+1:end) EoutAux(1:DiffPos)];   %Shift based on sampling sliding
                %             else%... but if the difference is negative, right-shift
                %                 EoutAux = ifft(fft(EoutAux).*exp(1j*2*pi*f*(DiffPos*Ta))); %Shift based on time change
                % %                 EoutAux = [EoutAux(end+DiffPos+1:end) EoutAux(1:end + ...
                % %                                                                 DiffPos)];%Shift based on sampling sliding
                %             end
                %
                %                 %                   Synchronizing
                %             %Because of reasons, sometimes it may be required to make a
                %             %synchronization process with the  end of the data stream as
                %             %well. This following verification check if the user set (or
                %             %not) a second synchronization process to be done.
                %             if SencondAdjust
                %                 PhaDel          = 0;
                %                 TimDel          = T;
                %                 [ESync1,ESync2] = DelayInterf(t,TimDel,PhaDel,EoutAux);
                %
                %                 ESync1 = ESync1.*conj(ESync1);
                %                 ESync2 = ESync2.*conj(ESync2);
                %                 Esync  = ESync2-ESync1;
                %
                %                 ESync1 = ESync1./max(abs(ESync1));
                %                 ESync2 = ESync2./max(abs(ESync2));
                %                 Esync  = Esync./max(abs(Esync));
                %
                %
                %                 %This synchronization process is based on the mean expected
                %                 %value. That means, the information mean value should be
                %                 %within one period of symbol. Thus, the mean value of the
                %                 %received signal is acquired and compare of the known
                %                 %sync-word to verify if this mean value is at the right
                %                 %possition.
                %
                %                 SyncAuxEnd   = Esync(end-SyncPos+1:end-IniSyncPos+1);
                %                 SyncedAuxEnd = SyncSymbEnd(end-SyncPos+1:end-IniSyncPos+1);
                %                 SyncAuxEnd(SyncAuxEnd<0)                 = 0;              %To keep the mean value above zero anything under is neglected
                %                 SyncAuxEnd(SyncAuxEnd>=mean(SyncAuxEnd)) = 1;              %Adding a flag to the first sample of the received mean value
                %                 SyncAuxEnd(SyncAuxEnd<mean(SyncAuxEnd))  = -1;             %All the others samples at set to the lowest level
                %
                %                 PosToSynEnd  = find(ismember(SyncAuxEnd,1));               %Finding where is the location of the first sample to synchronize
                %                 PosSyncedEnd = find(ismember(SyncedAuxEnd,1));             %Finding where is the location of the first sample to synchronize
                %
                %                 DiffPosEnd = PosToSynEnd(round(end/2))-PosSyncedEnd(...
                %                                                              round(end/2));
                %                 if DiffPosEnd>=0%If positive difference, left-shift...
                %                     EoutAux = ifft(fft(EoutAux).*exp(1j*2*pi*f*(...
                %                                                           DiffPosEnd*Ta)));%Shift based on time change
                % %                     EoutAux = [EoutAux(DiffPosEnd+1:end) EoutAux(1:...
                % %                                                             DiffPosEnd)];%Shift based on sampling sliding
                %                 else%... but if the difference is negative, right-shift
                %                     EoutAux = ifft(fft(EoutAux).*exp(1j*2*pi*f*(...
                %                                                           DiffPosEnd*Ta)));%Shift based on time change
                % %                     EoutAux = [EoutAux(end+DiffPosEnd+1:end) EoutAux(...
                % %                                                       1:end+DiffPosEnd)];%Shift based on sampling sliding
                %                 end
                %             end

                %%                   Plot for Qualitative analizes
                %             PrintInfo(Ploting*23,t(IniSyncPos:SyncPos),Esync(IniSyncPos:...
                %                 SyncPos),SyncSymb(IniSyncPos:SyncPos),ESync1(IniSyncPos:...
                %                                       SyncPos),ESync2(IniSyncPos:SyncPos));
                %%               Recovering the Information
                %Once the Income signal was synchronized it is possible to
                %split its components in phase and in quadrature. It is very
                %important to understand here that those components will carry
                %the actual data transmited. From it we can't recover the  I
                %and Q eletrical signals that were used to modulate our
                %carrier. Those signals get mingled together to composed the
                %actual data. Therefore, when we pass the income signal
                %throughout the receptor we will recover the TxData generated
                %and transmited. It is important that it get very clear here
                %because I toke two months to understand it, no one could
                %explain it well to me and those articles about it make this
                %information unclear. Thus, to sheed a light on this concept,
                %it is important to know that the Data is encoded in the I and
                %Q eletrical components that will get mixed together to
                %modulate the optical carrier. When this optical signal passes
                %throughout the MZ-Interferometer the acutal data (TX) will be
                %recover, not the I and Q component. But the data at the odd
                %possition (encoded in I) will be at the real component (in
                %phase), whereas the data at the even possition (encoded in Q)
                %will be at the imaginary component(in quadrature). The
                %MZ-Interferometer will be resposible to take apart the real
                %and imaginary components of the income optical field.

                E_rec3 = EoutAux./max(abs(EoutAux));                           %Normalizing income signal
                %For the interferometric process  take in account just the real
                %component it is needed a phase delay of 45° degrees;
                PhaDel = 1*pi/4;
                %Remember that at DQPSK modulation the information is stored at
                %the difference of phase between the current and previous
                %symbol hence this time delay of one symbol period is needed.
                TimDel = T;
                [EoutA,EoutB] = DelayInterf(t,TimDel,PhaDel,EoutAux);          %Coverting phase shift to amplitude variation
                %For the interferometric process  take in account just the real
                %component it is needed a phase delay of -45° degrees;
                PhaDel = -1*pi/4;
                %Remember that at DQPSK modulation the information is stored at
                %the difference of phase between the current and previous
                %symbol hence this time delay of one symbol period is needed.
                TimDel = T;
                [EoutC,EoutD] = DelayInterf(t,TimDel,PhaDel,EoutAux);          %Coverting phase shift to amplitude variation

                %The second moment is to transfer this information from the
                %optical domain to the eletrical domain for an eletronic
                %processing. It is done with the help of an photo diode.
                EoutA = EoutA.*conj(EoutA);
                EoutB = EoutB.*conj(EoutB);
                EoutC = EoutC.*conj(EoutC);
                EoutD = EoutD.*conj(EoutD);

                %The process with the photo diode is self-coherent, which means
                %the rusult will be a component of the signal centered in f=0
                %and another component centered at f=2*fc (frenquecy central).
                %Therefore, to remove the higher order component a low pass
                %filter will be used.
                %%           Creating the Reception Filter
                [BitFilt,~] = FiltroGaussiano(f,BWD2,CenFeq2,FiltOrd2);        %Creating filter for selection of the received signal
                BitFilt = fftshift(BitFilt);                                   %Shifting the filter for matching the received signal
                %%
                
                if RecFilBanPas
                    EoutA = ifft(fft(EoutA).*BitFilt);
                    EoutB = ifft(fft(EoutB).*BitFilt);
                    EoutC = ifft(fft(EoutC).*BitFilt);
                    EoutD = ifft(fft(EoutD).*BitFilt);
                end

                % The configuration here used is an balanced receiver as the
                %output of the Delay Interferometer has two signals resulting
                %from the constructive and destructive signal interaction.
                %%         Adding Noise to Received Signal
                %Just after the photo diod is where the noise must be added as
                %it is the main constrain of the system. Once the signal is
                %converted from optical to electrical and it digital values
                %were recevered, there is no point to add noise, because the
                %information was already received, it just needs decoding.
                %The noise here added will be gaussian basically it represents
                %the reception sensibility. In another words, how much the
                %signal must be about the noise for an acceptable recovering.

                if ReceptorNoise                                               %Verify whether the noise is to be added or not
                    [~,SigPowerA] = MeasPower(EoutA);                          %it own received signal or
                    [~,SigPowerB] = MeasPower(EoutB);
                    [~,SigPowerC] = MeasPower(EoutC);                          %it own received signal or
                    [~,SigPowerD] = MeasPower(EoutD);
                    SNR2 = CarSNR + 10*log10(2) - 10*log10(Nsamp) + 10*log10(10^0.3);
                    SNR = CarSNR + 10*log10(2) - 10*log10((t(NPPB))/(2*(t(NPPB))/(NPPB)));
                    SigPowerA2 = SigPowerA-30;%20*log10(SigPowerI);
                    SigPowerB2 = SigPowerB-30;%20*log10(SigPowerQ);
                    SigPowerC2 = SigPowerC-30;%20*log10(SigPowerI);
                    SigPowerD2 = SigPowerD-30;%20*log10(SigPowerQ);
                    EoutA = awgn(EoutA,SNR,SigPowerA2);
                    EoutB = awgn(EoutB,SNR,SigPowerB2);
                    EoutC = awgn(EoutC,SNR,SigPowerC2);
                    EoutD = awgn(EoutD,SNR,SigPowerD2);
                end

                
                if ~RecFilBanPas
                    EoutA = ifft(fft(EoutA).*BitFilt);
                    EoutB = ifft(fft(EoutB).*BitFilt);
                    EoutC = ifft(fft(EoutC).*BitFilt);
                    EoutD = ifft(fft(EoutD).*BitFilt);
                end
                EoutI = (EoutB - EoutA);
                EoutQ = (EoutD - EoutC);
                EoutI = EoutI-mean(EoutI);
                EoutQ = EoutQ-mean(EoutQ);
                EoutI = EoutI./max(abs(EoutI));                                %Normalizing the signal
                EoutQ = EoutQ./max(abs(EoutQ));                                %Normalizing the signal
                VetElecPowerI(CurrentTest,ThisCarr)= MeasPower(EoutI);
                VetElecPowerQ(CurrentTest,ThisCarr)= MeasPower(EoutQ);
                EoutAuxF = 20*log10(abs(fftshift(fft(EoutAux)./length(...
                    EoutAux))));
                [VetOptiPower(CurrentTest,ThisCarr),~]= findpeaks(EoutAuxF,...
                    'SortStr','descend','NPeaks',1);


%                 SyncAuxI   = EoutI(IniSyncPos:SyncPos);                        %Selecting just the symbol to synchronize
%                 SyncAuxQ   = EoutQ(IniSyncPos:SyncPos);                        %Selecting just the symbol to synchronize
% %                             figure;plot(SyncAuxQ);set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%                 SyncedAux  = SyncSymb(IniSyncPos:SyncPos);                     %Selecting just the symbol to synchronize
%                 %%                   Synchronizing
%                 %This synchronization process is based on the mean expected
%                 %value. That means, the information mean value should be within
%                 %one period of symbol. Thus, the mean value of the received
%                 %signal is acquired and compare of the known sync-word to
%                 %verify if this mean value is at the right possition. Which is
%                 %the midel point (peak) of the highest level at the sync period
%                 SyncAuxI(EoutI(IniSyncPos:SyncPos)<0.5*max((EoutI(IniSyncPos:SyncPos))))               = 0;   %To keep the mean value above zero anything under is neglected
%                 SyncAuxQ(EoutQ(IniSyncPos:SyncPos)<0.5*max((EoutQ(IniSyncPos:SyncPos))))               = 0;   %To keep the mean value above zero anything under is neglected
%                 SyncAuxI(SyncAuxI>=mean(SyncAuxI)) = 1;   %Adding a flag to the first sample of the received mean value
%                 SyncAuxQ(SyncAuxQ>=mean(SyncAuxQ)) = 1;   %Adding a flag to the first sample of the received mean value
%                 SyncAuxI(SyncAuxI<mean(SyncAuxI))  = -1;  %All the others samples at set to the lowest level
%                 SyncAuxQ(SyncAuxQ<mean(SyncAuxQ))  = -1;  %All the others samples at set to the lowest level
% 
%                 PosToSynI  = find(ismember(SyncAuxI,1));                       %Finding where is the location of the samples to synchronize
%                 PosToSynQ  = find(ismember(SyncAuxQ,1));                       %Finding where is the location of the samples to synchronize
%                 PosSynced = find(ismember(SyncedAux,1));                       %Finding where is the location of the samples to synchronize
% 
%                 DiffPosI = ExtDel*(PosToSynI(round(end/2)) - PosSynced(round(end/2)));  %Accounting the peak (midel point) displacement
%                 DiffPosQ = ExtDel*(PosToSynQ(round(end/2)) - PosSynced(round(end/2)));  %Accounting the peak (midel point) displacement
%                 if DiffPosI>=0%If the difference is positive, left-shift...
%                     %                 EoutI = ifft(fft(EoutI).*exp(1j*2*pi*f*(DiffPosI*Ta))); %Shift based on time change
%                     EoutI = [EoutI(DiffPosI+1:end) EoutI(1:DiffPosI)];   %Shift based on sampling sliding
%                     EoutQ = [EoutQ(DiffPosI+1:end) EoutQ(1:DiffPosI)];   %Shift based on sampling sliding
%                 else%... but if the difference is negative, right-shift
%                     %                 EoutI = ifft(fft(EoutI).*exp(-1j*2*pi*f*(DiffPosI*Ta))); %Shift based on time change
%                     EoutI = [EoutI(end+DiffPosI+1:end) EoutI(1:end+DiffPosI)]; %Shift based on sampling sliding
%                     EoutQ = [EoutQ(end+DiffPosI+1:end) EoutQ(1:end+DiffPosI)]; %Shift based on sampling sliding
%                 end
                %             if DiffPosQ>=0%If the difference is positive, left-shift...
                % %                 EoutQ = (ifft(fft(EoutQ).*exp(1j*2*pi*f*(DiffPosQ*Ta)))); %Shift based on time change
                %                 EoutQ = [EoutQ(DiffPosQ+1:end) EoutQ(1:DiffPosQ)];   %Shift based on sampling sliding
                %             else%... but if the difference is negative, right-shift
                % %                 EoutQ = (ifft(fft(EoutQ).*exp(-1j*2*pi*f*(DiffPosQ*Ta)))); %Shift based on time change
                %                 EoutQ = [EoutQ(end+DiffPosQ+1:end) EoutQ(1:end+DiffPosQ)]; %Shift based on sampling sliding
                %             end
                %%                   Plot for Qualitative analizes
%                              PrintInfo(PlotingThis*22,t(IniSyncPos:SyncPos),EoutI(IniSyncPos...
%                              :SyncPos),SyncSymb(IniSyncPos:SyncPos),EoutA(IniSyncPos:...
%                                                       SyncPos),EoutB(IniSyncPos:SyncPos));
%                              PrintInfo(PlotingThis*22,t(IniSyncPos:SyncPos),EoutQ(IniSyncPos...
%                              :SyncPos),SyncSymb(IniSyncPos:SyncPos),EoutC(IniSyncPos:...
%                                                       SyncPos),EoutD(IniSyncPos:SyncPos));
                %% Removing CP
                if AddCP
                    IxAux1  = EoutI(1:end - StuffSampels);
                    IxAux1  = reshape(IxAux1,(2*NumAmosCP+NPPB),NbDQPSK/2);
                    IxAux1  = IxAux1(1+NumAmosCP:end-NumAmosCP,:);
                    IxAux1  = reshape(IxAux1,1,NPPB*NbDQPSK/2);
                    EoutI  = IxAux1;

                    IxAux2  = EoutQ(1:end - StuffSampels);
                    IxAux2  = reshape(IxAux2,(2*NumAmosCP+NPPB),NbDQPSK/2);
                    IxAux2  = IxAux2(1+NumAmosCP:end-NumAmosCP,:);
                    IxAux2  = reshape(IxAux2,1,NPPB*NbDQPSK/2);
                    EoutQ  = IxAux2;

                    IxAux3  = E_rec3(1:end - StuffSampels);
                    IxAux3  = reshape(IxAux3,(2*NumAmosCP+NPPB),NbDQPSK/2);
                    IxAux3  = IxAux3(1+NumAmosCP:end-NumAmosCP,:);
                    IxAux3  = reshape(IxAux3,1,NPPB*NbDQPSK/2);
                    E_rec3  = IxAux3;
                end
                %% Taking the sampling the EVM meassurement
                clear IxAux;
                PosAuxEout1 = NPPB/2:NPPB:length(EoutI);%Varriable respossible to take just the samples at the middle of the symbol
                PosAuxEout2 = ((NPPB/2)+(NPPB/4)):NPPB:length(EoutI);%Varriable respossible to take just the samples at the middle of the symbol
                PosAuxEout3 = ((NPPB/2)-(NPPB/4)):NPPB:length(EoutI);%Varriable respossible to take just the samples at the middle of the symbol
                IxAux1      = EoutI(PosAuxEout1) + 1j.*EoutQ(PosAuxEout1);    %Normalizing the reference
                IxAux2      = EoutI(PosAuxEout2) + 1j.*EoutQ(PosAuxEout2);    %Normalizing the reference
                IxAux3      = EoutI(PosAuxEout3) + 1j.*EoutQ(PosAuxEout3);    %Normalizing the reference
                a=0;
                EvmMatRecA(1,:) = IxAux1;                       %Taking just the middle samples as references
                EvmMatRecA(2,:) = IxAux2;                       %Taking just the middle samples as references
                EvmMatRecA(3,:) = IxAux3;                       %Taking just the middle samples as references
                EvmMatRecA(4,:) = IxAux1;                       %Taking just the middle samples as references
                [EvmDBA(1), EvmPerA(1), EvmRmsA(1) ] = EvmCalc( EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux1 );
                [EvmDBA(2), EvmPerA(2), EvmRmsA(2) ] = EvmCalc( EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux2 );
                [EvmDBA(3), EvmPerA(3), EvmRmsA(3) ] = EvmCalc( EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux3 );
                [EvmDBA(4), EvmPerA(4), EvmRmsA(4) ] = EvmCalc( EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux1 );
                [EvmDBJA(1),EvmPerJA(1),EvmRmsJA(1)] = evm1(4,'dpsk',EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux1);
                [EvmDBJA(2),EvmPerJA(2),EvmRmsJA(2)] = evm1(4,'dpsk',EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux2);
                [EvmDBJA(3),EvmPerJA(3),EvmRmsJA(3)] = evm1(4,'dpsk',EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux3);
                [EvmDBJA(4),EvmPerJA(4),EvmRmsJA(4)] = evm1(4,'dpsk',EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux1);
%                 if ThisCarr==1
%                     figure;hold all;grid on;
%                     plot(real(IxAux1),imag(IxAux1),'rx');
%                     plot(real(IxAux2),imag(IxAux2),'yx');
%                     plot(real(IxAux3),imag(IxAux3),'gx');
%                     set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%                     a=0;
%                 end
                %%
%                 figure;plot(t(1:length(E_rec3)),abs(E_rec3),t(1:length(EoutQ)),EoutQ,t(1:length(EoutI)),EoutI);set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                %One important step on this project was the confirmation of our
                %results from a closed and theoretical equation that relates
                %the income optical field witht its respective data (as
                %descrived above in the receiving process). The result from
                %this equation was further compared with the result from the
                %MZ-Interferometer as an proof of concept. This equation can be
                %found at the book of Optical Fiber Telecommunications V B,
                %which one of the authors is Ivan P. Kaminow at the page 144.
                taux = t(1:length(E_rec3));
                faux = time2freq(taux);
                Ui = real(exp(1j*pi/4).*E_rec3.*conj(ifft(fft(E_rec3).*exp(-...
                    1j*2*pi*faux*T))));%The data at odd position
                Ui = Ui./max(abs(Ui));                                         %Normalizing the signal
                Uq = imag(exp(1j*pi/4).*E_rec3.*conj(ifft(fft(E_rec3).*exp(-...
                    1j*2*pi*faux*T))));%The data at the even position
                Uq = Uq./max(abs(Uq));                                         %Normalizing the signal
%##########################################################################
                PosIx = NPPB/2:NPPB:length(EoutI);                                %Possition of the central samples - but the number of samples per symbol is even ????
                IxAux = EoutI(PosIx);                                             %From the main received signal just a few samples are taken for further evaluation
                n     = 100;                                                   %The number of boxes to be filled up on the histogram process
                Perce = 0.7;
                %At this process each eye will be evaluated separetelly. It is
                %important to remember that the PAM4 has 4 levels, which means
                %three level of decissions that are exaclty the center of the
                %eye diagram.
                IxAuxAB = EoutI(PosIx);       %Taking just those values relative to the uper eye
                InterAB = linspace(Perce*(min(EoutI)),Perce*(max(EoutI)),n);                     %Building the histogram boxes
                EyeAB = hist(IxAuxAB,InterAB);                                 %filling up the boxes with samples that fit on them.
                
                QIxAuxAB = EoutQ(PosIx);       %Taking just those values relative to the uper eye
                QInterAB = linspace(Perce*(min(EoutQ)),Perce*(max(EoutQ)),n);                     %Building the histogram boxes
                QEyeAB = hist(QIxAuxAB,QInterAB);                                 %filling up the boxes with samples that fit on them.
                %The same process described for the uper level will be done at
                %the middle and lower eyes levels.
                IxAuxCD = EoutI(PosIx-NPPB/4);
                InterCD = linspace(Perce*(min(EoutI)),Perce*(max(EoutI)),n);%NPPB*2^n);
                EyeCD = hist(IxAuxCD,InterCD);
                
                QIxAuxCD = EoutQ(PosIx-NPPB/4);
                QInterCD = linspace(Perce*(min(EoutQ)),Perce*(max(EoutQ)),n);%NPPB*2^n);
                QEyeCD = hist(QIxAuxCD,QInterCD);
                
                IxAuxEF = EoutI(PosIx+NPPB/4);
                InterEF = linspace(Perce*(min(EoutI)),Perce*(max(EoutI)),n);%NPPB*2^n);
                EyeEF = hist(IxAuxEF,InterEF);
                
                QIxAuxEF = EoutQ(PosIx+NPPB/4);
                QInterEF = linspace(Perce*(min(EoutQ)),Perce*(max(EoutQ)),n);%NPPB*2^n);
                QEyeEF = hist(QIxAuxEF,QInterEF);
                
                %What we are looking for here are not where exist the
                %occurrences of values. The eye is where there is not samples,
                %this means, the decission level is a reagion between actual
                %levels thus this reagion does not contain any signal.
                EyeAB = ~EyeAB;                                                %Changing zeros to one - Zeros compose the eye region
                EyeCD = ~EyeCD;
                EyeEF = ~EyeEF;
                
                QEyeAB = ~QEyeAB;                                                %Changing zeros to one - Zeros compose the eye region
                QEyeCD = ~QEyeCD;
                QEyeEF = ~QEyeEF;
                %Starting variables that will be used to identify the eye
                %diagram.
                CountAB   = 1;
                CountCD   = 1;
                CountEF   = 1;
                
                QCountAB   = 1;
                QCountCD   = 1;
                QCountEF   = 1;
                
                SeqOnesAB = 0;
                SeqOnesCD = 0;
                SeqOnesEF = 0;
                
                QSeqOnesAB = 0;
                QSeqOnesCD = 0;
                QSeqOnesEF = 0;
                
                SeqFinAB  = 0;
                SeqFinCD  = 0;
                SeqFinEF  = 0;
                
                QSeqFinAB  = 0;
                QSeqFinCD  = 0;
                QSeqFinEF  = 0;
                
                SeqIniAB  = 1;
                SeqIniCD  = 1;
                SeqIniEF  = 1;
                
                QSeqIniAB  = 1;
                QSeqIniCD  = 1;
                QSeqIniEF  = 1;
                %The for loop will take account of every box with ones. It is
                %important to take note that the not operator was used in this
                %vector, therefore ones means zeros (the eye diagram -
                %possibly) and zeros means values abouve zeroa (not the eye).
                for kk=1:length(EyeAB)                                         %For every box
                    if EyeAB(kk)                                               %if it contains "1"
                        SeqOnesAB(SeqIniAB)=CountAB;                           %count this element as part of a consecutive sequency
                        CountAB = CountAB + 1;                                 %adds one to the counter of consecutive elements "1"
                        if kk==length(EyeAB)                                   %if the current box is the last box we got to an end
                            SeqFinAB(SeqIniAB) = kk;                           %The final sequency element is equal to its possition (kk)
                        end
                    else                                                       %else if the current box contains "0"
                        SeqFinAB(SeqIniAB) = kk-1;                             %the previous element was the last element of a consecutive sequency of ones
                        SeqIniAB = SeqIniAB + 1;                               %adds one to the counter of consecutive number (take in account how many sequencies there are)
                        CountAB = 1;                                           %reset the counter of consecutive elements to mark the start of a new consecutive sequence
                    end
                    
                    if QEyeAB(kk)                                               %if it contains "1"
                        QSeqOnesAB(QSeqIniAB)=QCountAB;                           %count this element as part of a consecutive sequency
                        QCountAB = QCountAB + 1;                                 %adds one to the counter of consecutive elements "1"
                        if kk==length(QEyeAB)                                   %if the current box is the last box we got to an end
                            QSeqFinAB(QSeqIniAB) = kk;                           %The final sequency element is equal to its possition (kk)
                        end
                    else                                                       %else if the current box contains "0"
                        QSeqFinAB(QSeqIniAB) = kk-1;                             %the previous element was the last element of a consecutive sequency of ones
                        QSeqIniAB = QSeqIniAB + 1;                               %adds one to the counter of consecutive number (take in account how many sequencies there are)
                        QCountAB = 1;                                           %reset the counter of consecutive elements to mark the start of a new consecutive sequence
                    end
                    
                    if EyeCD(kk)
                        SeqOnesCD(SeqIniCD)=CountCD;
                        CountCD = CountCD + 1;
                        if kk==length(EyeCD)
                            SeqFinCD(SeqIniCD) = kk;
                        end
                    else
                        SeqFinCD(SeqIniCD) = kk-1;
                        SeqIniCD = SeqIniCD + 1;
                        CountCD = 1;
                    end
                    
                    if QEyeCD(kk)
                        QSeqOnesCD(QSeqIniCD)=QCountCD;
                        QCountCD = QCountCD + 1;
                        if kk==length(QEyeCD)
                            QSeqFinCD(QSeqIniCD) = kk;
                        end
                    else
                        QSeqFinCD(QSeqIniCD) = kk-1;
                        QSeqIniCD = QSeqIniCD + 1;
                        QCountCD = 1;
                    end
                    
                    if EyeEF(kk)
                        SeqOnesEF(SeqIniEF)=CountEF;
                        CountEF = CountEF + 1;
                        if kk==length(EyeEF)
                            SeqFinEF(SeqIniEF) = kk;
                        end
                    else
                        SeqFinEF(SeqIniEF) = kk-1;
                        SeqIniEF = SeqIniEF + 1;
                        CountEF = 1;
                    end
                    
                    if QEyeEF(kk)
                        QSeqOnesEF(QSeqIniEF)=QCountEF;
                        QCountEF = QCountEF + 1;
                        if kk==length(QEyeEF)
                            QSeqFinEF(QSeqIniEF) = kk;
                        end
                    else
                        QSeqFinEF(QSeqIniEF) = kk-1;
                        QSeqIniEF = QSeqIniEF + 1;
                        QCountEF = 1;
                    end
                end
                
                
                %If the eye is open, which means there is a clear difference
                %between adjacent levels, the eye diagram will be the longest
                %sequence of ones.
                [MaxValAB,LocMaxAB]=max(SeqOnesAB);
                if LocMaxAB<2 || MaxValAB<2                                    %if any sequency was found or there is just one sequency it is a error thus
                    LevDec3 = 0.1;                                            %the decission level will be by default 0.67. Also, the other variables will
                    LocMaxAB = 1;                                              %will be set with values to not cause errors in the future
                    SeqFinAB(1)=2;
                    MaxValAB = 0;
                    InterAB(1)=LevDec3;
                else                                                           %if a sequency was found the middle of the eye will be the middle of the sequency
                    if (SeqFinAB(LocMaxAB)-MaxValAB/2)<1                       %if for some reason the final element of the sequency minus the half of its
                        LevDec3 = 0.1;                                         %length results in a negative value, something went very wrong, and by
                    else                                                       %default it will be set to 0.7
                        LevDec3 = InterAB(round(SeqFinAB(LocMaxAB)-MaxValAB/...
                            2));%Otherwise, the decission level is the middle of the sequency
                    end
                end
                [QMaxValAB,QLocMaxAB]=max(QSeqOnesAB);
                if QLocMaxAB<2 || QMaxValAB<2                                    %if any sequency was found or there is just one sequency it is a error thus
                    QLevDec3 = 0.1;                                            %the decission level will be by default 0.67. Also, the other variables will
                    QLocMaxAB = 1;                                              %will be set with values to not cause errors in the future
                    QSeqFinAB(1)=2;
                    QMaxValAB = 0;
                    QInterAB(1)=QLevDec3;
                else                                                           %if a sequency was found the middle of the eye will be the middle of the sequency
                    if (QSeqFinAB(QLocMaxAB)-QMaxValAB/2)<1                       %if for some reason the final element of the sequency minus the half of its
                        QLevDec3 = 0.1;                                         %length results in a negative value, something went very wrong, and by
                    else                                                       %default it will be set to 0.7
                        QLevDec3 = QInterAB(round(QSeqFinAB(QLocMaxAB)-QMaxValAB/...
                            2));%Otherwise, the decission level is the middle of the sequency
                    end
                end
                
                [MaxValCD,LocMaxCD]=max(SeqOnesCD);
                if LocMaxCD<2 || MaxValCD<2
                    LevDec2 = 0.1;
                    LocMaxCD = 1;
                    SeqFinCD(1)=2;
                    MaxValCD = 0;
                    InterCD(1)=LevDec2;
                else
                    if (SeqFinCD(LocMaxCD)-MaxValCD/2)<1
                        LevDec2 = 0.1;
                    else
                        LevDec2 = InterCD(round(SeqFinCD(LocMaxCD)-MaxValCD/2));
                    end
                end
                [QMaxValCD,QLocMaxCD]=max(QSeqOnesCD);
                if QLocMaxCD<2 || QMaxValCD<2
                    QLevDec2 = 0.1;
                    QLocMaxCD = 1;
                    QSeqFinCD(1)=2;
                    QMaxValCD = 0;
                    QInterCD(1)=LevDec2;
                else
                    if (QSeqFinCD(QLocMaxCD)-QMaxValCD/2)<1
                        QLevDec2 = 0.1;
                    else
                        QLevDec2 = QInterCD(round(QSeqFinCD(QLocMaxCD)-QMaxValCD/2));
                    end
                end
                
                [MaxValEF,LocMaxEF]=max(SeqOnesEF);
                if LocMaxEF<2 || MaxValEF<2
                    LevDec1 = 0.1;
                    LocMaxEF = 1;
                    SeqFinEF(1)=2;
                    MaxValEF = 0;
                    InterEF(1)=LevDec1;
                else
                    if (SeqFinEF(LocMaxEF)-MaxValEF/2)<1
                        LevDec1 = 0.1;
                    else
                        LevDec1 = InterEF(round(SeqFinEF(LocMaxEF)-MaxValEF/2));
                    end
                    
                end
                [QMaxValEF,QLocMaxEF]=max(QSeqOnesEF);
                if QLocMaxEF<2 || QMaxValEF<2
                    QLevDec1 = 0.1;
                    QLocMaxEF = 1;
                    QSeqFinEF(1)=2;
                    QMaxValEF = 0;
                    QInterEF(1)=QLevDec1;
                else
                    if (QSeqFinEF(QLocMaxEF)-QMaxValEF/2)<1
                        QLevDec1 = 0.1;
                    else
                        QLevDec1 = QInterEF(round(QSeqFinEF(QLocMaxEF)-QMaxValEF/2));
                    end
                    
                end
                
                %another way to measure the eye opening is the get all the
                %boxes and find all peaks on it, that will be a plato created
                %by the sequences of ones (which was zeros). From thos peaks,
                %the eye diagram will be the longer of them hence it will take
                %the most part of the vector that store the findpeaks result.
                %Thus, the middle of the eye will be basically the middle of
                %the peaks vector.
                LocAB = find(EyeAB);
                if isempty(LocAB)                                              %if for some reason there are no peaks, something went wrong.
                    LocAB = 1;
                    LevelDec3 = 0.1;%mean(Levels(3:4));                       %by default the decission level will be set to 0.65
                else
                    LevelDec3 = LevDec3;
                end
                QLocAB = find(QEyeAB);
                if isempty(QLocAB)                                              %if for some reason there are no peaks, something went wrong.
                    QLocAB = 1;
                    QLevelDec3 = 0.1;%mean(Levels(3:4));                       %by default the decission level will be set to 0.65
                else
                    QLevelDec3 = QLevDec3;
                end
                
                LocCD = find(EyeCD);
                if isempty(LocCD)
                    LocCD = 1;
                    LevelDec2 = 0.1;%mean(Levels(2:3));
                else
                    LevelDec2 = LevDec2;
                end
                QLocCD = find(QEyeCD);
                if isempty(QLocCD)
                    QLocCD = 1;
                    QLevelDec2 = 0.1;%mean(Levels(2:3));
                else
                    QLevelDec2 = QLevDec2;
                end
                
                LocEF = find(EyeEF);
                if isempty(LocEF)
                    LocEF = 1;
                    LevelDec1 = 0.1;%mean(Levels(1:2));
                else
                    LevelDec1 = LevDec1;
                end
                QLocEF = find(QEyeEF);
                if isempty(QLocEF)
                    QLocEF = 1;
                    QLevelDec1 = 0.1;%mean(Levels(1:2));
                else
                    QLevelDec1 = QLevDec1;
                end
%##########################################################################
                if CalcS
                    [~,~,EyeOpenI,EyeOpenHighI,EyeOpenLowI] = Olho(EoutI,T,NPPB,0);
                    [~,~,EyeOpenQ,EyeOpenHighQ,EyeOpenLowQ] = Olho(EoutQ,T,NPPB,0);
                    AberLevIS(CurrentTest,ThisCarr)  = EyeOpenI;
                    ValsLevIS(CurrentTest,ThisCarr)  = EyeOpenLowI + EyeOpenI/2;
                    AberLevQS(CurrentTest,ThisCarr)  = EyeOpenQ;
                    ValsLevQS(CurrentTest,ThisCarr)  = EyeOpenLowQ + EyeOpenQ/2;
                end
                
                AberLevAuxI(1)= InterAB(SeqFinAB(LocMaxAB))-InterAB(round(SeqFinAB(LocMaxAB)-MaxValAB));
                AberLevAuxI(3)= InterCD(SeqFinCD(LocMaxCD))-InterCD(round(SeqFinCD(LocMaxCD)-MaxValCD));
                AberLevAuxI(2)= InterEF(SeqFinEF(LocMaxEF))-InterEF(round(SeqFinEF(LocMaxEF)-MaxValEF));
                ValsLevAuxI(1)= LevDec3;
                ValsLevAuxI(3)= LevDec2;
                ValsLevAuxI(2)= LevDec1;
                AberLevAuxQ(1)= QInterAB(QSeqFinAB(QLocMaxAB))-QInterAB(round(QSeqFinAB(QLocMaxAB)-QMaxValAB));
                AberLevAuxQ(3)= QInterCD(QSeqFinCD(QLocMaxCD))-QInterCD(round(QSeqFinCD(QLocMaxCD)-QMaxValCD));
                AberLevAuxQ(2)= QInterEF(QSeqFinEF(QLocMaxEF))-QInterEF(round(QSeqFinEF(QLocMaxEF)-QMaxValEF));
                ValsLevAuxQ(1)= QLevDec3;
                ValsLevAuxQ(3)= QLevDec2;
                ValsLevAuxQ(2)= QLevDec1;

                %% Ploting the result for qualitative analizes
                
                if PrintinEye
                    PrintInfo(PlotingThis*24,EoutI,T,NPPB);
                    hold on;
                    plot(t((NPPB)/2),InterAB(round(SeqFinAB(LocMaxAB)-MaxValAB)),'kx');
                    plot(t((NPPB)/2),InterAB(SeqFinAB(LocMaxAB)),'ko');
                    plot(t((NPPB)/2),InterAB(round(SeqFinAB(LocMaxAB)-MaxValAB/2)),'kd');
                    plot(t((NPPB)/4),InterCD(round(SeqFinCD(LocMaxCD)-MaxValCD)),'gx');
                    plot(t((NPPB)/4),InterCD(SeqFinCD(LocMaxCD)),'go');
                    plot(t((NPPB)/4),InterCD(round(SeqFinCD(LocMaxCD)-MaxValCD/2)),'gd');
                    plot(t(3*(NPPB)/4),InterEF(round(SeqFinEF(LocMaxEF)-MaxValEF)),'bx');
                    plot(t(3*(NPPB)/4),InterEF(SeqFinEF(LocMaxEF)),'bo');
                    plot(t(3*(NPPB)/4),InterEF(round(SeqFinEF(LocMaxEF)-MaxValEF/2)),'bd');
                    %PrintInfo(Ploting*25,Ui,T,NPPB);
                    PrintInfo(PlotingThis*26,EoutQ,T,NPPB);
                    hold on;
                    plot(t((NPPB)/2),QInterAB(round(QSeqFinAB(QLocMaxAB)-QMaxValAB)),'kx');
                    plot(t((NPPB)/2),QInterAB(QSeqFinAB(QLocMaxAB)),'ko');
                    plot(t((NPPB)/2),QInterAB(round(QSeqFinAB(QLocMaxAB)-QMaxValAB/2)),'kd');
                    plot(t((NPPB)/4),QInterCD(round(QSeqFinCD(QLocMaxCD)-QMaxValCD)),'gx');
                    plot(t((NPPB)/4),QInterCD(QSeqFinCD(QLocMaxCD)),'go');
                    plot(t((NPPB)/4),QInterCD(round(QSeqFinCD(QLocMaxCD)-QMaxValCD/2)),'gd');
                    plot(t(3*(NPPB)/4),QInterEF(round(QSeqFinEF(QLocMaxEF)-QMaxValEF)),'bx');
                    plot(t(3*(NPPB)/4),QInterEF(QSeqFinEF(QLocMaxEF)),'bo');
                    plot(t(3*(NPPB)/4),QInterEF(round(QSeqFinEF(QLocMaxEF)-QMaxValEF/2)),'bd');
                    if CalcS
                        plot(t((NPPB)/2),EyeOpenLowI,'mx');
                        plot(t((NPPB)/2),EyeOpenHighI,'mo');
                        plot(t((NPPB)/2),EyeOpenLowI + EyeOpenI/2,'md');
                        plot(t((NPPB)/2),EyeOpenLowQ,'cx');
                        plot(t((NPPB)/2),EyeOpenHighQ,'co');
                        plot(t((NPPB)/2),EyeOpenLowQ + EyeOpenQ/2,'cd');
                    end
                    a=0;
                end
                %PrintInfo(Ploting*27,Uq,T,NPPB);
                %%

%                 Txaux1 = rectpulse(TxDataMat(ThisCarr,:),NPPB);
%                 Txaux2 = rectpulse(TxDataMat(ThisCarr+NumCarr,:),NPPB);
%                 Txaux1(Txaux1==0) = -1;
%                 Txaux2(Txaux2==0) = -1;
                %%  Ploting some results for qualitative analizes
                %             PrintInfo(PlotingThis*28,t(1:length(EoutI)),Txaux1,Txaux2,EoutI,...
                %                           EoutA,EoutB,EoutQ,EoutC,EoutD,real(Ui),real(Uq));
                %                       a=1;
                %%               Recovering the Information
                %After passing the optical signal to the eletrical domain, for
                %actually detect the data withing the signal the following
                %steps are needed.
                %
                %Finding Decission Levels:
                %The process for decoding the income signal will be based on
                %eletronic comparators. Inasmuch as the right decission level
                %must be acquired for accurately decide, within a symbol
                %periode, what that current leavel means (ones or zeros).
                %
                %The process hereafter of chosing the  decission levels is not
                %deterministic rather it is a statistic process. The main idea
                %is to take the decission level from the histogram generated
                %from the income signal stream.
                %
                %This process is realized inside the function Olho.
                %
                %Basicaly the decission level will be the minimal value of the
                %currente eye under evaluation plus the half of the its eye
                %opening.The following ilustration better describe this process
                %
                %Eye Limits:   + (1/2 * Eye Opening:)  =    Comparisson Limit:
                %
                %UperLevel  ______     ______________     _____
                %                 \   /      |       \   /
                %                  \ /       |        \ /
                %                   \ Half Eye Opening /     Decission Level
                %                  / \       |        / \
                %LowerLevel ______/   \______|_______/   \_____
                %
                %
                %Actualy Receiving Data:
                %Once the signal was processed the next step is through a
                %comparator decide the actual information received.
                %
                %
                %For accurately read the income signal is necessary more than
                %just have the decision level, it is also needed to know where
                %it occurs. For instance, if we decide to use just to take the
                %decision level and measure the mean value across and a portion
                %of the time period, which portion should we take? It would be
                %logical to take the central portion, but the bit may be
                %deformed in such way that the information gets concentrated
                %out of the middle part, which means the symbol is not
                %symmetric. The symmetry information can be acquired from the
                %eye diagram by measuring the longitudinal opening. The
                %following sketch better describes this process:
                %
                %                   Point of Symmetry
                %         _____     _______|_______     _____
                %              \   /               \   /
                %               \ /  Longitudinal   \ /
                %                \ ----------------- /
                %               / \    Opening      / \
                %         _____/   \_______________/   \_____
                %
                %With those two pieces of information, decision level and point
                %of symmetry, we have the X, Y coordinates for the centre of
                %the Eye Diagram. Therefore, as long as there is an opening on
                %it it will be possible to recover the transmitted information
                %without error... theoretically.
                %
                %As this process is also statistical, first we reshape the
                %income vector to analyze all periods at the same time.
%                 EyeSymMatI = reshape(EoutI(1+SyncPeriod*NPPB:end-SyncPeriod*...
%                     NPPB),NPPB,NbDQPSK/2-2*SyncPeriod);
%                 %Then we take the values that compose the decision level
%                 %because they will mark the point of symmetry.
%                 %
%                 %Firstly it was set the interval in which the histogram will be
%                 %build. It is based on the number of samples per bit period.
%                 IntervalI = linspace(min(EoutI(1+SyncPeriod*NPPB:end-...
%                     SyncPeriod*NPPB)),max(EoutI(1+SyncPeriod*NPPB:end-...
%                     SyncPeriod*NPPB)),2*NPPB);
%                 %Therefore, the MATLAB hist function returns the number of
%                 %occurrence of each interval.
%                 EyeMaxI = hist(EoutI,IntervalI);
%                 EyeMaxauxI = [0 EyeMaxI 0];                                    %Zeros are added at the EyeMax to auxiliate the finding peaks process
%                 [~,EyeLocI] = findpeaks(EyeMaxauxI,'MinPeakDistance',NPPB/4,...
%                     'SortStr','descend','NPeaks',2);%The peaks on the Eye profile will be the levels at the Eyes limit
%                 %From the location of the max values that occured, which means
%                 %the uper and lower level of the eye diagram it needs to take
%                 %the actual value that those occurences represent that is
%                 %withing the the Interval variable.
%                 ValToSeekI = IntervalI(EyeLocI-1);
%                 ValToSeekI = sort(ValToSeekI,'ascend');
%                 %The number of ocurrences is a statical measure therefore one
%                 %does not have control which interval will have the highest
%                 %peak, thus it is important to ordenate the values to be seek
%                 %from the lower part of the eye diagram to the uper part of the
%                 %eye diagram.
%                 OccuCountI = zeros(1,size(EyeSymMatI,1));                      %Auxiliar Variable for accounting.
%                 for kk=1:size(EyeSymMatI,1)                                    %For every sample within a symbol period
% 
%                     OccuCountI(kk) = OccuCountI(kk)+sum((EyeSymMatI(kk,:)>=...
%                         min(EoutI))&(EyeSymMatI(kk,:)<=UpeSymPer*EyeOpenLowI)); %Account all occurencies of the valeu 1
%                     OccuCountI(kk) = OccuCountI(kk)+sum((EyeSymMatI(kk,:)>=...
%                         LowSymPer*EyeOpenHighI)&(EyeSymMatI(kk,:)<=max(EoutI))); %Account all occurencies of the valeu 2
%                 end
                %The point of symmetry of the eye diagram will be where the
                %maximum number of occurrences were measured inasmuch as those
                %are the points where all the bits go to the center of the
                %symbol.  From the maximum number of occurrences, it can happen
                %for more than one sample within one symbol period, in the best
                %case, all samples would have the same accounting as it is
                %shown the ilustration above hence the symmetry will be at the
                %middle sample of this group of maximum occurrences. This value
                %can be found by the mean of the samples positions within the
                %symbol period. The problem with this approach is that the
                %signal must be synchronized with the maximum displacement of
                %a symbol period minus 25% of the eye Longitudinal opening if
                %the displacement is higher than that the point of symmetry
                %will be wrongly measured.
%                 [SymLocI] = round(mean(find(ismember(OccuCountI,max(...
%                     OccuCountI)))));
                %             [~,SymLocI] = findpeaks(OccuCountI,'SortStr','descend');       %The peak on the Eye profile will be the Symmetry level
                DataOdd = [];                                                  %Initialization of the vector that will store the income data
                DataOddM = [];                                                  %Initialization of the vector that will store the income data
                DataOddL = [];                                                  %Initialization of the vector that will store the income data
                DataOddU = [];                                                  %Initialization of the vector that will store the income data
                DataOddU2 = [];                                                  %Initialization of the vector that will store the income data
                DataOddU3 = [];                                                  %Initialization of the vector that will store the income data
                DataOddS = [];                                                  %Initialization of the vector that will store the income data
                for kk=1:NPPB:length(EoutI)                                    %The comparison process will be made for each symbol period
                    %                 MeanOfData = mean(EoutI((kk-1)+SymLocI(1)));
                    MeanOfData = mean(EoutI((kk-1)+NPPB/2));
                    if MeanOfData > LevelDec3%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                        DataOdd = [DataOdd 1];                                 %is 1
                    else                                                       %If it is the lowest level the incoming data
                        DataOdd = [DataOdd 0];                                 %is 0
                    end
                    MeanOfData = mean(EoutI((kk-1)+3*NPPB/4));
                    if MeanOfData > LevelDec1%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                        DataOddM = [DataOddM 1];                                 %is 1
                    else                                                       %If it is the lowest level the incoming data
                        DataOddM = [DataOddM 0];                                 %is 0
                    end
                    MeanOfData = mean(EoutI((kk-1)+NPPB/4));
                    if MeanOfData > LevelDec2%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                        DataOddL = [DataOddL 1];                                 %is 1
                    else                                                       %If it is the lowest level the incoming data
                        DataOddL = [DataOddL 0];                                 %is 0
                    end
                    MeanOfData = mean(EoutI((kk-1)+(NPPB/2)));
                    if MeanOfData > LevDefDpqsk%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                        DataOddU = [DataOddU 1];                                 %is 1
                    else                                                       %If it is the lowest level the incoming data
                        DataOddU = [DataOddU 0];                                 %is 0
                    end
                    MeanOfData = mean(EoutI((kk-1)+(NPPB/4)));
                    if MeanOfData > LevDefDpqsk%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                        DataOddU2 = [DataOddU2 1];                                 %is 1
                    else                                                       %If it is the lowest level the incoming data
                        DataOddU2 = [DataOddU2 0];                                 %is 0
                    end
                    MeanOfData = mean(EoutI((kk-1)+(3*NPPB/4)));
                    if MeanOfData > LevDefDpqsk%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                        DataOddU3 = [DataOddU3 1];                                 %is 1
                    else                                                       %If it is the lowest level the incoming data
                        DataOddU3 = [DataOddU3 0];                                 %is 0
                    end
                    if CalcS
                        MeanOfData = mean(EoutI((kk-1)+NPPB/2));
                        if MeanOfData > EyeOpenLowI + EyeOpenI/2%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                            DataOddS = [DataOddS 1];                                 %is 1
                        else                                                       %If it is the lowest level the incoming data
                            DataOddS = [DataOddS 0];                                 %is 0
                        end
                    end
                end
                %The identical process just described above will also be used
                %to recover the data at the even positions.
                %
                %As this process is also statistical, first we reshape the
                %income vector to analyze all periods at the same time.
%                 EyeSymMatQ = reshape(EoutQ(1+SyncPeriod*NPPB:end-SyncPeriod*...
%                     NPPB),NPPB,NbDQPSK/2-2*SyncPeriod);
%                 %Then we take the values that compose the decision level
%                 %because they will mark the point of symmetry.
%                 %
%                 %Firstly it was set the interval in which the histogram will be
%                 %build. It is based on the number of samples per bit period.
%                 IntervalQ = linspace(min(EoutQ(1+SyncPeriod*NPPB:end-...
%                     SyncPeriod*NPPB)),max(EoutQ(1+SyncPeriod*NPPB:end-...
%                     SyncPeriod*NPPB)),2*NPPB);
%                 %Therefore, the MATLAB hist function returns the number of
%                 %occurrence of each interval.
%                 EyeMaxQ = hist(EoutQ,IntervalQ);
%                 EyeMaxauxQ = [0 EyeMaxI 0];                                    %Zeros are added at the EyeMax to auxiliate the finding peaks process
%                 [~,EyeLocQ] = findpeaks(EyeMaxauxQ,'MinPeakDistance',NPPB/4,...
%                     'SortStr','descend','NPeaks',2);%The peaks on the Eye profile will be the levels at the Eyes limit
%                 ValToSeekQ = IntervalQ(EyeLocQ-1);
%                 ValToSeekQ = sort(ValToSeekQ,'ascend');
%                 OccuCountQ = zeros(1,size(EyeSymMatQ,1));                      %Auxiliar Variable for accounting.
%                 for kk=1:size(EyeSymMatQ,1)                                    %For every sample within a symbol period
%                     OccuCountQ(kk) = OccuCountQ(kk)+sum((EyeSymMatI(kk,:)>=...
%                         LowSymPer*ValToSeekQ(1))&(EyeSymMatQ(kk,:)<=UpeSymPer*...
%                         ValToSeekQ(1)));%Account all occurencies of the valeu 1
%                     OccuCountQ(kk) = OccuCountQ(kk)+sum((EyeSymMatQ(kk,:)>=...
%                         LowSymPer*ValToSeekQ(2))&(EyeSymMatQ(kk,:)<=UpeSymPer*...
%                         ValToSeekQ(2)));%Account all occurencies of the valeu 2
%                 end
%                 %             [~,SymLocQ] = findpeaks(OccuCountQ,'SortStr','descend');       %The peak on the Eye profile will be the Symmetry level
%                 [SymLocQ] = round(mean(find(ismember(OccuCountQ,max(...
%                     OccuCountQ)))));
                %##############################################################
                %######################Important###############################
                %The ber results for the Data Even is not as good as the
                %results of the Data Odd. One possible reason is the decision
                %point of symmetry hence, for testing, we change SymLocQ to
                %SymLocI for evaluation of improvement. It is expected as bouth
                %signal will the same point of symetry will perform in an equal
                %way. If it is confirmed the creation of SymLocQ will be
                %erased.
                %##############################################################
                %##############################################################
                DataEven = [];                                                 %Initialization of the vector that will store the income data
                DataEvenM = [];                                                 %Initialization of the vector that will store the income data
                DataEvenL = [];                                                 %Initialization of the vector that will store the income data
                DataEvenU = [];                                                 %Initialization of the vector that will store the income data
                DataEvenU2 = [];                                                 %Initialization of the vector that will store the income data
                DataEvenU3 = [];                                                 %Initialization of the vector that will store the income data
                DataEvenS = [];                                                 %Initialization of the vector that will store the income data
                for kk=1:NPPB:length(EoutQ)                                    %The comparison process will be made for each symbol period
                    %                 MeanOfData = mean(EoutQ((kk-1)+SymLocI(1)));
                    MeanOfData = mean(EoutQ((kk-1)+NPPB/2));
                    if MeanOfData > QLevelDec3%EyeOpenLowQ+EyeOpenQ/2                     %If it is the uper level the incoming data
                        DataEven = [DataEven 1];                               %is 1
                    else                                                       %If it is the lowest level the incoming data
                        DataEven = [DataEven 0];                               %is 0
                    end
                    MeanOfData = mean(EoutQ((kk-1)+3*NPPB/4));
                    if MeanOfData > QLevelDec1%EyeOpenLowQ+EyeOpenQ/2                     %If it is the uper level the incoming data
                        DataEvenM = [DataEvenM 1];                               %is 1
                    else                                                       %If it is the lowest level the incoming data
                        DataEvenM = [DataEvenM 0];                               %is 0
                    end
                    MeanOfData = mean(EoutQ((kk-1)+NPPB/4));
                    if MeanOfData > QLevelDec3%EyeOpenLowQ+EyeOpenQ/2                     %If it is the uper level the incoming data
                        DataEvenL = [DataEvenL 1];                               %is 1
                    else                                                       %If it is the lowest level the incoming data
                        DataEvenL = [DataEvenL 0];                               %is 0
                    end
                    MeanOfData = mean(EoutQ((kk-1)+(NPPB/2)));
                    if MeanOfData > LevDefDpqsk%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                        DataEvenU = [DataEvenU 1];                                 %is 1
                    else                                                       %If it is the lowest level the incoming data
                        DataEvenU = [DataEvenU 0];                                 %is 0
                    end
                    MeanOfData = mean(EoutQ((kk-1)+(NPPB/4)));
                    if MeanOfData > LevDefDpqsk%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                        DataEvenU2 = [DataEvenU2 1];                                 %is 1
                    else                                                       %If it is the lowest level the incoming data
                        DataEvenU2 = [DataEvenU2 0];                                 %is 0
                    end
                    MeanOfData = mean(EoutQ((kk-1)+(3*NPPB/4)));
                    if MeanOfData > LevDefDpqsk%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                        DataEvenU3 = [DataEvenU3 1];                                 %is 1
                    else                                                       %If it is the lowest level the incoming data
                        DataEvenU3 = [DataEvenU3 0];                                 %is 0
                    end
                    if CalcS
                        MeanOfData = mean(EoutI((kk-1)+NPPB/2));
                        if MeanOfData > EyeOpenLowQ + EyeOpenQ/2%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                            DataEvenS = [DataEvenS 1];                                 %is 1
                        else                                                       %If it is the lowest level the incoming data
                            DataEvenS = [DataEvenS 0];                                 %is 0
                        end
                    end
                end
                %%       Calculating the Bit Error Ratio (BER)
                %The final process here is to count the number of wrongdoings
                %of this whole process upon the transmited data for
                %quantitative analizes
                TxDataOdd  = TxDataMat(ThisCarr,1+SyncPeriod:end-SyncPeriod);
                TxDataEven = TxDataMat(ThisCarr+NumCarr,1+SyncPeriod:end-...
                    SyncPeriod);
                if CalcS
                    DataOddS  = DataOddS(1+SyncPeriod:end-SyncPeriod);
                    DataEvenS = DataEvenS(1+SyncPeriod:end-SyncPeriod);
                    BitErrOddS    = sum(xor(TxDataOdd,DataOddS));                      %Comparison between the Transmited and received and counting the differences
                    BitErrEvenS   = sum(xor(TxDataEven,DataEvenS));                    %Comparison between the Transmited and received and counting the differences
                    BerDQPSKS(CurrentTest,ThisCarr) = (BitErrOddS+BitErrEvenS)/...
                        ((NbDQPSK)-(2*SyncPeriod));%Calculating the ration of wrong bits and the total number of bits transmited
                end
                DataOdd  = DataOdd(1+SyncPeriod:end-SyncPeriod);
                DataEven = DataEven(1+SyncPeriod:end-SyncPeriod);
                DataOddM  = DataOddM(1+SyncPeriod:end-SyncPeriod);
                DataEvenM = DataEvenM(1+SyncPeriod:end-SyncPeriod);
                DataOddL  = DataOddL(1+SyncPeriod:end-SyncPeriod);
                DataEvenL = DataEvenL(1+SyncPeriod:end-SyncPeriod);
                DataOddU  = DataOddU(1+SyncPeriod:end-SyncPeriod);
                DataEvenU = DataEvenU(1+SyncPeriod:end-SyncPeriod);
                DataOddU2  = DataOddU2(1+SyncPeriod:end-SyncPeriod);
                DataEvenU2 = DataEvenU2(1+SyncPeriod:end-SyncPeriod);
                DataOddU3  = DataOddU3(1+SyncPeriod:end-SyncPeriod);
                DataEvenU3 = DataEvenU3(1+SyncPeriod:end-SyncPeriod);
                AberLevAuxI(4) = 0;
                ValsLevAuxI(4) = 0.01;
                AberLevAuxQ(4) = 0;
                ValsLevAuxQ(4) = 0.01;
                AberLevAuxI(5) = 0;
                ValsLevAuxI(5) = 0.01;
                AberLevAuxQ(5) = 0;
                ValsLevAuxQ(5) = 0.01;
                AberLevAuxI(6) = 0;
                ValsLevAuxI(6) = 0.01;
                AberLevAuxQ(6) = 0;
                ValsLevAuxQ(6) = 0.01;

                BitErrOdd(1)  = sum(xor(TxDataOdd,DataOdd));                      %Comparison between the Transmited and received and counting the differences
                BitErrEven(1) = sum(xor(TxDataEven,DataEven));                    %Comparison between the Transmited and received and counting the differences
                BitErrOdd(2)  = sum(xor(TxDataOdd,DataOddM));                      %Comparison between the Transmited and received and counting the differences
                BitErrEven(2) = sum(xor(TxDataEven,DataEvenM));                    %Comparison between the Transmited and received and counting the differences
                BitErrOdd(3)  = sum(xor(TxDataOdd,DataOddL));                      %Comparison between the Transmited and received and counting the differences
                BitErrEven(3) = sum(xor(TxDataEven,DataEvenL));                    %Comparison between the Transmited and received and counting the differences
                BitErrOdd(4)  = sum(xor(TxDataOdd,DataOddU));                      %Comparison between the Transmited and received and counting the differences
                BitErrEven(4) = sum(xor(TxDataEven,DataEvenU));                    %Comparison between the Transmited and received and counting the differences
                BitErrOdd(5)  = sum(xor(TxDataOdd,DataOddU2));                      %Comparison between the Transmited and received and counting the differences
                BitErrEven(5) = sum(xor(TxDataEven,DataEvenU2));                    %Comparison between the Transmited and received and counting the differences
                BitErrOdd(6)  = sum(xor(TxDataOdd,DataOddU3));                      %Comparison between the Transmited and received and counting the differences
                BitErrEven(6) = sum(xor(TxDataEven,DataEvenU3));                    %Comparison between the Transmited and received and counting the differences
                BerDQPSK(CurrentTest,ThisCarr) = (min(BitErrOdd)+min(BitErrEven))/...
                    ((NbDQPSK)-(2*SyncPeriod));%Calculating the ration of wrong bits and the total number of bits transmited
                
                AberIAux = max(AberLevAuxI(BitErrOdd<=min(BitErrOdd)));
                AberQAux = max(AberLevAuxQ(BitErrEven<=min(BitErrEven)));
                AberLevI(CurrentTest,ThisCarr)  = AberIAux;
                ValsLevI(CurrentTest,ThisCarr)  = max(ValsLevAuxI(AberLevAuxI==AberIAux));

                AberLevQ(CurrentTest,ThisCarr)  = AberQAux;
                ValsLevQ(CurrentTest,ThisCarr)  = max(ValsLevAuxQ(AberLevAuxQ==AberQAux));
                
                [~,LocAux] = max(AberLevAuxI==AberIAux);
                EvmDB(CurrentTest,ThisCarr)   = EvmDBA(LocAux);
                EvmPer(CurrentTest,ThisCarr)  = EvmPerA(LocAux);
                EvmRms(CurrentTest,ThisCarr)  = EvmRmsA(LocAux);
                EvmDBJ(CurrentTest,ThisCarr)  = EvmDBJA(LocAux);
                EvmPerJ(CurrentTest,ThisCarr) = EvmPerJA(LocAux);
                EvmRmsJ(CurrentTest,ThisCarr) = EvmRmsJA(LocAux);
                
                EvmMatRec(ObsCarrPos==ThisCarr,:) = EvmMatRecA(LocAux,:);                       %Taking just the middle samples as references
                %% Ploting the result for qualitative analizes
%                 PrintInfo(PlotingThis*29,TxDataOdd,DataOdd);
%                 PrintInfo(PlotingThis*30,TxDataEven,DataEven);
                %%
                %             berpos = 1:2:size(BerDQPSK,2);
                %             BerDQPSK(size(BerDQPSK,1),berpos)
                %             a=6;
                close all;
%                 BerDQPSK(CurrentTest,ThisCarr+1) = 1;
            end
        end
    case '4PAM'
        
        %%              Receiver 4PAM
        for ThisCarr=InitCarr:CarrPass:NumCarr                                           %For each carrier the same process of reception will be used.
            if CarrUsedDo(ThisCarr)
                %%  Reception Process: Handling the income optical field
                %At the first step the income field will be selected from the
                %optical FFT Output with the help of the maping vector
                %(previously described).
                EoutAux = EoutAux1(VetThisCarr==ObsCarrUsed(ThisCarr),:);
                if AddingNoiseP
    %                 figure;plot(f,20*log10(abs(fftshift(fft(EoutAux)./length(EoutAux)))));
    %                 [~,sigpower] = MeasPower(EoutAux);
                    [p,l] = findpeaks(20*log10(abs(fftshift(fft(EoutAux)./length(EoutAux)))),'SortStr','descend','NPeaks',1);
                    EoutAuxA = fftshift(fft(EoutAux));
                    EoutAuxA(f==f(l)) = EoutAuxA(f==f(l-3));
    %                 [testfilt,~] = FiltroGaussiano(f,12.5e9,f(l),5);           %Creating filter for selection of the received signal
    %                 testfilt = fftshift(testfilt);                                   %Shifting the filter for matching the received signal
    %                 EoutAux = ifft(fftshift(EoutAuxA.*testfilt));
                    EoutAuxA = ifft(fftshift(EoutAuxA));
    %                 figure;plot(f,20*log10(abs(fftshift(fft(EoutAuxA)./length(EoutAuxA)))));
    %                 [~,sigpower] = MeasPower(EoutAuxA);
                    [~,sigpower] = MeasPower(EoutAuxA);
                    EoutAux      = awgn(EoutAux,osnrp,sigpower-30);
                    [~,sigpower] = MeasPower(EoutAux);
                end
                %%            Fiber Time Delay Compensation
                switch Medium
                    case 'Fiber'
                        EoutAux = ifft(fft(EoutAux).*exp(1j*2*pi*f*(...
                            FiberDelay(ThisCarr)*Ta)));
                    otherwise
                end
                %The current incoming signal is them converted from the optical
                %domain to the eletrical domain with the help of an photo
                %detector.
                Ix = EoutAux.*conj(EoutAux);
                %The process with the photo diode is self-coherent, which means
                %the rusult will be a component of the signal centered in f=0
                %and another component centered at f=2*fc (frenquecy central).
                %Therefore, to remove the higher order component a low pass
                %filter will be used.
                %             switch Medium
                %                 case 'Fiber'
                %%           Creating the Reception Filter
                %             taux = t(1:length(Ix));
                %             faux = time2freq(taux);
                [BitFilt,~] = FiltroGaussiano(f,BWD2,CenFeq2,FiltOrd2);        %Creating filter for selection of the received signal
                BitFilt = fftshift(BitFilt);                                   %Shifting the filter for matching the received signal
                if RecFilBanPas
                    Ix = ifft(fft(Ix).*BitFilt);
                    Ix(1:4*(2*NumAmosCP+NPPB))         = Ix(6*(2*NumAmosCP+NPPB));
                    Ix(end+1-4*(2*NumAmosCP+NPPB):end) = Ix(6*(2*NumAmosCP+NPPB));
                    Ix = Ix - min(Ix);                                         %Removing the DC component from them Eletrical signal received

                    Ix = Ix./max(abs(Ix));                                     %Normalizing the eletrical signal (amplifying what is needed)
                end
                %                 otherwise
                %             end
                
                
                %%         Adding Noise to Received Signal
                %Just after the photo diod is where the noise must be added as
                %it is the main constrain of the system. Once the signal is
                %converted from optical to electrical and it digital values
                %were recevered, there is no point to add noise, because the
                %information was already received, it just needs decoding.
                %The noise here added will be gaussian basically it represents
                %the reception sensibility. In another words, how much the
                %signal must be about the noise for an acceptable recovering.
                if ReceptorNoise                                               %Verify whether the noise is to be added or not
                    [~,SigPower] = MeasPower(Ix);                              %it own received signal or
                    SigPower2 = SigPower-30;%10*log10(SigPower);
                    %SNR = CarSNR - 10*log10(2) + 10*log10(Nsamp) + 10*log10(1/2);
                    SNR2 = CarSNR + 10*log10(2) - 10*log10(Nsamp) - 10*log10(1/2) + 10*log10(10^0.6);
                    SNR = CarSNR + 10*log10(2) - 10*log10((t(NPPB))/(6*(t(2*NumAmosCP+NPPB))/(NPPB)));
                                    
                    Ix = awgn(Ix,SNR,SigPower2);
                end
                
                if ~RecFilBanPas
                    
%                     [BitFilt,~] = FiltroGaussiano(f,fc,CenFeq2,FiltOrd2);        %Creating filter for selection of the received signal
%                     BitFilt = fftshift(BitFilt);                                   %Shifting the filter for matching the received signal
                    Ix = ifft(fft(Ix).*BitFilt);
                    Ix(1:4*(2*NumAmosCP+NPPB))         = Ix(6*(2*NumAmosCP+NPPB));
                    Ix(end+1-4*(2*NumAmosCP+NPPB):end) = Ix(6*(2*NumAmosCP+NPPB));
                    Ix = Ix - min(Ix);                                         %Removing the DC component from them Eletrical signal received

                    Ix = Ix./max(abs(Ix));                                     %Normalizing the eletrical signal (amplifying what is needed)
                end
                
                VetElecPower(CurrentTest,ThisCarr)= MeasPower(Ix);
                EoutAuxF = 20*log10(abs(fftshift(fft(EoutAux)./length(...
                    EoutAux))));
                [VetOptiPower(CurrentTest,ThisCarr),~]= findpeaks(EoutAuxF,...
                    'SortStr','descend','NPeaks',1);
                %% Ploting the result for qualitative analizes
                %             PrintInfo(Ploting*31,f,20.*log10(abs(fftshift(fft(Ix)./...
                %                                                             length(Ix)))));
                %             PrintInfo(PlotingThis*32,t,Ix);
                %             PrintInfo(PlotingThis*33,Ix,T,NPPB);
                %             PrintInfo(Ploting*34,t(IniSyncPos:SyncPos),Ix(IniSyncPos:...
                %                                     SyncPos),SyncSymb(IniSyncPos:SyncPos));
                %%        synchronizing
                %For the reception process work properly, it is needed to
                %sincronized the recieved signal or the sampling process will
                %not work.
                AuxSync = (Ix(IniSyncPos:SyncPos));                            %Selecting the sync-word within the received signal
                AuxSync1 = AuxSync;
                %This synchronization process is based on the mean expected
                %value. That means, the information mean value should be within
                %one period of symbol. Thus, the mean value of the received
                %signal is acquired and compare of the known sync-word to
                %verify if this mean value is at the right possition
                AuxSync1(AuxSync1<0)               = 0;                        %To keep the mean value above zero anything under is neglected
                AuxSync1(AuxSync1>=mean(AuxSync1)) = 1;                        %Adding a flag to the first sample of the received mean value
                AuxSync1(AuxSync1<mean(AuxSync1))  = -1;                       %All the others samples at set to the lowest level
                AuxSync2                           = SyncSymb(IniSyncPos:...
                    SyncPos);%Selecting the sync-word within the known signal
                AuxSync2(AuxSync2>=mean(AuxSync2)) = 1;                        %Adding a flag to the first sample of the known mean value
                AuxSync2(AuxSync2<mean(AuxSync2))  = -1;                       %All the others samples at set to the lowest level
                
                PosToSyn  = find(ismember(AuxSync1,1));                        %Finding where is the location of the first sample to synchronize
                PosSyn = find(ismember(AuxSync2,1));                           %Finding where is the location of the first sample to synchronize
                
                AuxSyncCorr = PosToSyn(round(end/2))-PosSyn(round(end/2));
                
                %The difference between the PossitionTosynchronize and
                %Possitionsynchronized will be used to correct the time
                %shifting on the transmition and reception process.
                if AuxSyncCorr>=0%If the difference is positive, left-shift...
                    %                 Ix = ifft(fft(Ix).*exp(1j*2*pi*f*AuxSyncCorr));          %Shift based on time change
                    Ix = [Ix(AuxSyncCorr+1:end) Ix(1:AuxSyncCorr)];            %Shift based on sampling sliding
                else%... but if the difference is negative, right-shift
                    %                 Ix = ifft(fft(Ix).*exp(-1j*2*pi*f*AuxSyncCorr));         %Shift based on time change
                    Ix = [Ix(end+AuxSyncCorr+1:end) Ix(1:end + AuxSyncCorr)];  %Shift based on sampling sliding
                end
                if SencondAdjust
                    %For some reason that we could not understand sometimes the
                    %time (sampling) sliding of the signal is not equal
                    %throught the data stream. Thus, the second part of the
                    %synchronism process will be turn ON or OFF according to
                    %the user's will.
                    AuxSyncEnd     = (Ix(end-SyncPos+1:end-IniSyncPos-1));     %Selecting the sync-word within the received signal
                    SyncSymbEndAux = SyncSymbEnd(end-SyncPos+1:end-...
                        IniSyncPos-1);%Selecting the sync-word within the known signal
                    AuxSyncEnd1 = AuxSyncEnd;
                    AuxSyncEnd1(AuxSyncEnd1<0) = 0;                            %To keep the mean value above zero anything under is neglected
                    AuxSyncEnd1(AuxSyncEnd1>=mean(AuxSyncEnd1)) = 1;           %Adding a flag to the first sample of the received mean value
                    AuxSyncEnd1(AuxSyncEnd1<mean(AuxSyncEnd1)) = -1;           %All the others samples at set to the lowest level
                    AuxSyncEnd2 = SyncSymbEndAux;
                    AuxSyncEnd2(AuxSyncEnd2>=mean(AuxSyncEnd2)) = 1;           %Adding a flag to the first sample of the known mean value
                    AuxSyncEnd2(AuxSyncEnd2<mean(AuxSyncEnd2)) = -1;           %All the others samples at set to the lowest level
                    
                    
                    PosToSynEnd  = find(ismember(AuxSyncEnd1,1));              %Finding where is the location of the first sample to synchronize
                    PosSynEnd = find(ismember(AuxSyncEnd2,1));                 %Finding where is the location of the first sample to synchronize
                    
                    AuxSyncCorrEnd = PosToSynEnd(round(end/2))-PosSynEnd(...
                        round(end/2));
                    
                    %The difference between the PossitionTosynchronize and
                    %Possitionsynchronized will be used to correct the time
                    %shifting on the transmition and reception process.
                    if AuxSyncCorrEnd>=0%If possitive difference, left-shift...
                        %                 Ix = ifft(fft(Ix).*exp(1j*2*pi*f*AuxSyncCorrEnd));   %Shift based on time change
                        Ix = [Ix(AuxSyncCorrEnd+1:end) Ix(1:AuxSyncCorrEnd)];  %Shift based on sampling sliding
                    else%... but if the difference is negative, right-shift
                        %                 Ix = ifft(fft(Ix).*exp(-1j*2*pi*f*AuxSyncCorrEnd));  %Shift based on time change
                        Ix = [Ix(end-AuxSyncCorrEnd+1:end) Ix(1:end - ...
                            AuxSyncCorrEnd)];%Shift based on sampling sliding
                    end
                end
                
                %%          Ploting the result for qualitative analizes
%                             PrintInfo(PlotingThis*35,t(end-SyncPos+1:end-IniSyncPos+1),Ix(...
%                                 end-SyncPos+1:end-IniSyncPos+1),SyncSymb(end-SyncPos+1:...
%                                                                         end-IniSyncPos+1));
                %% Removing CP
                if AddCP
                    IxAux = Ix(1:end - StuffSampels);
                    IxAux = reshape(IxAux,(2*NumAmosCP+NPPB),Nb4Pam/2);
                    IxAux = IxAux(1+NumAmosCP:end-NumAmosCP,:);
                    IxAux = reshape(IxAux,1,NPPB*Nb4Pam/2);
                    Ix    = IxAux;
                end
                %% Taking the sampling the EVM meassurement
                PosAuxEout = NPPB/2:NPPB:length(Ix);%Varriable respossible to take just the samples at the middle of the symbol
                IxAux      = Ix(PosAuxEout);                           %Normalizing the reference
                EvmMatRec(ObsCarrPos==ThisCarr,:) = IxAux;                       %Taking just the middle samples as references
                [EvmDB(CurrentTest,ThisCarr), EvmPer(CurrentTest,ThisCarr), EvmRms(CurrentTest,ThisCarr) ] = EvmCalc( EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux );
                [EvmDBJ(CurrentTest,ThisCarr),EvmPerJ(CurrentTest,ThisCarr),EvmRmsJ(CurrentTest,ThisCarr)] = evm1(4,'pam',EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux);
                a=0;
%                 if ThisCarr==1
%                     figure;hold all;grid on;
%                     plot(IxAux,zeros(1,length(IxAux)),'rx');
%                     set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%                     a=0;
%                 end
%                 %%  Measuring the EVM
%                 IxPosAux = NPPB/2:NPPB:length(Ix);
%                 IxRec    = Ix(IxPosAux);
%                 [EvmRms(CurrentTest,ThisCarr),EvmDB(CurrentTest,ThisCarr),...
%                     EvmPer(CurrentTest,ThisCarr)] = EvmCalc( EvmMatRef(ThisCarr,...
%                     :),IxRec );
                %%          Ploting the result for qualitative analizes
                %%         Finding Decission Levels
                %The process for decoding the income signal will be based on
                %eletronic comparators. Inasmuch as the right decission level
                %must be acquired for accurately decide, within a symbol
                %periode, what that current leavel means (ones or zeros).
                %
                %The process hereafter of chosing the  decission levels is not
                %deterministic rather it is a statistic process. The main idea
                %is to take the decission level from the histogram generated
                %from the income signal stream.
                %
                %Firstly it was set the interval in which the histogram will be
                %build. It is based on the number of samples per bit period.
                %             IxToSeek = reshape(Ix(1+SyncPeriod*NPPB:end-SyncPeriod*NPPB...
                %                                        ),NPPB,(Nb4Pam/2) - (2*SyncPeriod));
                Interval = linspace(min(Ix(1+SyncPeriod:end-SyncPeriod)),max...
                    (Ix(1+SyncPeriod:end-SyncPeriod)),IntervalStep);
                %Therefore, the MATLAB hist function returns the number of
                %occurrence of each interval.
                EyeMax = hist(Ix,Interval);
                EyeMaxaux = [0 EyeMax 0];                                      %Zeros are added at the EyeMax to auxiliate the finding peaks process
                [~,EyeLoc] = findpeaks(EyeMaxaux,'MinPeakDistance',MinDist,...
                    'SortStr','descend','NPeaks',4);%,'MinPeakHeight',mean(EyeMax));%The peaks on the Eye profile will be the levels at the Eyes limit
                if length(EyeLoc)<4
                    [~,EyeLoc] = findpeaks(EyeMaxaux,'MinPeakDistance',...
                        MinDist*0.8,'SortStr','descend','NPeaks',4);%,'MinPeakHeight',mean(EyeMax)/4);%The peaks on the Eye profile will be the levels at the Eyes limit
                end
                %This variable will brings the eye profile of the input signal
                %that will be used to generate the decision level. Basicaly the
                %decission level will be the minimal value of the currente eye
                %under evaluation plus the half of the its eye opening. The
                %following ilustration better describe this process.
                %
                %Eye Limits:   + (1/2 * Eye Opening:)  =    Comparisson Limits:
                %
                %UperLevel  ______     ______________     _____
                %                 \   /      |       \   /
                %                  \ /       |        \ /
                %                   \ Half Eye Opening /     Decission Level 3
                %                  / \       |        / \
                %LowerLevel 3_____/   \______|_______/   \_____
                %                 \   /      |       \   /
                %                  \ /       |        \ /
                %                   \ Half Eye Opening /     Decission Level 2
                %                  / \       |        / \
                %LowerLevel 2_____/   \______|_______/   \_____
                %                 \   /      |       \   /
                %                  \ /       |        \ /
                %                   \ Half Eye Opening /     Decission Level 1
                %                  / \       |        / \
                %LowerLevel 1_____/   \______|_______/   \_____
                %
                %%           Ploting for Qualitative Analizes
                %             PrintInfo(Ploting*37,Interval,EyeMax);
                %%         Finding Decission Levels
                %It is not always possible to totaly recover the signal.
                %Depending on the configuration of the transmition and
                %reception system the eye diagram may be nonexistent. Which
                %means, there will not be a profile to be found therefore the
                %EyeLoc will not return the correct location. Inasmuch as the
                %detection process to works limts will be set accordling to the
                %amplitude of the received signal.
                if length(EyeLoc)<4%If it was not able to find the eye profile.
                    EyeLoc = [2 3 4 5];
                    Levels = [0 0 0.35 0.55 0.85];
                    %                 Levels = sort(Levels);
                    %                 LevelDec1 = mean(Ix) - 2*mean(Ix)/3 ;
                    %                 LevelDec2 = mean(Ix) ;
                    %                 LevelDec3 = mean(Ix) + 2*mean(Ix)/3 ;
                else%Whereas, if there is an profile the decission can be found
                    Levels = [Interval(EyeLoc(1)-1) Interval(EyeLoc(2)-1) ...
                        Interval(EyeLoc(3)-1) Interval(EyeLoc(4)-1)];
                    Levels = sort(Levels);
                    %                 LevelDec1 = Levels(1) + (Levels(2)-Levels(1))/2 ;
                    %                 LevelDec2 = Levels(2) + (Levels(3)-Levels(2))/2 ;
                    %                 LevelDec3 = Levels(3) + (Levels(4)-Levels(3))/2 ;
                end
                %##############################################################
                if CalcS
                    limiar1 = (1/2)*abs(max(Ix)-min(Ix))/3;
                    limiar2 = ((1/2)*abs(max(Ix)-min(Ix))/3) + abs(max(Ix)-min(Ix))/3;
                    limiar3 = ((1/2)*abs(max(Ix)-min(Ix))/3) + 2*(abs(max(Ix)-min(Ix))/3);
                    limiarPos1 = Interval>limiar1;
                    limiarPos2 = (Interval<=limiar1)&(Interval>limiar2);
                    limiarPos3 = (Interval<=limiar2)&(Interval>limiar3);
                    limiarPos4 = Interval<=limiar3;

                    EyeSymMa1 = reshape(Ix(1+SyncPeriod*NPPB:end-SyncPeriod*NPPB...
                                               ),NPPB,(Nb4Pam/2) - (2*SyncPeriod));
                    for kk = 1:size(EyeSymMa1,1)
                        EyeHi1    = find((EyeSymMa1(kk,:)<Levels(4))&(EyeSymMa1(kk,:)>limiar1));
                        EyeHi2    = find((EyeSymMa1(kk,:)<Levels(3))&(EyeSymMa1(kk,:)>limiar2));
                        EyeHi3    = find((EyeSymMa1(kk,:)<Levels(2))&(EyeSymMa1(kk,:)>limiar3));

                        EyeLo1    = find((EyeSymMa1(kk,:)<limiar1)&(EyeSymMa1(kk,:)>Levels(3)));
                        EyeLo2    = find((EyeSymMa1(kk,:)<limiar2)&(EyeSymMa1(kk,:)>Levels(2)));
                        EyeLo3    = find((EyeSymMa1(kk,:)<limiar3)&(EyeSymMa1(kk,:)>Levels(1)));

                        Hi1(kk)   = mean((EyeSymMa1(kk,EyeHi1)));
                        Hi2(kk)   = mean((EyeSymMa1(kk,EyeHi2)));
                        Hi3(kk)   = mean((EyeSymMa1(kk,EyeHi3)));

                        Lo1(kk)   = mean((EyeSymMa1(kk,EyeLo1)));
                        Lo2(kk)   = mean((EyeSymMa1(kk,EyeLo2)));
                        Lo3(kk)   = mean((EyeSymMa1(kk,EyeLo3)));

                        LevHi1(kk)= mean(Hi1) - std(EyeSymMa1(kk,EyeHi1));
                        LevHi2(kk)= mean(Hi2) - std(EyeSymMa1(kk,EyeHi2));
                        LevHi3(kk)= mean(Hi3) - std(EyeSymMa1(kk,EyeHi3));

                        LevLo1(kk)= mean(Lo1) + std(EyeSymMa1(kk,EyeLo1));
                        LevLo2(kk)= mean(Lo2) + std(EyeSymMa1(kk,EyeLo2));
                        LevLo3(kk)= mean(Lo3) + std(EyeSymMa1(kk,EyeLo3));

                        EyeAb1(kk)= LevHi1(kk) - LevLo1(kk);
                        EyeAb2(kk)= LevHi2(kk) - LevLo2(kk);
                        EyeAb3(kk)= LevHi3(kk) - LevLo3(kk);
                    end
                    EyeAbertura1 = mean(EyeAb1);
                    EyeAbertura2 = mean(EyeAb2);
                    EyeAbertura3 = mean(EyeAb3);
                    ThisDeciLev1 = mean(LevLo1) + EyeAbertura1/2;
                    ThisDeciLev2 = mean(LevLo2) + EyeAbertura2/2;
                    ThisDeciLev3 = mean(LevLo3) + EyeAbertura3/2;
                    AberLevS(1,ThisCarr,CurrentTest)  = EyeAbertura3;
                    AberLevS(2,ThisCarr,CurrentTest)  = EyeAbertura2;
                    AberLevS(3,ThisCarr,CurrentTest)  = EyeAbertura1;
                    ValsLevS(1,ThisCarr,CurrentTest)  = ThisDeciLev3;
                    ValsLevS(2,ThisCarr,CurrentTest)  = ThisDeciLev2;
                    ValsLevS(3,ThisCarr,CurrentTest)  = ThisDeciLev1;
                end
                
%                             for ii = 1:N1
%                               media(ii) = mean(olho(ii,:));
%                               nhi       = find(olho(ii,:) > media(ii));
%                               nlow      = find(olho(ii,:) < media(ii));
%                               hi(ii)    = mean(olho(ii,nhi));
%                               low(ii)   =  mean(olho(ii,nlow));
%                             %
%                               Lhi(ii)   = mean(hi)  - std(olho(ii,nhi));
%                               Llow(ii)  = mean(low) + std(olho(ii,nlow));
%                               Ab(ii)    = Lhi(ii) - Llow(ii);
%                             end
%                             Abertura   = mean(Ab);
                            a=1;
                %##############################################################
                %For accurately read the income signal is necessary more than
                %just have the decision level, it is also needed to know where
                %it occurs. For instance, if we decide to use just to take the
                %decision level and measure the mean value across and a portion
                %of the time period, which portion should we take? It would be
                %logical to take the central portion, but the bit may be
                %deformed in such way that the information gets concentrated
                %out of the middle part, which means the symbol is not
                %symmetric. The symmetry information can be acquired from the
                %eye diagram by measuring the longitudinal opening. The
                %following sketch better describes this process:
                %
                %                   Point of Symmetry
                %         _____     _______|_______     _____
                %              \   /               \   /
                %               \ /  Longitudinal   \ /
                %                \ ----------------- /
                %               / \    Opening      / \
                %         _____/   \_______________/   \_____
                %
                %With those two pieces of information, decision level and point
                %of symmetry, we have the X, Y coordinates for the centre of
                %the Eye Diagram. Therefore, as long as there is an opening on
                %it it will be possible to recover the transmitted information
                %without error, theoretically.
                
                %As this process is also statistical, first we reshape the
                %income vector to analyze all periods at the same time.
                %             EyeSymMat = reshape(Ix(1+SyncPeriod*NPPB:end-SyncPeriod*NPPB...
                %                 ),NPPB,(Nb4Pam/2) - (2*SyncPeriod));
                %             EyeSymMat = EyeSymMat(1+NPPB/4:end-NPPB/4,:);
                %Then we take the values that compose the decision level
                %because they will mark the point of symmetry.
                %From the location of the max values that occured, which means
                %the uper and lower level of the eye diagram it needs to take
                %the actual value that those occurences represent that is
                %withing the the Interval variable.
                %             ValToSeek = Interval(EyeLoc-1);
                %The number of ocurrences is a statical measure therefore one
                %does not have control which interval will have the highest
                %peak, thus it is important to ordenate the values to be seek
                %from the lower part of the eye diagram to the uper part of the
                %eye diagram.
                %             ValToSeek = sort(ValToSeek,'ascend');
                %             OccuCount = zeros(1,size(EyeSymMat,1));                        %Auxiliar Variable for accounting.
                %             for kk=1:size(EyeSymMat,1)                                     %For every sample within a symbol period
                %                 if (kk>=NPPB/4)&&(kk<=NPPB-NPPB/4)
                %                     OccuCount(kk) = OccuCount(kk) + sum((EyeSymMat(kk,:)...
                %                         >=LowSymPer*ValToSeek(1))&(EyeSymMat(kk,:)<=...
                %                         UpeSymPer*ValToSeek(1)));%Account all occurencies of the valeu 1
                %                     OccuCount(kk) = OccuCount(kk) + sum((EyeSymMat(kk,:)...
                %                         >=LowSymPer*ValToSeek(2))&(EyeSymMat(kk,:)<=...
                %                         UpeSymPer*ValToSeek(2)));%Account all occurencies of the valeu 2
                %                     OccuCount(kk) = OccuCount(kk) + sum((EyeSymMat(kk,:)...
                %                         >=LowSymPer*ValToSeek(3))&(EyeSymMat(kk,:)<=...
                %                         UpeSymPer*ValToSeek(3)));%Account all occurencies of the valeu 3
                %                     OccuCount(kk) = OccuCount(kk) + sum((EyeSymMat(kk,:)...
                %                         >=LowSymPer*ValToSeek(4))&(EyeSymMat(kk,:)<=...
                %                         UpeSymPer*ValToSeek(4)));%Account all occurencies of the valeu 4
                %                 end
                %             end
                
                %The point of symmetry of the eye diagram will be where the
                %maximum number of occurrences were measured inasmuch as those
                %are the points where all the bits go to the center of the
                %symbol.  From the maximum number of occurrences, it can happen
                %for more than one sample within one symbol period, in the best
                %case, all samples would have the same accounting as it is
                %shown the ilustration above hence the symmetry will be at the
                %middle sample of this group of maximum occurrences. This value
                %can be found at the highest peak founded.
                %             [~,SymLoc] = findpeaks(OccuCount,'SortStr','descend');       %The peak on the Eye profile will be the Symmetry level
                
                %             PosIx = SymLoc(1):NPPB:length(Ix);
                
                
                
                
                
                %To find the Eye opening this script will concider just the
                %samples at the center of a symbol. Thus, it is based on the
                %fact that the syncronism process was done succefully.
                %##############################################################
                PosIx = NPPB/2:NPPB:length(Ix);                                %Possition of the central samples - but the number of samples per symbol is even ????
                IxAux = Ix(PosIx);                                             %From the main received signal just a few samples are taken for further evaluation
                n     = 100;                                                   %The number of boxes to be filled up on the histogram process
                %At this process each eye will be evaluated separetelly. It is
                %important to remember that the PAM4 has 4 levels, which means
                %three level of decissions that are exaclty the center of the
                %eye diagram.
                IxAuxAB = IxAux((IxAux<=Levels(4))&(IxAux>=Levels(3)));        %Taking just those values relative to the uper eye
                InterAB = linspace(Levels(3),Levels(4),n);                     %Building the histogram boxes
                EyeAB = hist(IxAuxAB,InterAB);                                 %filling up the boxes with samples that fit on them.
                %The same process described for the uper level will be done at
                %the middle and lower eyes levels.
                IxAuxCD = IxAux((IxAux<=Levels(3))&(IxAux>=Levels(2)));
                InterCD = linspace(Levels(2),Levels(3),n);%NPPB*2^n);
                EyeCD = hist(IxAuxCD,InterCD);
                
                IxAuxEF = IxAux((IxAux<=Levels(2))&(IxAux>=Levels(1)));
                InterEF = linspace(Levels(1),Levels(2),n);%NPPB*2^n);
                EyeEF = hist(IxAuxEF,InterEF);
                
                %What we are looking for here are not where exist the
                %occurrences of values. The eye is where there is not samples,
                %this means, the decission level is a reagion between actual
                %levels thus this reagion does not contain any signal.
                EyeAB = ~EyeAB;                                                %Changing zeros to one - Zeros compose the eye region
                EyeCD = ~EyeCD;
                EyeEF = ~EyeEF;
                %Starting variables that will be used to identify the eye
                %diagram.
                CountAB   = 1;
                CountCD   = 1;
                CountEF   = 1;
                
                SeqOnesAB = 0;
                SeqOnesCD = 0;
                SeqOnesEF = 0;
                
                SeqFinAB  = 0;
                SeqFinCD  = 0;
                SeqFinEF  = 0;
                
                SeqIniAB  = 1;
                SeqIniCD  = 1;
                SeqIniEF  = 1;
                %The for loop will take account of every box with ones. It is
                %important to take note that the not operator was used in this
                %vector, therefore ones means zeros (the eye diagram -
                %possibly) and zeros means values abouve zeroa (not the eye).
                for kk=1:length(EyeAB)                                         %For every box
                    if EyeAB(kk)                                               %if it contains "1"
                        SeqOnesAB(SeqIniAB)=CountAB;                           %count this element as part of a consecutive sequency
                        CountAB = CountAB + 1;                                 %adds one to the counter of consecutive elements "1"
                        if kk==length(EyeAB)                                   %if the current box is the last box we got to an end
                            SeqFinAB(SeqIniAB) = kk;                           %The final sequency element is equal to its possition (kk)
                        end
                    else                                                       %else if the current box contains "0"
                        SeqFinAB(SeqIniAB) = kk-1;                             %the previous element was the last element of a consecutive sequency of ones
                        SeqIniAB = SeqIniAB + 1;                               %adds one to the counter of consecutive number (take in account how many sequencies there are)
                        CountAB = 1;                                           %reset the counter of consecutive elements to mark the start of a new consecutive sequence
                    end
                    
                    if EyeCD(kk)
                        SeqOnesCD(SeqIniCD)=CountCD;
                        CountCD = CountCD + 1;
                        if kk==length(EyeCD)
                            SeqFinCD(SeqIniCD) = kk;
                        end
                    else
                        SeqFinCD(SeqIniCD) = kk-1;
                        SeqIniCD = SeqIniCD + 1;
                        CountCD = 1;
                    end
                    
                    if EyeEF(kk)
                        SeqOnesEF(SeqIniEF)=CountEF;
                        CountEF = CountEF + 1;
                        if kk==length(EyeEF)
                            SeqFinEF(SeqIniEF) = kk;
                        end
                    else
                        SeqFinEF(SeqIniEF) = kk-1;
                        SeqIniEF = SeqIniEF + 1;
                        CountEF = 1;
                    end
                end
                
                
                %If the eye is open, which means there is a clear difference
                %between adjacent levels, the eye diagram will be the longest
                %sequence of ones.
                [MaxValAB,LocMaxAB]=max(SeqOnesAB);
                if LocMaxAB<1 || MaxValAB<2                                    %if any sequency was found or there is just one sequency it is a error thus
                    LevDec3 = 0.6953;                                            %the decission level will be by default 0.67. Also, the other variables will
                    LocMaxAB = 1;                                              %will be set with values to not cause errors in the future
                    SeqFinAB(1)=2;
                    MaxValAB = 0;
                    InterAB(1)=LevDec3;
                else                                                           %if a sequency was found the middle of the eye will be the middle of the sequency
                    if (SeqFinAB(LocMaxAB)-MaxValAB/2)<1                       %if for some reason the final element of the sequency minus the half of its
                        LevDec3 = 0.6953;                                         %length results in a negative value, something went very wrong, and by
                    else                                                       %default it will be set to 0.7
                        LevDec3 = InterAB(round(SeqFinAB(LocMaxAB)-MaxValAB/...
                            2));%Otherwise, the decission level is the middle of the sequency
                    end
                end
                [MaxValCD,LocMaxCD]=max(SeqOnesCD);
                if LocMaxCD<1 || MaxValCD<2
                    LevDec2 = 0.4013;
                    LocMaxCD = 1;
                    SeqFinCD(1)=2;
                    MaxValCD = 0;
                    InterCD(1)=LevDec2;
                else
                    if (SeqFinCD(LocMaxCD)-MaxValCD/2)<1
                        LevDec2 = 0.4013;
                    else
                        LevDec2 = InterCD(round(SeqFinCD(LocMaxCD)-MaxValCD/2));
                    end
                end
                
                [MaxValEF,LocMaxEF]=max(SeqOnesEF);
                if LocMaxEF<1 || MaxValEF<2
                    LevDec1 = 0.1877;
                    LocMaxEF = 1;
                    SeqFinEF(1)=2;
                    MaxValEF = 0;
                    InterEF(1)=LevDec1;
                else
                    if (SeqFinEF(LocMaxEF)-MaxValEF/2)<1
                        LevDec1 = 0.1877;
                    else
                        LevDec1 = InterEF(round(SeqFinEF(LocMaxEF)-MaxValEF/2));
                    end
                    
                end
                
                %another way to measure the eye opening is the get all the
                %boxes and find all peaks on it, that will be a plato created
                %by the sequences of ones (which was zeros). From thos peaks,
                %the eye diagram will be the longer of them hence it will take
                %the most part of the vector that store the findpeaks result.
                %Thus, the middle of the eye will be basically the middle of
                %the peaks vector.
                LocAB = find(EyeAB);
                if isempty(LocAB)                                              %if for some reason there are no peaks, something went wrong.
                    LocAB = 1;
                    LevelDec3 = 0.6953;%mean(Levels(3:4));                       %by default the decission level will be set to 0.65
                else
                    if DecMod%InterAB(LocAB(round(end/2)))<=LevDec3
                        LevelDec3 = LevDec3;
                    else
                        LevelDec3 = InterAB(LocAB(round(end/2)));
                    end
                end
                
                LocCD = find(EyeCD);
                if isempty(LocCD)
                    LocCD = 1;
                    LevelDec2 = 0.4013;%mean(Levels(2:3));
                else
                    if DecMod%InterCD(LocCD(round(end/2)))>=LevDec2
                        LevelDec2 = LevDec2;
                    else
                        LevelDec2 = InterCD(LocCD(round(end/2)));
                    end
                end
                
                LocEF = find(EyeEF);
                if isempty(LocEF)
                    LocEF = 1;
                    LevelDec1 = 0.1877;%mean(Levels(1:2));
                else
                    if DecMod%InterEF(LocEF(round(end/2)))<=LevDec1
                        LevelDec1 = LevDec1;
                    else
                        LevelDec1 = InterEF(LocEF(round(end/2)));
                    end
                end
                
                AberLev(1,ThisCarr,CurrentTest)  = abs(InterAB(SeqFinAB(LocMaxAB)) - InterAB(SeqFinAB(LocMaxAB)-MaxValAB+1));
                AberLev(2,ThisCarr,CurrentTest)  = abs(InterCD(SeqFinCD(LocMaxCD)) - InterCD(SeqFinCD(LocMaxCD)-MaxValCD+1));
                AberLev(3,ThisCarr,CurrentTest)  = abs(InterEF(SeqFinEF(LocMaxEF)) - InterEF(SeqFinEF(LocMaxEF)-MaxValEF+1));
                ValsLev(1,ThisCarr,CurrentTest)  = LevDec3;
                ValsLev(2,ThisCarr,CurrentTest)  = LevDec2;
                ValsLev(3,ThisCarr,CurrentTest)  = LevDec1;
                ValsLev2(1,ThisCarr,CurrentTest) = InterAB(LocAB(round(end/2)));
                ValsLev2(2,ThisCarr,CurrentTest) = InterCD(LocCD(round(end/2)));
                ValsLev2(3,ThisCarr,CurrentTest) = InterEF(LocEF(round(end/2)));
                %%           Ploting for Qualitative Analizes
                if PrintinEye
                    PrintInfo(PlotingThis*36,Ix,T,NPPB);
                    hold all;
                    plot(t(NPPB/2),LevDec1,'bd');plot(t(NPPB/2),InterEF(SeqFinEF(LocMaxEF)-1),'bo');plot(t(NPPB/2),InterEF(round(SeqFinEF(LocMaxEF)-MaxValEF+1)),'bx');
                    plot(t(NPPB/2),LevDec2,'gd');plot(t(NPPB/2),InterCD(SeqFinCD(LocMaxCD)-1),'go');plot(t(NPPB/2),InterCD(round(SeqFinCD(LocMaxCD)-MaxValCD+1)),'gx');
                    plot(t(NPPB/2),LevDec3,'kd');plot(t(NPPB/2),InterAB(SeqFinAB(LocMaxAB)-1),'ko');plot(t(NPPB/2),InterAB(round(SeqFinAB(LocMaxAB)-MaxValAB+1)),'kx');
                    if CalcS
                        plot(t(NPPB/2),ThisDeciLev3,'cd');plot(t(NPPB/2),mean(LevHi3),'co');plot(t(NPPB/2),mean(LevLo3),'cx');
                        plot(t(NPPB/2),ThisDeciLev2,'rd');plot(t(NPPB/2),mean(LevHi2),'ro');plot(t(NPPB/2),mean(LevLo2),'rx');
                        plot(t(NPPB/2),ThisDeciLev1,'md');plot(t(NPPB/2),mean(LevHi1),'mo');plot(t(NPPB/2),mean(LevLo1),'mx');
                    end
                    a=1;
                end
%                             plot(t(NPPB/2),InterAB(LocAB(round(end/2))),'kx');
%                             plot(t(NPPB/2),InterCD(LocCD(round(end/2))),'gx');
%                             plot(t(NPPB/2),InterEF(LocEF(round(end/2))),'bx');
%                             plot([t(SymLoc(1)) t(SymLoc(1))],[0 1]);
                %
                %             drawnow;
                %             PrintInfo(PlotingThis*38,EyeMax,EyeLoc-1,Interval);
                % %                         figure;plot(SeqOnesEF);figure;plot(SeqOnesCD);figure;plot(SeqOnesAB)
                %             LevDec3
                %             a1=InterAB(LocAB(round(end/2)))
                %             b1=ThisDeciLev1
                %             LevDec2
                %             a2=InterCD(LocCD(round(end/2)))
                %             b2=ThisDeciLev2
                %             LevDec1
                %             c1=InterEF(LocEF(round(end/2)))
                %             c2=ThisDeciLev3
                %             PrintInfo(Ploting*39,NPPB,OccuCount);
                %             PrintInfo(Ploting*40,OccuCount,SymLoc);
                %%      Actualy Receiving Data
                %Once the signal was processed the next step is through a
                %comparator decide the actual information received.
                IxRec    = [];                                                 %Initialization of the vector that will store the income data
                IxRecDef = [];                                                 %Initialization of the vector that will store the income data
                IxRecDeS = [];                                                 %Initialization of the vector that will store the income data
                for kk=1:NPPB:length(Ix)                                       %The comparison process will be made for each symbol period
                    %                 midaux = round(mean(SymLoc(1:round(end/2))));
                    midaux = NPPB/2;%SymLoc(1);
                    aux1 = Ix((kk-1)+midaux);     %An small portion of the income signal is take for evaluation
                    MeanRec = mean(aux1);                                      %Measuring the avarage value of the samples taken
                    %Verifying the interval for each symbol received.
                    if MeanRec <= LevelDec1                                    %If it is the lowest level the incoming data
                        IxRec = [IxRec 0 0];                                   %is 01 (1)
                    elseif (MeanRec <= LevelDec2)&&(MeanRec > LevelDec1)       %If it is the second level the incoming data
                        IxRec = [IxRec 0 1];                                   %is 00 (0)
                    elseif (MeanRec <= LevelDec3)&&(MeanRec > LevelDec2)       %If it is the tird level the incoming data
                        IxRec = [IxRec 1 1];                                   %is 10 (2)
                    elseif MeanRec > LevelDec3                                 %If it is the uper level the incoming data
                        IxRec = [IxRec 1 0];                                   %is 11 (3)
                    else                                                       %If for some misteriose reason neither of previous verification were sucedded
                        IxRec = [IxRec 0 0];                                   %by default the current data is set to be 00 (0)
                    end
                    
                    %Verifying the interval for each symbol received.
                    if MeanRec <= DecLevDef1                                   %If it is the lowest level the incoming data
                        IxRecDef = [IxRecDef 0 0];                             %is 01 (1)
                    elseif (MeanRec <= DecLevDef2)&&(MeanRec > DecLevDef1)     %If it is the second level the incoming data
                        IxRecDef = [IxRecDef 0 1];                             %is 00 (0)
                    elseif (MeanRec <= DecLevDef3)&&(MeanRec > DecLevDef2)     %If it is the tird level the incoming data
                        IxRecDef = [IxRecDef 1 1];                             %is 10 (2)
                    elseif MeanRec > DecLevDef3                                %If it is the uper level the incoming data
                        IxRecDef = [IxRecDef 1 0];                             %is 11 (3)
                    else                                                       %If for some misteriose reason neither of previous verification were sucedded
                        IxRecDef = [IxRecDef 0 0];                             %by default the current data is set to be 00 (0)
                    end
                    
                    if CalcS
                        if MeanRec <= ThisDeciLev1                                   %If it is the lowest level the incoming data
                            IxRecDeS = [IxRecDeS 0 0];                             %is 01 (1)
                        elseif (MeanRec <= ThisDeciLev2)&&(MeanRec > ThisDeciLev1)     %If it is the second level the incoming data
                            IxRecDeS = [IxRecDeS 0 1];                             %is 00 (0)
                        elseif (MeanRec <= ThisDeciLev3)&&(MeanRec > ThisDeciLev2)     %If it is the tird level the incoming data
                            IxRecDeS = [IxRecDeS 1 1];                             %is 10 (2)
                        elseif MeanRec > ThisDeciLev3                                %If it is the uper level the incoming data
                            IxRecDeS = [IxRecDeS 1 0];                             %is 11 (3)
                        else                                                       %If for some misteriose reason neither of previous verification were sucedded
                            IxRecDeS = [IxRecDeS 0 0];                             %by default the current data is set to be 00 (0)
                        end
                    end
                end
                %%           Ploting for Qualitative Analizes
                %             PrintInfo(Ploting*41,t(length(TxDataMat(ThisCarr,:))),Nb4Pam...
                %                                            /2,TxDataMat(ThisCarr,:),IxRec);
                %%       Calculating the Bit Error Ratio (BER)
                %The final process here is to count the number of wrongdoings
                %of this whole process upon the transmited data for
                %quantitative analizes
                BitErr = sum(xor(TxDataMat(ThisCarr,1+2*SyncPeriod:end-2*...
                    SyncPeriod),IxRec(1+2*SyncPeriod:end-2*SyncPeriod)));%Comparison between the Transmited and received and counting the differences
                BitErrAux1 = BitErr;
                BitErrAux2 = sum(xor(TxDataMat(ThisCarr,1+2*SyncPeriod:end-...
                    2*SyncPeriod),IxRecDef(1+2*SyncPeriod:end-2*SyncPeriod)));%Comparison between the Transmited and received and counting the differences
                if BitErr ~= 0
                    if BitErrAux2<BitErrAux1
                        BitErr = BitErrAux2;
                    end
                else
                    DecLevDef1 = LevelDec1;
                    DecLevDef2 = LevelDec2;
                    DecLevDef3 = LevelDec3;
                end
                if CalcS
                    BitErrS = sum(xor(TxDataMat(ThisCarr,1+2*SyncPeriod:end-2*...
                        SyncPeriod),IxRecDeS(1+2*SyncPeriod:end-2*SyncPeriod)));%Comparison between the Transmited and received and counting the differences
                    Ber4PAMS(CurrentTest,ThisCarr) = BitErrS/(Nb4Pam-(2*SyncPeriod));%Calculating the ration of wrong bits and the total number of bits transmited
                end
                Ber4PAM(CurrentTest,ThisCarr) = BitErr/(Nb4Pam-(2*SyncPeriod));%Calculating the ration of wrong bits and the total number of bits transmited
%                 Ber4PAM(CurrentTest,ThisCarr+1) = 1;%Calculating the ration of wrong bits and the total number of bits transmited
                %             if BitErr>0
                %                 berpos = 1:2:size(Ber4PAM,2);
                %                 Ber4PAM(size(Ber4PAM,1),berpos)
                %                 drawnow;%pause(1);
                %                 a=2;
                %             end
                %             berpos = 1:2:size(Ber4PAM,2);
                %             Ber4PAM(size(Ber4PAM,1),berpos)
                %             drawnow;%pause(1);
                %             a=3;
                %             drawnow;
                %             pause(0.5);
                            close all;
                
            end
        end
    otherwise
        %%              Receiver OOK
%         AberLev(CurrentTest,:) = zeros(1,NumCarr);
        for ThisCarr=InitCarr:CarrPass:NumCarr                                           %For each carrier the same process of reception will be used.
            if CarrUsedDo(ThisCarr)
                %%  Reception Process: Handling the income optical field
                %At the first step the income field will be selected from the
                %optical FFT Output with the help of the maping vector
                %(previously described).
                EoutAux = EoutAux1(VetThisCarr==ObsCarrUsed(ThisCarr),:);
                if AddingNoiseP
    %                 figure;plot(f,20*log10(abs(fftshift(fft(EoutAux)./length(EoutAux)))));
    %                 [~,sigpower] = MeasPower(EoutAux);
                    [p,l] = findpeaks(20*log10(abs(fftshift(fft(EoutAux)./length(EoutAux)))),'SortStr','descend','NPeaks',1);
                    EoutAuxA = fftshift(fft(EoutAux));
                    EoutAuxA(f==f(l)) = EoutAuxA(f==f(l-3));
    %                 [testfilt,~] = FiltroGaussiano(f,12.5e9,f(l),5);           %Creating filter for selection of the received signal
    %                 testfilt = fftshift(testfilt);                                   %Shifting the filter for matching the received signal
    %                 EoutAux = ifft(fftshift(EoutAuxA.*testfilt));
                    EoutAuxA = ifft(fftshift(EoutAuxA));
    %                 figure;plot(f,20*log10(abs(fftshift(fft(EoutAuxA)./length(EoutAuxA)))));
    %                 [~,sigpower] = MeasPower(EoutAux);
                    [~,sigpower] = MeasPower(EoutAuxA);
                    EoutAux      = awgn(EoutAux,osnrp,sigpower-30);
                    [~,sigpower] = MeasPower(EoutAux);
                end
                %%            Fiber Time Delay Compensation
                switch Medium
                    case 'Fiber'
                        EoutAux = ifft(fft(EoutAux).*exp(1j*2*pi*f*(...
                            FiberDelay((ThisCarr))*Ta)));
                    otherwise
                end
                
                %The current incoming signal is them converted from the optical
                %domain to the eletrical domain with the help of an photo
                %detector.
                Ix =EoutAux.*conj(EoutAux);
                %%           Creating the Reception Filter
                [BitFilt,~] = FiltroGaussiano(f,BWD,CenFeq,FiltOrd);           %Creating filter for selection of the received signal
                BitFilt = fftshift(BitFilt);                                   %Shifting the filter for matching the received signal
                %%
                if RecFilBanPas
                    Ix = ifft(fft(Ix).*BitFilt);
                    Ix = Ix - min(Ix);                                         %Removing the DC component from them Eletrical signal received
                    Ix = Ix./max(abs(Ix));                                     %Normalizing the eletrical signal (amplifying what is needed)
                end
                
                if ReceptorNoise
                    [~,SigPower] = MeasPower(Ix-mean(Ix));
                    SigPower2 = SigPower-30;%10*log10(SigPower);
                    SNR2 = CarSNR + 10*log10(1) - 10*log10(Nsamp) - 10*log10(1/2) + 10*log10(10^0.3);
                    SNR = CarSNR + 10*log10(1) - 10*log10((t(NPPB))/(3*(t(2*NumAmosCP+NPPB))/(NPPB)));
                    Ix = awgn(Ix,SNR,SigPower2);
                end
                if ~RecFilBanPas
                    Ix = ifft(fft(Ix).*BitFilt);
                    Ix = Ix - min(Ix);                                         %Removing the DC component from them Eletrical signal received
                    Ix = Ix./max(abs(Ix));                                     %Normalizing the eletrical signal (amplifying what is needed)
                end
                VetElecPower(CurrentTest,ThisCarr)= MeasPower(Ix);
                EoutAuxF = 20*log10(abs(fftshift(fft(EoutAux)./length(...
                    EoutAux))));
                [VetOptiPower(CurrentTest,ThisCarr),~]= findpeaks(EoutAuxF,...
                    'SortStr','descend','NPeaks',1);
                %% Ploting the result for qualitative analizes
                %             PrintInfo(Ploting*42,f,20.*log10(abs(fftshift(fft(Ix)./...
                %                                                             length(Ix)))));
                %             PrintInfo(Ploting*43,t,TxSigMat(ThisCarr,:)./max(TxSigMat(ThisCarr,:)),Ix);
                %             PrintInfo(Ploting*44,Ix,T,NPPB);
                %             PrintInfo(Ploting*45,t(IniSyncPos:SyncPos),Ix(IniSyncPos:...
                %                                     SyncPos),SyncSymb(IniSyncPos:SyncPos));
                %%        synchronizing
                %For the reception process work properly, it is needed to
                %sincronized the recieved signal or the sampling process will
                %not work.
                
                AuxSync = (Ix(IniSyncPos:SyncPos));                            %Selecting the sync-word within the received signal
                AuxSync1 = AuxSync;
                %This synchronization process is based on the mean expected
                %value. That means, the information mean value should be within
                %one period of symbol. Thus, the mean value of the received
                %signal is acquired and compare of the known sync-word to
                %verify if this mean value is at the right possition
                AuxSync1(AuxSync1<0) = 0;                                      %To keep the mean value above zero anything under is neglected
                AuxSync1(AuxSync1>=mean(AuxSync1)) = 1;                        %Adding a flag to the first sample of the received mean value
                AuxSync1(AuxSync1<mean(AuxSync1)) = -1;                        %All the others samples at set to the lowest level
                AuxSync2 = SyncSymb(IniSyncPos:SyncPos);                       %Selecting the sync-word within the known signal
                AuxSync2(AuxSync2>=mean(AuxSync2)) = 1;                        %Adding a flag to the first sample of the known mean value
                AuxSync2(AuxSync2<mean(AuxSync2)) = -1;                        %All the others samples at set to the lowest level
                
                PosToSyn  = find(ismember(AuxSync1,1));
                PosSyn    = find(ismember(AuxSync2,1));
                
                AuxSyncCorr = PosToSyn(round(end/2))-PosSyn(round(end/2));
                %The difference between the PossitionTosynchronize and
                %Possitionsynchronized will be used to correct the time
                %shifting on the transmition and reception process.
                
                if AuxSyncCorr>=0%If the difference is positive, left-shift...
                    %                 Ix = ifft(fft(Ix).*exp(1j*2*pi*f*(AuxSyncCorr*Ta)));     %Shift based on time change
                    Ix = [Ix(AuxSyncCorr+1:end) Ix(1:AuxSyncCorr)];            %Shift based on sampling sliding
                else%... but if the difference is negative, right-shift
                    %                 Ix = ifft(fft(Ix).*exp(-1j*2*pi*f*(AuxSyncCorr*Ta)));    %Shift based on time change
                    Ix = [Ix(end+AuxSyncCorr+1:end) Ix(1:end + AuxSyncCorr)];  %Shift based on sampling sliding
                end
                
                
                %%        synchronizing
                if SencondAdjust
                    %For some reason that we could not understand sometimes the
                    %time (sampling) sliding of the signal is not equal
                    %throught the data stream. Thus, the second part of the
                    %synchronism process will be turn ON or OFF according to the
                    %user's will.
                    AuxSyncEnd = (Ix(end-(SyncPos+1*abs(AuxSyncCorr))+1:end-...
                        (2*NPPB+1*abs(AuxSyncCorr))));%Selecting the sync-word within the received signal
                    SyncSymbEndAux = SyncSymbEnd(end-SyncPos+1:end-2*NPPB);    %Selecting the sync-word within the known signal
                    AuxSyncEnd1 = AuxSyncEnd;
                    AuxSyncEnd1(AuxSyncEnd1<0) = 0;                            %To keep the mean value above zero anything under is neglected
                    AuxSyncEnd1(AuxSyncEnd1>=mean(AuxSyncEnd1)) = 1;           %Adding a flag to the first sample of the received mean value
                    AuxSyncEnd1(AuxSyncEnd1<mean(AuxSyncEnd1)) = -1;           %All the others samples at set to the lowest level
                    AuxSyncEnd2 = SyncSymbEndAux;
                    AuxSyncEnd2(AuxSyncEnd2>=mean(AuxSyncEnd2)) = 1;           %Adding a flag to the first sample of the known mean value
                    AuxSyncEnd2(AuxSyncEnd2<mean(AuxSyncEnd2)) = -1;           %All the others samples at set to the lowest level
                    
                    PosToSynEnd  = find(ismember(AuxSyncEnd1,1));
                    PosSynEnd = find(ismember(AuxSyncEnd2,1));
                    
                    AuxSyncCorrEnd = PosToSynEnd(round(end/2))-PosSynEnd(...
                        round(end/2));
                    
                    %The difference between the PossitionTosynchronize and
                    %Possitionsynchronized will be used to correct the time
                    %shifting on the transmition and reception process.
                    
                    if AuxSyncCorrEnd>=0%If possitive difference, left-shift...
                        %                     Ix = ifft(fft(Ix).*exp(1j*2*pi*f*(AuxSyncCorrEnd*...
                        %                                                                    Ta)));%Shift based on time change
                        Ix = [Ix(AuxSyncCorrEnd+1:end) Ix(1:AuxSyncCorrEnd)];
                    else%... but if the difference is negative, right-shift
                        %                     Ix = ifft(fft(Ix).*exp(-1j*2*pi*f*(AuxSyncCorrEnd*...
                        %                                                                    Ta)));%Shift based on time change
                        Ix = [Ix(end-AuxSyncCorrEnd+1:end) Ix(1:end - ...
                            AuxSyncCorrEnd)];%Shift based on sampling sliding
                    end
                end
                
                %% Ploting the result for qualitative analizes
                
%                             PrintInfo(PlotingThis*45,t(IniSyncPos:SyncPos),Ix(IniSyncPos:SyncPos),SyncSymb(IniSyncPos:SyncPos));
%                             PrintInfo(PlotingThis*46,t(1:length(Ix)),Ix);
                %%                  Removing CP
                if AddCP
                    IxAux = Ix(1:end - StuffSampels);
                    IxAux = reshape(IxAux,(2*NumAmosCP+NPPB),Nb-NumBitDesc);
                    IxAux = IxAux(1+NumAmosCP:end-NumAmosCP,:);
                    IxAux = reshape(IxAux,1,NPPB*(Nb-NumBitDesc));
                    Ix    = IxAux;
                end
                %% Taking the sampling the EVM meassurement
%                 evm = comm.EVM('MaximumEVMOutputPort',true,...
%                 'XPercentileEVMOutputPort',true, 'XPercentileValue',90,...
%                 'SymbolCountOutputPort',true);
                PosAuxEout = NPPB/2:NPPB:length(Ix);                       %Varriable respossible to take just the samples at the middle of the symbol
                IxAux      = Ix(PosAuxEout);                               %Normalizing the reference
                a=0;
                EvmMatRec(ObsCarrPos==ThisCarr,:) = IxAux;                       %Taking just the middle samples as references
                [EvmDB(CurrentTest,ThisCarr), EvmPer(CurrentTest,ThisCarr), EvmRms(CurrentTest,ThisCarr) ] = EvmCalc( EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux );
                [EvmDBJ(CurrentTest,ThisCarr),EvmPerJ(CurrentTest,ThisCarr),EvmRmsJ(CurrentTest,ThisCarr)] = evm1(2,'pam',EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux);
%                 IxRef = EvmMatRef(ObsCarrPos==ThisCarr,:);
%                 [MatEvem,maxEVM,pctEVM,numSym] = evm(abs(IxRef),abs(IxAux));
%                 EvmRmsM(CurrentTest,ThisCarr) = MatEvem;
%                 EvmPerM(CurrentTest,ThisCarr) = MatEvem*100;
%                 EvmDBM(CurrentTest,ThisCarr) = 20*log10(MatEvem);
%                 if ThisCarr==1
%                     figure;hold all;grid on;
%                     plot(IxAux,zeros(1,length(IxAux)),'rx');
%                     set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%                     a=0;
%                 end
                %%  Measuring the EVM
%                 IxPosAux = NPPB/2:NPPB:length(Ix);
%                 IxRec    = Ix(IxPosAux);
%                 [EvmRms(CurrentTest,ThisCarr),EvmDB(CurrentTest,ThisCarr),...
%                     EvmPer(CurrentTest,ThisCarr)] = EvmCalc( EvmMatRef(ThisCarr,...
%                     :),IxRec );
%##########################################################################
                PosIx = NPPB/2:NPPB:length(Ix);                                %Possition of the central samples - but the number of samples per symbol is even ????
                IxAux = Ix(PosIx);                                             %From the main received signal just a few samples are taken for further evaluation
                n     = 100;                                                   %The number of boxes to be filled up on the histogram process
                MinDist = n/3;
                Perce = 0.5;
                %At this process each eye will be evaluated separetelly. It is
                %important to remember that the PAM4 has 4 levels, which means
                %three level of decissions that are exaclty the center of the
                %eye diagram.
                Interval = linspace(Perce*(min(Ix)),Perce*(max(Ix)),n);                     %Building the histogram boxes
                
                %Therefore, the MATLAB hist function returns the number of
                %occurrence of each interval.
                EyeMax = hist(Ix,Interval);
                EyeMaxaux = [0 EyeMax 0];                                      %Zeros are added at the EyeMax to auxiliate the finding peaks process
                [~,EyeLoc] = findpeaks(EyeMaxaux,'MinPeakDistance',MinDist,'SortStr','descend','NPeaks',2);%The peaks on the Eye profile will be the levels at the Eyes limit
                
                if length(EyeLoc)<2%If it was not able to find the eye profile.
                    [~,EyeLoc] = findpeaks(EyeMaxaux,'MinPeakDistance',MinDist/2,'SortStr','descend','NPeaks',2);%The peaks on the Eye profile will be the levels at the Eyes limit
                end
                
                if length(EyeLoc)<2%If it was not able to find the eye profile.
                    EyeLoc = [2 3];
                    Levels = [0 -0.2 0.35];
                else%Whereas, if there is an profile the decission can be found
                    Levels = [Interval(EyeLoc(1)-1) Interval(EyeLoc(2)-1)];
                    Levels = sort(Levels);
                    %                 LevelDec1 = Levels(1) + (Levels(2)-Levels(1))/2 ;
                    %                 LevelDec2 = Levels(2) + (Levels(3)-Levels(2))/2 ;
                    %                 LevelDec3 = Levels(3) + (Levels(4)-Levels(3))/2 ;
                end
                IxAuxAB = IxAux((IxAux<=Levels(2))&(IxAux>=Levels(1)));
                Inter = linspace(Levels(1),Levels(2),n);                     %Building the histogram boxes
                Inter2 = linspace(min(IxAux),max(IxAux),n);                     %Building the histogram boxes
                EyeAB = hist(IxAuxAB,Inter);                                 %filling up the boxes with samples that fit on them.
                EyeCD = hist(IxAux,Inter2);                                 %filling up the boxes with samples that fit on them.
                
                %What we are looking for here are not where exist the
                %occurrences of values. The eye is where there is not samples,
                %this means, the decission level is a reagion between actual
                %levels thus this reagion does not contain any signal.
                EyeAB = ~EyeAB;                                                %Changing zeros to one - Zeros compose the eye region
                EyeCD = ~EyeCD;                                                %Changing zeros to one - Zeros compose the eye region
                %Starting variables that will be used to identify the eye
                %diagram.
                Count   = 1;
                Count2   = 1;
                
                SeqOnes = 0;
                SeqOnes2 = 0;
                
                SeqIni  = 1;
                SeqIni2  = 1;
                
                SeqFin  = 0;
                SeqFin2  = 0;
                %The for loop will take account of every box with ones. It is
                %important to take note that the not operator was used in this
                %vector, therefore ones means zeros (the eye diagram -
                %possibly) and zeros means values abouve zeroa (not the eye).
                for kk=1:length(EyeAB)                                         %For every box
                    if EyeAB(kk)                                               %if it contains "1"
                        SeqOnes(SeqIni)=Count;                           %count this element as part of a consecutive sequency
                        Count = Count + 1;                                 %adds one to the counter of consecutive elements "1"
                        if kk==length(EyeAB)                                   %if the current box is the last box we got to an end
                            SeqFin(SeqIni) = kk;                           %The final sequency element is equal to its possition (kk)
                        end
                    else                                                       %else if the current box contains "0"
                        SeqFin(SeqIni) = kk-1;                             %the previous element was the last element of a consecutive sequency of ones
                        SeqIni = SeqIni + 1;                               %adds one to the counter of consecutive number (take in account how many sequencies there are)
                        Count = 1;                                           %reset the counter of consecutive elements to mark the start of a new consecutive sequence
                    end
                    
                    if EyeCD(kk)                                               %if it contains "1"
                        SeqOnes2(SeqIni2)=Count2;                           %count this element as part of a consecutive sequency
                        Count2 = Count2 + 1;                                 %adds one to the counter of consecutive elements "1"
                        if kk==length(EyeCD)                                   %if the current box is the last box we got to an end
                            SeqFin2(SeqIni2) = kk;                           %The final sequency element is equal to its possition (kk)
                        end
                    else                                                       %else if the current box contains "0"
                        SeqFin2(SeqIni2) = kk-1;                             %the previous element was the last element of a consecutive sequency of ones
                        SeqIni2 = SeqIni2 + 1;                               %adds one to the counter of consecutive number (take in account how many sequencies there are)
                        Count2 = 1;                                           %reset the counter of consecutive elements to mark the start of a new consecutive sequence
                    end
                end
                
                
                %If the eye is open, which means there is a clear difference
                %between adjacent levels, the eye diagram will be the longest
                %sequence of ones.
                [MaxVal,LocMax]=max(SeqOnes);
                if LocMax<2 || MaxVal<2                                    %if any sequency was found or there is just one sequency it is a error thus
                    LevDec = 0.0;                                            %the decission level will be by default 0.67. Also, the other variables will
                    LocMax = 1;                                              %will be set with values to not cause errors in the future
                    SeqFin(1)=2;
                    MaxVal = 0;
                    Inter(1)=LevDec;
                else                                                           %if a sequency was found the middle of the eye will be the middle of the sequency
                    if (SeqFin(LocMax)-MaxVal/2)<1                       %if for some reason the final element of the sequency minus the half of its
                        LevDec = 0.0;                                         %length results in a negative value, something went very wrong, and by
                    else                                                       %default it will be set to 0.7
                        LevDec = Inter(round(SeqFin(LocMax)-MaxVal/2));%Otherwise, the decission level is the middle of the sequency
                    end
                end
                [MaxVal2,LocMax2]=max(SeqOnes2);
                if LocMax2<2 || MaxVal2<2                                    %if any sequency was found or there is just one sequency it is a error thus
                    LevDec2 = 0.0;                                            %the decission level will be by default 0.67. Also, the other variables will
                    LocMax2 = 1;                                              %will be set with values to not cause errors in the future
                    SeqFin2(1)=2;
                    MaxVal2 = 0;
                    Inter2(1)=LevDec2;
                else                                                           %if a sequency was found the middle of the eye will be the middle of the sequency
                    if (SeqFin2(LocMax2)-MaxVal2/2)<1                       %if for some reason the final element of the sequency minus the half of its
                        LevDec2 = 0.0;                                         %length results in a negative value, something went very wrong, and by
                    else                                                       %default it will be set to 0.7
                        LevDec2 = Inter2(round(SeqFin2(LocMax2)-MaxVal2/2));%Otherwise, the decission level is the middle of the sequency
                    end
                end
                
                %another way to measure the eye opening is the get all the
                %boxes and find all peaks on it, that will be a plato created
                %by the sequences of ones (which was zeros). From thos peaks,
                %the eye diagram will be the longer of them hence it will take
                %the most part of the vector that store the findpeaks result.
                %Thus, the middle of the eye will be basically the middle of
                %the peaks vector.
                Loc = find(EyeAB);
                if isempty(Loc)                                              %if for some reason there are no peaks, something went wrong.
                    Loc = 1;
                    LevelDec = 0.0;%mean(Levels(3:4));                       %by default the decission level will be set to 0.65
                else
                    LevelDec = LevDec;
                end
                Loc2 = find(EyeCD);
                if isempty(Loc2)                                              %if for some reason there are no peaks, something went wrong.
                    Loc2 = 1;
                    LevelDec2 = 0.0;%mean(Levels(3:4));                       %by default the decission level will be set to 0.65
                else
                    LevelDec2 = LevDec2;
                end
%##########################################################################

                if CalcS
                    [~,~,EyeOpen,EyeOpenHigh,EyeOpenLow] = Olho((Ix),T,NPPB,0);
                    AberLevS(CurrentTest,ThisCarr)= EyeOpen;
                    ValsLevS(CurrentTest,ThisCarr)= EyeOpenLow + EyeOpen/2;
                else
                    [~,~,EyeOpen,EyeOpenHigh,EyeOpenLow] = Olho((Ix),T,NPPB,0);
                end
                %% Ploting the result for qualitative analizes
                
                if PrintinEye
                    PrintInfo(PlotingThis*47,(Ix),T,NPPB);
                    hold on;
                    plot(t((NPPB)/2),Inter(round(SeqFin(LocMax)-MaxVal)),'kx');
                    plot(t((NPPB)/2),Inter(SeqFin(LocMax)),'ko');
                    plot(t((NPPB)/2),Inter(round(SeqFin(LocMax)-MaxVal/2)),'kd');
                    plot(t((NPPB)/2),Inter2(round(SeqFin2(LocMax2)-MaxVal2)),'bx');
                    plot(t((NPPB)/2),Inter2(SeqFin2(LocMax2)),'bo');
                    plot(t((NPPB)/2),Inter2(round(SeqFin2(LocMax2)-MaxVal2/2)),'bd');
                    if CalcS
                        plot(t((NPPB)/2),EyeOpenLow,'mx');
                        plot(t((NPPB)/2),EyeOpenHigh,'mo');
                        plot(t((NPPB)/2),EyeOpenLow + EyeOpen/2,'md');
                    end
                    a=0;
                end
                AberLev(CurrentTest,ThisCarr)= Inter(SeqFin(LocMax))-Inter(round(SeqFin(LocMax)-MaxVal));
                ValsLev(CurrentTest,ThisCarr)= LevDec;
                %%         Finding Decission Levels
                %The process for decoding the income signal will be based on
                %eletronic comparators. Inasmuch as the right decission level
                %must be acquired for accurately decide, within a symbol
                %periode, what that current leavel means (ones or zeros).
                %
                %The process hereafter of chosing the  decission levels is not
                %deterministic rather it is a statistic process. The main idea
                %is to take the decission level from the histogram generated
                %from the income signal stream.
                %
                %This process is realized inside the function Olho.
                %
                %Basicaly the decission level will be the minimal value of the
                %currente eye under evaluation plus the half of the its eye
                %opening.The following ilustration better describe this process
                %
                %Eye Limits:   + (1/2 * Eye Opening:)  =    Comparisson Limit:
                %
                %UperLevel  ______     ______________     _____
                %                 \   /      |       \   /
                %                  \ /       |        \ /
                %                   \ Half Eye Opening /     Decission Level
                %                  / \       |        / \
                %LowerLevel ______/   \______|_______/   \_____
                %
                %
                %
                %For accurately read the income signal is necessary more than
                %just have the decision level, it is also needed to know where
                %it occurs. For instance, if we decide to use just to take the
                %decision level and measure the mean value across and a portion
                %of the time period, which portion should we take? It would be
                %logical to take the central portion, but the bit may be
                %deformed in such way that the information gets concentrated
                %out of the middle part, which means the symbol is not
                %symmetric. The symmetry information can be acquired from the
                %eye diagram by measuring the longitudinal opening. The
                %following sketch better describes this process:
                %
                %                   Point of Symmetry
                %         _____     _______|_______     _____
                %              \   /               \   /
                %               \ /  Longitudinal   \ /
                %                \ ----------------- /
                %               / \    Opening      / \
                %         _____/   \_______________/   \_____
                %
                %With those two pieces of information, decision level and point
                %of symmetry, we have the X, Y coordinates for the centre of
                %the Eye Diagram. Therefore, as long as there is an opening on
                %it it will be possible to recover the transmitted information
                %without error... theoretically.
                %
                %As this process is also statistical, first we reshape the
                %income vector to analyze all periods at the same time.
                
%                 EyeSymMat = reshape(Ix(1+SyncPeriod*NPPB:end-SyncPeriod*...
%                     NPPB),NPPB,(Nb-NumBitDesc)-2*SyncPeriod);
                %Then we take the values that compose the decision level
                %because they will mark the point of symmetry.
                %
                %Firstly it was set the interval in which the histogram will be
                %build. It is based on the number of samples per bit period.
%                 Interval = linspace(min(Ix),max(Ix),2*NPPB);
                %Therefore, the MATLAB hist function returns the number of
                %occurrence of each interval.
%                 EyeMax = hist(Ix,Interval);
%                 EyeMaxaux = [0 EyeMax 0];                                      %Zeros are added at the EyeMax to auxiliate the finding peaks process
%                 [~,EyeLoc] = findpeaks(EyeMaxaux,'MinPeakDistance',NPPB/4,...
%                     'SortStr','descend','NPeaks',2);%The peaks on the Eye profile will be the levels at the Eyes limit
                
                %             figure;findpeaks(EyeMaxaux,'MinPeakDistance',NPPB/4,'SortStr','descend','NPeaks',2);
                
%                 ValToSeek = Interval(EyeLoc-1);
                %From the location of the max values that occured, which means
                %the uper and lower level of the eye diagram it needs to take
                %the actual value that those occurences represent that is
                %withing the the Interval variable.
%                 ValToSeek = sort(ValToSeek,'ascend');
                %The number of ocurrences is a statical measure therefore one
                %does not have control which interval will have the highest
                %peak, thus it is important to ordenate the values to be seek
                %from the lower part of the eye diagram to the uper part of the
                %eye diagram.
%                 OccuCount = zeros(1,size(EyeSymMat,1));                        %Auxiliar Variable for accounting.
%                 for kk=1:size(EyeSymMat,1)                                     %For every sample within a symbol period
%                     OccuCount(kk) = OccuCount(kk)+sum((EyeSymMat(kk,:)>=min(...
%                         Ix))&(EyeSymMat(kk,:)<=UpeSymPer*EyeOpenLow));%Account all occurencies of the value 1
%                     OccuCount(kk) = OccuCount(kk)+sum((EyeSymMat(kk,:)>=...
%                         LowSymPer*EyeOpenHigh)&(EyeSymMat(kk,:)<=max(Ix)));%Account all occurencies of the value 2
%                 end
                %The point of symmetry of the eye diagram will be where the
                %maximum number of occurrences were measured inasmuch as those
                %are the points where all the bits go to the center of the
                %symbol.  From the maximum number of occurrences, it can happen
                %for more than one sample within one symbol period, in the best
                %case, all samples would have the same accounting as it is
                %shown the ilustration above hence the symmetry will be at the
                %middle sample of this group of maximum occurrences. This value
                %can be found by the mean of the samples positions within the
                %symbol period. The problem with this approach is that the
                %signal must be synchronized with the maximum displacement of
                %a symbol period minus 25% of the eye Longitudinal opening if
                %the displacement is higher than that the point of symmetry
                %will be wrongly measured.
%                 [SymLoc] = round(mean(find(ismember(OccuCount,max(OccuCount)...
%                     ))));%The peak on the Eye profile will be the Symmetry level
                %%           Ploting for Qualitative Analizes
                %             PrintInfo(Ploting*48,NPPB,OccuCount);
                %             PrintInfo(Ploting*49,OccuCount,SymLoc);
                
                %%      Actualy Receiving Data
                %Once the signal was processed the next step is through a
                %comparator decide the actual information received.
                EoutCorr = [];                                                 %Initialization of the vector that will store the income data
                EoutCorrD = [];                                                 %Initialization of the vector that will store the income data
                EoutCorr2 = [];                                                 %Initialization of the vector that will store the income data
                EoutCorrS = [];                                                 %Initialization of the vector that will store the income data
                for kk=1:NPPB:length(Ix)                                       %The comparison process will be made for each symbol period
                    %An small portion of the income signal is take for
                    %evaluation by measuring the avarage value of the samples
                    %taken
                    %                 CalcMean = mean((Ix((kk-1)+SymLoc(1))));
                    CalcMean = mean((Ix((kk-1)+NPPB/2)));
                    %Verifying the interval for each symbol received.
                    if CalcMean >= LevDec%EyeOpenLow+EyeOpen/2                        %If it is the uper level the incoming data
                        EoutCorr = [EoutCorr 1];                               %is 1
                    else                                                       %If it is the lowest level the incoming data
                        EoutCorr = [EoutCorr 0];                               %is 0
                    end
                    CalcMean2 = mean((Ix((kk-1)+NPPB/2)));
                    %Verifying the interval for each symbol received.
                    if CalcMean2 >= LevDec2%EyeOpenLow+EyeOpen/2                        %If it is the uper level the incoming data
                        EoutCorr2 = [EoutCorr2 1];                               %is 1
                    else                                                       %If it is the lowest level the incoming data
                        EoutCorr2 = [EoutCorr2 0];                               %is 0
                    end
                    if CalcMean2 >= EyeOpenLow+EyeOpen/2                        %If it is the uper level the incoming data
                        EoutCorrD = [EoutCorrD 1];                               %is 1
                    else                                                       %If it is the lowest level the incoming data
                        EoutCorrD = [EoutCorrD 0];                               %is 0
                    end
                    if CalcS
                        CalcMeanS = mean((Ix((kk-1)+NPPB/2)));
                        %Verifying the interval for each symbol received.
                        if CalcMeanS >= EyeOpenLow + EyeOpen/2%EyeOpenLow+EyeOpen/2                        %If it is the uper level the incoming data
                            EoutCorrS = [EoutCorrS 1];                               %is 1
                        else                                                       %If it is the lowest level the incoming data
                            EoutCorrS = [EoutCorrS 0];                               %is 0
                        end
                    end
                end
                
                %PrintInfo(Ploting*50,linspace(0,t(end),length(EoutCorr))...
                %                                           ,TxDataMat(ThisCarr,:),EoutCorr);
                %%       Calculating the Bit Error Ratio (BER)
                %The final process here is to count the number of wrongdoings
                %of this whole process upon the transmited data for
                %quantitative analizes
                if CalcS
                    BitErrS = sum(xor(TxDataMat(ThisCarr,1+SyncPeriod:end-...
                        SyncPeriod),EoutCorrS(1+SyncPeriod:end-SyncPeriod)));%Comparison between the Transmited and received and counting the differences
                    BerOOKS(CurrentTest,ThisCarr) = BitErrS/((Nb-NumBitDesc)-...
                        2*SyncPeriod);%Calculating the ration of wrong bits and the total number of bits transmited
                end
                BitErr = sum(xor(TxDataMat(ThisCarr,1+SyncPeriod:end-...
                    SyncPeriod),EoutCorr(1+SyncPeriod:end-SyncPeriod)));%Comparison between the Transmited and received and counting the differences
                BitErr2 = sum(xor(TxDataMat(ThisCarr,1+SyncPeriod:end-...
                    SyncPeriod),EoutCorr2(1+SyncPeriod:end-SyncPeriod)));%Comparison between the Transmited and received and counting the differences
                BitErrD = sum(xor(TxDataMat(ThisCarr,1+SyncPeriod:end-...
                    SyncPeriod),EoutCorrD(1+SyncPeriod:end-SyncPeriod)));%Comparison between the Transmited and received and counting the differences
                if BitErr2<=BitErr
                    if BitErr2<=BitErrD
                        BerOOK(CurrentTest,ThisCarr) = BitErr2/((Nb-NumBitDesc)-...
                        2*SyncPeriod);%Calculating the ration of wrong bits and the total number of bits transmited
                        AberLev(CurrentTest,ThisCarr)= Inter2(SeqFin2(LocMax2))-Inter2(round(SeqFin2(LocMax2)-MaxVal2));
                        ValsLev(CurrentTest,ThisCarr)= LevDec2;
                    else
                        BerOOK(CurrentTest,ThisCarr) = BitErrD/((Nb-NumBitDesc)-...
                        2*SyncPeriod);%Calculating the ration of wrong bits and the total number of bits transmited
                        AberLev(CurrentTest,ThisCarr)= EyeOpen;
                        ValsLev(CurrentTest,ThisCarr)= EyeOpenLow+EyeOpen/2;
                    end
                else
                    if BitErr<=BitErrD
                        BerOOK(CurrentTest,ThisCarr) = BitErr/((Nb-NumBitDesc)-...
                        2*SyncPeriod);%Calculating the ration of wrong bits and the total number of bits transmited
                    else
                        BerOOK(CurrentTest,ThisCarr) = BitErrD/((Nb-NumBitDesc)-...
                        2*SyncPeriod);%Calculating the ration of wrong bits and the total number of bits transmited
                        AberLev(CurrentTest,ThisCarr)= EyeOpen;
                        ValsLev(CurrentTest,ThisCarr)= EyeOpenLow+EyeOpen/2;
                    end
                end
%                             berpos = 1:2:size(BerOOK,2);
%                             BerOOK(size(BerOOK,1),berpos)
%                             a=1;
                            close all;
%                             BerOOK(CurrentTest,ThisCarr+1) = 1;
                
            end
            %         end
        end
end
if exist('BerDPSK','var')
    if size(BerDPSK,2)<NumCarr
        BerDPSK(1,end+1) = 0;
    end
elseif exist('BerDQPSK','var')
    if size(BerDQPSK,2)<NumCarr
        BerDQPSK(1,end+1) = 0;
    end
elseif exist('Ber4PAM','var')
    if size(Ber4PAM,2)<NumCarr
        Ber4PAM(1,end+1) = 0;
    end
elseif exist('BerOOK','var')
    if size(BerOOK,2)<NumCarr
        BerOOK(1,end+1) = 0;
    end
end
a=0;

% clearvars -except UsedModula TestBundle ThisModula CurrentTest TxDataMat...
% CurrentModula ThisTest ThisModula OfcName BerOOK Ber4PAM BerDQPSK OSNRPC...
% MaxNumStag EoutUpRec BerDPSK CurrentMedium OcsToTest CurrentOCS t f T...
%                                                                   PulseResp

