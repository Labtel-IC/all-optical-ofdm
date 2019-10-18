%c
%c                                                       ..'Ž`'..'Ž`..'Ž`..                                                   
%c       File: DatumReception File For the All Optical OFDM
%c(Main fail to call each idividual parameters)
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
PARDatumUpstrRecInputData;

%%            Recovering the signal
% At this point the FFT will be inplemented. The received signal need to be
% given as an parameter. This following function recursively implement the
% FFT Operation.
% [EoutAux1,~] = OpticalFFT(t,T,MaxNumStag,Count,EoutRec,E0,ActStag);
if FFTSplit
    [EoutAux1,~,VetThisCarr]=OpticalFFTN(t,T,MaxNumStag,EoutRec);
else
    VetThisCarr = 1:NumCarr;
    EoutAux1    = PARSelectEachCarrier(EoutRec,NumCarr,f,fin,12.5e9,5);
end
%%  Ploting some results for qualitative analizes
%PrintInfo(Ploting*12,f,db(abs((fftshift(fft(EoutMod))))/length(EoutMod))...
%                                                              ,EoutAux1,fc);
%%      Recovering the Data
%This is basicaly the final step, so far, in this process. Here, the
%transmited signal will be received and processed to recover the actual
%data within it.
switch Modulation
    case 'DPSK'
        %%                   Receiver DPSK
        parfor ThisCarr=InitCarr:NumCarr                                           %For each carrier the same process of reception will be used.
            if ~mod(ThisCarr,2)
                %%  Reception Process: Handling the income optical field
                %At the first step the income field will be selected from 
                %the optical FFT Output with the help of the maping vector
                %(previously described).
                EoutAux = EoutAux1(VetThisCarr==ThisCarr,:);
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
                %sincronized the recieved signal or the sampling process 
                %will not work.

                %At this first moment a small part of the income signal
                %needs be analized for the sincronism process. As we are 
                %looking for the actual synchronis symbol, which is the 
                %bigest phase shift of the income signal. The next delay 
                %interferometer convert this phase shift to a amplitude 
                %variation.

                %The phase delay is not needed because the optical field 
                %will be analized as a whole as it is composed only by the 
                %real part
                PhaDel          = 0;                                           
                %Remember that at DPSK modulation the information is stored
                %at the phase difference of the current symbol with the 
                %previous one hence this time delay is needed to analyse 
                %each symbol by analizes of its interaction of the current 
                %symbel with the previous one.
                TimDel          = T;
                [ESync1,ESync2] = DelayInterf(t,TimDel,PhaDel,EoutAux);    %Coverting phase shift to amplitude variation

                %The second moment is to transfer this information from the
                %optical domain to the eletrical domain for an eletronic
                %processing. It is done with the help of an photo diode. 
                %The configuration here used is an balanced receiver as the 
                %output of the Delay Interferometer has two signals 
                %resulting from the constructive and destructive signal 
                %interaction. 
                ESync1 = ESync1.*conj(ESync1); 
                ESync2 = ESync2.*conj(ESync2);

                %%           Creating the Reception Filter
                [BitFilt,~] = FiltroGaussiano(f,BWD2,CenFeq2,FiltOrd2);    %Creating filter for selection of the received signal
                BitFilt = fftshift(BitFilt);                               %Shifting the filter for matching the received signal 
                %%
                Esync  = ESync2-ESync1;
                %%         Adding Noise to Received Signal
                %Just after the photo diod is where the noise must be added 
                %as it is the main constrain of the system. Once the signal 
                %is converted from optical to electrical and it digital 
                %values were recevered, there is no point to add noise, 
                %because the information was already received, it just 
                %needs decoding. The noise here added will be gaussian 
                %basically it represents the reception sensibility. In 
                %another words, how much the signal must be about the noise 
                %for an acceptable recovering.
                if ReceptorNoise                                           %Verify whether the noise is to be added or not
                    if SnrRef                                              %Verifing if the referency for the noise is
                        SigPower = MeasPower(Esync);                       %it own received signal or
                    else
                        SigPower = MeasPower(PowRef(ThisCarr))*1e6;        %The pilot carrier ussed to test the channel
                    end
                    SigPower2 = 20*log10(SigPower);
                    Esync = awgn(Esync,CarSNR,SigPower2);
                end

                VetElecPowerF(ThisCarr)= MeasPower(Esync);
                EoutAuxF = 20*log10(abs(fftshift(fft(EoutAux)./length(...
                                                               EoutAux))));
                [VetOptiPowerF(ThisCarr),~]= findpeaks(EoutAuxF,...
                                           'SortStr','descend','NPeaks',1);
                                       
                Esync = ifft(fft(Esync).*BitFilt);                         %Filter is used to remove higher order components

                ESync1 = ESync1./max(abs(ESync1));                         %Normalizing the signal
                ESync2 = ESync2./max(abs(ESync2));                         %Normalizing the signal
                Esync  = Esync./max(abs(Esync));                           %Normalizing the signal
                
                SyncAux   = Esync(IniSyncPos:SyncPos);                     %Selecting just the symbol to synchronize
                SyncedAux = SyncSymb(IniSyncPos:SyncPos);                  %Selecting just the symbol to synchronize
                %%                   Plot for Qualitative analizes
    %             PrintInfo(Ploting*15,t(IniSyncPos:SyncPos),Esync(IniSyncPos:...
    %             SyncPos),SyncSymb(IniSyncPos:SyncPos),ESync1(IniSyncPos:...
    %                                       SyncPos),ESync2(IniSyncPos:SyncPos));

                %%                   Synchronizing
                %This synchronization process is based on the mean expected
                %value. That means, the information mean value should be 
                %within one period of symbol. Thus, the mean value of the 
                %received signal is acquired and compare of the known 
                %sync-word to verify if this mean value is at the right 
                %possition. Which is the midel point (peak) of the highest 
                %level at the sync period
                SyncAux(SyncAux<0)              = 0;                       %To keep the mean value above zero anything under is neglected
                SyncAux(SyncAux>=mean(SyncAux)) = 1;                       %Adding a flag to the first sample of the received mean value
                SyncAux(SyncAux<mean(SyncAux))  = -1;                      %All the others samples at set to the lowest level

                PosToSyn  = find(ismember(SyncAux,1));                     %Finding where is the location of the samples to synchronize
                PosSynced = find(ismember(SyncedAux,1));                   %Finding where is the location of the samples to synchronize

                DiffPos = PosToSyn(round(end/2)) - PosSynced(round(end/2));%Accounting the peak (midel point) displacement
                if DiffPos~=0
                    if DiffPos>0%If the difference is positive, left-shift...
    %                     Esync = ifft(fft(Esync).*exp(1j*2*pi*f*(DiffPos*Ta))); %Shift based on time change     
                        Esync = [Esync(DiffPos+1:end) Esync(1:DiffPos)];   %Shift based on sampling sliding
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
                    %This synchronization process is based on the mean 
                    %expected value. That means, the information mean value 
                    %should be within one period of symbol. Thus, the mean 
                    %value of the received signal is acquired and compare 
                    %of the known sync-word to verify if this mean value is 
                    %at the right possition.

                    SyncAuxEnd   = Esync(end-SyncPos+1:end-IniSyncPos+1);
                    SyncAuxEnd(SyncAuxEnd<0)                 = 0;          %To keep the mean value above zero anything under is neglected
                    SyncAuxEnd(SyncAuxEnd>=mean(SyncAuxEnd)) = 1;          %Adding a flag to the first sample of the received mean value
                    SyncAuxEnd(SyncAuxEnd<mean(SyncAuxEnd))  = -1;         %All the others samples at set to the lowest level

                    PosToSynEnd  = find(ismember(SyncAuxEnd,1));           %Finding where is the location of the first sample to synchronize
                    PosSyncedEnd = find(ismember(SyncedAuxEnd,1));         %Finding where is the location of the first sample to synchronize

                    DiffPosEnd = PosToSynEnd(round(end/2))-PosSyncedEnd(...
                                                             round(end/2));
                     if DiffPosEnd~=0
                        if DiffPosEnd>0%If positive difference, left-shift
    %                         Esync = ifft(fft(Esync).*exp(1j*2*pi*f*(...
    %                                                     DiffPosEnd*Ta)));%Shift based on time change     
                            Esync = [Esync(DiffPosEnd+1:end) Esync(1:...
                                                              DiffPosEnd)];%Shift based on sampling sliding
                        else%but if the difference is negative, right-shift
    %                         Esync = ifft(fft(Esync).*exp(1j*2*pi*f*(...
    %                                                     DiffPosEnd*Ta)));%Shift based on time change     
                            Esync = [Esync(end+DiffPosEnd+1:end) Esync(...
                                                        1:end+DiffPosEnd)];%Shift based on sampling sliding
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
            [~,~,EyeOpen,EyeOpenHigh,EyeOpenLow] = Olho(Esync,T,NPPB,0);     
            
            AberLev1(ThisCarr)  = EyeOpen;
            ValsLev1(ThisCarr)  = EyeOpenLow + EyeOpen/2; 
            
            Txaux1 = rectpulse(TxDataMat(ThisCarr,:),NPPB);
            Txaux1(Txaux1==0) = -1;
           %% Ploting the result for qualitative analizes
            %PrintInfo(Ploting*17,Esync,T,NPPB);
            %PrintInfo(Ploting*18,t(1:length(Esync)),Txaux1,Esync);
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
            EyeSymMat = reshape(Esync(1+SyncPeriod*NPPB:end-SyncPeriod*...
                                           NPPB),NPPB,NbDPSK-2*SyncPeriod);  
            %Then we take the values that compose the decision level 
            %because they will mark the point of symmetry.
            %
            %Firstly it was set the interval in which the histogram will be
            %build. It is based on the number of samples per bit period.
            Interval = linspace(min(Esync(1+SyncPeriod*NPPB:end-...
            SyncPeriod*NPPB)),max(Esync(1+SyncPeriod*NPPB:end-SyncPeriod...
                                                           *NPPB)),2*NPPB);
            %Therefore, the MATLAB hist function returns the number of
            %occurrence of each interval.
            EyeMax = hist(Esync,Interval);
            EyeMaxaux = [0 EyeMax 0];                                      %Zeros are added at the EyeMax to auxiliate the finding peaks process
            [~,EyeLoc] = findpeaks(EyeMaxaux,'MinPeakDistance',NPPB/4,...
                                           'SortStr','descend','NPeaks',2);%The peaks on the Eye profile will be the levels at the Eyes limit
            
%             figure;findpeaks(EyeMaxaux,'MinPeakDistance',NPPB/4,'SortStr','descend','NPeaks',2);
            
            %From the location of the max values that occured, which means
            %the uper and lower level of the eye diagram it needs to take
            %the actual value that those occurences represent that is
            %withing the the Interval variable.
            ValToSeek = Interval(EyeLoc-1);                                
            %The number of ocurrences is a statical measure therefore one
            %does not have control which interval will have the highest
            %peak, thus it is important to ordenate the values to be seek
            %from the lower part of the eye diagram to the uper part of the
            %eye diagram.
            ValToSeek = sort(ValToSeek,'ascend');                          
            OccuCount = zeros(1,size(EyeSymMat,1));                        %Auxiliar Variable for accounting.
            for kk=1:size(EyeSymMat,1)                                     %For every sample within a symbol period
                OccuCount(kk) = OccuCount(kk)+sum((EyeSymMat(kk,:)>=min(...
                          Esync))&(EyeSymMat(kk,:)<=UpeSymPer*EyeOpenLow));%Account all occurencies of the value 1
                OccuCount(kk) = OccuCount(kk)+sum((EyeSymMat(kk,:)>=...
                     LowSymPer*EyeOpenHigh)&(EyeSymMat(kk,:)<=max(Esync)));%Account all occurencies of the value 2
            end
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
            [SymLoc] = round(mean(find(ismember(OccuCount,max(OccuCount)...
                                                                      )))); 
%             figure;findpeaks(OccuCount,'SortStr','descend');
            
            %% Actualy Receiving Data:
            Data = [];                                                     %Initialization of the vector that will store the income data
            for kk=1:NPPB:length(Esync)                                    %The comparison process will be made for each symbol period
%                 MeanOfData = mean(Esync((kk-1)+SymLoc(1)));
                MeanOfData = mean(Esync((kk-1)+NPPB/2));
                if MeanOfData > EyeOpenLow+EyeOpen/2                       %If it is the uper level the incoming data
                    Data = [Data 1];                                       %is 1
                else                                                       %If it is the lowest level the incoming data
                    Data = [Data 0];                                       %is 0
                end
            end
            %%       Calculating the Bit Error Ratio (BER)
            %The final process here is to count the number of wrongdoings
            %of this whole process upon the transmited data for 
            %quantitative analizes
            TxData  = TxDataMat(ThisCarr,1+SyncPeriod:end-SyncPeriod);
            Data    = Data(1+SyncPeriod:end-SyncPeriod);
            
            BitErr = sum(xor(TxData,Data));                                %Comparison between the Transmited and received and counting the differences
            BerDPSK(ThisCarr) = (BitErr)/(NbDPSK-(2*SyncPeriod));          %Calculating the ration of wrong bits and the total number of bits transmited
            
            %% Ploting the result for qualitative analizes
            %PrintInfo(Ploting*19,TxData,Data);
            %%
%             berpos = 1:2:size(BerDPSK,2);
%             BerDPSK(size(BerDPSK,1),berpos)
%             a=6;
%             close all;
            end
        end
    case 'DQPSK'
        %%                   Receiver DQPSK
        parfor ThisCarr=InitCarr:NumCarr                                           %For each carrier the same process of reception will be used.
            if ~mod(ThisCarr,2)
                %%  Reception Process: Handling the income optical field
                %At the first step the income field will be selected from 
                %the optical FFT Output with the help of the maping vector
                %(previously described).
                EoutAux = EoutAux1(VetThisCarr==ThisCarr,:);
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
                %Thus, for the reception process work properly, it is 
                %needed to sincronized the recieved signal or the sampling 
                %process will not work.

                %At this first moment a small part of the income signal 
                %needs be analized for the sincronism process. As we are 
                %not looking for the actual data rather for the synchronis 
                %symbol, which is the bigest phase shift of the income 
                %signal. The next delay interferometer convert this phase 
                %shift to a amplitude variation.

%                 %The phase delay is not needed because the optical field 
%                 %will be analized as a whole.
%                 PhaDel          = 0;                                           
%                 %Remember that at DQPSK modulation the information is 
%                 %stored at the difference of phase between the current and 
%                 %previous symbol hence this time delay is needed.
%                 TimDel          = T;
%                 [ESync1,ESync2] = DelayInterf(t,TimDel,PhaDel,EoutAux);    %Coverting phase shift to amplitude variation
% 
%                 %The second moment is to transfer this information from the
%                 %optical domain to the eletrical domain for an eletronic
%                 %processing. It is done with the help of an photo diode. 
%                 %The configuration here used is an balanced receiver as the 
%                 %output of the Delay Interferometer has two signals 
%                 %resulting from the constructive and destructive signal 
%                 %interaction. 
%                 ESync1 = ESync1.*conj(ESync1);                                 
%                 ESync2 = ESync2.*conj(ESync2);
%                 Esync  = ESync2-ESync1;
% 
%                 ESync1 = ESync1./max(abs(ESync1));                         %Normalizing the signal
%                 ESync2 = ESync2./max(abs(ESync2));                         %Normalizing the signal
%                 Esync  = Esync./max(abs(Esync));                           %Normalizing the signal
%             
%                 if ReceptorNoise
%                     if SnrRef
%                         SigPower = MeasPower(Esync);
%                     else
%                         SigPower = MeasPower(PowRef(ThisCarr))*1e6;
%                     end
%                     SigPower2 = 20*log10(SigPower);
%                     Esync = awgn(Esync,CarSNR,SigPower2);
%                 end
% 
%                 SyncAux   = Esync(IniSyncPos:SyncPos);                     %Selecting just the symbol to synchronize
%                 SyncedAux = SyncSymb(IniSyncPos:SyncPos);                  %Selecting just the symbol to synchronize
%                 %%                   Plot for Qualitative analizes
%     %              PrintInfo(Ploting*22,t(IniSyncPos:SyncPos),Esync(IniSyncPos...
%     %              :SyncPos),SyncSymb(IniSyncPos:SyncPos),ESync1(IniSyncPos:...
%     %                                       SyncPos),ESync2(IniSyncPos:SyncPos));
%                 %%                   Synchronizing
%                 %This synchronization process is based on the mean expected
%                 %value. That means, the information mean value should be 
%                 %within one period of symbol. Thus, the mean value of the 
%                 %received signal is acquired and compare of the known 
%                 %sync-word to verify if this mean value is at the right 
%                 %possition. Which is the midel point (peak) of the highest 
%                 %level at the sync period
%                 SyncAux(SyncAux<0)              = 0;                       %To keep the mean value above zero anything under is neglected
%                 SyncAux(SyncAux>=mean(SyncAux)) = 1;                       %Adding a flag to the first sample of the received mean value
%                 SyncAux(SyncAux<mean(SyncAux))  = -1;                      %All the others samples at set to the lowest level
% 
%                 PosToSyn  = find(ismember(SyncAux,1));                     %Finding where is the location of the samples to synchronize
%                 PosSynced = find(ismember(SyncedAux,1));                   %Finding where is the location of the samples to synchronize
% 
%                 DiffPos = PosToSyn(round(end/2)) - PosSynced(round(end/2));%Accounting the peak (midel point) displacement
%                 if DiffPos>=0%If the difference is positive, left-shift...
%                     EoutAux = ifft(fft(EoutAux).*exp(1j*2*pi*f*(DiffPos*...
%                                                                      Ta)));%Shift based on time change     
%     %                 EoutAux = [EoutAux(DiffPos+1:end) EoutAux(1:DiffPos)];   %Shift based on sampling sliding
%                 else%... but if the difference is negative, right-shift
%                     EoutAux = ifft(fft(EoutAux).*exp(1j*2*pi*f*(DiffPos*...
%                                                                      Ta)));%Shift based on time change     
%     %                 EoutAux = [EoutAux(end+DiffPos+1:end) EoutAux(1:end + ...
%     %                                                                 DiffPos)];%Shift based on sampling sliding
%                 end
% 
%                     %                   Synchronizing
%                 %Because of reasons, sometimes it may be required to make a
%                 %synchronization process with the  end of the data stream 
%                 %as well. This following verification check if the user set 
%                 %(or not) a second synchronization process to be done.
%                 if SencondAdjust
%                     PhaDel          = 0;
%                     TimDel          = T;
%                     [ESync1,ESync2] = DelayInterf(t,TimDel,PhaDel,EoutAux);
% 
%                     ESync1 = ESync1.*conj(ESync1);
%                     ESync2 = ESync2.*conj(ESync2);
%                     Esync  = ESync2-ESync1;
% 
%                     ESync1 = ESync1./max(abs(ESync1)); 
%                     ESync2 = ESync2./max(abs(ESync2)); 
%                     Esync  = Esync./max(abs(Esync)); 
% 
% 
%                     %This synchronization process is based on the mean 
%                     %expected value. That means, the information mean value 
%                     %should be within one period of symbol. Thus, the mean 
%                     %value of the received signal is acquired and compare 
%                     %of the known sync-word to verify if this mean value is 
%                     %at the right possition.
% 
%                     SyncAuxEnd   = Esync(end-SyncPos+1:end-IniSyncPos+1);
% %                     SyncedAuxEnd = SyncSymbEnd(end-SyncPos+1:end-IniSyncPos+1);
%                     SyncAuxEnd(SyncAuxEnd<0)                 = 0;          %To keep the mean value above zero anything under is neglected
%                     SyncAuxEnd(SyncAuxEnd>=mean(SyncAuxEnd)) = 1;          %Adding a flag to the first sample of the received mean value
%                     SyncAuxEnd(SyncAuxEnd<mean(SyncAuxEnd))  = -1;         %All the others samples at set to the lowest level
% 
%                     PosToSynEnd  = find(ismember(SyncAuxEnd,1));           %Finding where is the location of the first sample to synchronize
%                     PosSyncedEnd = find(ismember(SyncedAuxEnd,1));         %Finding where is the location of the first sample to synchronize
% 
%                     DiffPosEnd = PosToSynEnd(round(end/2))-PosSyncedEnd(...
%                                                              round(end/2));
%                     if DiffPosEnd>=0%If positive difference, left-shift...
%                         EoutAux = ifft(fft(EoutAux).*exp(1j*2*pi*f*(...
%                                                           DiffPosEnd*Ta)));%Shift based on time change     
%     %                     EoutAux = [EoutAux(DiffPosEnd+1:end) EoutAux(1:...
%     %                                                             DiffPosEnd)];%Shift based on sampling sliding
%                     else%... but if the difference is negative, right-shift
%                         EoutAux = ifft(fft(EoutAux).*exp(1j*2*pi*f*(...
%                                                           DiffPosEnd*Ta)));%Shift based on time change     
%     %                     EoutAux = [EoutAux(end+DiffPosEnd+1:end) EoutAux(...
%     %                                                       1:end+DiffPosEnd)];%Shift based on sampling sliding
%                     end
%                 end

                %%                   Plot for Qualitative analizes
    %             PrintInfo(Ploting*23,t(IniSyncPos:SyncPos),Esync(IniSyncPos:...
    %                 SyncPos),SyncSymb(IniSyncPos:SyncPos),ESync1(IniSyncPos:...
    %                                       SyncPos),ESync2(IniSyncPos:SyncPos));
                %%               Recovering the Information
                %Once the Income signal was synchronized it is possible to 
                %split its components in phase and in quadrature. It is very
                %important to understand here that those components will 
                %carry the actual data transmited. From it we can't recover 
                %the  I and Q eletrical signals that were used to modulate 
                %our carrier. Those signals get mingled together to 
                %composed the actual data. Therefore, when we pass the 
                %income signalthroughout the receptor we will recover the 
                %TxData generated and transmited. It is important that it 
                %get very clear here because I toke two months to 
                %understand it, no one could explain it well to me and 
                %those articles about it make thisinformation unclear. Thus
                %, to sheed a light on this concept, it is important to 
                %know that the Data is encoded in the I and Q eletrical 
                %components that will get mixed together to modulate the 
                %optical carrier. When this optical signal passes 
                %throughout the MZ-Interferometer the acutal data (TX) will 
                %be recover, not the I and Q component. But the data at the 
                %odd possition (encoded in I) will be at the real component 
                %(in phase), whereas the data at the even possition (
                %encoded in Q) will be at the imaginary component(in
                %quadrature). The MZ-Interferometer will be resposible to 
                %take apart the realand imaginary components of the income 
                %optical field.

                E_rec3 = EoutAux./max(abs(EoutAux));                       %Normalizing income signal
                %For the interferometric process  take in account just the 
                %real component it is needed a phase delay of 45° degrees;
                PhaDel = 1*pi/4;                                
                %Remember that at DQPSK modulation the information is 
                %stored at the difference of phase between the current and 
                %previous symbol hence this time delay of one symbol period 
                %is needed.
                TimDel = T;
                [EoutA,EoutB] = DelayInterf(t,TimDel,PhaDel,EoutAux);      %Coverting phase shift to amplitude variation
                %For the interferometric process  take in account just the 
                %real component it is needed a phase delay of -45° degrees;
                PhaDel = -1*pi/4;
                %Remember that at DQPSK modulation the information is 
                %stored at the difference of phase between the current and 
                %previous symbol hence this time delay of one symbol period 
                %is needed.
                TimDel = T;
                [EoutC,EoutD] = DelayInterf(t,TimDel,PhaDel,EoutAux);      %Coverting phase shift to amplitude variation

                %The second moment is to transfer this information from the
                %optical domain to the eletrical domain for an eletronic
                %processing. It is done with the help of an photo diode.
                EoutA = EoutA.*conj(EoutA);
                EoutB = EoutB.*conj(EoutB);
                EoutC = EoutC.*conj(EoutC);
                EoutD = EoutD.*conj(EoutD);

                %The process with the photo diode is self-coherent, which
                %means the rusult will be a component of the signal 
                %centered in f=0and another component centered at f=2*fc 
                %(frenquecy central). Therefore, to remove the higher order 
                %component a low pass filter will be used.
                %%           Creating the Reception Filter
                [BitFilt,~] = FiltroGaussiano(f,BWD2,CenFeq2,FiltOrd2);    %Creating filter for selection of the received signal
                BitFilt = fftshift(BitFilt);                               %Shifting the filter for matching the received signal 
                %%
                EoutA = ifft(fft(EoutA).*BitFilt);
                EoutB = ifft(fft(EoutB).*BitFilt);
                EoutC = ifft(fft(EoutC).*BitFilt);
                EoutD = ifft(fft(EoutD).*BitFilt);

                %The configuration here used is an balanced receiver as 
                %the output of the Delay Interferometer has two signals 
                %resulting from the constructive and destructive signal 
                %interaction. 
                EoutI = (EoutB - EoutA);
                EoutQ = (EoutD - EoutC);
                %%         Adding Noise to Received Signal
                %Just after the photo diod is where the noise must be added 
                %as it is the main constrain of the system. Once the signal 
                %is converted from optical to electrical and it digital 
                %values were recevered, there is no point to add noise, 
                %because the information was already received, it just 
                %needs decoding. The noise here added will be gaussian 
                %basically it represents the reception sensibility. In 
                %another words, how much the signal must be about the noise 
                %for an acceptable recovering.

                if ReceptorNoise                                           %Verify whether the noise is to be added or not
                    if SnrRef                                              %Verifing if the referency for the noise is
                        SigPowerI = MeasPower(EoutI);                      %it own received signal or
                        SigPowerQ = MeasPower(EoutQ);
                    else
                        SigPowerI = MeasPower(PowRef(ThisCarr))*1e6;       %The pilot carrier ussed to test the channel
                        SigPowerQ = MeasPower(PowRef(ThisCarr))*1e6;
                    end
                    SigPowerI2 = 20*log10(SigPowerI);
                    SigPowerQ2 = 20*log10(SigPowerQ);
                    EoutI = awgn(EoutI,CarSNR,SigPowerI2);
                    EoutQ = awgn(EoutQ,CarSNR,SigPowerQ2);
                end

                VetElecPowerIF(ThisCarr) = MeasPower(EoutI);
                VetElecPowerQF(ThisCarr) = MeasPower(EoutQ);
                EoutAuxF = 20*log10(abs(fftshift(fft(EoutAux)./length(...
                                                               EoutAux))));
                [VetOptiPowerF(ThisCarr),~] = findpeaks(EoutAuxF,...
                                           'SortStr','descend','NPeaks',1);

                EoutI = EoutI./max(abs(EoutI));                            %Normalizing the signal
                EoutQ = EoutQ./max(abs(EoutQ));                            %Normalizing the signal
            
                SyncAuxI   = EoutI(IniSyncPos:SyncPos);                        %Selecting just the symbol to synchronize
                SyncAuxQ   = EoutQ(IniSyncPos:SyncPos);                        %Selecting just the symbol to synchronize
    %             figure;plot(SyncAuxQ);set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                SyncedAux  = SyncSymb(IniSyncPos:SyncPos);                     %Selecting just the symbol to synchronize
                %%                   Synchronizing
                %This synchronization process is based on the mean expected
                %value. That means, the information mean value should be within
                %one period of symbol. Thus, the mean value of the received
                %signal is acquired and compare of the known sync-word to
                %verify if this mean value is at the right possition. Which is
                %the midel point (peak) of the highest level at the sync period
                SyncAuxI(SyncAuxI<0)               = 0;                        %To keep the mean value above zero anything under is neglected
                SyncAuxQ(SyncAuxQ<0)               = 0;                        %To keep the mean value above zero anything under is neglected
                SyncAuxI(SyncAuxI>=mean(SyncAuxI)) = 1;                        %Adding a flag to the first sample of the received mean value
                SyncAuxQ(SyncAuxQ>=mean(SyncAuxQ)) = 1;                        %Adding a flag to the first sample of the received mean value
                SyncAuxI(SyncAuxI<mean(SyncAuxI))  = -1;                       %All the others samples at set to the lowest level
                SyncAuxQ(SyncAuxQ<mean(SyncAuxQ))  = -1;                       %All the others samples at set to the lowest level

                PosToSynI  = find(ismember(SyncAuxI,1));                       %Finding where is the location of the samples to synchronize
                PosToSynQ  = find(ismember(SyncAuxQ,1));                       %Finding where is the location of the samples to synchronize
                PosSynced = find(ismember(SyncedAux,1));                       %Finding where is the location of the samples to synchronize

                DiffPosI = PosToSynI(round(end/2)) - PosSynced(round(end/2));  %Accounting the peak (midel point) displacement
                DiffPosQ = PosToSynQ(round(end/2)) - PosSynced(round(end/2));  %Accounting the peak (midel point) displacement
                if DiffPosI>=0%If the difference is positive, left-shift...
    %                 EoutI = ifft(fft(EoutI).*exp(1j*2*pi*f*(DiffPosI*Ta))); %Shift based on time change     
                    EoutI = [EoutI(DiffPosI+1:end) EoutI(1:DiffPosI)];   %Shift based on sampling sliding    
                    EoutQ = [EoutQ(DiffPosI+1:end) EoutQ(1:DiffPosI)];   %Shift based on sampling sliding
                else%... but if the difference is negative, right-shift
    %                 EoutI = ifft(fft(EoutI).*exp(-1j*2*pi*f*(DiffPosI*Ta))); %Shift based on time change     
                    EoutI = [EoutI(end+DiffPosI+1:end) EoutI(1:end+DiffPosI)]; %Shift based on sampling sliding
                    EoutQ = [EoutQ(end+DiffPosI+1:end) EoutQ(1:end+DiffPosI)]; %Shift based on sampling sliding
                end
    %             if DiffPosQ>=0%If the difference is positive, left-shift...
    % %                 EoutQ = (ifft(fft(EoutQ).*exp(1j*2*pi*f*(DiffPosQ*Ta)))); %Shift based on time change     
    %                 EoutQ = [EoutQ(DiffPosQ+1:end) EoutQ(1:DiffPosQ)];   %Shift based on sampling sliding
    %             else%... but if the difference is negative, right-shift
    % %                 EoutQ = (ifft(fft(EoutQ).*exp(-1j*2*pi*f*(DiffPosQ*Ta)))); %Shift based on time change     
    %                 EoutQ = [EoutQ(end+DiffPosQ+1:end) EoutQ(1:end+DiffPosQ)]; %Shift based on sampling sliding
    %             end
                %%                   Plot for Qualitative analizes
    %              PrintInfo(PlotingThis*22,t(IniSyncPos:SyncPos),EoutI(IniSyncPos...
    %              :SyncPos),SyncSymb(IniSyncPos:SyncPos),EoutA(IniSyncPos:...
    %                                       SyncPos),EoutB(IniSyncPos:SyncPos));
    %              PrintInfo(PlotingThis*22,t(IniSyncPos:SyncPos),EoutQ(IniSyncPos...
    %              :SyncPos),SyncSymb(IniSyncPos:SyncPos),EoutC(IniSyncPos:...
    %                                       SyncPos),EoutD(IniSyncPos:SyncPos));


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
                %One important step on this project was the confirmation of 
                %our results from a closed and theoretical equation that 
                %relates the income optical field witht its respective data 
                %(as descrived above in the receiving process). The result 
                %from this equation was further compared with the result 
                %from the MZ-Interferometer as an proof of concept. This 
                %equation can be found at the book of Optical Fiber 
                %Telecommunications V B, which one of the authors is Ivan 
                %P. Kaminow at the page 144.
                taux = t(1:length(E_rec3));
                faux = time2freq(taux);
                Ui = real(exp(1j*pi/4).*E_rec3.*conj(ifft(fft(E_rec3).*...
                                                   exp(-1j*2*pi*faux*T))));%The data at odd position
                Ui = Ui./max(abs(Ui));                                     %Normalizing the signal
                Uq = imag(exp(1j*pi/4).*E_rec3.*conj(ifft(fft(E_rec3).*...
                                                   exp(-1j*2*pi*faux*T))));%The data at the even position
                Uq = Uq./max(abs(Uq));                                     %Normalizing the signal

                [~,~,EyeOpenI,EyeOpenHighI,EyeOpenLowI] = Olho(EoutI,T,...
                                                                   NPPB,0);
                [~,~,EyeOpenQ,~,EyeOpenLowQ]            = Olho(EoutQ,T,...
                                                                   NPPB,0);

                AberLev1(ThisCarr)  = EyeOpenI;
                ValsLev1(ThisCarr)  = EyeOpenLowI + EyeOpenI/2;
            
                AberLev2(ThisCarr)  = EyeOpenQ;
                ValsLev02(ThisCarr) = EyeOpenLowQ + EyeOpenQ/2;  
                %% Ploting the result for qualitative analizes
                %PrintInfo(Ploting*24,EoutI,T,NPPB);
                %PrintInfo(Ploting*25,Ui,T,NPPB);
                %PrintInfo(Ploting*26,EoutQ,T,NPPB);
                %PrintInfo(Ploting*27,Uq,T,NPPB);
                %%

                Txaux1 = rectpulse(TxDataMat1(ThisCarr,:),NPPB);
                Txaux2 = rectpulse(TxDataMat2(ThisCarr,:),NPPB);
                Txaux1(Txaux1==0) = -1;
                Txaux2(Txaux2==0) = -1;
                %%  Ploting some results for qualitative analizes
    %             PrintInfo(Ploting*28,t(1:length(EoutI)),Txaux1,Txaux2,EoutI,...
    %                           EoutA,EoutB,EoutQ,EoutC,EoutD,real(Ui),real(Uq));
                %%               Recovering the Information
                %After passing the optical signal to the eletrical domain, 
                %for actually detect the data withing the signal the 
                %following steps are needed.
                %
                %Finding Decission Levels:
                %The process for decoding the income signal will be based 
                %on eletronic comparators. Inasmuch as the right decission 
                %level must be acquired for accurately decide, within a 
                %symboleriode, what that current leavel means (ones or 
                %zeros).
                %
                %The process hereafter of chosing the  decission levels is 
                %not deterministic rather it is a statistic process. The 
                %main idea is to take the decission level from the 
                %histogram generated from the income signal stream.
                %
                %This process is realized inside the function Olho.
                %
                %Basicaly the decission level will be the minimal value of 
                %the currente eye under evaluation plus the half of the its 
                %eye opening.The following ilustration better describe this 
                %process
                %
                %Eye Limits:   + (1/2 * Eye Opening:)    Comparisson Limit:
                %
                %UperLevel  ______     ______________     _____
                %                 \   /      |       \   /                    
                %                  \ /       |        \ /                     
                %                   \ Half Eye Opening /   Decission Level
                %                  / \       |        / \                     
                %LowerLevel ______/   \______|_______/   \_____               
                %
                %
                %Actualy Receiving Data:
                %Once the signal was processed the next step is through a
                %comparator decide the actual information received.
                %
                %
                %For accurately read the income signal is necessary more 
                %than just have the decision level, it is also needed to 
                %know where it occurs. For instance, if we decide to use 
                %just to take the decision level and measure the mean value 
                %across and a portion of the time period, which portion 
                %should we take? It would be logical to take the central 
                %portion, but the bit may be deformed in such way that the 
                %information gets concentrated out of the middle part, 
                %which means the symbol is not symmetric. The symmetry 
                %information can be acquired from the eye diagram by 
                %measuring the longitudinal opening. The following sketch 
                %better describes this process:
                %
                %                   Point of Symmetry          
                %         _____     _______|_______     _____
                %              \   /               \   /                    
                %               \ /  Longitudinal   \ /                     
                %                \ ----------------- /     
                %               / \    Opening      / \                     
                %         _____/   \_______________/   \_____               
                %
                %With those two pieces of information, decision level and 
                %point of symmetry, we have the X, Y coordinates for the 
                %centre of the Eye Diagram. Therefore, as long as there is
                %an opening on it it will be possible to recover the 
                %transmitted information without error... theoretically.
                %
                %As this process is also statistical, first we reshape the 
                %income vector to analyze all periods at the same time.
                EyeSymMatI = reshape(EoutI(1+SyncPeriod*NPPB:end-...
                             SyncPeriod*NPPB),NPPB,NbDQPSK/2-2*SyncPeriod);  
                %Then we take the values that compose the decision level 
                %because they will mark the point of symmetry.
                %
                %Firstly it was set the interval in which the histogram 
                %will be build. It is based on the number of samples per 
                %bit period.
                IntervalI = linspace(min(EoutI(1+SyncPeriod*NPPB:end-...
                SyncPeriod*NPPB)),max(EoutI(1+SyncPeriod*NPPB:end-...
                                                 SyncPeriod*NPPB)),2*NPPB);
                %Therefore, the MATLAB hist function returns the number of
                %occurrence of each interval.
                EyeMaxI = hist(EoutI,IntervalI);
                EyeMaxauxI = [0 EyeMaxI 0];                                %Zeros are added at the EyeMax to auxiliate the finding peaks process
                [~,EyeLocI] = findpeaks(EyeMaxauxI,'MinPeakDistance',...
                                    NPPB/4,'SortStr','descend','NPeaks',2);%The peaks on the Eye profile will be the levels at the Eyes limit
                %From the location of the max values that occured, which 
                %means the uper and lower level of the eye diagram it needs 
                %to take the actual value that those occurences represent 
                %that is withing the the Interval variable.
                ValToSeekI = IntervalI(EyeLocI-1);                                
                ValToSeekI = sort(ValToSeekI,'ascend');
                %The number of ocurrences is a statical measure therefore 
                %one does not have control which interval will have the 
                %highest peak, thus it is important to ordenate the values 
                %to be seek from the lower part of the eye diagram to the 
                %uper part of the eye diagram.
                OccuCountI = zeros(1,size(EyeSymMatI,1));                  %Auxiliar Variable for accounting.
                for kk=1:size(EyeSymMatI,1)                                %For every sample within a symbol period

                    OccuCountI(kk) = OccuCountI(kk)+sum((EyeSymMatI(kk,:...
                    )>=min(EoutI))&(EyeSymMatI(kk,:)<=UpeSymPer*...
                                                             EyeOpenLowI));%Account all occurencies of the valeu 1
                    OccuCountI(kk) = OccuCountI(kk)+sum((EyeSymMatI(kk,:...
                    )>=LowSymPer*EyeOpenHighI)&(EyeSymMatI(kk,:)<=max(...
                                                                  EoutI)));%Account all occurencies of the valeu 2
                end
                %The point of symmetry of the eye diagram will be where the 
                %maximum number of occurrences were measured inasmuch as 
                %those  are the points where all the bits go to the center 
                %of the symbol.  From the maximum number of occurrences, it 
                %can happen for more than one sample within one symbol 
                %period, in the best case, all samples would have the same 
                %accounting as it is shown the ilustration above hence the 
                %symmetry will be at the middle sample of this group of 
                %maximum occurrences. This value can be found by the mean 
                %of the samples positions within the symbol period. The 
                %problem with this approach is that the signal must be 
                %synchronized with the maximum displacement of a symbol 
                %period minus 25% of the eye Longitudinal opening if the 
                %displacement is higher than that the point of symmetry 
                %will be wrongly measured.
                [SymLocI] = round(mean(find(ismember(OccuCountI,max(...
                                                           OccuCountI)))));             
    %             [~,SymLocI] = findpeaks(OccuCountI,'SortStr','descend'); %The peak on the Eye profile will be the Symmetry level
                DataOdd = [];                                              %Initialization of the vector that will store the income data
                for kk=1:NPPB:length(EoutI)                                %The comparison process will be made for each symbol period
%                     MeanOfData = mean(EoutI((kk-1)+SymLocI(1)));
                    MeanOfData = mean(EoutI((kk-1)+NPPB/2));
                    if MeanOfData > EyeOpenLowI+EyeOpenI/2                 %If it is the uper level the incoming data
                        DataOdd = [DataOdd 1];                             %is 1
                    else                                                   %If it is the lowest level the incoming data
                        DataOdd = [DataOdd 0];                             %is 0
                    end
                end
                %The identical process just described above will also be used
                %to recover the data at the even positions.
                %
                %As this process is also statistical, first we reshape the 
                %income vector to analyze all periods at the same time.
                EyeSymMatQ = reshape(EoutQ(1+SyncPeriod*NPPB:end-...
                             SyncPeriod*NPPB),NPPB,NbDQPSK/2-2*SyncPeriod);  
                %Then we take the values that compose the decision level 
                %because they will mark the point of symmetry.
                %
                %Firstly it was set the interval in which the histogram will be
                %build. It is based on the number of samples per bit period.
                IntervalQ = linspace(min(EoutQ(1+SyncPeriod*NPPB:end-...
                SyncPeriod*NPPB)),max(EoutQ(1+SyncPeriod*NPPB:end-...
                                                 SyncPeriod*NPPB)),2*NPPB);
                %Therefore, the MATLAB hist function returns the number of
                %occurrence of each interval.
                EyeMaxQ = hist(EoutQ,IntervalQ);
                EyeMaxauxQ = [0 EyeMaxI 0];                                %Zeros are added at the EyeMax to auxiliate the finding peaks process
                [~,EyeLocQ] = findpeaks(EyeMaxauxQ,'MinPeakDistance',...
                                    NPPB/4,'SortStr','descend','NPeaks',2);%The peaks on the Eye profile will be the levels at the Eyes limit
                ValToSeekQ = IntervalQ(EyeLocQ-1);       
                ValToSeekQ = sort(ValToSeekQ,'ascend');                         
                OccuCountQ = zeros(1,size(EyeSymMatQ,1));                  %Auxiliar Variable for accounting.
                for kk=1:size(EyeSymMatQ,1)                                %For every sample within a symbol period
                    OccuCountQ(kk) = OccuCountQ(kk)+sum((EyeSymMatI(kk,...
                    :)>=LowSymPer*ValToSeekQ(1))&(EyeSymMatQ(kk,:)<=...
                                                 UpeSymPer*ValToSeekQ(1)));%Account all occurencies of the valeu 1
                    OccuCountQ(kk) = OccuCountQ(kk)+sum((EyeSymMatQ(kk,...
                    :)>=LowSymPer*ValToSeekQ(2))&(EyeSymMatQ(kk,:)<=...
                                                 UpeSymPer*ValToSeekQ(2)));%Account all occurencies of the valeu 2
                end
    %             [~,SymLocQ] = findpeaks(OccuCountQ,'SortStr','descend'); %The peak on the Eye profile will be the Symmetry level
                [SymLocQ] = round(mean(find(ismember(OccuCountQ,max(...
                                                           OccuCountQ)))));
                %##########################################################
                %######################Important###########################
                %The ber results for the Data Even is not as good as the
                %results of the Data Odd. One possible reason is the 
                %decision oint of symmetry hence, for testing, we change 
                %SymLocQ to SymLocI for evaluation of improvement. It is 
                %expected as bouth signal will the same point of symetry 
                %will perform in an equal way. If it is confirmed the 
                %creation of SymLocQ will be erased.
                %##########################################################
                %##########################################################
                DataEven = [];                                             %Initialization of the vector that will store the income data
                for kk=1:NPPB:length(EoutQ)                                %The comparison process will be made for each symbol period
%                     MeanOfData = mean(EoutQ((kk-1)+SymLocI(1)));
                    MeanOfData = mean(EoutQ((kk-1)+NPPB/2));
                    if MeanOfData > EyeOpenLowQ+EyeOpenQ/2                 %If it is the uper level the incoming data
                        DataEven = [DataEven 1];                           %is 1
                    else                                                   %If it is the lowest level the incoming data
                        DataEven = [DataEven 0];                           %is 0
                    end
                end
                %%       Calculating the Bit Error Ratio (BER)
                %The final process here is to count the number of 
                %wrongdoings of this whole process upon the transmited data 
                %for quantitative analizes
                TxDataOdd  = TxDataMat1(ThisCarr,1+SyncPeriod:length(...
                    TxDataMat1)-SyncPeriod);
                TxDataEven = TxDataMat2(ThisCarr,1+SyncPeriod:length(...
                    TxDataMat2)-SyncPeriod);
                DataOdd  = DataOdd(1+SyncPeriod:length(DataOdd)-...
                                                             SyncPeriod);
                DataEven = DataEven(1+SyncPeriod:length(DataEven)-...
                                                               SyncPeriod);

                BitErrOdd          = sum(xor(TxDataOdd,DataOdd));          %Comparison between the Transmited and received and counting the differences
                BitErrEven         = sum(xor(TxDataEven,DataEven));        %Comparison between the Transmited and received and counting the differences
                BerDQPSK(ThisCarr) = (BitErrOdd+BitErrEven)/((NbDQPSK)-(...
                                                            2*SyncPeriod));%Calculating the ration of wrong bits and the total number of bits transmited
                %% Ploting the result for qualitative analizes
                %PrintInfo(Ploting*29,TxDataOdd,DataOdd);
                %PrintInfo(Ploting*30,TxDataEven,DataEven);
                %%
%                 berpos = 1:2:size(BerDQPSK,2);
%                 BerDQPSK(size(BerDQPSK,1),berpos)
                a=6;
                close all;
            end
        end
    case '4PAM'
        %%              Receiver 4PAM
        for ThisCarr=InitCarr:NumCarr                                           %For each carrier the same process of reception will be used.
            if ~mod(ThisCarr,2)
                %%  Reception Process: Handling the income optical field
                %At the first step the income field will be selected from 
                %the optical FFT Output with the help of the maping vector
                %(previously described).
                EoutAux = EoutAux1(VetThisCarr==ThisCarr,:);     
                %%            Fiber Time Delay Compensation
                switch Medium
                    case 'Fiber'
                        EoutAux = ifft(fft(EoutAux).*exp(1j*2*pi*f*(...
                                                FiberDelay(ThisCarr)*Ta)));
                    otherwise
                end
    %             EoutAux = EoutRec;              
                %The current incoming signal is them converted from the 
                %optical domain to the eletrical domain with the help of an 
                %photo detector.
    %             EoutAux = EoutAux./max(abs(EoutAux));
    %             EoutAux = EoutAux + 1;
                Ix = EoutAux.*conj(EoutAux);

                if ReceptorNoise
                    if SnrRef
                        SigPower = MeasPower(Ix);
                    else
                        SigPower = MeasPower(PowRef(ThisCarr))*1e6;
                    end
                    SigPower2 = 20*log10(SigPower);
                    Ix = awgn(Ix,CarSNR,SigPower2);
                end
                
                VetElecPowerF(ThisCarr)= MeasPower(Ix);
                EoutAuxF = 20*log10(abs(fftshift(fft(EoutAux)./length(...
                                                               EoutAux))));
                [VetOptiPowerF(ThisCarr),~]= findpeaks(EoutAuxF,...
                                           'SortStr','descend','NPeaks',1);
                %The process with the photo diode is self-coherent, which 
                %means the rusult will be a component of the signal 
                %centered in f=0 and another component centered at f=2*fc 
                %(frenquecy central). Therefore, to remove the higher order 
                %component a low pass filter will be used.
                %%           Creating the Reception Filter
                
                [BitFilt,~] = FiltroGaussiano(f,BWD2,CenFeq2,FiltOrd2);    %Creating filter for selection of the received signal
                BitFilt = fftshift(BitFilt);                               %Shifting the filter for matching the received signal 
                
                Ix = ifft(fft(Ix).*BitFilt);
                
                Ix(1:4*(2*NumAmosCP+NPPB)) = Ix(6*(2*NumAmosCP+NPPB));
                
                Ix(end + 1 - 4*(2*NumAmosCP+NPPB):end) = Ix(6*(2*...
                                                          NumAmosCP+NPPB));
                Ix = Ix - min(Ix);                                         %Removing the DC component from them Eletrical signal received

                Ix = 1.*Ix./max(abs(Ix));                                  %Normalizing the eletrical signal (amplifying what is needed)
            
            
                %% Ploting the result for qualitative analizes
    %             PrintInfo(Ploting*31,f,20.*log10(abs(fftshift(fft(Ix)./...
    %                                                             length(Ix)))));
    %             PrintInfo(PlotingThis*32,t,Ix);
    %             PrintInfo(PlotingThis*33,Ix,T,NPPB);
    %             PrintInfo(Ploting*34,t(IniSyncPos:SyncPos),Ix(IniSyncPos:...
    %                                     SyncPos),SyncSymb(IniSyncPos:SyncPos));
                %%        synchronizing
                %For the reception process work properly, it is needed to
                %sincronized the recieved signal or the sampling process 
                %will not work.
                AuxSync = (Ix(IniSyncPos:SyncPos));                        %Selecting the sync-word within the received signal
                AuxSync1 = AuxSync;
                %This synchronization process is based on the mean expected
                %value. That means, the information mean value should be 
                %within one period of symbol. Thus, the mean value of the 
                %received signal is acquired and compare of the known 
                %sync-word to verify if this mean value is at the right 
                %possition
                AuxSync1(AuxSync1<0) = 0;                                  %To keep the mean value above zero anything under is neglected
                AuxSync1(AuxSync1>=mean(AuxSync1)) = 1;                    %Adding a flag to the first sample of the received mean value
                AuxSync1(AuxSync1<mean(AuxSync1)) = -1;                    %All the others samples at set to the lowest level
                AuxSync2 = SyncSymb(IniSyncPos:SyncPos);                   %Selecting the sync-word within the known signal
                AuxSync2(AuxSync2>=mean(AuxSync2)) = 1;                    %Adding a flag to the first sample of the known mean value
                AuxSync2(AuxSync2<mean(AuxSync2)) = -1;                    %All the others samples at set to the lowest level

                PosToSyn  = find(ismember(AuxSync1,1));                    %Finding where is the location of the first sample to synchronize
                PosSyn = find(ismember(AuxSync2,1));                       %Finding where is the location of the first sample to synchronize

                AuxSyncCorr = PosToSyn(round(end/2))-PosSyn(round(end/2));

                %The difference between the PossitionTosynchronize and 
                %Possitionsynchronized will be used to correct the time 
                %shifting on the transmition and reception process.                     
                if AuxSyncCorr>=0%If the difference is positive, left-shift...
    %                 Ix = ifft(fft(Ix).*exp(1j*2*pi*f*AuxSyncCorr));      %Shift based on time change     
                    Ix = [Ix(AuxSyncCorr+1:end) Ix(1:AuxSyncCorr)];        %Shift based on sampling sliding
                else%... but if the difference is negative, right-shift
    %                 Ix = ifft(fft(Ix).*exp(-1j*2*pi*f*AuxSyncCorr));     %Shift based on time change     
                    Ix = [Ix(end+AuxSyncCorr+1:end) Ix(1:end + ...
                                                             AuxSyncCorr)];%Shift based on sampling sliding
                end
                if SencondAdjust
                    %For some reason that we could not understand sometimes 
                    %time (sampling) sliding of the signal is not equal
                    %throught the data stream. Thus, the second part of the
                    %synchronism process will be turn ON or OFF according to 
                    %the user's will.
                    AuxSyncEnd     = (Ix(end-SyncPos+1:end-IniSyncPos-1)); %Selecting the sync-word within the received signal
%                     SyncSymbEndAux = SyncSymbEnd(end-SyncPos+1:end-IniSyncPos-1);%Selecting the sync-word within the known signal
                    AuxSyncEnd1 = AuxSyncEnd;
                    AuxSyncEnd1(AuxSyncEnd1<0) = 0;                        %To keep the mean value above zero anything under is neglected
                    AuxSyncEnd1(AuxSyncEnd1>=mean(AuxSyncEnd1)) = 1;       %Adding a flag to the first sample of the received mean value
                    AuxSyncEnd1(AuxSyncEnd1<mean(AuxSyncEnd1)) = -1;       %All the others samples at set to the lowest level
                    AuxSyncEnd2 = SyncSymbEndAux;
                    AuxSyncEnd2(AuxSyncEnd2>=mean(AuxSyncEnd2)) = 1;       %Adding a flag to the first sample of the known mean value
                    AuxSyncEnd2(AuxSyncEnd2<mean(AuxSyncEnd2)) = -1;       %All the others samples at set to the lowest level


                    PosToSynEnd  = find(ismember(AuxSyncEnd1,1));          %Finding where is the location of the first sample to synchronize
                    PosSynEnd = find(ismember(AuxSyncEnd2,1));             %Finding where is the location of the first sample to synchronize

                    AuxSyncCorrEnd = PosToSynEnd(round(end/2))-PosSynEnd...
                                                            (round(end/2));

                    %The difference between the PossitionTosynchronize and 
                    %Possitionsynchronized will be used to correct the time 
                    %shifting on the transmition and reception process.
                    if AuxSyncCorrEnd>=0%If possitive difference,left-shift
        %                Ix = ifft(fft(Ix).*exp(1j*2*pi*f*AuxSyncCorrEnd));%Shift based on time change     
                        Ix = [Ix(AuxSyncCorrEnd+1:end) Ix(1:...
                                                          AuxSyncCorrEnd)];%Shift based on sampling sliding
                    else%... but if the difference is negative, right-shift
        %               Ix = ifft(fft(Ix).*exp(-1j*2*pi*f*AuxSyncCorrEnd));%Shift based on time change     
                        Ix = [Ix(end+AuxSyncCorrEnd+1:end) Ix(1:end + ...
                                                          AuxSyncCorrEnd)];%Shift based on sampling sliding
                    end
                end

                %% Removing CP
                if AddCP
                    IxAux = Ix(1:end - StuffSampels);
                    IxAux = reshape(IxAux,(2*NumAmosCP+NPPB),Nb4Pam/2);
                    IxAux = IxAux(1+NumAmosCP:end-NumAmosCP,:);
                    IxAux = reshape(IxAux,1,NPPB*Nb4Pam/2);
                    Ix    = IxAux;
                end
                %%  Measuring the EVM
                IxPosAux = NPPB/2:NPPB:length(Ix);
                IxRec    = Ix(IxPosAux);
                [EvmRmsF(ThisCarr),EvmDBF(ThisCarr),EvmPerF(ThisCarr)] =...
                       EvmCalc(EvmMatRef(ThisCarr,1:end-NumBitDesc),IxRec);
                %%          Ploting the result for qualitative analizes
    %             PrintInfo(Ploting*35,t(end-SyncPos+1:end-IniSyncPos+1),Ix(...
    %                 end-SyncPos+1:end-IniSyncPos+1),SyncSymb(end-SyncPos+1:...
    %                                                         end-IniSyncPos+1));
    %             PrintInfo(Ploting*36,Ix,T,NPPB);
                %%         Finding Decission Levels
                %The process for decoding the income signal will be based 
                %oneletronic comparators. Inasmuch as the right decission 
                %level must be acquired for accurately decide, within a 
                %symbol periode, what that current leavel means (ones or 
                %zeros).
                %
                %The process hereafter of chosing the  decission levels is 
                %not deterministic rather it is a statistic process. The 
                %main idea is to take the decission level from the 
                %histogram generated from the income signal stream.
                %
                %Firstly it was set the interval in which the histogram 
                %will be build. It is based on the number of samples per 
                %bit period.
    %             IxToSeek = reshape(Ix(1+SyncPeriod*NPPB:end-SyncPeriod*NPPB...
    %                                        ),NPPB,(Nb4Pam/2) - (2*SyncPeriod));
                Interval = linspace(min(Ix(1+SyncPeriod:end-SyncPeriod))...
                       ,max(Ix(1+SyncPeriod:end-SyncPeriod)),IntervalStep);
                %Therefore, the MATLAB hist function returns the number of
                %occurrence of each interval.
                EyeMax = hist(Ix,Interval);
                EyeMaxaux = [0 EyeMax 0];                                  %Zeros are added at the EyeMax to auxiliate the finding peaks process
                [~,EyeLoc] = findpeaks(EyeMaxaux,'MinPeakDistance',...
                                   MinDist,'SortStr','descend','NPeaks',4);%The peaks on the Eye profile will be the levels at the Eyes limit
                if length(EyeLoc)<4
                    [~,EyeLoc] = findpeaks(EyeMaxaux,'MinPeakDistance',...
                    MinDist/2,'SortStr','descend','NPeaks',4,...
                                             'MinPeakHeight',mean(EyeMax));%The peaks on the Eye profile will be the levels at the Eyes limit
                end
                %This variable will brings the eye profile of the input 
                %signal that will be used to generate the decision level. 
                %Basicaly the decission level will be the minimal value of 
                %the currente eye under evaluation plus the half of the its 
                %eye opening. The following ilustration better describe 
                %this process.
                %
                %Eye Limits:   + (1/2 * Eye Opening:)  =Comparisson Limits:
                %
                %UperLevel  ______     ______________     _____
                %                 \   /      |       \   /                    
                %                  \ /       |        \ /                     
                %                   \ Half Eye Opening /  Decission Level 3
                %                  / \       |        / \                     
                %LowerLevel 3_____/   \______|_______/   \_____               
                %                 \   /      |       \   /                    
                %                  \ /       |        \ /                     
                %                   \ Half Eye Opening /  Decission Level 2
                %                  / \       |        / \                     
                %LowerLevel 2_____/   \______|_______/   \_____              
                %                 \   /      |       \   /                    
                %                  \ /       |        \ /                     
                %                   \ Half Eye Opening /  Decission Level 1
                %                  / \       |        / \                             
                %LowerLevel 1_____/   \______|_______/   \_____             
                %
                %%           Ploting for Qualitative Analizes
    %             PrintInfo(Ploting*37,Interval,EyeMax);
                %%         Finding Decission Levels
                %It is not always possible to totaly recover the signal.
                %Depending on the configuration of the transmition and
                %reception system the eye diagram may be nonexistent. Which
                %means, there will not be a profile to be found therefore 
                %the EyeLoc will not return the correct location. Inasmuch 
                %as the detection process to works limts will be set 
                %accordling to the amplitude of the received signal.
                if length(EyeLoc)<4                                        %If it was not able to find the eye profile.
                    EyeLoc = [2 3 4 5];
                    Levels = [0 0 0.35 0.55 0.85];
                else                                                       %Whereas, if there is an profile the decission can be found
                    Levels = [Interval(EyeLoc(1)-1) Interval(EyeLoc(2)-...
                           1) Interval(EyeLoc(3)-1) Interval(EyeLoc(4)-1)];
                    Levels = sort(Levels);
                end
                %##########################################################
%                 limiar1 = Levels(3)+(Levels(4)-Levels(3))/2;
%                 limiar2 = Levels(2)+(Levels(3)-Levels(2))/2;
%                 limiar3 = Levels(1)+(Levels(2)-Levels(1))/2;
%                 limiarPos1 = Interval>limiar1;
%                 limiarPos2 = (Interval<=limiar1)&(Interval>limiar2);
%                 limiarPos3 = (Interval<=limiar2)&(Interval>limiar3);
%                 limiarPos4 = Interval<=limiar3;
% 
%                 EyeSymMa1 = reshape(Ix(1+SyncPeriod*NPPB:end-SyncPeriod*NPPB...
%                                            ),NPPB,(Nb4Pam/2) - (2*SyncPeriod));
%                 for kk = 1:size(EyeSymMa1,1)
%                     EyeHi1    = find((EyeSymMa1(kk,:)<Levels(4))&(EyeSymMa1(kk,:)>limiar1));
%                     EyeHi2    = find((EyeSymMa1(kk,:)<Levels(3))&(EyeSymMa1(kk,:)>limiar2));
%                     EyeHi3    = find((EyeSymMa1(kk,:)<Levels(2))&(EyeSymMa1(kk,:)>limiar3));
% 
%                     EyeLo1    = find((EyeSymMa1(kk,:)<limiar1)&(EyeSymMa1(kk,:)>Levels(3)));
%                     EyeLo2    = find((EyeSymMa1(kk,:)<limiar2)&(EyeSymMa1(kk,:)>Levels(2)));
%                     EyeLo3    = find((EyeSymMa1(kk,:)<limiar3)&(EyeSymMa1(kk,:)>Levels(1)));
% 
%                     Hi1(kk)   = mean((EyeSymMa1(kk,EyeHi1)));
%                     Hi2(kk)   = mean((EyeSymMa1(kk,EyeHi2)));
%                     Hi3(kk)   = mean((EyeSymMa1(kk,EyeHi3)));
% 
%                     Lo1(kk)   = mean((EyeSymMa1(kk,EyeLo1)));
%                     Lo2(kk)   = mean((EyeSymMa1(kk,EyeLo2)));
%                     Lo3(kk)   = mean((EyeSymMa1(kk,EyeLo3)));
% 
%                     LevHi1(kk)= mean(Hi1) - std(EyeSymMa1(kk,EyeHi1));
%                     LevHi2(kk)= mean(Hi2) - std(EyeSymMa1(kk,EyeHi2));
%                     LevHi3(kk)= mean(Hi3) - std(EyeSymMa1(kk,EyeHi3));
% 
%                     LevLo1(kk)= mean(Lo1) + std(EyeSymMa1(kk,EyeLo1));
%                     LevLo2(kk)= mean(Lo2) + std(EyeSymMa1(kk,EyeLo2));
%                     LevLo3(kk)= mean(Lo3) + std(EyeSymMa1(kk,EyeLo3));
% 
%                     EyeAb1(kk)= LevHi1(kk) - LevLo1(kk);
%                     EyeAb2(kk)= LevHi2(kk) - LevLo2(kk);
%                     EyeAb3(kk)= LevHi3(kk) - LevLo3(kk);
%                 end
%                 EyeAbertura1 = mean(EyeAb1);
%                 EyeAbertura2 = mean(EyeAb2);
%                 EyeAbertura3 = mean(EyeAb3);
%                 ThisDeciLev1 = mean(LevLo1) + EyeAbertura1/2;
%                 ThisDeciLev2 = mean(LevLo2) + EyeAbertura2/2;
%                 ThisDeciLev3 = mean(LevLo3) + EyeAbertura3/2;
%                 a=1;
                %##########################################################
                %For accurately read the income signal is necessary more 
                %than  just have the decision level, it is also needed to 
                %know where it occurs. For instance, if we decide to use 
                %just to take the decision level and measure the mean value 
                %across and a portion of the time period, which portion 
                %should we take? It would be logical to take the central
                %portion, but the bit may be deformed in such way that the 
                %information gets concentrated out of the middle part, 
                %which means the symbol is not symmetric. The symmetry 
                %information can be acquired from the eye diagram by 
                %measuring the longitudinal opening. The following sketch
                % better describes this process:
                %
                %                   Point of Symmetry          
                %         _____     _______|_______     _____
                %              \   /               \   /                    
                %               \ /  Longitudinal   \ /                     
                %                \ ----------------- /     
                %               / \    Opening      / \                     
                %         _____/   \_______________/   \_____               
                %
                %With those two pieces of information, decision level and  
                %point of symmetry, we have the X, Y coordinates for the 
                %centre of the Eye Diagram. Therefore, as long as there is 
                %an opening on it it will be possible to recover the 
                %transmitted information without error, theoretically.

                %As this process is also statistical, first we reshape the 
                %income vector to analyze all periods at the same time.
                EyeSymMat = reshape(Ix(1+SyncPeriod*NPPB:end-SyncPeriod*...
                                   NPPB),NPPB,(Nb4Pam/2) - (2*SyncPeriod));
    %             EyeSymMat = EyeSymMat(1+NPPB/4:end-NPPB/4,:);
                %Then we take the values that compose the decision level 
                %because they will mark the point of symmetry.
                %From the location of the max values that occured, which 
                %means the uper and lower level of the eye diagram it needs 
                %to take the actual value that those occurences represent 
                %that is withing the the Interval variable.
                ValToSeek = Interval(EyeLoc-1);     
                %The number of ocurrences is a statical measure therefore 
                %one does not have control which interval will have the 
                %highest peak, thus it is important to ordenate the values 
                %to be seek from the lower part of the eye diagram to the 
                %uper part of the eye diagram.
                ValToSeek = sort(ValToSeek,'ascend');                           
                OccuCount = zeros(1,size(EyeSymMat,1));                    %Auxiliar Variable for accounting.
                for kk=1:size(EyeSymMat,1)                                 %For every sample within a symbol period
                    if (kk>=NPPB/4)&&(kk<=NPPB-NPPB/4)
                        OccuCount(kk) = OccuCount(kk) + sum((EyeSymMat(...
                        kk,:)>=LowSymPer*ValToSeek(1))&(EyeSymMat(kk,:)...
                                                <=UpeSymPer*ValToSeek(1)));%Account all occurencies of the valeu 1
                        OccuCount(kk) = OccuCount(kk) + sum((EyeSymMat(...
                        kk,:)>=LowSymPer*ValToSeek(2))&(EyeSymMat(kk,:)...
                                                <=UpeSymPer*ValToSeek(2)));%Account all occurencies of the valeu 2
                        OccuCount(kk) = OccuCount(kk) + sum((EyeSymMat(...
                        kk,:)>=LowSymPer*ValToSeek(3))&(EyeSymMat(kk,:)...
                                                <=UpeSymPer*ValToSeek(3)));%Account all occurencies of the valeu 3
                        OccuCount(kk) = OccuCount(kk) + sum((EyeSymMat(...
                        kk,:)>=LowSymPer*ValToSeek(4))&(EyeSymMat(kk,:)...
                                                <=UpeSymPer*ValToSeek(4)));%Account all occurencies of the valeu 4
                    end
                end

                %The point of symmetry of the eye diagram will be where the 
                %maximum number of occurrences were measured inasmuch as 
                %those are the points where all the bits go to the center 
                %of the symbol.  From the maximum number of occurrences, it 
                %can happen for more than one sample within one symbol 
                %period, in the best case, all samples would have the same
                %accounting as it is shown the ilustration above hence the 
                %symmetry will be at the middle sample of this group of 
                %maximum occurrences. This value can be found at the 
                %highest peak founded.
                [~,SymLoc] = findpeaks(OccuCount,'SortStr','descend');     %The peak on the Eye profile will be the Symmetry level

    %             PosIx = SymLoc(1):NPPB:length(Ix);
                PosIx = NPPB/2:NPPB:length(Ix);
                IxAux = Ix(PosIx);
                n=100;

                IxAuxAB = IxAux((IxAux<=Levels(4))&(IxAux>=Levels(3)));
                InterAB = linspace(Levels(3),Levels(4),n);
                EyeAB = hist(IxAuxAB,InterAB);
                EyeAB = ~EyeAB;
                CountAB=1;
                SeqOnesAB=0;
                SeqFinAB = 0;
                SeqIniAB=1;
                for kk=1:length(EyeAB)
                    if EyeAB(kk)
                        SeqOnesAB(SeqIniAB)=CountAB;
                        CountAB = CountAB + 1;
                        if kk==length(EyeAB)
                            SeqFinAB(end + 1) = kk;
                        end
                    else
                        SeqFinAB(SeqIniAB) = kk-1;
                        SeqIniAB = SeqIniAB + 1;
                        CountAB = 1;
                    end
                end
                [MaxValAB,LocMaxAB]=max(SeqOnesAB);
                if LocMaxAB<2 || MaxValAB<2
                    LevDec3 = 0.67;
                    LocMaxAB = 1;
                    SeqFinAB(1)=2;
                    MaxValAB = 0;
                    InterAB(1)=LevDec3;
                else
                    if (SeqFinAB(LocMaxAB)-MaxValAB/2)<1
                        LevDec3 = 0.7;
                    else
                        LevDec3 = InterAB(round(SeqFinAB(LocMaxAB)-...
                                                              MaxValAB/2));
                    end
                end

                LocAB = find(EyeAB);
                if isempty(LocAB)
                    LocAB = 1;
                    LevelDec3 = 0.65;%mean(Levels(3:4));
                else
                    if DecMod%InterAB(LocAB(round(end/2)))<=LevDec3
                        LevelDec3 = LevDec3;
                    else
                        LevelDec3 = InterAB(LocAB(round(end/2)));
                    end
                end

                IxAuxCD = IxAux((IxAux<=Levels(3))&(IxAux>=Levels(2)));
                InterCD = linspace(Levels(2),Levels(3),n);%NPPB*2^n);
                EyeCD = hist(IxAuxCD,InterCD);
                EyeCD = ~EyeCD;
                CountCD=1;
                SeqOnesCD=0;
                SeqFinCD = 0;
                SeqIniCD=1;
                for kk=1:length(EyeCD)
                    if EyeCD(kk)
                        SeqOnesCD(SeqIniCD)=CountCD;
                        CountCD = CountCD + 1;
                        if kk==length(EyeCD)
                            SeqFinCD(end + 1) = kk;
                        end
                    else
                        SeqFinCD(SeqIniCD) = kk-1;
                        SeqIniCD = SeqIniCD + 1;
                        CountCD = 1;
                    end
                end
                [MaxValCD,LocMaxCD]=max(SeqOnesCD);
                if LocMaxCD<2 || MaxValCD<2
                    LevDec2 = 0.45;
                    LocMaxCD = 1;
                    SeqFinCD(1)=2;
                    MaxValCD = 0;
                    InterCD(1)=LevDec2;
                else
                    if (SeqFinCD(LocMaxCD)-MaxValCD/2)<1
                        LevDec2 = 0.5;
                    else
                        LevDec2 = InterCD(round(SeqFinCD(LocMaxCD)-...
                                                              MaxValCD/2));
                    end
                end

    %             EyeCD = ~EyeCD;
                LocCD = find(EyeCD);
                if isempty(LocCD)
                    LocCD = 1;
                    LevelDec2 = 0.35;%mean(Levels(2:3));
                else
                    if DecMod%InterCD(LocCD(round(end/2)))>=LevDec2
                        LevelDec2 = LevDec2;
                    else
                        LevelDec2 = InterCD(LocCD(round(end/2)));
                    end
                end

                IxAuxEF = IxAux((IxAux<=Levels(2))&(IxAux>=Levels(1)));
                InterEF = linspace(Levels(1),Levels(2),n);%NPPB*2^n);
                EyeEF = hist(IxAuxEF,InterEF);
                EyeEF = ~EyeEF;
                CountEF=1;
                SeqOnesEF=0;
                SeqFinEF = 0;
                SeqIniEF=1;
                for kk=1:length(EyeEF)
                    if EyeEF(kk)
                        SeqOnesEF(SeqIniEF)=CountEF;
                        CountEF = CountEF + 1;
                        if kk==length(EyeEF)
                            SeqFinEF(end + 1) = kk;
                        end
                    else
                        SeqFinEF(SeqIniEF) = kk-1;
                        SeqIniEF = SeqIniEF + 1;
                        CountEF = 1;
                    end
                end
                [MaxValEF,LocMaxEF]=max(SeqOnesEF);
                if LocMaxEF<2 || MaxValEF<2
                    LevDec1 = 0.2;
                    LocMaxEF = 1;
                    SeqFinEF(1)=2;
                    MaxValEF = 0;
                    InterEF(1)=LevDec1;
                else
                    if (SeqFinEF(LocMaxEF)-MaxValEF/2)<1
                        LevDec1 = 0.3;
                    else
                        LevDec1 = InterEF(round(SeqFinEF(LocMaxEF)-...
                                                              MaxValEF/2));
                    end
                end

    %             EyeEF = ~EyeEF;
                LocEF = find(EyeEF);
                if isempty(LocEF)
                    LocEF = 1;
                    LevelDec1 = 0.12;%mean(Levels(1:2));
                else
                    if DecMod%InterEF(LocEF(round(end/2)))<=LevDec1
                        LevelDec1 = LevDec1;
                    else
                        LevelDec1 = InterEF(LocEF(round(end/2)));
                    end
                end
                AberLev1(ThisCarr) = abs(InterAB(SeqFinAB(LocMaxAB)-1) - InterAB(SeqFinAB(LocMaxAB)-MaxValAB+1));
                AberLev2(ThisCarr) = abs(InterCD(SeqFinCD(LocMaxCD)-1) - InterCD(SeqFinCD(LocMaxCD)-MaxValCD+1));
                AberLev3(ThisCarr) = abs(InterEF(SeqFinEF(LocMaxEF)-1) - InterEF(SeqFinEF(LocMaxEF)-MaxValEF+1));
                ValsLev1(ThisCarr)  = LevDec3;
                ValsLev02(ThisCarr)  = LevDec2;
                ValsLev3(ThisCarr)  = LevDec1;
                ValsLev21(ThisCarr) = InterAB(LocAB(round(end/2)));
                ValsLev22(ThisCarr) = InterCD(LocCD(round(end/2)));
                ValsLev23(ThisCarr) = InterEF(LocEF(round(end/2)));
                %%           Ploting for Qualitative Analizes
    %             hold all;
    %             plot(t(NPPB/2),LevDec1,'bd');%plot(t(NPPB/2),InterEF(SeqFin(LocMaxEF)-4),'bo');plot(t(NPPB/2),InterEF(round(SeqFin(LocMaxEF)-MaxValEF+4)),'bx');
    %             plot(t(NPPB/2),LevDec2,'gd');%plot(t(NPPB/2),InterCD(SeqFin(LocMaxCD)-4),'go');plot(t(NPPB/2),InterCD(abs(round(SeqFin(LocMaxCD)-MaxValCD+4))),'gx');
    %             plot(t(NPPB/2),LevDec3,'kd');%plot(t(NPPB/2),InterAB(SeqFin(LocMaxAB)-4),'ko');plot(t(NPPB/2),InterAB(round(SeqFin(LocMaxAB)-MaxValAB+4)),'kx');
    %             plot(t(NPPB/2),InterAB(LocAB(round(end/2))),'kx');
    %             plot(t(NPPB/2),InterCD(LocCD(round(end/2))),'gx');
    %             plot(t(NPPB/2),InterEF(LocEF(round(end/2))),'bx');
    %             plot([t(SymLoc(1)) t(SymLoc(1))],[0 1]);
    %             
    %             drawnow;PrintInfo(PlotingThis*38,EyeMax,EyeLoc-1,Interval);
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
                IxRec = [];                                                %Initialization of the vector that will store the income data
                for kk=1:NPPB:length(Ix)                                   %The comparison process will be made for each symbol period
    %                 midaux = round(mean(SymLoc(1:round(end/2))));
                    midaux = NPPB/2;%SymLoc(1);
                    aux1 = Ix((kk-1)+midaux);     %An small portion of the income signal is take for evaluation
                    MeanRec = mean(aux1);                                  %Measuring the avarage value of the samples taken
                    %Verifying the interval for each symbol received. 
                    if MeanRec <= LevelDec1                                %If it is the lowest level the incoming data
                        IxRec = [IxRec 0 0];                               %is 01 (1)
                    elseif (MeanRec <= LevelDec2)&&(MeanRec > LevelDec1)   %If it is the second level the incoming data 
                        IxRec = [IxRec 0 1];                               %is 00 (0)
                    elseif (MeanRec <= LevelDec3)&&(MeanRec > LevelDec2)   %If it is the tird level the incoming data
                        IxRec = [IxRec 1 1];                               %is 10 (2)
                    elseif MeanRec > LevelDec3                             %If it is the uper level the incoming data
                        IxRec = [IxRec 1 0];                               %is 11 (3)
                    else                                                   %If for some misteriose reason neither of previous verification were sucedded
                        IxRec = [IxRec 0 0];                               %by default the current data is set to be 00 (0)
                    end
                end
                %%           Ploting for Qualitative Analizes
    %             PrintInfo(Ploting*41,t(length(TxDataMat(ThisCarr,:))),Nb4Pam...
    %                                            /2,TxDataMat(ThisCarr,:),IxRec);
                %%       Calculating the Bit Error Ratio (BER)
                %The final process here is to count the number of wrongdoings
                %of this whole process upon the transmited data for 
                %quantitative analizes
                BitErr = sum(xor(TxDataMat(ThisCarr,1+2*SyncPeriod:end-...
                    2*SyncPeriod),IxRec(1+2*SyncPeriod:end-2*SyncPeriod)));%Comparison between the Transmited and received and counting the differences
                BitErrAux1 = BitErr;
                if BitErr ~= 0
                    IxRec = [];                                                    %Initialization of the vector that will store the income data
                    for kk=1:NPPB:length(Ix)                                       %The comparison process will be made for each symbol period
                        %                 midaux = round(mean(SymLoc(1:round(end/2))));
                        midaux = NPPB/2;%SymLoc(1);
                        aux1 = Ix((kk-1)+midaux);     %An small portion of the income signal is take for evaluation
                        MeanRec = mean(aux1);                                      %Measuring the avarage value of the samples taken
                        %Verifying the interval for each symbol received.
                        if MeanRec <= DecLevDef1                                    %If it is the lowest level the incoming data
                            IxRec = [IxRec 0 0];                                   %is 01 (1)
                        elseif (MeanRec <= DecLevDef2)&&(MeanRec > DecLevDef1)       %If it is the second level the incoming data
                            IxRec = [IxRec 0 1];                                   %is 00 (0)
                        elseif (MeanRec <= DecLevDef3)&&(MeanRec > DecLevDef2)       %If it is the tird level the incoming data
                            IxRec = [IxRec 1 1];                                   %is 10 (2)
                        elseif MeanRec > DecLevDef3                                 %If it is the uper level the incoming data
                            IxRec = [IxRec 1 0];                                   %is 11 (3)
                        else                                                       %If for some misteriose reason neither of previous verification were sucedded
                            IxRec = [IxRec 0 0];                                   %by default the current data is set to be 00 (0)
                        end
                    end
                    BitErr = [];
                    BitErr = sum(xor(TxDataMat(ThisCarr,1+2*SyncPeriod:end-2*...
                        SyncPeriod),IxRec(1+2*SyncPeriod:end-2*SyncPeriod)));%Comparison between the Transmited and received and counting the differences
                    BitErrAux2 = BitErr;
                    if BitErr ~= 0
                        if BitErrAux1<BitErrAux2
                            BitErr = BitErrAux1;
%                             DecLevDef1 = (DecLevDef1+LevelDec1)/2;
%                             DecLevDef2 = (DecLevDef2+LevelDec2)/2;
%                             DecLevDef3 = (DecLevDef3+LevelDec3)/2;
                        end
                    end
%                 else
%                     DecLevDef1 = LevelDec1;
%                     DecLevDef2 = LevelDec2;
%                     DecLevDef3 = LevelDec3;
                end
                Ber4PAM(ThisCarr) = BitErr/(Nb4Pam-(2*SyncPeriod));%Calculating the ration of wrong bits and the total number of bits transmited
                
%                 close all;
            end
        end
    otherwise
        %%              Receiver OOK
        parfor ThisCarr=InitCarr:NumCarr                                           %For each carrier the same process of reception will be used.
            if ~mod(ThisCarr,2)
                %%  Reception Process: Handling the income optical field
                %At the first step the income field will be selected from 
                %the optical FFT Output with the help of the maping vector
                %(previously described).
                EoutAux = EoutAux1(VetThisCarr==ThisCarr,:);
                %%            Fiber Time Delay Compensation
                switch Medium
                    case 'Fiber'
                        EoutAux = ifft(fft(EoutAux).*exp(1j*2*pi*f*(...
                                                FiberDelay(ThisCarr)*Ta)));
                    otherwise
                end

                %The current incoming signal is them converted from the 
                %optical domain to the eletrical domain with the help of an 
                %photo detector.
                Ix =EoutAux.*conj(EoutAux);
                if ReceptorNoise
                    if SnrRef
                        SigPower = MeasPower(Ix);
                    else
                        SigPower = MeasPower(PowRef(ThisCarr))*1e6;
                    end
                    SigPower2 = 20*log10(SigPower);
                    Ix = awgn(Ix,CarSNR,SigPower2);
                end
                
                VetElecPowerF(ThisCarr)= MeasPower(Ix);
                EoutAuxF = 20*log10(abs(fftshift(fft(EoutAux)./length(...
                                                               EoutAux))));
                [VetOptiPowerF(ThisCarr),~]= findpeaks(EoutAuxF,...
                                           'SortStr','descend','NPeaks',1);
                %%           Creating the Reception Filter
                [BitFilt,~] = FiltroGaussiano(f,BWD,CenFeq,FiltOrd);       %Creating filter for selection of the received signal
                BitFilt = fftshift(BitFilt);                               %Shifting the filter for matching the received signal 

                Ix = ifft(fft(Ix).*BitFilt);                               %Filtering the signal for better detection
                if NRZPolarity
                    Ix = Ix - mean(Ix);                                    %Removing any off-set that may exist
                else
                    Ix = Ix - min(Ix);                                     %Removing any off-set that may exist
                end 
                Ix = Ix./max(abs(Ix));                                     %Normalizing the eletrical signal (amplifying what is needed)
            
                %% Ploting the result for qualitative analizes
    %             PrintInfo(Ploting*42,f,20.*log10(abs(fftshift(fft(Ix)./...
    %                                                             length(Ix)))));
    %             PrintInfo(Ploting*43,t,TxSigMat(ThisCarr,:)./max(TxSigMat(ThisCarr,:)),Ix);
    %             PrintInfo(Ploting*44,Ix,T,NPPB);
    %             PrintInfo(Ploting*45,t(IniSyncPos:SyncPos),Ix(IniSyncPos:...
    %                                     SyncPos),SyncSymb(IniSyncPos:SyncPos));
                %%        synchronizing
                %For the reception process work properly, it is needed to
                %sincronized the recieved signal or the sampling process 
                %will not work.

                AuxSync = (Ix(IniSyncPos:SyncPos));                        %Selecting the sync-word within the received signal
                AuxSync1 = AuxSync;
                %This synchronization process is based on the mean expected
                %value. That means, the information mean value should be 
                %within one period of symbol. Thus, the mean value of the 
                %received signal is acquired and compare of the known 
                %sync-word to verify if this mean value is at the right 
                %possition
                AuxSync1(AuxSync1<0) = 0;                                  %To keep the mean value above zero anything under is neglected
                AuxSync1(AuxSync1>=mean(AuxSync1)) = 1;                    %Adding a flag to the first sample of the received mean value
                AuxSync1(AuxSync1<mean(AuxSync1)) = -1;                    %All the others samples at set to the lowest level
                AuxSync2 = SyncSymb(IniSyncPos:SyncPos);                   %Selecting the sync-word within the known signal
                AuxSync2(AuxSync2>=mean(AuxSync2)) = 1;                    %Adding a flag to the first sample of the known mean value
                AuxSync2(AuxSync2<mean(AuxSync2)) = -1;                    %All the others samples at set to the lowest level

                PosToSyn  = find(ismember(AuxSync1,1));                  
                PosSyn    = find(ismember(AuxSync2,1));               

                AuxSyncCorr = PosToSyn(round(end/2))-PosSyn(round(end/2));
                %The difference between the PossitionTosynchronize and 
                %Possitionsynchronized will be used to correct the time 
                %shifting on the transmition and reception process.

                if AuxSyncCorr>=0%If the difference is positive, left-shift...
    %                 Ix = ifft(fft(Ix).*exp(1j*2*pi*f*(AuxSyncCorr*Ta))); %Shift based on time change  
                    Ix = [Ix(AuxSyncCorr+1:end) Ix(1:AuxSyncCorr)];        %Shift based on sampling sliding
                else%... but if the difference is negative, right-shift
    %                 Ix = ifft(fft(Ix).*exp(-1j*2*pi*f*(AuxSyncCorr*Ta)));%Shift based on time change  
                    Ix = [Ix(end+AuxSyncCorr+1:end) Ix(1:end + ...
                                                             AuxSyncCorr)];%Shift based on sampling sliding
                end


                %%        synchronizing
                if SencondAdjust
                    %For some reason that we could not understand sometimes 
                    %the time (sampling) sliding of the signal is not equal
                    %throught the data stream. Thus, the second part of the
                    %synchronism process will be turn ON or OFF according 
                    %to the user's will.
                    AuxSyncEnd = (Ix(end-(SyncPos+1*abs(AuxSyncCorr))+1:...
                                         end-(2*NPPB+1*abs(AuxSyncCorr))));%Selecting the sync-word within the received signal
                    AuxSyncEnd1 = AuxSyncEnd;
                    AuxSyncEnd1(AuxSyncEnd1<0) = 0;                        %To keep the mean value above zero anything under is neglected
                    AuxSyncEnd1(AuxSyncEnd1>=mean(AuxSyncEnd1)) = 1;       %Adding a flag to the first sample of the received mean value
                    AuxSyncEnd1(AuxSyncEnd1<mean(AuxSyncEnd1)) = -1;       %All the others samples at set to the lowest level
                    AuxSyncEnd2 = SyncSymbEndAux;
                    AuxSyncEnd2(AuxSyncEnd2>=mean(AuxSyncEnd2)) = 1;       %Adding a flag to the first sample of the known mean value
                    AuxSyncEnd2(AuxSyncEnd2<mean(AuxSyncEnd2)) = -1;       %All the others samples at set to the lowest level

                    PosToSynEnd  = find(ismember(AuxSyncEnd1,1));              
                    PosSynEnd = find(ismember(AuxSyncEnd2,1));                

                    AuxSyncCorrEnd = PosToSynEnd(round(end/2))-PosSynEnd...
                                                            (round(end/2));

                    %The difference between the PossitionTosynchronize and 
                    %Possitionsynchronized will be used to correct the time 
                    %shifting on the transmition and reception process.

                    if AuxSyncCorrEnd>=0%If possitive difference,left-shift
    %                     Ix = ifft(fft(Ix).*exp(1j*2*pi*f*(AuxSyncCorrEnd*...
    %                                                                Ta)));%Shift based on time change  
                        Ix = [Ix(AuxSyncCorrEnd+1:end) Ix(1:...
                                                          AuxSyncCorrEnd)];
                    else%... but if the difference is negative, right-shift
    %                     Ix = ifft(fft(Ix).*exp(-1j*2*pi*f*(AuxSyncCorrEnd*...
    %                                                                Ta)));%Shift based on time change  
                        Ix = [Ix(end+AuxSyncCorrEnd+1:end) Ix(1:end + ...
                                                          AuxSyncCorrEnd)];%Shift based on sampling sliding
                    end
                end

                %% Ploting the result for qualitative analizes
    %             PrintInfo(Ploting*46,t(1:length(Ix)),TxSigMat(ThisCarr,:)./max(TxSigMat(ThisCarr,:)),Ix);
                %%                  Removing CP
                if AddCP
                    IxAux = Ix(1:end - StuffSampels);
                    IxAux = reshape(IxAux,(2*NumAmosCP+NPPB),Nb-...
                                                               NumBitDesc);
                    IxAux = IxAux(1+NumAmosCP:end-NumAmosCP,:);
                    IxAux = reshape(IxAux,1,NPPB*(Nb-NumBitDesc));
                    Ix    = IxAux;
                end
                %%  Measuring the EVM
                IxPosAux = NPPB/2:NPPB:length(Ix);
                IxRec    = Ix(IxPosAux);
                [EvmRmsF(ThisCarr),EvmDBF(ThisCarr),EvmPerF(ThisCarr)] =...
                     EvmCalc(EvmMatRef(ThisCarr,1:end - NumBitDesc),IxRec);
                %% Ploting the result for qualitative analizes

    %             PrintInfo(Ploting*47,(Ix),T,NPPB);
                [~,~,EyeOpen,EyeOpenHigh,EyeOpenLow] = Olho((Ix),T,NPPB,0);
                AberLev01(ThisCarr) = EyeOpen;
                ValsLev01(ThisCarr)  = EyeOpenLow + EyeOpen/2;
                %%         Finding Decission Levels
                %The process for decoding the income signal will be based 
                %on eletronic comparators. Inasmuch as the right decission 
                %level must be acquired for accurately decide, within a 
                %symbol periode, what that current leavel means (ones or 
                %zeros).
                %
                %The process hereafter of chosing the  decission levels is 
                %not deterministic rather it is a statistic process. The 
                %main idea is to take the decission level from the 
                %histogram generated from the income signal stream.
                %
                %This process is realized inside the function Olho.
                %
                %Basicaly the decission level will be the minimal value of 
                %the currente eye under evaluation plus the half of the its 
                %eye opening.The following ilustration better describe this 
                %process
                %
                %Eye Limits:   + (1/2 * Eye Opening:)  =Comparisson Limit:
                %
                %UperLevel  ______     ______________     _____
                %                 \   /      |       \   /                    
                %                  \ /       |        \ /                     
                %                   \ Half Eye Opening /  Decission Level
                %                  / \       |        / \                     
                %LowerLevel ______/   \______|_______/   \_____               
                %
                %
                %
                %For accurately read the income signal is necessary more 
                %than just have the decision level, it is also needed to 
                %know where it occurs. For instance, if we decide to use 
                %just to take the decision level and measure the mean value 
                %across and a portion of the time period, which portion 
                %should we take? It would be logical to take the central 
                %portion, but the bit may be deformed in such way that the 
                %information gets concentrated out of the middle part, 
                %which means the symbol is not symmetric. The symmetry 
                %information can be acquired from the eye diagram by 
                %measuring the longitudinal opening. The following sketch 
                %better describes this process:
                %
                %                   Point of Symmetry          
                %         _____     _______|_______     _____
                %              \   /               \   /                    
                %               \ /  Longitudinal   \ /                     
                %                \ ----------------- /     
                %               / \    Opening      / \                     
                %         _____/   \_______________/   \_____               
                %
                %With those two pieces of information, decision level and 
                %point of symmetry, we have the X, Y coordinates for the 
                %centre of the Eye Diagram. Therefore, as long as there is 
                %an opening on it it will be possible to recover the 
                %transmitted information without error... theoretically.
                %
                %As this process is also statistical, first we reshape the 
                %income vector to analyze all periods at the same time.

                EyeSymMat = reshape(Ix(1+SyncPeriod*NPPB:end-SyncPeriod*...
                                      NPPB),NPPB,(Nb-NumBitDesc)-2*SyncPeriod);  
                %Then we take the values that compose the decision level 
                %because they will mark the point of symmetry.
                %
                %Firstly it was set the interval in which the histogram 
                %will be build. It is based on the number of samples per 
                %bit period.
                Interval = linspace(min(Ix),max(Ix),2*NPPB);
                %Therefore, the MATLAB hist function returns the number of
                %occurrence of each interval.
                EyeMax = hist(Ix,Interval);
                EyeMaxaux = [0 EyeMax 0];                                  %Zeros are added at the EyeMax to auxiliate the finding peaks process
                [~,EyeLoc] = findpeaks(EyeMaxaux,'MinPeakDistance',NPPB/...
                                         4,'SortStr','descend','NPeaks',2);%The peaks on the Eye profile will be the levels at the Eyes limit


                ValToSeek = Interval(EyeLoc-1);             
                %From the location of the max values that occured, which 
                %means the uper and lower level of the eye diagram it needs 
                %to take the actual value that those occurences represent 
                %that is withing the the Interval variable.                   
                ValToSeek = sort(ValToSeek,'ascend');
                %The number of ocurrences is a statical measure therefore 
                %one does not have control which interval will have the 
                %highest peak, thus it is important to ordenate the values 
                %to be seek from the lower part of the eye diagram to the 
                %uper part of the eye diagram.
                OccuCount = zeros(1,size(EyeSymMat,1));                    %Auxiliar Variable for accounting.
                for kk=1:size(EyeSymMat,1)                                 %For every sample within a symbol period
                    OccuCount(kk) = OccuCount(kk)+sum((EyeSymMat(kk,:)>=...
                         min(Ix))&(EyeSymMat(kk,:)<=UpeSymPer*EyeOpenLow));%Account all occurencies of the value 1
                    OccuCount(kk) = OccuCount(kk)+sum((EyeSymMat(kk,:)>=...
                        LowSymPer*EyeOpenHigh)&(EyeSymMat(kk,:)<=max(Ix)));%Account all occurencies of the value 2
                end
                %The point of symmetry of the eye diagram will be where the 
                %maximum number of occurrences were measured inasmuch as 
                %those are the points where all the bits go to the center 
                %of the symbol.  From the maximum number of occurrences, it 
                %can happen for more than one sample within one symbol 
                %period, in the best case, all samples would have the same 
                %accounting as it is shown the ilustration above hence the 
                %symmetry will be at the middle sample of this group of 
                %maximum occurrences. This value can be found by the mean 
                %of the samples positions within the symbol period. The 
                %problem with this approach is that the signal must be 
                %synchronized with the maximum displacement of a symbol 
                %period minus 25% of the eye Longitudinal opening if the 
                %displacement is higher than that the point of symmetry 
                %will be wrongly measured.
                [SymLoc] = round(mean(find(ismember(OccuCount,max(...
                                                            OccuCount)))));%The peak on the Eye profile will be the Symmetry level
                %%           Ploting for Qualitative Analizes
    %             PrintInfo(Ploting*48,NPPB,OccuCount);
    %             PrintInfo(Ploting*49,OccuCount,SymLoc);

                %%      Actualy Receiving Data
                %Once the signal was processed the next step is through a
                %comparator decide the actual information received.
                EoutCorr = [];                                             %Initialization of the vector that will store the income data
                for kk=1:NPPB:length(Ix)                                   %The comparison process will be made for each symbol period
                    %An small portion of the income signal is take for 
                    %evaluation by measuring the avarage value of the 
                    %samples taken
%                     CalcMean = mean((Ix((kk-1)+SymLoc(1))));
                    CalcMean = mean((Ix((kk-1)+ NPPB/2)));
                    %Verifying the interval for each symbol received. 
                    if CalcMean >= EyeOpenLow+EyeOpen/2                    %If it is the uper level the incoming data
                        EoutCorr = [EoutCorr 1];                           %is 1
                    else                                                   %If it is the lowest level the incoming data
                        EoutCorr = [EoutCorr 0];                           %is 0
                    end
                end

                %PrintInfo(Ploting*50,linspace(0,t(end),length(EoutCorr))...
    %                                           ,TxDataMat(ThisCarr,:),EoutCorr);
                %%       Calculating the Bit Error Ratio (BER)
                %The final process here is to count the number of 
                %wrongdoings of this whole process upon the transmited data 
                %for quantitative analizes
                switch UsedModula
                    case 1
                        BitErr = sum(xor(TxDataMat(ThisCarr,1+SyncPeriod:end-...
                               SyncPeriod),EoutCorr(1+SyncPeriod:end-SyncPeriod)));%Comparison between the Transmited and received and counting the differences
                        BerOOK(ThisCarr) = BitErr/((Nb-NumBitDesc)-2*SyncPeriod);  %Calculating the ration of wrong bits and the total number of bits transmited
                    case 2
                        BitErr = sum(xor(TxDataMat(ThisCarr,1+SyncPeriod:end-...
                               SyncPeriod),EoutCorr(1+SyncPeriod:end-SyncPeriod)));%Comparison between the Transmited and received and counting the differences
                        Ber4PAM(ThisCarr) = BitErr/((Nb-NumBitDesc)-2*SyncPeriod);  %Calculating the ration of wrong bits and the total number of bits transmited
                    case 3
                        BitErr = sum(xor(TxDataMat(ThisCarr,1+SyncPeriod:end-...
                               SyncPeriod),EoutCorr(1+SyncPeriod:end-SyncPeriod)));%Comparison between the Transmited and received and counting the differences
                        BerDQPSK(ThisCarr) = BitErr/((Nb-NumBitDesc)-2*SyncPeriod);  %Calculating the ration of wrong bits and the total number of bits transmited
                    case 4
                        BitErr = sum(xor(TxDataMat(ThisCarr,1+SyncPeriod:end-...
                               SyncPeriod),EoutCorr(1+SyncPeriod:end-SyncPeriod)));%Comparison between the Transmited and received and counting the differences
                        BerDPSK(ThisCarr) = BitErr/((Nb-NumBitDesc)-2*SyncPeriod);  %Calculating the ration of wrong bits and the total number of bits transmited
                    otherwise
                        BitErr = sum(xor(TxDataMat(ThisCarr,1+SyncPeriod:end-...
                               SyncPeriod),EoutCorr(1+SyncPeriod:end-SyncPeriod)));%Comparison between the Transmited and received and counting the differences
                        BerOOK(ThisCarr) = BitErr/((Nb-NumBitDesc)-2*SyncPeriod);  %Calculating the ration of wrong bits and the total number of bits transmited
                end
                
%                 berpos = 1:2:size(BerOOK,2);
%                 BerOOK(size(BerOOK,1),berpos)
%                 a=1;
                close all;
            end
        end
end
a=0;

clear EoutAux1
% clearvars -except UsedModula TestBundle ThisModula CurrentTest TxDataMat...
% CurrentModula ThisTest ThisModula OfcName BerOOK Ber4PAM BerDQPSK OSNRPC...
% MaxNumStag EoutUpRec BerDPSK CurrentMedium OcsToTest CurrentOCS t f T...
%                                                                   PulseResp

