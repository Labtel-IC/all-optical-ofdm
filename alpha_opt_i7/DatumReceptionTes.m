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
DatumReceptionInputData;                                                   % Configuration File of DatumReception.

%%            Recovering the signal
% At this point the FFT will be inplemented. The received signal need to be
% given as an parameter. This following function recursively implement the
% FFT Operation.
% [EoutAux1,~] = OpticalFFT(t,T,MaxNumStag,Count,EoutRec,E0,ActStag);
[EoutAux1,~,VetThisCarr]=OpticalFFTN(t,T,MaxNumStag,EoutRec);
%%  Ploting some results for qualitative analizes
PrintInfo(Ploting*12,f,db(abs((fftshift(fft(EoutMod))))/length(EoutMod))...
                                                             ,EoutAux1,fc);
%%      Recovering the Data
%This is basicaly the final step, so far, in this process. Here, the
%transmited signal will be received and processed to recover the actual
%data within it.
switch Modulation
    case 'DPSK'
        if exist('CurrentOCS','var')
            BerDPSK(1,:,CurrentOCS) = OSNRPC(1:NumCarr);
        end
        %%                   Receiver DPSK
        for ThisCarr=1:2:NumCarr                                           %For each carrier the same process of reception will be used.
            %%           Creating the Reception Filter
            [BitFilt,~] = FiltroGaussiano(f,BWD2,CenFeq2,FiltOrd2);        %Creating filter for selection of the received signal
            BitFilt = fftshift(BitFilt);                                   %Shifting the filter for matching the received signal 
            %%  Reception Process: Handling the income optical field
            %At the first step the income field will be selected from the
            %optical FFT Output with the help of the maping vector
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
            PrintInfo(Ploting*13,t,abs(EoutAux));
            PrintInfo(Ploting*14,f,20*log10(abs(fftshift(fft(EoutAux)./...
                                                       length(EoutAux)))));
            
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
            Esync  = ESync2-ESync1;
            Esync = ifft(fft(Esync).*BitFilt);                             %Filter is used to remove higher order components
            
            ESync1 = ESync1./max(abs(ESync1));                             %Normalizing the signal
            ESync2 = ESync2./max(abs(ESync2));                             %Normalizing the signal
            Esync  = Esync./max(abs(Esync));                               %Normalizing the signal
            
            SyncAux   = Esync(IniSyncPos:SyncPos);                         %Selecting just the symbol to synchronize
            SyncedAux = SyncSymb(IniSyncPos:SyncPos);                      %Selecting just the symbol to synchronize
            %%                   Plot for Qualitative analizes
            PrintInfo(PlotingThis*15,t(IniSyncPos:SyncPos),Esync(IniSyncPos:...
            SyncPos),SyncSymb(IniSyncPos:SyncPos),ESync1(IniSyncPos:...
                                      SyncPos),ESync2(IniSyncPos:SyncPos));

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
            PrintInfo(Ploting*16,t(end-SyncPos+1:end-IniSyncPos+1),Esync...
            (end-SyncPos+1:end-IniSyncPos+1),SyncSymb(end-SyncPos+1:end-...
            IniSyncPos+1),ESync1(end-SyncPos+1:end-IniSyncPos+1),ESync2(...
                                          end-SyncPos+1:end-IniSyncPos+1));
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
                        Esync = [Esync(end+DiffPosEnd+1:end) Esync(...
                                                        1:end+DiffPosEnd)];%Shift based on sampling sliding
                    end
                 end
            end
           %% Ploting the result for qualitative analizes
            PrintInfo(PlotingThis*17,Esync,T,NPPB);
            
            [~,~,EyeOpen,EyeOpenHigh,EyeOpenLow] = Olho(Esync,T,NPPB,0);            
            Txaux1 = rectpulse(TxDataMat(ThisCarr,:),NPPB);
            Txaux1(Txaux1==0) = -1;
            
            PrintInfo(Ploting*18,t,Txaux1,Esync);
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
                                               NPPB),NPPB,Nb-2*SyncPeriod);  
            %Then we take the values that compose the decision level 
            %because they will mark the point of symmetry.
            %
            %Firstly it was set the interval in which the histogram will be
            %build. It is based on the number of samples per bit period.
            Interval = linspace(min(Esync),max(Esync),2*NPPB);
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
                          Esync))&(EyeSymMat(kk,:)<=UpeSymPer*EyeOpenLow)); %Account all occurencies of the value 1
                OccuCount(kk) = OccuCount(kk)+sum((EyeSymMat(kk,:)>=...
                     LowSymPer*EyeOpenHigh)&(EyeSymMat(kk,:)<=max(Esync))); %Account all occurencies of the value 2
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
                MeanOfData = mean(Esync((kk-1)+SymLoc(1)));
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
            if exist('CurrentOCS','var')
                BerDPSK(CurrentTest+1,ThisCarr,CurrentOCS) = (BitErr)/((...
                                                       Nb)-(2*SyncPeriod));
            else
                BerDPSK(CurrentTest,ThisCarr) = (BitErr)/((Nb)-(2*...
                                                              SyncPeriod));%Calculating the ration of wrong bits and the total number of bits transmited
            end
            %% Ploting the result for qualitative analizes
            PrintInfo(Ploting*19,TxData,Data);
            %%
            a=6;
            close all;
            
        end
    case 'DQPSK'
        
        if exist('CurrentOCS','var')
            BerDQPSK(1,:,CurrentOCS) = OSNRPC(1:NumCarr);
        end
        %%                   Receiver DQPSK
        for ThisCarr=1:2:NumCarr                                           %For each carrier the same process of reception will be used.
            %%           Creating the Reception Filter
            [BitFilt,~] = FiltroGaussiano(f,BWD2,CenFeq2,FiltOrd2);        %Creating filter for selection of the received signal
            BitFilt = fftshift(BitFilt);                                   %Shifting the filter for matching the received signal 
            %%  Reception Process: Handling the income optical field
            %At the first step the income field will be selected from the
            %optical FFT Output with the help of the maping vector
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
            PrintInfo(Ploting*20,t,abs(EoutAux));
            PrintInfo(Ploting*21,f,20*log10(abs(fftshift(fft(EoutAux)./...
                                                       length(EoutAux)))));
            
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
            
            %The phase delay is not needed because the optical field will 
            %be analized as a whole.
            PhaDel          = 0;                                           
            %Remember that at DQPSK modulation the information is stored at
            %the difference of phase between the current and previous
            %symbol hence this time delay is needed.
            TimDel          = T;
            [ESync1,ESync2] = DelayInterf(t,TimDel,PhaDel,EoutAux);        %Coverting phase shift to amplitude variation
            
            %The second moment is to transfer this information from the
            %optical domain to the eletrical domain for an eletronic
            %processing. It is done with the help of an photo diode. The
            %configuration here used is an balanced receiver as the output
            %of the Delay Interferometer has two signals resulting from
            %the constructive and destructive signal interaction. 
            ESync1 = ESync1.*conj(ESync1);                                 
            ESync2 = ESync2.*conj(ESync2);
            Esync  = ESync2-ESync1;
            
            ESync1 = ESync1./max(abs(ESync1));                             %Normalizing the signal
            ESync2 = ESync2./max(abs(ESync2));                             %Normalizing the signal
            Esync  = Esync./max(abs(Esync));                               %Normalizing the signal
            
            SyncAux   = Esync(IniSyncPos:SyncPos);                         %Selecting just the symbol to synchronize
            SyncedAux = SyncSymb(IniSyncPos:SyncPos);                      %Selecting just the symbol to synchronize
            %%                   Plot for Qualitative analizes
             PrintInfo(Ploting*22,t(IniSyncPos:SyncPos),Esync(IniSyncPos...
             :SyncPos),SyncSymb(IniSyncPos:SyncPos),ESync1(IniSyncPos:...
                                      SyncPos),ESync2(IniSyncPos:SyncPos));
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
            if DiffPos>=0%If the difference is positive, left-shift...
                EoutAux = ifft(fft(EoutAux).*exp(1j*2*pi*f*(DiffPos*Ta))); %Shift based on time change     
%                 EoutAux = [EoutAux(DiffPos+1:end) EoutAux(1:DiffPos)];   %Shift based on sampling sliding
            else%... but if the difference is negative, right-shift
                EoutAux = ifft(fft(EoutAux).*exp(1j*2*pi*f*(DiffPos*Ta))); %Shift based on time change     
%                 EoutAux = [EoutAux(end+DiffPos+1:end) EoutAux(1:end + ...
%                                                                 DiffPos)];%Shift based on sampling sliding
            end
            
                %                   Synchronizing
            %Because of reasons, sometimes it may be required to make a
            %synchronization process with the  end of the data stream as
            %well. This following verification check if the user set (or
            %not) a second synchronization process to be done.
            if SencondAdjust
                PhaDel          = 0;
                TimDel          = T;
                [ESync1,ESync2] = DelayInterf(t,TimDel,PhaDel,EoutAux);

                ESync1 = ESync1.*conj(ESync1);
                ESync2 = ESync2.*conj(ESync2);
                Esync  = ESync2-ESync1;

                ESync1 = ESync1./max(abs(ESync1)); 
                ESync2 = ESync2./max(abs(ESync2)); 
                Esync  = Esync./max(abs(Esync)); 
                

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
                if DiffPosEnd>=0%If positive difference, left-shift...
                    EoutAux = ifft(fft(EoutAux).*exp(1j*2*pi*f*(...
                                                          DiffPosEnd*Ta)));%Shift based on time change     
%                     EoutAux = [EoutAux(DiffPosEnd+1:end) EoutAux(1:...
%                                                             DiffPosEnd)];%Shift based on sampling sliding
                else%... but if the difference is negative, right-shift
                    EoutAux = ifft(fft(EoutAux).*exp(1j*2*pi*f*(...
                                                          DiffPosEnd*Ta)));%Shift based on time change     
%                     EoutAux = [EoutAux(end+DiffPosEnd+1:end) EoutAux(...
%                                                       1:end+DiffPosEnd)];%Shift based on sampling sliding
                end
            end
            %%                   Plot for Qualitative analizes
            PrintInfo(Ploting*23,t(IniSyncPos:SyncPos),Esync(IniSyncPos:...
                SyncPos),SyncSymb(IniSyncPos:SyncPos),ESync1(IniSyncPos:...
                                      SyncPos),ESync2(IniSyncPos:SyncPos));
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
            EoutA = ifft(fft(EoutA).*BitFilt);
            EoutB = ifft(fft(EoutB).*BitFilt);
            EoutC = ifft(fft(EoutC).*BitFilt);
            EoutD = ifft(fft(EoutD).*BitFilt);
            
            % The configuration here used is an balanced receiver as the 
            %output of the Delay Interferometer has two signals resulting 
            %from the constructive and destructive signal interaction. 
            EoutI = (EoutB - EoutA);
            EoutQ = (EoutD - EoutC);
            
            EoutI = EoutI./max(abs(EoutI));                                %Normalizing the signal
            EoutQ = EoutQ./max(abs(EoutQ));                                %Normalizing the signal
            
            %One important step on this project was the confirmation of our
            %results from a closed and theoretical equation that relates
            %the income optical field witht its respective data (as
            %descrived above in the receiving process). The result from
            %this equation was further compared with the result from the
            %MZ-Interferometer as an proof of concept. This equation can be
            %found at the book of Optical Fiber Telecommunications V B,
            %which one of the authors is Ivan P. Kaminow at the page 144.
            Ui = real(exp(1j*pi/4).*E_rec3.*conj(ifft(fft(E_rec3).*exp(-...
                                                           1j*2*pi*f*T))));%The data at odd position
            Ui = Ui./max(abs(Ui));                                         %Normalizing the signal
            Uq = imag(exp(1j*pi/4).*E_rec3.*conj(ifft(fft(E_rec3).*exp(-...
                                                           1j*2*pi*f*T))));%The data at the even position
            Uq = Uq./max(abs(Uq));                                         %Normalizing the signal
            %% Ploting the result for qualitative analizes
            PrintInfo(Ploting*24,EoutI,T,NPPB);
            %%
            [~,~,EyeOpenI,EyeOpenHighI,EyeOpenLowI] = Olho(EoutI,T,NPPB,0);
            %%  Ploting some results for qualitative analizes
            PrintInfo(Ploting*25,Ui,T,NPPB);
            %%
            Olho(Ui,T,NPPB,0);
            %%  Ploting some results for qualitative analizes
            PrintInfo(Ploting*26,EoutQ,T,NPPB);
            %%
            [~,~,EyeOpenQ,~,EyeOpenLowQ] = Olho(EoutQ,T,NPPB,Ploting);
            %%  Ploting some results for qualitative analizes
            PrintInfo(Ploting*27,Uq,T,NPPB);
            %%
            Olho(Uq,T,NPPB,Ploting);
                        
            Txaux1 = rectpulse(TxDataMat(ThisCarr,:),NPPB);
            Txaux2 = rectpulse(TxDataMat(ThisCarr+NumCarr,:),NPPB);
            Txaux1(Txaux1==0) = -1;
            Txaux2(Txaux2==0) = -1;
            %%  Ploting some results for qualitative analizes
            PrintInfo(Ploting*28,t,Txaux1,Txaux2,EoutI,EoutA,EoutB,EoutQ...
                                           ,EoutC,EoutD,real(Ui),real(Uq));
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
            EyeSymMatI = reshape(EoutI(1+SyncPeriod*NPPB:end-SyncPeriod*...
                                               NPPB),NPPB,Nb-2*SyncPeriod);  
            %Then we take the values that compose the decision level 
            %because they will mark the point of symmetry.
            %
            %Firstly it was set the interval in which the histogram will be
            %build. It is based on the number of samples per bit period.
            IntervalI = linspace(min(EoutI),max(EoutI),2*NPPB);
            %Therefore, the MATLAB hist function returns the number of
            %occurrence of each interval.
            EyeMaxI = hist(EoutI,IntervalI);
            EyeMaxauxI = [0 EyeMaxI 0];                                    %Zeros are added at the EyeMax to auxiliate the finding peaks process
            [~,EyeLocI] = findpeaks(EyeMaxauxI,'MinPeakDistance',NPPB/4,...
                                           'SortStr','descend','NPeaks',2);%The peaks on the Eye profile will be the levels at the Eyes limit
            %From the location of the max values that occured, which means
            %the uper and lower level of the eye diagram it needs to take
            %the actual value that those occurences represent that is
            %withing the the Interval variable.
            ValToSeekI = IntervalI(EyeLocI-1);                                
            ValToSeekI = sort(ValToSeekI,'ascend');
            %The number of ocurrences is a statical measure therefore one
            %does not have control which interval will have the highest
            %peak, thus it is important to ordenate the values to be seek
            %from the lower part of the eye diagram to the uper part of the
            %eye diagram.
            OccuCountI = zeros(1,size(EyeSymMatI,1));                      %Auxiliar Variable for accounting.
            for kk=1:size(EyeSymMatI,1)                                    %For every sample within a symbol period
                
                OccuCountI(kk) = OccuCountI(kk)+sum((EyeSymMatI(kk,:)>=...
                    min(EoutI))&(EyeSymMatI(kk,:)<=UpeSymPer*EyeOpenLowI)); %Account all occurencies of the valeu 1
                OccuCountI(kk) = OccuCountI(kk)+sum((EyeSymMatI(kk,:)>=...
                   LowSymPer*EyeOpenHighI)&(EyeSymMatI(kk,:)<=max(EoutI))); %Account all occurencies of the valeu 2
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
            [SymLocI] = round(mean(find(ismember(OccuCountI,max(...
                                                           OccuCountI)))));             
%             [~,SymLocI] = findpeaks(OccuCountI,'SortStr','descend');       %The peak on the Eye profile will be the Symmetry level
            DataOdd = [];                                                  %Initialization of the vector that will store the income data
            for kk=1:NPPB:length(EoutI)                                    %The comparison process will be made for each symbol period
                MeanOfData = mean(EoutI((kk-1)+SymLocI(1)));
                if MeanOfData > EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
                    DataOdd = [DataOdd 1];                                 %is 1
                else                                                       %If it is the lowest level the incoming data
                    DataOdd = [DataOdd 0];                                 %is 0
                end
            end
            %The identical process just described above will also be used
            %to recover the data at the even positions.
            %
            %As this process is also statistical, first we reshape the 
            %income vector to analyze all periods at the same time.
            EyeSymMatQ = reshape(EoutQ(1+SyncPeriod*NPPB:end-SyncPeriod*...
                                               NPPB),NPPB,Nb-2*SyncPeriod);  
            %Then we take the values that compose the decision level 
            %because they will mark the point of symmetry.
            %
            %Firstly it was set the interval in which the histogram will be
            %build. It is based on the number of samples per bit period.
            IntervalQ = linspace(min(EoutQ),max(EoutQ),2*NPPB);
            %Therefore, the MATLAB hist function returns the number of
            %occurrence of each interval.
            EyeMaxQ = hist(EoutQ,IntervalQ);
            EyeMaxauxQ = [0 EyeMaxI 0];                                    %Zeros are added at the EyeMax to auxiliate the finding peaks process
            [~,EyeLocQ] = findpeaks(EyeMaxauxQ,'MinPeakDistance',NPPB/4,...
                                           'SortStr','descend','NPeaks',2);%The peaks on the Eye profile will be the levels at the Eyes limit
            ValToSeekQ = IntervalQ(EyeLocQ-1);       
            ValToSeekQ = sort(ValToSeekQ,'ascend');                         
            OccuCountQ = zeros(1,size(EyeSymMatQ,1));                      %Auxiliar Variable for accounting.
            for kk=1:size(EyeSymMatQ,1)                                    %For every sample within a symbol period
                OccuCountQ(kk) = OccuCountQ(kk)+sum((EyeSymMatI(kk,:)>=...
                LowSymPer*ValToSeekQ(1))&(EyeSymMatQ(kk,:)<=UpeSymPer*...
                                                           ValToSeekQ(1)));%Account all occurencies of the valeu 1
                OccuCountQ(kk) = OccuCountQ(kk)+sum((EyeSymMatQ(kk,:)>=...
                LowSymPer*ValToSeekQ(2))&(EyeSymMatQ(kk,:)<=UpeSymPer*...
                                                           ValToSeekQ(2)));%Account all occurencies of the valeu 2
            end
%             [~,SymLocQ] = findpeaks(OccuCountQ,'SortStr','descend');       %The peak on the Eye profile will be the Symmetry level
            [SymLocQ] = round(mean(find(ismember(OccuCountQ,max(...
                                                           OccuCountQ)))));
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
            for kk=1:NPPB:length(EoutQ)                                    %The comparison process will be made for each symbol period
                MeanOfData = mean(EoutQ((kk-1)+SymLocI(1)));
                if MeanOfData > EyeOpenLowQ+EyeOpenQ/2                     %If it is the uper level the incoming data
                    DataEven = [DataEven 1];                               %is 1
                else                                                       %If it is the lowest level the incoming data
                    DataEven = [DataEven 0];                               %is 0
                end
            end
            %%       Calculating the Bit Error Ratio (BER)
            %The final process here is to count the number of wrongdoings
            %of this whole process upon the transmited data for 
            %quantitative analizes
            TxDataOdd  = TxDataMat(ThisCarr,1+SyncPeriod:end-SyncPeriod);
            TxDataEven = TxDataMat(ThisCarr+NumCarr,1+SyncPeriod:end-...
                                                             SyncPeriod);
            DataOdd  = DataOdd(1+SyncPeriod:end-SyncPeriod);
            DataEven = DataEven(1+SyncPeriod:end-SyncPeriod);
            
            BitErrOdd  = sum(xor(TxDataOdd,DataOdd));                      %Comparison between the Transmited and received and counting the differences
            BitErrEven = sum(xor(TxDataEven,DataEven));                    %Comparison between the Transmited and received and counting the differences
            if exist('CurrentOCS','var')
                BerDQPSK(CurrentTest+1,ThisCarr,CurrentOCS) = (BitErrOdd+...
                                       BitErrEven)/((2*Nb)-(2*SyncPeriod));%Calculating the ration of wrong bits and the total number of bits transmited
            else
                BerDQPSK(CurrentTest,ThisCarr) = (BitErrOdd+BitErrEven)/...
                                                   ((2*Nb)-(2*SyncPeriod));%Calculating the ration of wrong bits and the total number of bits transmited
            end
            %% Ploting the result for qualitative analizes
            PrintInfo(Ploting*29,TxDataOdd,DataOdd);
            PrintInfo(Ploting*30,TxDataEven,DataEven);
            %%
            a=6;
            close all;
            
        end
    case '4PAM'
        
        if exist('CurrentOCS','var')
            Ber4PAM(1,:,CurrentOCS) = OSNRPC(1:NumCarr);
        end
        %%              Receiver 4PAM
        for ThisCarr=1:2:NumCarr                                           %For each carrier the same process of reception will be used.
            %%           Creating the Reception Filter
            [BitFilt,~] = FiltroGaussiano(f,BWD2,CenFeq2,FiltOrd2);        %Creating filter for selection of the received signal
            BitFilt = fftshift(BitFilt);                                   %Shifting the filter for matching the received signal 
            %%  Reception Process: Handling the income optical field
            %At the first step the income field will be selected from the
            %optical FFT Output with the help of the maping vector
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
            %The current incoming signal is them converted from the optical
            %domain to the eletrical domain with the help of an photo
            %detector.
%             EoutAux = EoutAux./max(abs(EoutAux));
%             EoutAux = EoutAux + 1;
            Ix = EoutAux.*conj(EoutAux);
            %The process with the photo diode is self-coherent, which means
            %the rusult will be a component of the signal centered in f=0
            %and another component centered at f=2*fc (frenquecy central).
            %Therefore, to remove the higher order component a low pass
            %filter will be used.
%             switch Medium
%                 case 'Fiber'
%                     Ix = ifft(fft(Ix).*BitFilt);
%                 otherwise
%             end
            
            Ix = Ix - min(Ix);                                             %Removing the DC component from them Eletrical signal received
            
            Ix = 1.*Ix./max(abs(Ix));                                      %Normalizing the eletrical signal (amplifying what is needed)
            %% Ploting the result for qualitative analizes
            PrintInfo(Ploting*31,f,20.*log10(abs(fftshift(fft(Ix)./...
                                                            length(Ix)))));
            PrintInfo(Ploting*32,t,Ix);
            PrintInfo(Ploting*33,Ix,T,NPPB);
            PrintInfo(Ploting*34,t(IniSyncPos:SyncPos),Ix(IniSyncPos:...
                                    SyncPos),SyncSymb(IniSyncPos:SyncPos));
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
                    Ix = [Ix(end+AuxSyncCorrEnd+1:end) Ix(1:end + ...
                                                          AuxSyncCorrEnd)];%Shift based on sampling sliding
                end
            end
            
            %%          Ploting the result for qualitative analizes
            PrintInfo(Ploting*35,t(end-SyncPos+1:end-IniSyncPos+1),Ix(...
                end-SyncPos+1:end-IniSyncPos+1),SyncSymb(end-SyncPos+1:...
                                                        end-IniSyncPos+1));
            PrintInfo(PlotingThis*36,Ix,T,NPPB);
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
            Interval = linspace(min(Ix(1+SyncPeriod:end-SyncPeriod)),max...
                           (Ix(1+SyncPeriod:end-SyncPeriod)),IntervalStep);
            %Therefore, the MATLAB hist function returns the number of
            %occurrence of each interval.
            EyeMax = hist(Ix,Interval);
            EyeMaxaux = [0 EyeMax 0];                                      %Zeros are added at the EyeMax to auxiliate the finding peaks process
            [~,EyeLoc] = findpeaks(EyeMaxaux,'MinPeakDistance',MinDist,...
                                           'SortStr','descend','NPeaks',4);%The peaks on the Eye profile will be the levels at the Eyes limit
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
            PrintInfo(Ploting*37,Interval,EyeMax);
            PrintInfo(Ploting*38,EyeMax,EyeLoc-1);
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
                Levels = [0 1/3 2/3 1];
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
            EyeSymMat = reshape(Ix,NPPB,Nb);  
            %Then we take the values that compose the decision level 
            %because they will mark the point of symmetry.
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
                OccuCount(kk) = OccuCount(kk) + sum((EyeSymMat(kk,:)>=...
                LowSymPer*ValToSeek(1))&(EyeSymMat(kk,:)<=UpeSymPer*...
                                                            ValToSeek(1)));%Account all occurencies of the valeu 1
                OccuCount(kk) = OccuCount(kk) + sum((EyeSymMat(kk,:)>=...
                LowSymPer*ValToSeek(2))&(EyeSymMat(kk,:)<=UpeSymPer*...
                                                            ValToSeek(2)));%Account all occurencies of the valeu 2
                OccuCount(kk) = OccuCount(kk) + sum((EyeSymMat(kk,:)>=...
                LowSymPer*ValToSeek(3))&(EyeSymMat(kk,:)<=UpeSymPer*...
                                                            ValToSeek(3)));%Account all occurencies of the valeu 3
                OccuCount(kk) = OccuCount(kk) + sum((EyeSymMat(kk,:)>=...
                LowSymPer*ValToSeek(4))&(EyeSymMat(kk,:)<=UpeSymPer*...
                                                            ValToSeek(4)));%Account all occurencies of the valeu 4
            end

            %The point of symmetry of the eye diagram will be where the 
            %maximum number of occurrences were measured inasmuch as those 
            %are the points where all the bits go to the center of the 
            %symbol.  From the maximum number of occurrences, it can happen 
            %for more than one sample within one symbol period, in the best 
            %case, all samples would have the same accounting as it is 
            %shown the ilustration above hence the symmetry will be at the 
            %middle sample of this group of maximum occurrences. This value 
            %can be found at the highest peak founded.
            [~,SymLoc] = findpeaks(OccuCount,'SortStr','descend');         %The peak on the Eye profile will be the Symmetry level
            
            PosIx = SymLoc(1):NPPB:length(Ix);
            IxAux = Ix(PosIx);
            n=0;
            
            InterAB = linspace(Levels(3),Levels(4),NPPB*2^n);
            EyeAB = hist(IxAux,InterAB);
            EyeAB = ~EyeAB;
            LocAB = find(EyeAB);
            if isempty(LocAB)
                LevelDec3 = mean(Levels(3:4));
            else
                LevelDec3 = InterAB(LocAB(round(end/2)));
            end
            
            InterCD = linspace(Levels(2),Levels(3),NPPB*2^n);
            EyeCD = hist(IxAux,InterCD);
            EyeCD = ~EyeCD;
            LocCD = find(EyeCD);
            if isempty(LocCD)
                LevelDec2 = mean(Levels(2:3));
            else
                LevelDec2 = InterCD(LocCD(round(end/2)));
            end
            
            InterEF = linspace(Levels(1),Levels(2),NPPB*2^n);
            EyeEF = hist(IxAux,InterEF);
            EyeEF = ~EyeEF;
            LocEF = find(EyeEF);
            if isempty(LocEF)
                LevelDec1 = mean(Levels(1:2));
            else
                LevelDec1 = InterEF(LocEF(round(end/2)));
            end
            %%           Ploting for Qualitative Analizes
            PrintInfo(Ploting*39,NPPB,OccuCount);
            PrintInfo(Ploting*40,OccuCount,SymLoc);
            %%      Actualy Receiving Data
            %Once the signal was processed the next step is through a
            %comparator decide the actual information received.
            IxRec = [];                                                    %Initialization of the vector that will store the income data
            for kk=1:NPPB:length(Ix)                                       %The comparison process will be made for each symbol period
                aux1 = Ix((kk-1)+SymLoc(1));                               %An small portion of the income signal is take for evaluation
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
            end
            %%           Ploting for Qualitative Analizes
            PrintInfo(Ploting*41,FinalTime,Nb,TxDataMat(ThisCarr,:),IxRec);
            %%       Calculating the Bit Error Ratio (BER)
            %The final process here is to count the number of wrongdoings
            %of this whole process upon the transmited data for 
            %quantitative analizes
            BitErr = sum(xor(TxDataMat(ThisCarr,1+2*SyncPeriod:end-2*...
                      SyncPeriod),IxRec(1+2*SyncPeriod:end-2*SyncPeriod)));%Comparison between the Transmited and received and counting the differences
            if exist('CurrentOCS','var')
                Ber4PAM(CurrentTest+1,ThisCarr,CurrentOCS) = BitErr/(2*Nb-...
                                                           (2*SyncPeriod));
            else
                Ber4PAM(CurrentTest,ThisCarr) = BitErr/(2*Nb-(2*...
                                                              SyncPeriod));%Calculating the ration of wrong bits and the total number of bits transmited
            end
            a=3;
%             drawnow;
%             pause(0.5);
            close all;
        end
    otherwise
        if exist('CurrentOCS','var')
            BerOOK(1,:,CurrentOCS) = OSNRPC(1:NumCarr);
        end
        %%              Receiver OOK
        for ThisCarr=1:2:NumCarr                                           %For each carrier the same process of reception will be used.
            %%           Creating the Reception Filter
            [BitFilt,~] = FiltroGaussiano(f,BWD,CenFeq,FiltOrd);           %Creating filter for selection of the received signal
            BitFilt = fftshift(BitFilt);                                   %Shifting the filter for matching the received signal 
            %%  Reception Process: Handling the income optical field
            %At the first step the income field will be selected from the
            %optical FFT Output with the help of the maping vector
            %(previously described).
            EoutDemd = EoutAux1(VetThisCarr==ThisCarr,:);
            %%            Fiber Time Delay Compensation
            switch Medium
                case 'Fiber'
                    EoutDemd = ifft(fft(EoutDemd).*exp(1j*2*pi*f*(...
                                                FiberDelay(ThisCarr)*Ta)));
                otherwise
            end
                  
            %The current incoming signal is them converted from the optical
            %domain to the eletrical domain with the help of an photo
            %detector.
            Ix =EoutDemd.*conj(EoutDemd);
            Ix = ifft(fft(Ix).*BitFilt);                                   %Filtering the signal for better detection
            Ix = Ix - min(Ix);                                             %Removing any off-set that may exist 
            Ix = Ix./max(abs(Ix));                                         %Normalizing the eletrical signal (amplifying what is needed)
            %% Ploting the result for qualitative analizes
            PrintInfo(Ploting*42,f,20.*log10(abs(fftshift(fft(Ix)./...
                                                            length(Ix)))));
            PrintInfo(Ploting*43,t,TxSigMat(ThisCarr,:)./max(TxSigMat(...
                                                          ThisCarr,:)),Ix);
            PrintInfo(Ploting*44,Ix,T,NPPB);
            PrintInfo(Ploting*45,t(IniSyncPos:SyncPos),Ix(IniSyncPos:...
                                    SyncPos),SyncSymb(IniSyncPos:SyncPos));
            %%        synchronizing
            %For the reception process work properly, it is needed to
            %sincronized the recieved signal or the sampling process will
            %not work.
            
            AuxSync = (Ix(IniSyncPos:SyncPos));                                %Selecting the sync-word within the received signal
            AuxSync1 = AuxSync;
            %This synchronization process is based on the mean expected
            %value. That means, the information mean value should be within
            %one period of symbol. Thus, the mean value of the received
            %signal is acquired and compare of the known sync-word to
            %verify if this mean value is at the right possition
            AuxSync1(AuxSync1<0) = 0;                                      %To keep the mean value above zero anything under is neglected
            AuxSync1(AuxSync1>=mean(AuxSync1)) = 1;                        %Adding a flag to the first sample of the received mean value
            AuxSync1(AuxSync1<mean(AuxSync1)) = -1;                        %All the others samples at set to the lowest level
            AuxSync2 = SyncSymb(IniSyncPos:SyncPos);                           %Selecting the sync-word within the known signal
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
                SyncSymbEndAux = SyncSymbEnd(end-SyncPos+1:end-2*NPPB);      %Selecting the sync-word within the known signal
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
                    Ix = [Ix(end+AuxSyncCorrEnd+1:end) Ix(1:end + ...
                                                          AuxSyncCorrEnd)];%Shift based on sampling sliding
                end
            end
            %% Ploting the result for qualitative analizes
            PrintInfo(Ploting*46,t,TxSigMat(ThisCarr,:)./max(TxSigMat(...
                                                          ThisCarr,:)),Ix);
          
            PrintInfo(Ploting*47,abs(Ix),T,NPPB);
            [~,~,EyeOpen,EyeOpenHigh,EyeOpenLow] = Olho(abs(Ix),T,NPPB,0);
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

            EyeSymMat = reshape(Ix(1+SyncPeriod*NPPB:end-SyncPeriod*...
                                               NPPB),NPPB,Nb-2*SyncPeriod);  
            %Then we take the values that compose the decision level 
            %because they will mark the point of symmetry.
            %
            %Firstly it was set the interval in which the histogram will be
            %build. It is based on the number of samples per bit period.
            Interval = linspace(min(Ix),max(Ix),2*NPPB);
            %Therefore, the MATLAB hist function returns the number of
            %occurrence of each interval.
            EyeMax = hist(Ix,Interval);
            EyeMaxaux = [0 EyeMax 0];                                      %Zeros are added at the EyeMax to auxiliate the finding peaks process
            [~,EyeLoc] = findpeaks(EyeMaxaux,'MinPeakDistance',NPPB/4,...
                                           'SortStr','descend','NPeaks',2);%The peaks on the Eye profile will be the levels at the Eyes limit
            
%             figure;findpeaks(EyeMaxaux,'MinPeakDistance',NPPB/4,'SortStr','descend','NPeaks',2);
            
            ValToSeek = Interval(EyeLoc-1);             
            %From the location of the max values that occured, which means
            %the uper and lower level of the eye diagram it needs to take
            %the actual value that those occurences represent that is
            %withing the the Interval variable.                   
            ValToSeek = sort(ValToSeek,'ascend');
            %The number of ocurrences is a statical measure therefore one
            %does not have control which interval will have the highest
            %peak, thus it is important to ordenate the values to be seek
            %from the lower part of the eye diagram to the uper part of the
            %eye diagram.
            OccuCount = zeros(1,size(EyeSymMat,1));                        %Auxiliar Variable for accounting.
            for kk=1:size(EyeSymMat,1)                                     %For every sample within a symbol period
                OccuCount(kk) = OccuCount(kk)+sum((EyeSymMat(kk,:)>=min(...
                             Ix))&(EyeSymMat(kk,:)<=UpeSymPer*EyeOpenLow));%Account all occurencies of the value 1
                OccuCount(kk) = OccuCount(kk)+sum((EyeSymMat(kk,:)>=...
                        LowSymPer*EyeOpenHigh)&(EyeSymMat(kk,:)<=max(Ix)));%Account all occurencies of the value 2
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
                                                                      ))));%The peak on the Eye profile will be the Symmetry level
            %%           Ploting for Qualitative Analizes
            PrintInfo(Ploting*48,NPPB,OccuCount);
            PrintInfo(Ploting*49,OccuCount,SymLoc);

            %%      Actualy Receiving Data
            %Once the signal was processed the next step is through a
            %comparator decide the actual information received.
            EoutCorr = [];                                                 %Initialization of the vector that will store the income data
            for kk=1:NPPB:length(Ix)                                       %The comparison process will be made for each symbol period
                %An small portion of the income signal is take for 
                %evaluation by measuring the avarage value of the samples 
                %taken
                CalcMean = mean(abs(Ix((kk-1)+SymLoc(1))));
                %Verifying the interval for each symbol received. 
                if CalcMean >= EyeOpenLow+EyeOpen/2                        %If it is the uper level the incoming data
                    EoutCorr = [EoutCorr 1];                               %is 1
                else                                                       %If it is the lowest level the incoming data
                    EoutCorr = [EoutCorr 0];                               %is 0
                end
            end

            PrintInfo(Ploting*50,t,rectpulse(TxDataMat(ThisCarr,:),NPPB)...
                                                ,rectpulse(EoutCorr,NPPB));
            %%       Calculating the Bit Error Ratio (BER)
            %The final process here is to count the number of wrongdoings
            %of this whole process upon the transmited data for 
            %quantitative analizes
            BitErr = sum(xor(TxDataMat(ThisCarr,1+SyncPeriod*NPPB:end-...
            SyncPeriod*NPPB),EoutCorr(1+SyncPeriod*NPPB:end-SyncPeriod*...
                                                                   NPPB)));%Comparison between the Transmited and received and counting the differences
            if exist('CurrentOCS','var')
                BerOOK(CurrentTest+1,ThisCarr,CurrentOCS) = BitErr/(Nb-2*...
                                                               SyncPeriod);
            else
                BerOOK(CurrentTest,ThisCarr) = BitErr/(Nb-2*SyncPeriod);       %Calculating the ration of wrong bits and the total number of bits transmited
            end
            a=1;
            close all;
        end
%         end
end

a=0;


%% ########################################################################
%##########################################################################
%##########################################################################
%                 This was the end for the reception of data.
%   The folowing lines are exemplifications for the optical FFT processs.
%##########################################################################
%##########################################################################
%##########################################################################


%%
if Elucidation
    %% Simple example of how the Optical FFT was implemented
    % Them main code that implements an all optical FFT was previouly described
    % at those following lines is described step by step this process.
    %% First Stage
    % At this point the scrip is implementing an FFT for 2 subcaries. The basic
    % idea is to use an interferometre Delay to create an controled
    % desconstrutive interferency between the input signal;
     %split
     % Initially the signal was splited in two different path. The signal that
     % passes within the uper arm will suffe no defformation.
    Eout11 = ExpCoef.*(cos(kl).*EoutRec + MutCoef*1j*sin(kl).*0);
    Eout12 = ExpCoef.*(cos(kl).*0 + MutCoef*1j*sin(kl).*EoutRec);
     %Dellay and change phase
     %Mean while the signal that passes through the Lower arm will be dellayed
     %and its phase will be changed.
    aux1 = Eout12(t<=Del0);
    Eout12d = [Eout12(end-(length(aux1)-1):end) Eout12(1:end-(length(aux1)))];
    Eout12d = Eout12d.*exp(-1j*Pha0);
    %Ploting results for qualitative analizes
    if Ploting
        figure;
        plot(t,abs(Eout12),t,abs(Eout12d))
    end
    % Thus, the signal will be coupled again. This process was implemented with
    % an teoretical device.
    Eout1a = ExpCoef.*(cos(kl).*Eout11 + MutCoef*1j*sin(kl).*Eout12d);
    Eout1b = ExpCoef.*(cos(kl).*Eout12d + MutCoef*1j*sin(kl).*Eout11);

    %Ploting results for qualitative analizes
    if Ploting
        figure;
        hold all;
        plot(f,db(abs((fftshift(fft(EoutT))))/length(EoutT)));
        plot(f,db(abs((fftshift(fft(Eout1a))))/length(Eout1a)));
        plot(f,db(abs((fftshift(fft(Eout1b))))/length(Eout1b)));
        axis([1e11 3e11 -370 10]);
    end
    %% Second Stage
    % From the first stage two signals where generated. This stage will be
    % exactly as the first. With the difference that now, there will be two
    % delay interferometers. It is also important to mention that there will be
    % variation on the input parameters of the DI such as different time and
    % phase delay from the first stage
      %% Up arm
      % This is the first Delay Interferomentro within this stage.
     %split
     % Initially the signal was splited in two different path. The signal that
     % passes within the uper arm will suffe no defformation.
    Eout21a = ExpCoef.*(cos(kl).*Eout1a + MutCoef*1j*sin(kl).*0);
    Eout22a = ExpCoef.*(cos(kl).*0 + MutCoef*1j*sin(kl).*Eout1a);
     %Dellay and change phase
     %Mean while the signal that passes through the Lower arm will be dellayed
     %and its phase will be changed.

    %Ploting results for qualitative analizes
    if Ploting
        figure;
        hold all;
        plot(f,db(abs((fftshift(fft(EoutT))))/length(EoutT)));
        plot(f,db(abs((fftshift(fft(Eout21a))))/length(Eout21a)));
        plot(f,db(abs((fftshift(fft(Eout22a))))/length(Eout22a)));
        axis([1e11 3e11 -370 10]);
    end

    % Eout21a = (1/2).*Eout1a;
    % Eout22a = (1/2).*Eout1a;
     %Dellay and change phase
    aux2a = Eout22a(t<=Del1a);
    Eout22da = [Eout22a(end-(length(aux2a)-1):end) Eout22a(1:end-(length(aux2a)))];
    Eout22da = Eout22da.*exp(-1j*Pha1a);
    %Ploting results for qualitative analizes
    if Ploting
        figure;
        plot(t,abs(Eout22a),t,abs(Eout22da))
    end
    %Ploting results for qualitative analizes
    if Ploting
        figure;
        hold all;
        plot(f,db(abs((fftshift(fft(EoutT))))/length(EoutT)));
        plot(f,db(abs((fftshift(fft(Eout21a))))/length(Eout21a)));
        plot(f,db(abs((fftshift(fft(Eout22da))))/length(Eout22da)));
        axis([1e11 3e11 -370 10]);
    end
    % Eout2a = Eout21a + Eout22da;
    % Eout2aa = (1/2).*Eout2a;
    % Eout2ab = (1/2).*Eout2a;
    % kl =(90)*pi/180;
    Eout2aa = ExpCoef.*(cos(kl).*Eout21a + MutCoef*1j*sin(kl).*Eout22da);
    Eout2ab = ExpCoef.*(cos(kl).*Eout22da + MutCoef*1j*sin(kl).*Eout21a);

    %Ploting results for qualitative analizes
    if Ploting
        figure;
        hold all;
        plot(f,db(abs((fftshift(fft(EoutT))))/length(EoutT)));
        plot(f,db(abs((fftshift(fft(Eout2aa))))/length(Eout2aa)));
        plot(f,db(abs((fftshift(fft(Eout2ab))))/length(Eout2ab)));
        axis([1e11 3e11 -370 10]);
    end
      %% Down arm
      %The second signal generate by the first stage will enter here
     %split
     %As previouly mentioned. The frist part is to split the input signal in
     %two parts that will travel through different paths.
    Eout21b = ExpCoef.*(cos(kl).*Eout1b + MutCoef*1j*sin(kl).*0);
    Eout22b = ExpCoef.*(cos(kl).*0 + MutCoef*1j*sin(kl).*Eout1b);
    % Eout21b = (1/2).*Eout1b;
    % Eout22b = (1/2).*Eout1b;
     %Dellay and change phase
    aux2b = Eout22b(t<=Del1b);
    Eout22db = [Eout22b(end-(length(aux2b)-1):end) Eout22b(1:end-(length(aux2b)))];
    Eout22db = Eout22db.*exp(-1j*Pha1b);
    %Ploting results for qualitative analizes
    if Ploting
        figure;
        plot(t,abs(Eout22b),t,abs(Eout22db))
    end
    % Eout2b = Eout21b + Eout22db;
    % Eout2ba = (1/2).*Eout2b;
    % Eout2bb = (1/2).*Eout2b;
    Eout2ba = ExpCoef.*(cos(kl).*Eout21b + MutCoef*1j*sin(kl).*Eout22db);
    Eout2bb = ExpCoef.*(cos(kl).*Eout22db + MutCoef*1j*sin(kl).*Eout21b);

    if Ploting
        figure;
        hold all;
        plot(f,db(abs((fftshift(fft(EoutT))))/length(EoutT)));
        plot(f,db(abs((fftshift(fft(Eout2ba))))/length(Eout2ba)));
        plot(f,db(abs((fftshift(fft(Eout2bb))))/length(Eout2bb)));
        axis([1e11 3e11 -370 10]);
    end
    %% Third Stage
    % Finaly after the second stage, now we have 4 signals that was generated
    % each one would filter one sub carrier.
      %% First Arm
     %split
    Eout31aa = ExpCoef.*(cos(kl).*Eout2aa + MutCoef*1j*sin(kl).*0);
    Eout32aa = ExpCoef.*(cos(kl).*0 + MutCoef*1j*sin(kl).*Eout2aa);
    if Ploting
        figure;
        hold all;
        plot(f,db(abs((fftshift(fft(EoutT))))/length(EoutT)));
        plot(f,db(abs((fftshift(fft(Eout31aa))))/length(Eout31aa)));
        plot(f,db(abs((fftshift(fft(Eout32aa))))/length(Eout32aa)));
        axis([1e11 3e11 -370 10]);
    end
    % Eout31aa = (1/2).*Eout2aa;
    % Eout32aa = (1/2).*Eout2aa;
     %Dellay and change phase
    aux3aa = Eout32aa(t<=Del2aa);
    Eout32daa = [Eout32aa(end-(length(aux3aa)-1):end) Eout32aa(1:end-(length(aux3aa)))];
    Eout32daa = Eout32daa.*exp(-1j*Pha2aa);
    if Ploting
        figure;
        plot(t,abs(Eout32aa),t,abs(Eout32daa))
    end
    % Eout3aa = Eout31aa + Eout32daa;
    % Eout3aaa = (1/2).*Eout3aa;
    % Eout3aab = (1/2).*Eout3aa;
    Eout3aaa = ExpCoef.*(cos(kl).*Eout31aa + MutCoef*1j*sin(kl).*Eout32daa);
    Eout3aab = ExpCoef.*(cos(kl).*Eout32daa + MutCoef*1j*sin(kl).*Eout31aa);

    if Ploting
        figure;
        hold all;
        plot(f,db(abs((fftshift(fft(EoutT))))/length(EoutT)));
        plot(f,db(abs((fftshift(fft(Eout3aaa))))/length(Eout3aaa)));
        plot(f,db(abs((fftshift(fft(Eout3aab))))/length(Eout3aab)));
        axis([1e11 3e11 -370 10]);
    end
      %% Second Arm
     %split
    Eout31ab = ExpCoef.*(cos(kl).*Eout2ab + MutCoef*1j*sin(kl).*0);
    Eout32ab = ExpCoef.*(cos(kl).*0 + MutCoef*1j*sin(kl).*Eout2ab);
    % Eout31ab = (1/2).*Eout2ab;
    % Eout32ab = (1/2).*Eout2ab;
     %Dellay and change phase
    aux3ab = Eout32ab(t<=Del2ab);
    Eout32dab = [Eout32ab(end-(length(aux3ab)-1):end) Eout32ab(1:end-(length(aux3ab)))];
    Eout32dab = Eout32dab.*exp(-1j*Pha2ab);
    if Ploting
        figure;
        plot(t,abs(Eout32ab),t,abs(Eout32dab))
    end
    % Eout3ab = Eout31ab + Eout32dab;
    % Eout3aba = (1/2).*Eout3ab;
    % Eout3abb = (1/2).*Eout3ab;
    Eout3aba = ExpCoef.*(cos(kl).*Eout31ab + MutCoef*1j*sin(kl).*Eout32dab);
    Eout3abb = ExpCoef.*(cos(kl).*Eout32dab + MutCoef*1j*sin(kl).*Eout31ab);

    if Ploting
        figure;
        hold all;
        plot(f,db(abs((fftshift(fft(EoutT))))/length(EoutT)));
        plot(f,db(abs((fftshift(fft(Eout3aba))))/length(Eout3aba)));
        plot(f,db(abs((fftshift(fft(Eout3abb))))/length(Eout3abb)));
        axis([1e11 3e11 -370 10]);
    end
      %% Third Arm
     %split
    Eout31ba = ExpCoef.*(cos(kl).*Eout2ba + MutCoef*1j*sin(kl).*0);
    Eout32ba = ExpCoef.*(cos(kl).*0 + MutCoef*1j*sin(kl).*Eout2ba);
    % Eout31ba = (1/2).*Eout2ba;
    % Eout32ba = (1/2).*Eout2ba;
     %Dellay and change phase
    aux3ba = Eout32ba(t<=Del2ba);
    Eout32dba = [Eout32ba(end-(length(aux3ba)-1):end) Eout32ba(1:end-(length(aux3ba)))];
    Eout32dba = Eout32dba.*exp(-1j*Pha2ba);
    if Ploting
        figure;
        plot(t,abs(Eout32ba),t,abs(Eout32dba))
    end
    % Eout3ba = Eout31ba + Eout32dba;
    % Eout3baa = (1/2).*Eout3ba;
    % Eout3bab = (1/2).*Eout3ba;
    Eout3baa = ExpCoef.*(cos(kl).*Eout31ba + MutCoef*1j*sin(kl).*Eout32dba);
    Eout3bab = ExpCoef.*(cos(kl).*Eout32dba + MutCoef*1j*sin(kl).*Eout31ba);

    if Ploting
        figure;
        hold all;
        plot(f,db(abs((fftshift(fft(EoutT))))/length(EoutT)));
        plot(f,db(abs((fftshift(fft(Eout3baa))))/length(Eout3baa)));
        plot(f,db(abs((fftshift(fft(Eout3bab))))/length(Eout3bab)));
        axis([1e11 3e11 -370 10]);
    end
      %% Forth Arm
     %split
    Eout31bb = ExpCoef.*(cos(kl).*Eout2bb + MutCoef*1j*sin(kl).*0);
    Eout32bb = ExpCoef.*(cos(kl).*0 + MutCoef*1j*sin(kl).*Eout2bb);
    % Eout31bb = (1/2).*Eout2bb;
    % Eout32bb = (1/2).*Eout2bb;
     %Dellay and change phase
    aux3bb = Eout32bb(t<=Del2bb);
    Eout32dbb = [Eout32bb(end-(length(aux3bb)-1):end) Eout32bb(1:end-(length(aux3bb)))];
    Eout32dbb = Eout32dbb.*exp(-1j*Pha2bb);
    if Ploting
        figure;
        plot(t,abs(Eout32bb),t,abs(Eout32dbb))
    end
    % Eout3bb = Eout31bb + Eout32dbb;
    % Eout3bba = (1/2).*Eout3bb;
    % Eout3bbb = (1/2).*Eout3bb;
    Eout3bba = ExpCoef.*(cos(kl).*Eout31bb + MutCoef*1j*sin(kl).*Eout32dbb);
    Eout3bbb = ExpCoef.*(cos(kl).*Eout32dbb + MutCoef*1j*sin(kl).*Eout31bb);

    if Ploting
        figure;
        hold all;
        plot(f,db(abs((fftshift(fft(EoutT))))/length(EoutT)));
        plot(f,db(abs((fftshift(fft(Eout3bba))))/length(Eout3bba)));
        plot(f,db(abs((fftshift(fft(Eout3bbb))))/length(Eout3bbb)));
        axis([1e11 3e11 -370 10]);
    end
    % 
    % aux = Eout3bbb(t<=GrupDel);
    % Eout3bbb = [Eout3bbb(end-(length(aux3bb)-1):end) Eout3bbb(1:end-(length(aux3bb)))];
    clear Eout11 Eout12 Eout12d Eout1 Eout21a Eout22a Eout22da Eout2a ...
    Eout21b Eout22b Eout22db Eout2b Eout31aa Eout32aa Eout32daa Eout3aa ...
    Eout31ab Eout32ab Eout32dab Eout3ab Eout31ba Eout32ba Eout32dba Eout3ba...
    Eout31bb Eout32bb Eout32dbb Eout3bb;
    %%
    if Ploting
        figure;
        hold all;
        plot(f,db(abs((fftshift(fft(EoutT))))/length(EoutT)));
        plot(f,db(abs((fftshift(fft(Eout3aaa))))/length(Eout3aaa)));
        plot(f,db(abs((fftshift(fft(Eout3aab))))/length(Eout3aab)));
        plot(f,db(abs((fftshift(fft(Eout3aba))))/length(Eout3aba)));
        plot(f,db(abs((fftshift(fft(Eout3abb))))/length(Eout3abb)));
        plot(f,db(abs((fftshift(fft(Eout3baa))))/length(Eout3baa)));
        plot(f,db(abs((fftshift(fft(Eout3bab))))/length(Eout3bab)));
        plot(f,db(abs((fftshift(fft(Eout3bba))))/length(Eout3bba)));
        plot(f,db(abs((fftshift(fft(Eout3bbb))))/length(Eout3bbb)));
        axis([1e11 3e11 -370 10]);
    end
    Eout1 = ifft(fft(Eout3aaa).*SelecFilt3aaa);
    Eout1Demd = abs(Eout1).^2;
    % Eout1Demd = Eout1Demd-2;
    % Eout1Demd = Eout1Demd-min(abs(Eout1Demd));
    Eout1Demd = Eout1Demd./max(abs(Eout1Demd));
    Olho( Eout1Demd,T,NPPB );
    if Ploting
        figure;
        hold all;
        plot(fb,db(abs((fftshift(fft(EoutT))))/length(EoutT)));
        plot(fb,db(abs((fftshift(fft(Eout3aaa))))/length(Eout3aaa)));
        plot(fb,db(abs((fftshift(fft(Eout1))))/length(Eout1)));
        plot(fb,db(abs(fftshift(SelecFilt3aaa))));
        axis([-4e10 6e10 -500 40]);
        a=1;
    end
    if Ploting
        figure;
        hold all;
        plot(tb,TxDataNrz);
        plot(tb,abs(Eout1Demd));
        a=1;
    end

    Eout4 = ifft(fft(Eout3bba).*SelecFilt3bba);
    Eout4Demd = abs(Eout4).^2;
    % Eout4Demd = Eout4Demd-min(abs(Eout4Demd));
    Eout4Demd = Eout4Demd./max(abs(Eout4Demd));
    Olho( Eout4Demd,T,NPPB );
    if Ploting
        figure;
        hold all;
        plot(fb,db(abs((fftshift(fft(EoutT))))/length(EoutT)));
        plot(fb,db(abs((fftshift(fft(Eout3bba))))/length(Eout3bba)));
        plot(fb,db(abs((fftshift(fft(Eout4))))/length(Eout4)));
        plot(fb,db(abs(fftshift(SelecFilt3bba))));
        axis([1e10 9e10 -500 40]);
        a=1;
    end
    if Ploting
        figure;
        hold all;
        plot(tb,TxDataNrz);
        plot(tb,abs(Eout4Demd));
        a=1;
    end
    a=1;
end


clearvars -except UsedModula TestBundle ThisModula CurrentTest ...
CurrentModula ThisTest ThisModula OfcName BerOOK Ber4PAM BerDQPSK ...
                BerDPSK CurrentMedium OSNRPC OcsToTest CurrentOCS PulseResp

