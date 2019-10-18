function [BerDPSK,ValsLev,AberLev,Esync,InterAB,LocMaxAB,MaxValAB,SeqFinAB,TxData,Data] = RecDowDepsk(MaxNumStag,IfftOrSum,T,Ta,Esync1,AddCP,StuffSampels,NumAmosCP,NPPB,NbDPSK,CurTesSiz,SyncPeriod,ThisCarr,TxDataMat,LevDefDpqsk)
    %SyncAux   = Esync(IniSyncPos:SyncPos,:);%Selecting just the symbol to synchronize
    %SyncedAux = SyncSymb(IniSyncPos:SyncPos,:);%Selecting just the symbol to synchronize

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
    %SyncAux(SyncAux<0)              = 0;                           %To keep the mean value above zero anything under is neglected
    %SyncAux(SyncAux>=mean(SyncAux)) = 1;                           %Adding a flag to the first sample of the received mean value
    %SyncAux(SyncAux<mean(SyncAux))  = -1;                          %All the others samples at set to the lowest level

    %PosToSyn  = find(ismember(SyncAux,1));                         %Finding where is the location of the samples to synchronize
    %PosSynced = find(ismember(SyncedAux,1));                       %Finding where is the location of the samples to synchronize

    %DiffPos = PosToSyn(round(end/2)) - PosSynced(round(end/2));    %Accounting the peak (midel point) displacement
    sn=((1/2)*(((1/2)^MaxNumStag)-1))/((1/2)-1);
    if IfftOrSum==1
        DiffPos = round((sn/(3/2.55^IfftOrSum))*T/Ta);
    else
        DiffPos = round((sn/(3/2.55^IfftOrSum))*T/Ta);
    end
    Esync1 = [Esync1(DiffPos+1:end,:);Esync1(1:DiffPos,:)];%Shift based on sampling sliding
    %if DiffPos~=0
    %if DiffPos>0%If the difference is positive, left-shift...
    %Esync = ifft(fft(Esync).*exp(1j*2*pi*f*(DiffPos*Ta))); %Shift based on time change
    %Esync = [Esync(DiffPos+1:end) Esync(1:DiffPos)];       %Shift based on sampling sliding
    %else%... but if the difference is negative, right-shift
    %Esync = ifft(fft(Esync).*exp(1j*2*pi*f*(DiffPos*Ta))); %Shift based on time change
    %Esync = [Esync(end+DiffPos+1:end) Esync(1:end + ...
    %DiffPos)];%Shift based on sampling sliding
    %end
    %end

    %Because of reasons, sometimes it may be required to make a
    %synchronization process with the  end of the data stream as
    %well. This following verification check if the user set (or
    %not) a second synchronization process to be done.

    %%          Ploting the result for qualitative analizes
    %             PrintInfo(Ploting*16,t(end-SyncPos+1:end-IniSyncPos+1),Esync...
    %             (end-SyncPos+1:end-IniSyncPos+1),SyncSymb(end-SyncPos+1:end-...
    %             IniSyncPos+1),ESync1(end-SyncPos+1:end-IniSyncPos+1),ESync2(...
    %                                           end-SyncPos+1:end-IniSyncPos+1));
    %if SencondAdjust
    %%                   Synchronizing
    %This synchronization process is based on the mean expected
    %value. That means, the information mean value should be
    %within one period of symbol. Thus, the mean value of the
    %received signal is acquired and compare of the known
    %sync-word to verify if this mean value is at the right
    %possition.

    %SyncAuxEnd   = Esync(end-SyncPos+1:end-IniSyncPos+1);
    %SyncedAuxEnd = SyncSymbEnd(end-SyncPos+1:end-IniSyncPos+1);
    %SyncAuxEnd(SyncAuxEnd<0)                 = 0;              %To keep the mean value above zero anything under is neglected
    %SyncAuxEnd(SyncAuxEnd>=mean(SyncAuxEnd)) = 1;              %Adding a flag to the first sample of the received mean value
    %SyncAuxEnd(SyncAuxEnd<mean(SyncAuxEnd))  = -1;             %All the others samples at set to the lowest level

    %PosToSynEnd  = find(ismember(SyncAuxEnd,1));               %Finding where is the location of the first sample to synchronize
    %PosSyncedEnd = find(ismember(SyncedAuxEnd,1));             %Finding where is the location of the first sample to synchronize

    %DiffPosEnd = PosToSynEnd(round(end/2))-PosSyncedEnd(...
    %round(end/2));
    %if DiffPosEnd~=0
    %if DiffPosEnd>0%If positive difference, left-shift...
    %Esync = ifft(fft(Esync).*exp(1j*2*pi*f*(...
    %DiffPosEnd*Ta)));%Shift based on time change
    %Esync = [Esync(DiffPosEnd+1:end) Esync(1:...
    %DiffPosEnd)];%Shift based on sampling sliding
    %else%... but if the difference is negative, right-shift
    %Esync = ifft(fft(Esync).*exp(1j*2*pi*f*(...
    %DiffPosEnd*Ta)));%Shift based on time change
    %Esync = [Esync(end-DiffPosEnd+1:end) Esync(...
    %1:end-DiffPosEnd)];%Shift based on sampling sliding
    %end
    %end
    %end

%     if ThisCarr==125
%         IxAux  = Esync1(1:end - StuffSampels,:);
%         Olho((IxAux(:)),t(2*NumAmosCP+NPPB),(2*NumAmosCP+NPPB),1,5);
%         EyeToPlot(CurrentTest,1:length(IxAux(:))) = IxAux(:);
%         save(['savingforplotingeye' num2str(VetSnr)],'EyeToPlot','VetSnr','IfftOrSum','UsedModula','T','NPPB');
%     end
    %if ThisCarr==126
    %EyeToPlot(CurrentTest,1:length(Esync(:))) = Esync(:);
    %save(['savingforplotingeye' num2str(VetSnr)],'EyeToPlot','VetSnr','IfftOrSum','UsedModula','T','NPPB');
    %end
    %% Removing CP
    if AddCP
        IxAux  = Esync1(1:end - StuffSampels,:);
        IxAux  = reshape(IxAux,(2*NumAmosCP+NPPB),NbDPSK*CurTesSiz);
        IxAux  = IxAux(1+NumAmosCP:end-NumAmosCP,:);
        IxAux  = reshape(IxAux,NPPB*NbDPSK,CurTesSiz);
        Esync1  = IxAux;
    end

    %% Taking the sampling the EVM meassurement
    %clear IxAux;
    %PosAuxEout1 = NPPB/2:NPPB:length(Esync);                   %Varriable respossible to take just the samples at the middle of the symbol
    %PosAuxEout2 = ((NPPB/2)+(NPPB/16)):NPPB:length(Esync);      %Varriable respossible to take just the samples at the middle of the symbol
    %PosAuxEout3 = ((NPPB/2)-(NPPB/16)):NPPB:length(Esync);      %Varriable respossible to take just the samples at the middle of the symbol
    %IxAux1      = Esync(PosAuxEout1);                          %Normalizing the reference
    %IxAux2      = Esync(PosAuxEout2);    %Normalizing the reference
    %IxAux3      = Esync(PosAuxEout3);    %Normalizing the reference
    %a=a+0;
    %EvmMatRecA(1,:) = IxAux1;                       %Taking just the middle samples as references
    %EvmMatRecA(2,:) = IxAux2;                       %Taking just the middle samples as references
    %EvmMatRecA(3,:) = IxAux3;                       %Taking just the middle samples as references
    %EvmMatRecA(4,:) = IxAux1;                       %Taking just the middle samples as references
    %[EvmDBA(1), EvmPerA(1), EvmRmsA(1) ] = EvmCalc( EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux1 );
    %[EvmDBA(2), EvmPerA(2), EvmRmsA(2) ] = EvmCalc( EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux2 );
    %[EvmDBA(3), EvmPerA(3), EvmRmsA(3) ] = EvmCalc( EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux3 );
    %[EvmDBA(4), EvmPerA(4), EvmRmsA(4) ] = EvmCalc( EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux1 );
    %[EvmDBJA(1),EvmPerJA(1),EvmRmsJA(1)] = evm1(2,'dpsk',EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux1);
    %[EvmDBJA(2),EvmPerJA(2),EvmRmsJA(2)] = evm1(2,'dpsk',EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux2);
    %[EvmDBJA(3),EvmPerJA(3),EvmRmsJA(3)] = evm1(2,'dpsk',EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux3);
    %[EvmDBJA(4),EvmPerJA(4),EvmRmsJA(4)] = evm1(2,'dpsk',EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux1);
    %##########################################################################
    Esync = Esync1(1+SyncPeriod*NPPB:end-SyncPeriod*NPPB,:);
    Esync = reshape(Esync,1,size(Esync,1)*size(Esync,2));
    PosIx = NPPB/2:NPPB:size(Esync,2);%Possition of the central samples - but the number of samples per symbol is even ????
    %IxAux = Esync(PosIx,:);%From the main received signal just a few samples are taken for further evaluation
    n     = 100;%The number of boxes to be filled up on the histogram process
    Perce = 0.8;
    %At this process each eye will be evaluated separetelly. It is
    %important to remember that the PAM4 has 4 levels, which means
    %three level of decissions that are exaclty the center of the
    %eye diagram.

    %BerDPSKAux = zeros(1,size(Esync,2));
    %AberLevAux = zeros(1,size(Esync,2));
    %ValsLevAux = zeros(1,size(Esync,2));
    %for jj=1:size(Esync,2)
    IxAuxAB = Esync(PosIx);%Taking just those values relative to the uper eye
    InterAB = linspace(Perce*(min(Esync)),Perce*(max(Esync)),n);                     %Building the histogram boxes
    EyeAB = hist(IxAuxAB,InterAB);                                 %filling up the boxes with samples that fit on them.

    %What we are looking for here are not where exist the
    %occurrences of values. The eye is where there is not samples,
    %this means, the decission level is a reagion between actual
    %levels thus this reagion does not contain any signal.
    EyeAB = ~EyeAB;                                                %Changing zeros to one - Zeros compose the eye region
    %Starting variables that will be used to identify the eye
    %diagram.
    CountAB   = 1;
    SeqOnesAB = zeros(1,length(EyeAB));
    SeqFinAB  = zeros(1,length(EyeAB));
    SeqIniAB  = 1;
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
    %another way to measure the eye opening is the get all the
    %boxes and find all peaks on it, that will be a plato created
    %by the sequences of ones (which was zeros). From thos peaks,
    %the eye diagram will be the longer of them hence it will take
    %the most part of the vector that store the findpeaks result.
    %Thus, the middle of the eye will be basically the middle of
    %the peaks vector.
    LocAB = find(EyeAB);
    if ~any(LocAB)                                              %if for some reason there are no peaks, something went wrong.
        %LocAB = 1;
        LevelDec3 = 0.0;%mean(Levels(3:4));                       %by default the decission level will be set to 0.65
    else
        LevelDec3 = LevDec3;
    end
    %##########################################################################

    %[~,~,EyeOpenQ,EyeOpenHighQ,EyeOpenLowQ] = Olho(EoutQ,T,NPPB,0);
    AberLev = InterAB(SeqFinAB(LocMaxAB))-InterAB(round(SeqFinAB(LocMaxAB)-MaxValAB));
    ValsLev = LevDec3;
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
    ThisDataSize = NPPB/2:NPPB:size(Esync,2);
    %ThisDataPos  = 1:NPPB:size(Esync,2);
    Data  = zeros(1,length(ThisDataSize));                                                  %Initialization of the vector that will store the income data
    DataU = zeros(1,length(ThisDataSize));                                                  %Initialization of the vector that will store the income data
    % DataS = zeros(1,length(ThisDataSize));                                                  %Initialization of the vector that will store the income data
    for kk=1:length(ThisDataSize)%length(Esync(ThisDataSize))                                    %The comparison process will be made for each symbol period
        %MeanOfData = mean(EoutI((kk-1)+SymLocI(1)));
        MeanOfData = mean(Esync((kk-1)*NPPB+NPPB/2));
        if MeanOfData > LevelDec3%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
            Data(kk) = 1;                                 %is 1
        else                                                       %If it is the lowest level the incoming data
            Data(kk) = 0;                                 %is 0
        end
        %MeanOfData = mean(Esync((kk-1)*NPPB+(NPPB/2)));
        if MeanOfData > LevDefDpqsk%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
            DataU(kk) = 1;                                 %is 1
        else                                                       %If it is the lowest level the incoming data
            DataU(kk) = 0;                                 %is 0
        end
    %     if CalcS==1
    %         %MeanOfData = mean(Esync((kk-1)*NPPB+(NPPB/2)));
    %         if MeanOfData > EyeOpenLowI + EyeOpenI/2%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
    %             DataS(kk) = 1;                                 %is 1
    %         else                                                       %If it is the lowest level the incoming data
    %             DataS(kk) = 0;                                 %is 0
    %         end
    %     end
    end
    %%       Calculating the Bit Error Ratio (BER)
    %The final process here is to count the number of wrongdoings
    %of this whole process upon the transmited data for
    %quantitative analizes

    TxDataA = TxDataMat(ThisCarr,:);
    TxDataB = reshape(TxDataA,NbDPSK,CurTesSiz);
    TxDataC = TxDataB(1+SyncPeriod:end-SyncPeriod,:);
    TxData  = reshape(TxDataC,1,size(TxDataC,1)*size(TxDataC,2));
    % if CalcS
    %     %DataS = DataS(1+SyncPeriod:end-SyncPeriod);
    %     BitErrS = sum(xor(TxData,DataS));%Comparison between the Transmited and received and counting the differences
    %     BerDPSKS = BitErrS/((NbDPSK-(2*SyncPeriod))*CurTesSiz);%Calculating the ration of wrong bits and the total number of bits transmited
    % end
    %Data  = Data(1+SyncPeriod:end-SyncPeriod);
    %DataM  = DataM(1+SyncPeriod:end-SyncPeriod);
    %DataL  = DataL(1+SyncPeriod:end-SyncPeriod);
    %DataU  = DataU(1+SyncPeriod:end-SyncPeriod);
    %AberLevAuxI(4) = 0;
    %ValsLevAuxI(4) = 0.01;

    BitErr1  = sum(xor(TxData,Data));                        %Comparison between the Transmited and received and counting the differences
    %BitErr(2)  = sum(xor(TxData,DataM));                       %Comparison between the Transmited and received and counting the differences
    %BitErr(3)  = sum(xor(TxData,DataL));                       %Comparison between the Transmited and received and counting the differences
    BitErr4  = sum(xor(TxData,DataU));                       %Comparison between the Transmited and received and counting the differences
    BerDPSK = BitErr1/((NbDPSK-(2*SyncPeriod))*CurTesSiz);%Calculating the ration of wrong bits and the total number of bits transmited
    if BitErr4<BitErr1
        BerDPSK = BitErr4/((NbDPSK-(2*SyncPeriod))*CurTesSiz);
    end
end