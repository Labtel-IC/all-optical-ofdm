function [AberLevI,ValsLevI,AberLevQ,ValsLevQ,BerDQPSK,EoutI,EoutQ,TxDataOdd,TxDataEven,DataOdd,DataEven,LocMaxAB,MaxValAB,SeqFinAB,QMaxValAB,QLocMaxAB,QSeqFinAB,QInterAB,InterAB,IxAuxAB,QIxAuxAB]=RecDowDqpsk(EoutA,EoutB,EoutC,EoutD,T,Ta,IfftOrSum,MaxNumStag,StuffSampels,NbDQPSK,CurTesSiz,NumAmosCP,NPPB,SyncPeriod,TxDataMat,ThisCarr,LevDefDpqsk,AddCP,NumCarr)
    
    EoutI = (EoutB - EoutA);
    EoutQ = (EoutD - EoutC);
    EmeaI = mean(EoutI);
    EmeaQ = mean(EoutQ);
    EmeaI = repmat(EmeaI,size(EoutI,1),1);
    EmeaQ = repmat(EmeaQ,size(EoutQ,1),1);
    EoutI = EoutI-EmeaI;
    EoutQ = EoutQ-EmeaQ;
    EmaxI = max(EoutI);
    EmaxQ = max(EoutQ);
    EmaxI = repmat(EmaxI,size(EoutI,1),1);
    EmaxQ = repmat(EmaxQ,size(EoutQ,1),1);
    EoutI = EoutI./EmaxI;                                %Normalizing the signal
    EoutQ = EoutQ./EmaxQ;                                %Normalizing the signal
    %EoutAuxF = 20*log10(abs(fftshift(fft(Ix)./length(Ix))));
    %for jj=1:size(Ix,2)
    %VetElecPowerI(CurrentTest,ThisCarr,jj)= MeasPower(EoutI(:,jj));
    %VetElecPowerQ(CurrentTest,ThisCarr,jj)= MeasPower(EoutQ(:,jj));
    %VetElecPower(CurrentTest,ThisCarr,jj)= MeasPower(Ix(:,jj));
    %[VetOptiPower(CurrentTest,ThisCarr,jj),~]= findpeaks(EoutAuxF(:,jj),'SortStr','descend','NPeaks',1);
    %end

    %SyncAuxI   = EoutI(IniSyncPos:SyncPos);                        %Selecting just the symbol to synchronize
    %SyncAuxQ   = EoutQ(IniSyncPos:SyncPos);                        %Selecting just the symbol to synchronize
    %                 PrintInfo(Ploting*22,t(IniSyncPos:SyncPos),EoutI(IniSyncPos:SyncPos),SyncSymb(IniSyncPos:SyncPos),EoutA(IniSyncPos:SyncPos),EoutB(IniSyncPos:SyncPos));

    %SyncedAux  = SyncSymb(IniSyncPos:SyncPos);                     %Selecting just the symbol to synchronize
    %                   Synchronizing
    %This synchronization process is based on the mean expected
    %value. That means, the information mean value should be within
    %one period of symbol. Thus, the mean value of the received
    %signal is acquired and compare of the known sync-word to
    %verify if this mean value is at the right possition. Which is
    %the midel point (peak) of the highest level at the sync period
    %SyncAuxI(EoutI(IniSyncPos:SyncPos)<0.5*max((EoutI(IniSyncPos:SyncPos))))               = 0;   %To keep the mean value above zero anything under is neglected
    %SyncAuxQ(EoutQ(IniSyncPos:SyncPos)<0.5*max((EoutQ(IniSyncPos:SyncPos))))               = 0;   %To keep the mean value above zero anything under is neglected
    %SyncAuxI(SyncAuxI>=mean(SyncAuxI)) = 1;   %Adding a flag to the first sample of the received mean value
    %SyncAuxQ(SyncAuxQ>=mean(SyncAuxQ)) = 1;   %Adding a flag to the first sample of the received mean value
    %SyncAuxI(SyncAuxI<mean(SyncAuxI))  = -1;  %All the others samples at set to the lowest level
    %SyncAuxQ(SyncAuxQ<mean(SyncAuxQ))  = -1;  %All the others samples at set to the lowest level

    %PosToSynI  = find(ismember(SyncAuxI,1));                       %Finding where is the location of the samples to synchronize
    %PosToSynQ  = find(ismember(SyncAuxQ,1));                       %Finding where is the location of the samples to synchronize
    %PosSynced = find(ismember(SyncedAux,1));                       %Finding where is the location of the samples to synchronize

    %DiffPosI = ExtDel*(PosToSynI(round(end/2)) - PosSynced(round(end/2)));  %Accounting the peak (midel point) displacement
    %DiffPosQ = ExtDel*(PosToSynQ(round(end/2)) - PosSynced(round(end/2)));  %Accounting the peak (midel point) displacement
    
    sn=((1/2)*(((1/2)^MaxNumStag)-1))/((1/2)-1);
    if IfftOrSum==1
        DiffPosI = round((sn/(3/2.55^IfftOrSum))*T/Ta);
        DiffPosQ = round((sn/(3/2.55^IfftOrSum))*T/Ta);
    else
        DiffPosI = round((sn/(3/2.55^IfftOrSum))*T/Ta);
        DiffPosQ = round((sn/(3/2.55^IfftOrSum))*T/Ta);
    end
%     sn=((1/2)*(((1/2)^MaxNumStag)-1))/((1/2)-1);
%     DiffPosI = IfftOrSum*round((sn/(3/2.55^IfftOrSum))*T/Ta);
%     DiffPosQ = IfftOrSum*round((sn/(3/2.55^IfftOrSum))*T/Ta);
    EoutI = [EoutI(DiffPosI+1:end,:);EoutI(1:DiffPosI,:)];   %Shift based on sampling sliding
    EoutQ = [EoutQ(DiffPosQ+1:end,:);EoutQ(1:DiffPosQ,:)];   %Shift based on sampling sliding
    %if DiffPosI>=0%If the difference is positive, left-shift...
    %EoutI = ifft(fft(EoutI).*exp(1j*2*pi*f*(DiffPosI*Ta))); %Shift based on time change
    %EoutI = [EoutI(DiffPosI+1:end) EoutI(1:DiffPosI)];   %Shift based on sampling sliding
    %EoutQ = [EoutQ(DiffPosI+1:end) EoutQ(1:DiffPosI)];   %Shift based on sampling sliding
    %else%... but if the difference is negative, right-shift
    %EoutI = ifft(fft(EoutI).*exp(-1j*2*pi*f*(DiffPosI*Ta))); %Shift based on time change
    %EoutI = [EoutI(end+DiffPosI+1:end) EoutI(1:end+DiffPosI)]; %Shift based on sampling sliding
    %EoutQ = [EoutQ(end+DiffPosI+1:end) EoutQ(1:end+DiffPosI)]; %Shift based on sampling sliding
    %end
    %%                   Plot for Qualitative analizes
    %                              PrintInfo(Ploting*22,t(IniSyncPos:SyncPos),EoutI(IniSyncPos...
    %                              :SyncPos),SyncSymb(IniSyncPos:SyncPos),EoutA(IniSyncPos:...
    %                                                       SyncPos),EoutB(IniSyncPos:SyncPos));
    %                              PrintInfo(Ploting*22,t(IniSyncPos:SyncPos),EoutQ(IniSyncPos...
    %                              :SyncPos),SyncSymb(IniSyncPos:SyncPos),EoutC(IniSyncPos:...
    %                                                       SyncPos),EoutD(IniSyncPos:SyncPos));
    %% Removing CP
    if AddCP
        IxAux1  = EoutI(1:end - StuffSampels,:);
        IxAux1  = reshape(IxAux1,(2*NumAmosCP+NPPB),(NbDQPSK/2)*CurTesSiz);
        IxAux1  = IxAux1(1+NumAmosCP:end-NumAmosCP,:);
        IxAux1  = reshape(IxAux1,NPPB*(NbDQPSK/2),CurTesSiz);
        EoutI  = IxAux1;

        IxAux2  = EoutQ(1:end - StuffSampels,:);
        IxAux2  = reshape(IxAux2,(2*NumAmosCP+NPPB),(NbDQPSK/2)*CurTesSiz);
        IxAux2  = IxAux2(1+NumAmosCP:end-NumAmosCP,:);
        IxAux2  = reshape(IxAux2,NPPB*(NbDQPSK/2),CurTesSiz);
        EoutQ  = IxAux2;
    end
    %% Taking the sampling the EVM meassurement
    %clear IxAux;
    %PosAuxEout1 = NPPB/2:NPPB:length(EoutI);%Varriable respossible to take just the samples at the middle of the symbol
    %PosAuxEout2 = ((NPPB/2)+(NPPB/4)):NPPB:length(EoutI);%Varriable respossible to take just the samples at the middle of the symbol
    %PosAuxEout3 = ((NPPB/2)-(NPPB/4)):NPPB:length(EoutI);%Varriable respossible to take just the samples at the middle of the symbol
    %IxAux1      = EoutI(PosAuxEout1) + 1j.*EoutQ(PosAuxEout1);    %Normalizing the reference
    %IxAux2      = EoutI(PosAuxEout2) + 1j.*EoutQ(PosAuxEout2);    %Normalizing the reference
    %IxAux3      = EoutI(PosAuxEout3) + 1j.*EoutQ(PosAuxEout3);    %Normalizing the reference
    %a=a+0;
    %EvmMatRecA(1,:) = IxAux1;                       %Taking just the middle samples as references
    %EvmMatRecA(2,:) = IxAux2;                       %Taking just the middle samples as references
    %EvmMatRecA(3,:) = IxAux3;                       %Taking just the middle samples as references
    %EvmMatRecA(4,:) = IxAux1;                       %Taking just the middle samples as references
    %[EvmDBA(1), EvmPerA(1), EvmRmsA(1) ] = EvmCalc( EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux1 );
    %[EvmDBA(2), EvmPerA(2), EvmRmsA(2) ] = EvmCalc( EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux2 );
    %[EvmDBA(3), EvmPerA(3), EvmRmsA(3) ] = EvmCalc( EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux3 );
    %[EvmDBA(4), EvmPerA(4), EvmRmsA(4) ] = EvmCalc( EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux1 );
    %[EvmDBJA(1),EvmPerJA(1),EvmRmsJA(1)] = evm1(4,'dpsk',EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux1);
    %[EvmDBJA(2),EvmPerJA(2),EvmRmsJA(2)] = evm1(4,'dpsk',EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux2);
    %[EvmDBJA(3),EvmPerJA(3),EvmRmsJA(3)] = evm1(4,'dpsk',EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux3);
    %[EvmDBJA(4),EvmPerJA(4),EvmRmsJA(4)] = evm1(4,'dpsk',EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux1);
    %%
    %                 figure;plot(t(1:length(E_rec3)),abs(E_rec3),t(1:length(EoutQ)),EoutQ,t(1:length(EoutI)),EoutI);set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
    %One important step on this project was the confirmation of our
    %results from a closed and theoretical equation that relates
    %the income optical field witht its respective data (as
    %descrived above in the receiving process). The result from
    %this equation was further compared with the result from the
    %MZ-Interferometer as an proof of concept. This equation can be
    %found at the book of Optical Fiber Telecommunications V B,
    %which one of the authors is Ivan P. Kaminow at the page 144.
    %                                 taux = t(1:length(E_rec3));
    %                                 faux = time2freq(taux);
    %                                 Ui = real(exp(1j*pi/4).*E_rec3.*conj(ifft(fft(E_rec3).*exp(-...
    %                                     1j*2*pi*faux*T))));%The data at odd position
    %                                 Ui = Ui./max(abs(Ui));                                         %Normalizing the signal
    %                                 Uq = imag(exp(1j*pi/4).*E_rec3.*conj(ifft(fft(E_rec3).*exp(-...
    %                                     1j*2*pi*faux*T))));%The data at the even position
    %                                 Uq = Uq./max(abs(Uq));                                         %Normalizing the signal
    %##########################################################################
    EoutI = EoutI(1+SyncPeriod*NPPB:end-SyncPeriod*NPPB,:);
    EoutI = reshape(EoutI,1,size(EoutI,1)*size(EoutI,2));
    EoutQ = EoutQ(1+SyncPeriod*NPPB:end-SyncPeriod*NPPB,:);
    EoutQ = reshape(EoutQ,1,size(EoutQ,1)*size(EoutQ,2));
    PosIx = NPPB/2:NPPB:length(EoutI);                                %Possition of the central samples - but the number of samples per symbol is even ????
    
    %IxAux = EoutI(PosIx);                                             %From the main received signal just a few samples are taken for further evaluation
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
    %What we are looking for here are not where exist the
    %occurrences of values. The eye is where there is not samples,
    %this means, the decission level is a reagion between actual
    %levels thus this reagion does not contain any signal.
    EyeAB = ~EyeAB;                                                %Changing zeros to one - Zeros compose the eye region
    QEyeAB = ~QEyeAB;                                                %Changing zeros to one - Zeros compose the eye region
    %Starting variables that will be used to identify the eye
    %diagram.
    CountAB   = 1;
    QCountAB   = 1;
    SeqOnesAB = zeros(1,length(EyeAB));
    QSeqOnesAB = zeros(1,length(EyeAB));
    SeqFinAB  = zeros(1,length(EyeAB));
    QSeqFinAB  = zeros(1,length(EyeAB));
    SeqIniAB  = 1;
    QSeqIniAB  = 1;
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
        LevelDec3 = 0.1;%mean(Levels(3:4));                       %by default the decission level will be set to 0.65
    else
        LevelDec3 = LevDec3;
    end
    QLocAB = find(QEyeAB);
    if ~any(QLocAB)                                              %if for some reason there are no peaks, something went wrong.
        %QLocAB = 1;
        QLevelDec3 = 0.1;%mean(Levels(3:4));                       %by default the decission level will be set to 0.65
    else
        QLevelDec3 = QLevDec3;
    end
    %##########################################################################

    AberLevI = InterAB(SeqFinAB(LocMaxAB))-InterAB(round(SeqFinAB(LocMaxAB)-MaxValAB));
    ValsLevI = LevDec3;
    AberLevQ = QInterAB(QSeqFinAB(QLocMaxAB))-QInterAB(round(QSeqFinAB(QLocMaxAB)-QMaxValAB));
    ValsLevQ = QLevDec3;

    %% Ploting the result for qualitative analizes

    %if ThisCarr==126
    %EyeToPlot(CurrentTest,1:length([EoutI EoutQ])) = [EoutI EoutQ];
    %save(['savingforplotingeye' num2str(VetSnr)],'EyeToPlot','VetSnr','IfftOrSum','UsedModula','T','NPPB');
    %end
    TxDataOddA = TxDataMat(ThisCarr,:);
    TxDataOddB = reshape(TxDataOddA,NbDQPSK/2,CurTesSiz);
    TxDataOddC = TxDataOddB(1+SyncPeriod:end-SyncPeriod,:);
    TxDataOdd  = reshape(TxDataOddC,1,size(TxDataOddC,1)*size(TxDataOddC,2));
    TxDataEvenA = TxDataMat(ThisCarr+NumCarr,:);
    TxDataEvenB = reshape(TxDataEvenA,NbDQPSK/2,CurTesSiz);
    TxDataEvenC = TxDataEvenB(1+SyncPeriod:end-SyncPeriod,:);
    TxDataEven  = reshape(TxDataEvenC,1,size(TxDataEvenC,1)*size(TxDataEvenC,2));
    %%  Ploting some results for qualitative analizes
    %             PrintInfo(Ploting*28,t(1:length(EoutI)),Txaux1,Txaux2,EoutI,...
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
    %This process is realized inside the function Olho_mex.
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
    ThisDataSize = NPPB/2:NPPB:length(EoutI);
    %ThisDataPos  = 1:NPPB:length(EoutI);
    DataOdd = zeros(1,length(ThisDataSize));                                                  %Initialization of the vector that will store the income data
    DataOddU = zeros(1,length(ThisDataSize));                                                  %Initialization of the vector that will store the income data
    for kk=1:length(ThisDataSize)%length(EoutI(ThisDataSize))                                    %The comparison process will be made for each symbol period
        %MeanOfData = mean(EoutI((kk-1)+SymLocI(1)));
        MeanOfData = mean(EoutI((kk-1)*NPPB+NPPB/2));
        if MeanOfData > LevelDec3%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
            DataOdd(kk) = 1;                                 %is 1
        else                                                       %If it is the lowest level the incoming data
            DataOdd(kk) = 0;                                 %is 0
        end
        %MeanOfData = mean(EoutI((kk-1)*NPPB+(NPPB/2)));
        if MeanOfData > LevDefDpqsk%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
            DataOddU(kk) = 1;                                 %is 1
        else                                                       %If it is the lowest level the incoming data
            DataOddU(kk) = 0;                                 %is 0
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
    ThisDataSize = NPPB/2:NPPB:length(EoutQ);
    %ThisDataPos  = 1:NPPB:length(EoutQ);
    DataEven = zeros(1,length(ThisDataSize));                                                 %Initialization of the vector that will store the income data
    DataEvenU = zeros(1,length(ThisDataSize));                                                 %Initialization of the vector that will store the income data
    for kk=1:length(ThisDataSize)%length(EoutQ(ThisDataSize))                                    %The comparison process will be made for each symbol period
        %                 MeanOfData = mean(EoutQ((kk-1)+SymLocI(1)));
        MeanOfData = mean(EoutQ((kk-1)*NPPB+NPPB/2));
        if MeanOfData > QLevelDec3%EyeOpenLowQ+EyeOpenQ/2                     %If it is the uper level the incoming data
            DataEven(kk) = 1;                               %is 1
        else                                                       %If it is the lowest level the incoming data
            DataEven(kk) = 0;                               %is 0
        end
        %MeanOfData = mean(EoutQ((kk-1)*NPPB+(NPPB/2)));
        if MeanOfData > LevDefDpqsk%EyeOpenLowI+EyeOpenI/2                     %If it is the uper level the incoming data
            DataEvenU(kk) = 1;                                 %is 1
        else                                                       %If it is the lowest level the incoming data
            DataEvenU(kk) = 0;                                 %is 0
        end
    end
    %%       Calculating the Bit Error Ratio (BER)
    %The final process here is to count the number of wrongdoings
    %of this whole process upon the transmited data for
    %quantitative analizes

    BitErrOdd1  = sum(xor(TxDataOdd,DataOdd));                      %Comparison between the Transmited and received and counting the differences
    BitErrEven1 = sum(xor(TxDataEven,DataEven));                    %Comparison between the Transmited and received and counting the differences
    BitErrOdd4  = sum(xor(TxDataOdd,DataOddU));                      %Comparison between the Transmited and received and counting the differences
    BitErrEven4 = sum(xor(TxDataEven,DataEvenU));                    %Comparison between the Transmited and received and counting the differences
    BerDQPSK = (BitErrOdd1+BitErrEven1)/(((NbDQPSK)-(2*SyncPeriod))*CurTesSiz);%Calculating the ration of wrong bits and the total number of bits transmited
    if (BitErrOdd4+ BitErrEven4) < (BitErrOdd1+ BitErrEven1)
        BerDQPSK = (BitErrOdd4+BitErrEven4)/(((NbDQPSK)-(2*SyncPeriod))*CurTesSiz);%Calculating the ration of wrong bits and the total number of bits transmited
    end
    
end