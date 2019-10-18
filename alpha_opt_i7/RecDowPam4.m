function [LevDec1,LevDec2,LevDec3,Ix,Ber4PAM,DecLevDef3,DecLevDef2,DecLevDef1,TxData,IxRecDef,IxRec,AberLev1,AberLev2,AberLev3,ValsLev1,ValsLev2,ValsLev3,ValsLev21,ValsLev22,ValsLev23,InterAB,InterCD,InterEF,SeqFinAB,SeqFinCD,SeqFinEF,LocMaxAB,LocMaxCD,LocMaxEF,MaxValAB,MaxValCD,MaxValEF,Levels] = RecDowPam4(Ix1,T,Ta,MaxNumStag,StuffSampels,NumAmosCP,NPPB,CurTesSiz,Nb4Pam,IntervalStep,MinDist,DecLevDef1,DecLevDef2,DecLevDef3,TxDataMat,ThisCarr,IfftOrSum,AddCP,SyncPeriod,DecMod,FFTSplit)

sn=((1/2)*(((1/2)^MaxNumStag)-1))/((1/2)-1);
AuxSyncCorr = IfftOrSum*round((sn/(2/2^FFTSplit))*T/Ta);
Ix1 = [Ix1(AuxSyncCorr+1:end,:);Ix1(1:AuxSyncCorr,:)];
%The difference between the PossitionTosynchronize and
%Possitionsynchronized will be used to correct the time
%shifting on the transmition and reception process.
%if AuxSyncCorr>=0%If the difference is positive, left-shift...
%                 Ix = ifft(fft(Ix).*exp(1j*2*pi*f*AuxSyncCorr));          %Shift based on time change
%Ix = [Ix(AuxSyncCorr+1:end) Ix(1:AuxSyncCorr)];            %Shift based on sampling sliding
%else%... but if the difference is negative, right-shift
%                 Ix = ifft(fft(Ix).*exp(-1j*2*pi*f*AuxSyncCorr));         %Shift based on time change
%Ix = [Ix(end+AuxSyncCorr+1:end) Ix(1:end + AuxSyncCorr)];  %Shift based on sampling sliding
%end
%if SencondAdjust
%For some reason that we could not understand sometimes the
%time (sampling) sliding of the signal is not equal
%throught the data stream. Thus, the second part of the
%synchronism process will be turn ON or OFF according to
%the user's will.
%AuxSyncEnd     = (Ix(end-SyncPos+1:end-IniSyncPos-1));     %Selecting the sync-word within the received signal
%SyncSymbEndAux = SyncSymbEnd(end-SyncPos+1:end-...
%IniSyncPos-1);%Selecting the sync-word within the known signal
%AuxSyncEnd1 = AuxSyncEnd;
%AuxSyncEnd1(AuxSyncEnd1<0) = 0;                            %To keep the mean value above zero anything under is neglected
%AuxSyncEnd1(AuxSyncEnd1>=mean(AuxSyncEnd1)) = 1;           %Adding a flag to the first sample of the received mean value
%AuxSyncEnd1(AuxSyncEnd1<mean(AuxSyncEnd1)) = -1;           %All the others samples at set to the lowest level
%AuxSyncEnd2 = SyncSymbEndAux;
%AuxSyncEnd2(AuxSyncEnd2>=mean(AuxSyncEnd2)) = 1;           %Adding a flag to the first sample of the known mean value
%AuxSyncEnd2(AuxSyncEnd2<mean(AuxSyncEnd2)) = -1;           %All the others samples at set to the lowest level


%PosToSynEnd  = find(ismember(AuxSyncEnd1,1));              %Finding where is the location of the first sample to synchronize
%PosSynEnd = find(ismember(AuxSyncEnd2,1));                 %Finding where is the location of the first sample to synchronize

%AuxSyncCorrEnd = PosToSynEnd(round(end/2))-PosSynEnd(...
%round(end/2));

%The difference between the PossitionTosynchronize and
%Possitionsynchronized will be used to correct the time
%shifting on the transmition and reception process.
%if AuxSyncCorrEnd>=0%If possitive difference, left-shift...
%                 Ix = ifft(fft(Ix).*exp(1j*2*pi*f*AuxSyncCorrEnd));   %Shift based on time change
%Ix = [Ix(AuxSyncCorrEnd+1:end) Ix(1:AuxSyncCorrEnd)];  %Shift based on sampling sliding
%else%... but if the difference is negative, right-shift
%Ix = ifft(fft(Ix).*exp(-1j*2*pi*f*AuxSyncCorrEnd));  %Shift based on time change
%Ix = [Ix(end-AuxSyncCorrEnd+1:end) Ix(1:end - ...
%AuxSyncCorrEnd)];%Shift based on sampling sliding
%end
%end

%%          Ploting the result for qualitative analizes
%                             PrintInfo(Ploting*35,t(end-SyncPos+1:end-IniSyncPos+1),Ix(...
%                                 end-SyncPos+1:end-IniSyncPos+1),SyncSymb(end-SyncPos+1:...
%                                                                         end-IniSyncPos+1));
%% Removing CP
if AddCP==1 
    IxAux = Ix1(1:end - StuffSampels,:);
    IxAux = reshape(IxAux,(2*NumAmosCP+NPPB),(Nb4Pam/2)*CurTesSiz);
    IxAux = IxAux(1+NumAmosCP:end-NumAmosCP,:);
    IxAux = reshape(IxAux,NPPB*Nb4Pam/2,CurTesSiz);
    Ix1    = IxAux;
end
%% Taking the sampling the EVM meassurement
Ix = Ix1(1+SyncPeriod*NPPB:end-SyncPeriod*NPPB,:);
Ix = reshape(Ix,1,size(Ix,1)*size(Ix,2));
%PosAuxEout = NPPB/2:NPPB:length(Ix);%Varriable respossible to take just the samples at the middle of the symbol
%IxAux      = Ix(PosAuxEout);                           %Normalizing the reference
%RxSymbAmos(CurrentTest,ObsCarrPos==ThisCarr,:) = IxAux;%RxSymbAmosUp = [];
%EvmMatRec(ObsCarrPos==ThisCarr,:) = IxAux;                       %Taking just the middle samples as references
%[EvmDB(CurrentTest,ThisCarr), EvmPer(CurrentTest,ThisCarr), EvmRms(CurrentTest,ThisCarr) ] = EvmCalc( EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux );
%[EvmDBJ(CurrentTest,ThisCarr),EvmPerJ(CurrentTest,ThisCarr),EvmRmsJ(CurrentTest,ThisCarr)] = evm1(4,'pam',EvmMatRef(ObsCarrPos==ThisCarr,:),IxAux);
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
Interval = linspace(min(Ix),max(Ix),IntervalStep);
%Therefore, the MATLAB hist function returns the number of
%occurrence of each interval.
EyeMax = hist(Ix,Interval);
EyeMaxaux = [0 EyeMax 0];                                      %Zeros are added at the EyeMax to auxiliate the finding peaks process
[~,EyeLoc] = findpeaks(EyeMaxaux,'MinPeakDistance',MinDist,'SortStr','descend','NPeaks',4);%,'MinPeakHeight',mean(EyeMax));%The peaks on the Eye profile will be the levels at the Eyes limit
if length(EyeLoc)<4
    [~,EyeLoc] = findpeaks(EyeMaxaux,'MinPeakDistance',MinDist*0.8,'SortStr','descend','NPeaks',4);%,'MinPeakHeight',mean(EyeMax)/4);%The peaks on the Eye profile will be the levels at the Eyes limit
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
    %EyeLoc = [2 3 4 5];
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

SeqOnesAB = zeros(1,length(EyeAB));
SeqOnesCD = zeros(1,length(EyeAB));
SeqOnesEF = zeros(1,length(EyeAB));

SeqFinAB  = zeros(1,length(EyeAB));
SeqFinCD  = zeros(1,length(EyeAB));
SeqFinEF  = zeros(1,length(EyeAB));

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
if ~any(LocAB)                                              %if for some reason there are no peaks, something went wrong.
    LocAB = 1;
    LevelDec3 = 0.6953;%mean(Levels(3:4));                       %by default the decission level will be set to 0.65
else
    if DecMod==1%InterAB(LocAB(round(end/2)))<=LevDec3
        LevelDec3 = LevDec3;
    else
        LevelDec3 = InterAB(LocAB(round(end/2)));
    end
end

LocCD = find(EyeCD);
if ~any(LocCD)
    LocCD = 1;
    LevelDec2 = 0.4013;%mean(Levels(2:3));
else
    if DecMod==1%InterCD(LocCD(round(end/2)))>=LevDec2
        LevelDec2 = LevDec2;
    else
        LevelDec2 = InterCD(LocCD(round(end/2)));
    end
end

LocEF = find(EyeEF);
if ~any(LocEF)
    LocEF = 1;
    LevelDec1 = 0.1877;%mean(Levels(1:2));
else
    if DecMod==1%InterEF(LocEF(round(end/2)))<=LevDec1
        LevelDec1 = LevDec1;
    else
        LevelDec1 = InterEF(LocEF(round(end/2)));
    end
end

AberLev1  = abs(InterAB(SeqFinAB(LocMaxAB)-1) - InterAB(SeqFinAB(LocMaxAB)-MaxValAB+1));
AberLev2  = abs(InterCD(SeqFinCD(LocMaxCD)-1) - InterCD(SeqFinCD(LocMaxCD)-MaxValCD+1));
AberLev3  = abs(InterEF(SeqFinEF(LocMaxEF)-1) - InterEF(SeqFinEF(LocMaxEF)-MaxValEF+1));
ValsLev1  = LevDec3;
ValsLev2  = LevDec2;
ValsLev3  = LevDec1;
ValsLev21 = InterAB(LocAB(round(length(LocAB)/2)));
ValsLev22 = InterCD(LocCD(round(length(LocCD)/2)));
ValsLev23 = InterEF(LocEF(round(length(LocEF)/2)));
%%      Actualy Receiving Data
%Once the signal was processed the next step is through a
%comparator decide the actual information received.
ThisDataSize = NPPB/2:NPPB:length(Ix);
IxRec    = zeros(1,2*size(ThisDataSize,2));%Initialization of the vector that will store the income data
IxRecDef = zeros(1,2*size(ThisDataSize,2));%Initialization of the vector that will store the income data
ContBit  = 1;
ContBit1  = 1;
for kk=1:length(ThisDataSize)%NPPB:length(Ix)                                       %The comparison process will be made for each symbol period
    %                 midaux = round(mean(SymLoc(1:round(end/2))));
    midaux = NPPB/2;%SymLoc(1);
    aux1 = Ix((kk-1)*NPPB+midaux+1);     %An small portion of the income signal is take for evaluation
    MeanRec = mean(aux1);                                      %Measuring the avarage value of the samples taken
    %Verifying the interval for each symbol received.
    if MeanRec <= LevelDec1                                    %If it is the lowest level the incoming data
        IxRec(ContBit) = 0;
        ContBit  = 1 + ContBit;
        IxRec(ContBit) = 0;
        ContBit  = 1 + ContBit;
        %IxRec = [IxRec 0 0];                                   %is 01 (1)
    elseif (MeanRec <= LevelDec2)&&(MeanRec > LevelDec1)       %If it is the second level the incoming data
        IxRec(ContBit) = 0;
        ContBit  = 1 + ContBit;
        IxRec(ContBit) = 1;
        ContBit  = 1 + ContBit;
        %IxRec = [IxRec 0 1];                                   %is 00 (0)
    elseif (MeanRec <= LevelDec3)&&(MeanRec > LevelDec2)       %If it is the tird level the incoming data
        IxRec(ContBit) = 1;
        ContBit  = 1 + ContBit;
        IxRec(ContBit) = 1;
        ContBit  = 1 + ContBit;
        %IxRec = [IxRec 1 1];                                   %is 10 (2)
    elseif MeanRec > LevelDec3                                 %If it is the uper level the incoming data
        IxRec(ContBit) = 1;
        ContBit  = 1 + ContBit;
        IxRec(ContBit) = 0;
        ContBit  = 1 + ContBit;
        %IxRec = [IxRec 1 0];                                   %is 11 (3)
    else                                                       %If for some misteriose reason neither of previous verification were sucedded
        IxRec(ContBit) = 0;
        ContBit  = 1 + ContBit;
        IxRec(ContBit) = 0;
        ContBit  = 1 + ContBit;
        %IxRec = [IxRec 0 0];                                   %by default the current data is set to be 00 (0)
    end
    
    %Verifying the interval for each symbol received.
    if MeanRec <= DecLevDef1                                   %If it is the lowest level the incoming data
        IxRecDef(ContBit1) = 0;
        ContBit1  = 1 + ContBit1;
        IxRecDef(ContBit1) = 0;
        ContBit1  = 1 + ContBit1;
        %IxRecDef = [IxRecDef 0 0];                             %is 01 (1)
    elseif (MeanRec <= DecLevDef2)&&(MeanRec > DecLevDef1)     %If it is the second level the incoming data
        IxRecDef(ContBit1) = 0;
        ContBit1  = 1 + ContBit1;
        IxRecDef(ContBit1) = 1;
        ContBit1  = 1 + ContBit1;
        %IxRecDef = [IxRecDef 0 1];                             %is 00 (0)
    elseif (MeanRec <= DecLevDef3)&&(MeanRec > DecLevDef2)     %If it is the tird level the incoming data
        IxRecDef(ContBit1) = 1;
        ContBit1  = 1 + ContBit1;
        IxRecDef(ContBit1) = 1;
        ContBit1  = 1 + ContBit1;
        %IxRecDef = [IxRecDef 1 1];                             %is 10 (2)
    elseif MeanRec > DecLevDef3                                %If it is the uper level the incoming data
        IxRecDef(ContBit1) = 1;
        ContBit1  = 1 + ContBit1;
        IxRecDef(ContBit1) = 0;
        ContBit1  = 1 + ContBit1;
        %IxRecDef = [IxRecDef 1 0];                             %is 11 (3)
    else                                                       %If for some misteriose reason neither of previous verification were sucedded
        IxRecDef(ContBit1) = 0;
        ContBit1  = 1 + ContBit1;
        IxRecDef(ContBit1) = 0;
        ContBit1  = 1 + ContBit1;
        %IxRecDef = [IxRecDef 0 0];                             %by default the current data is set to be 00 (0)
    end
end

TxDataA = TxDataMat(ThisCarr,:);
TxDataB = reshape(TxDataA,Nb4Pam,CurTesSiz);
TxDataC = TxDataB(1+2*SyncPeriod:end-2*SyncPeriod,:);
TxData  = reshape(TxDataC,1,size(TxDataC,1)*size(TxDataC,2));
BitErr = sum(xor(TxData,IxRec));%Comparison between the Transmited and received and counting the differences
BitErrAux1 = BitErr;
BitErrAux2 = sum(xor(TxData,IxRecDef));%Comparison between the Transmited and received and counting the differences
if BitErr ~= 0
    if BitErrAux2<BitErrAux1
        BitErr = BitErrAux2;
    end
else
    DecLevDef1 = LevelDec1;
    DecLevDef2 = LevelDec2;
    DecLevDef3 = LevelDec3;
end
Ber4PAM = BitErr/((Nb4Pam-(4*SyncPeriod))*CurTesSiz);

end

