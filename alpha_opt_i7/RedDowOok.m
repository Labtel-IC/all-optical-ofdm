function [AberLev,ValsLev,Inter,SeqFin,LocMax,MaxVal,TxData,EoutCorrD,EoutCorr2,EoutCorr,BerOOK,LocMax2,SeqFin2,MaxVal2,Inter2] = RedDowOok(Ix,TxDataMat,NPPB,ThisCarr,Nb,NumBitDesc,SyncPeriod,CurTesSiz,EyeOpenLow,EyeOpen)

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
    %EyeLoc = [2 3];
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

SeqOnes  = zeros(1,length(EyeAB));
SeqOnes2  = zeros(1,length(EyeAB));

SeqIni  = 1;
SeqIni2  = 1;

SeqFin  = zeros(1,length(EyeAB));
SeqFin2  = zeros(1,length(EyeAB));
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
    LevDec = 0.50;                                            %the decission level will be by default 0.67. Also, the other variables will
    LocMax = 1;                                              %will be set with values to not cause errors in the future
    SeqFin(1)=2;
    MaxVal = 0;
    Inter(1)=LevDec;
else                                                           %if a sequency was found the middle of the eye will be the middle of the sequency
    if (SeqFin(LocMax)-MaxVal/2)<1                       %if for some reason the final element of the sequency minus the half of its
        LevDec = 0.50;                                         %length results in a negative value, something went very wrong, and by
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
if ~any(Loc)                                              %if for some reason there are no peaks, something went wrong.
    Loc = 1;
    LevelDec = 0.0;%mean(Levels(3:4));                       %by default the decission level will be set to 0.65
else
    LevelDec = LevDec;
end
Loc2 = find(EyeCD);
if ~any(Loc2)                                              %if for some reason there are no peaks, something went wrong.
    Loc2 = 1;
    LevelDec2 = 0.0;%mean(Levels(3:4));                       %by default the decission level will be set to 0.65
else
    LevelDec2 = LevDec2;
end
%##########################################################################

AberLev = Inter(SeqFin(LocMax))-Inter(round(SeqFin(LocMax)-MaxVal));
ValsLev = LevDec;
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
%%      Actualy Receiving Data
%Once the signal was processed the next step is through a
%comparator decide the actual information received.
%ThisDataPos = 1:NPPB:length(Ix);
ThisDataSize = NPPB/2:NPPB:length(Ix);
EoutCorr = zeros(1,length(ThisDataSize));                                                 %Initialization of the vector that will store the income data
EoutCorrD = zeros(1,length(ThisDataSize));                                                 %Initialization of the vector that will store the income data
EoutCorr2 = zeros(1,length(ThisDataSize));                                                 %Initialization of the vector that will store the income data
for kk=1:length(ThisDataSize)%length(Ix(ThisDataSize))                                       %The comparison process will be made for each symbol period
    %An small portion of the income signal is take for
    %evaluation by measuring the avarage value of the samples
    %taken
    %                 CalcMean = mean((Ix((kk-1)+SymLoc(1))));
    CalcMean = mean((Ix((kk-1)*NPPB+NPPB/2)));
    %Verifying the interval for each symbol received.
    if CalcMean >= LevDec%EyeOpenLow+EyeOpen/2                        %If it is the uper level the incoming data
        EoutCorr(kk) = 1;                               %is 1
    else                                                       %If it is the lowest level the incoming data
        EoutCorr(kk) = 0;                               %is 0
    end
    CalcMean2 = mean((Ix((kk-1)*NPPB+NPPB/2)));
    %Verifying the interval for each symbol received.
    if CalcMean2 >= LevDec2%EyeOpenLow+EyeOpen/2                        %If it is the uper level the incoming data
        EoutCorr2(kk) = 1;                               %is 1
    else                                                       %If it is the lowest level the incoming data
        EoutCorr2(kk) = 0;                               %is 0
    end
    if CalcMean2 >= EyeOpenLow+EyeOpen/2                        %If it is the uper level the incoming data
        EoutCorrD(kk) = 1;                               %is 1
    else                                                       %If it is the lowest level the incoming data
        EoutCorrD(kk) = 0;                               %is 0
    end
end

TxDataA = TxDataMat(ThisCarr,:);
TxDataB = reshape(TxDataA,(Nb-NumBitDesc),CurTesSiz);
TxDataC = TxDataB(1+SyncPeriod:end-SyncPeriod,:);
TxData  = reshape(TxDataC,1,size(TxDataC,1)*size(TxDataC,2));
%%       Calculating the Bit Error Ratio (BER)
%The final process here is to count the number of wrongdoings
%of this whole process upon the transmited data for
%quantitative analizes
BitErr = sum(xor(TxData,EoutCorr));%Comparison between the Transmited and received and counting the differences
BitErr2 = sum(xor(TxData,EoutCorr2));%Comparison between the Transmited and received and counting the differences
BitErrD = sum(xor(TxData,EoutCorrD));%Comparison between the Transmited and received and counting the differences
if BitErr2<=BitErr
    if BitErr2<=BitErrD
        BerOOK = BitErr2/(((Nb-NumBitDesc)-2*SyncPeriod)*CurTesSiz);%Calculating the ration of wrong bits and the total number of bits transmited
        AberLev = Inter2(SeqFin2(LocMax2))-Inter2(round(SeqFin2(LocMax2)-MaxVal2));
        ValsLev = LevDec2;
    else
        BerOOK  = BitErrD/(((Nb-NumBitDesc)-2*SyncPeriod)*CurTesSiz);%Calculating the ration of wrong bits and the total number of bits transmited
        AberLev = EyeOpen;
        ValsLev = EyeOpenLow+EyeOpen/2;
    end
else
    if BitErr<=BitErrD
        BerOOK = BitErr/(((Nb-NumBitDesc)-2*SyncPeriod)*CurTesSiz);%Calculating the ration of wrong bits and the total number of bits transmited
    else
        BerOOK  = BitErrD/(((Nb-NumBitDesc)-2*SyncPeriod)*CurTesSiz);%Calculating the ration of wrong bits and the total number of bits transmited
        AberLev = EyeOpen;
        ValsLev = EyeOpenLow+EyeOpen/2;
    end
end
end