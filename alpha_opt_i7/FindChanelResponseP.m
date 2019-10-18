function [ ChanelEqualizer ] = FindChanelResponseP(OfdMod,M,NFFT,Ns,NumFra,MZ_Input_File,SelecFilt,Medium,DifNsN,fin,FBWD,Order,ExtSam,FiberDelay,equ1,equ2)
%%                  FindChanelResponse
    %This function is resposible to find the equalizer for each carrier
    if nargin<14
        FiberDelay = [];
    end
    VetSnr = 1e9;
    Selecting = 1;
%% Loading data
    MainInputData;
    load([OfcName '.mat']);
    InitCarr     = InitCarrDo;
    MaxNumStagTx = nextpow2(NumCarr+1);                                    %Next number, power of 2, by the given number of carriers
    % MaxNumCarrTx = 2^MaxNumStagTx;                                       %With the number of stages the actual nubmer of carriers can be calculated
    MaxNumStagT  = nextpow2(NumCarr);                                      %Next number, power of 2, by the given number of carriers
    MaxNumStag = nextpow2(NumCarr);                                        %Next number, power of 2, by the given number of carriers
%     MaxNumCarr = 2^MaxNumStag;                                           %With the number of stages the actual nubmer of carriers can be calculated
    if SelecSetUP                                                          %Vefify is the split process will be by filter...
        EoutTx=SelectEachCarrier(Eout,NumCarr,f,fin,FBWD,Order,fc);        %This function is resposible to split each carrier
        VetThisCarrTx = (RefCarr-1)+1:(RefCarr-1)+NumCarr;                 %Keeping track of the right carrier sequence
    else                                                                   %... or by the OFFT
%As the OFFT has a periodic response the OFC needs to be constrained other 
%whise carrier multiple carrier may interfir with other channels, This
%first selection was done with a broad pass-band filter. A higher OFFT
%order can also be used, although it may increase the computational time
%and my not be exactly feasible in the real world.
        if Selecting                                                    
            EoutT = ifft(fft(Eout).*SelecFilt);
%                 PrintInfo(PlotingThis*1,f,EoutT);
%                 a=0;
        else
            EoutT = Eout;
        end                                                                 
        [EoutTx,~,VetThisCarrTx]=OpticalFFTN(t,T,MaxNumStagTx,EoutT);      %This function is resposible to split each carrier
    end
%% Transmiting DownStream
    for kk=InitCarrDo:CarrPass:NumCarr                                     %Generating different data for each carrier
        TxData          = zeros(1,Nb);                                     %Generation random information
        TxData(end/2)   = 1;
        TxDataMat(kk,:) = TxData;                                          %Saving data for future comparison
%         TxDataPara      = TxData.';                                        %Converting from serial to paralel 
%         switch OfdMod                                                      %Sellect which modulation was used
%             case 'QAM'
%                 TxSymb  = qammod(TxDataPara,M);                            %Modulating information by QAM
%             otherwise
%                 TxSymb  = dpskmod(TxDataPara,M);                           %Modulating information DPSK as default
%         end
%         TxSymbConj      = conj(flipud(TxSymb));
%         TxSymbH = [zeros(1,NumFra);TxSymb;zeros(1,NumFra);TxSymbConj];     %Hermitian Symetri
%         TxSymbC1 = [TxSymb;TxSymb(1:DifNsN,:)];
%         TxSymbConjC      = conj(flipud(TxSymbC1));
% %         TxSymbC = [zeros(1,NumFra);TxSymbC1;zeros(1,NumFra);TxSymbConjC];
% %         TxSymbC = [TxSymbH(1:size(TxSymbH)/2,:);TxSymbH(2:DifNsN,:);0;0;TxSymbH(end-DifNsN+2:end,:);TxSymbH(1+size(TxSymbH)/2:end,:);];
%         TxSymbC = zeros(NFFT,NumFra);
% %     %                 figure;hold all;plot(TxSymbH);
%         TxSymbC(1:size(TxSymbH)/2,:) = TxSymbH(1:size(TxSymbH)/2,:);
%         TxSymbC(end-(size(TxSymbH)/2)+1:end,:) = TxSymbH(1+size(TxSymbH)/2:end,:);
%     %                 plot(TxSymbC);
%         TxSymbMod       = ifft(TxSymbC,NFFT);
%         AuxTxSymbMod(1,:)  = TxSymbMod(:); 
        SigTx  = rectpulse(TxData,NPPB);
        SigTx(end/2) = 1;
        TxSig = SigTx;
%         TxSig = SigTx.';
%Assigning the eletrical signal to one drive of the MZM - The aplitude of 
%the signal can be controlled by the variable DatGai, which can be 
%understood as an gain for the eletrical signal or an atenuation. The 
%second signal will be similar with the only difference a phase shift of pi
        TxSig = 0.95*(TxSig./max(abs(TxSig)));
        U.U1t = TxSig;
        U.U2t = exp(-1j*pi).*TxSig;
%As both signals will have the mostrly the same characteristics with the 
%only difference the phase shift of 180 degress. The MZM-I will be working 
%on the Push-Pull configuration. It is necessary to reduce the Chirp noise
%to zero.
        [EoutModAux,~]=Mach_Zehnder_Modulator_simplificado(t,EoutTx(...
                         VetThisCarrTx==(RefCarr-1+kk),:),U,MZ_Input_File);%Modulating individual carriers
        EoutModTem(ObsCarrPos(kk),1:length(EoutModAux)) = EoutModAux;      %Adding the current Modulated output to the final OutPut
%             if kk==1
%                 figure;hold all;grid on;
%                 plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
%                 set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%             end
        if (~mod(CarrPass,2))&&(SendingUp)
            clear EoutModAux;
            EoutModAux = EoutTx(VetThisCarrTx==(RefCarr-1+kk+1),:);
            EoutModTem(ObsCarrPos(kk+1),1:length(EoutModAux)) = EoutModAux;%Adding the current Modulated output to the final OutPut
        end
    end

    if IfftOrSum
        if size(EoutModTem,1)>1
            [EoutMod,~,~] = OpticalIFFTN(t,T,MaxNumStagT,EoutModTem);
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


    switch Medium
        case 'B2B'
            EoutRec = EoutMod;
        case 'Fiber'
            [EoutRec,~,~]=Fibra_Monomodo1(t,EoutMod,lambdac,T,FiberLength,0,f);
    %         PrintInfo(PlotingThis*11,EoutRec.*conj(EoutRec)./max(abs(EoutRec...
    %                                                  .*conj(EoutRec))),T,NPPB);
        otherwise
            EoutRec = EoutMod;
    end
%%  Receiving DownStream
    
    if FFTSplit
        [EoutAux1,~,VetThisCarr]=OpticalFFTN(t,T,MaxNumStag,EoutRec);
    else
        VetThisCarr = ObsCarrPos;
        if SendingUp&&(NumCarr==2)
            EoutAux1    = SelectEachCarrier(EoutRec,NumCarr,f,fin,1.0*fc,11,fc);
        else
            EoutAux1(1,:)    = EoutRec;
            EoutAux1(2,:)    = 0;
        end
    end

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
        [BitFilt,~] = FiltroGaussiano(f,fc,0,1);                           %Creating filter for selection of the received signal
        BitFilt = fftshift(BitFilt);                                       %Shifting the filter for matching the received signal
        Ix = ifft(fft(Ix).*BitFilt);                                     %Filtering the signal for better detection
        Ix = Ix - min(Ix);                                                 %Removing any off-set that may exist
        Ix = Ix./max(abs(Ix));                                             %Normalizing the eletrical signal (amplifying what is needed)
%         RecSig = reshape(Ix,NumFra,Ns);
        SigRecep = intdump(Ix,ExtSam*NPPB);
        TxSigToPlot = zeros(1,length(SigRecep));
        TxSigToPlot(end/2) = 1;
% %         SigRecep = reshape(SigRecep,NumFra,NFFT);
%         SigRecep1 = fft(SigRecep.',NFFT);
%         SigRecep2 = SigRecep1(2:Ns+1,:);
%         SigRecep2 = SigRecep2(:).';
%         switch OfdMod
%             case 'QAM'
%                 TxSigToPlot  = qammod(TxDataMat(ThisCarr,:),M);            %Modulating information
%             otherwise
% %                 TxSig = reshape(TxDataMat(ThisCarr,:),NumFra,NFFT);
%                 TxSigToPlot  = dpskmod(TxDataMat(ThisCarr,:),M);                           %Modulating information
% %                 TxSigToPlot  = TxSigToPlot(:);
%         end
%         if exist('equ1','var')
%             SigRecep2 = ((SigRecep2)./equ1(ObsCarrUsed(ThisCarr),:));
%         end
%         if exist('equ2','var')
%             SigRecep2 = ((SigRecep2)./equ2(ObsCarrUsed(ThisCarr),:));
%         end
        equ        = (SigRecep);                               % Coeficientes do Equalizador
        
%         figure;
%         hold all;
%         plot(TxSigToPlot,'*','LineWidth',2);
%         plot(SigRecep2,'o','LineWidth',2);
%         set(gcf,'units','normalized','outerposition',[0 0 1 1]);
        
        ChanelEqualizer(ObsCarrUsed(ThisCarr),:)        = equ;

    end
    clear EoutModTem EoutModAux EoutMod;
    %% Transmiting UpStream
    for kk=InitCarrUp:CarrPass:NumCarr                                     %Generating different data for each carrier
        TxData          = zeros(1,Nb);                        %Generation random information
        TxData(end/2)   = 1;
        TxDataMat(kk,:) = TxData;                                       %Saving data for future comparison
%         TxDataPara      = TxData.';                                        %Converting from serial to paralel 
%         switch OfdMod                                                      %Sellect which modulation was used
%             case 'QAM'
%                 TxSymb  = qammod(TxDataPara,M);                            %Modulating information by QAM
%             otherwise
%                 TxSymb  = dpskmod(TxDataPara,M);                           %Modulating information DPSK as default
%         end
%         TxSymbConj      = conj(flipud(TxSymb));
%         TxSymbH = [zeros(1,NumFra);TxSymb;zeros(1,NumFra);TxSymbConj];     %Hermitian Symetri
%         TxSymbC1 = [TxSymb;TxSymb(1:DifNsN,:)];
%         TxSymbConjC      = conj(flipud(TxSymbC1));
% %         TxSymbC = [zeros(1,NumFra);TxSymbC1;zeros(1,NumFra);TxSymbConjC];
% %         TxSymbC = [TxSymbH(1:size(TxSymbH)/2,:);TxSymbH(2:DifNsN,:);0;0;TxSymbH(end-DifNsN+2:end,:);TxSymbH(1+size(TxSymbH)/2:end,:);];
%         TxSymbC = zeros(NFFT,NumFra);
% %     %                 figure;hold all;plot(TxSymbH);
%         TxSymbC(1:size(TxSymbH)/2,:) = TxSymbH(1:size(TxSymbH)/2,:);
%         TxSymbC(end-(size(TxSymbH)/2)+1:end,:) = TxSymbH(1+size(TxSymbH)/2:end,:);
% %     %                 plot(TxSymbC);
%         TxSymbMod       = ifft(TxSymbC,NFFT);
%         AuxTxSymbMod(1,:)  = TxSymbMod(:); 
        SigTx  = rectpulse(TxData,NPPB);
        SigTx(end/2) = 1;
        TxSig = SigTx;
%         TxSig = SigTx.';
%Assigning the eletrical signal to one drive of the MZM - The aplitude of 
%the signal can be controlled by the variable DatGai, which can be 
%understood as an gain for the eletrical signal or an atenuation. The 
%second signal will be similar with the only difference a phase shift of pi
        TxSig = 0.95*(TxSig./max(abs(TxSig)));
        U.U1t = TxSig;
        U.U2t = exp(-1j*pi).*TxSig;
%As both signals will have the mostrly the same characteristics with the 
%only difference the phase shift of 180 degress. The MZM-I will be working 
%on the Push-Pull configuration. It is necessary to reduce the Chirp noise
%to zero.
        [EoutModAux,~]=Mach_Zehnder_Modulator_simplificado(t,...
                   EoutTx(VetThisCarrTx==(RefCarr-1+kk),:),U,MZ_Input_File);%Modulating individual carriers
        EoutModTem(ObsCarrPos(kk),1:length(EoutModAux)) = EoutModAux;      %Adding the current Modulated output to the final OutPut
%             if kk==1
%                 figure;hold all;grid on;
%                 plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
%                 set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%             end
    end
    if IfftOrSum
        if size(EoutModTem,1)>1
            [EoutMod,~,~] = OpticalIFFTN(t,T,MaxNumStagT,EoutModTem);
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
    switch Medium
        case 'B2B'
            EoutRec = EoutMod;
        case 'Fiber'
            [EoutRec,~,~]=Fibra_Monomodo1(t,EoutMod,lambdac,T,FiberLength,0,f);
    %         PrintInfo(PlotingThis*11,EoutRec.*conj(EoutRec)./max(abs(EoutRec...
    %                                                  .*conj(EoutRec))),T,NPPB);
        otherwise
            EoutRec = EoutMod;
    end
%% Receving UpStream
    [EoutAux1,~,VetThisCarr]=OpticalFFTN(t,T,MaxNumStag,EoutRec);
    InitCarr = InitCarrUp;
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
        [BitFilt,~] = FiltroGaussiano(f,fc,0,1);                           %Creating filter for selection of the received signal
        BitFilt = fftshift(BitFilt);                                       %Shifting the filter for matching the received signal
        Ix = ifft(fft(Ix).*BitFilt);                                     %Filtering the signal for better detection
        Ix = Ix - min(Ix);                                                 %Removing any off-set that may exist
        Ix = Ix./max(abs(Ix));                                             %Normalizing the eletrical signal (amplifying what is needed)
        SigRecep = intdump(Ix,ExtSam*NPPB);
        TxSigToPlot = zeros(1,length(SigRecep));
        TxSigToPlot(end/2) = 1;
% %         SigRecep = reshape(SigRecep,NumFra,NFFT);
%         SigRecep1 = fft(SigRecep.',NFFT);
%         SigRecep2 = SigRecep1(2:Ns+1,:);
%         SigRecep2 = SigRecep2(:).';
%         switch OfdMod
%             case 'QAM'
%                 TxSigToPlot  = qammod(TxDataMat(ThisCarr,:),M);            %Modulating information
%             otherwise
% %                 TxSig = reshape(TxDataMat(ThisCarr,:),NumFra,NFFT);
%                 TxSigToPlot  = dpskmod(TxDataMat(ThisCarr,:),M);                           %Modulating information
% %                 TxSigToPlot  = TxSigToPlot(:);
%         end
%         if exist('equ1','var')
%             SigRecep2 = ((SigRecep2)./equ1(ObsCarrUsed(ThisCarr),:));
%         end
%         if exist('equ2','var')
%             SigRecep2 = ((SigRecep2)./equ2(ObsCarrUsed(ThisCarr),:));
%         end
        equ        = (SigRecep);                               % Coeficientes do Equalizador

%         figure;
%         hold all;
%         plot(TxSigToPlot,'*','LineWidth',2);
%         plot(SigRecep2,'o','LineWidth',2);
%         set(gcf,'units','normalized','outerposition',[0 0 1 1]);

        ChanelEqualizer(ObsCarrUsed(ThisCarr),:)        = equ;
    end
end

