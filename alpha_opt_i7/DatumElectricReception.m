%%

Ix = EoutAux.*NormFact;

if Ploting
    figure;hold all;grid on;plot(f,20*log10(abs(fftshift(fft(Ix)./length(Ix)))));axis([-25e9 37.5e9 -200 0]);set(gcf,'units','normalized','outerposition',[0 0 1 1]);
end 
Ix = Ix(1+NPSTUf:end);
Ix = reshape(Ix,length(Ix)/NumFra,NumFra);
Ix = intdump(Ix,NPPOF);                                        %Downsampling the income signal
SigRecepA = Ix(1+NPOFEX:end,:);

switch SelModTp
    case 'BP'
        %                     SigRecep  = intdump(SigRecepA,OvSam);
        [BitFiltEle,~] = FiltroGaussiano(ff,OBw/2,0,1);        %Creating filter for selection of the received signal
        %                     [ BitFiltEle ] = Filtro_Retangular(OBw,0,ff);
        BitFiltEle = fftshift(BitFiltEle);
        for jj=1:NumFra
            ModFilta(:,jj) = BitFiltEle;
        end
        SigRecepB   = SigRecepA;
        if Ploting
            figure;
            hold all;
            plot(ff,20*log10(abs(fftshift(fft(SigRecepB)./length(SigRecepB)))));
            plot(ff,20*log10(abs(fftshift(BitFiltEle))));
        end
        if SelModFilt
            SigRecepB = ifft(fft(SigRecepB).*ModFilta);
        end
        SigRecep  = intdump(SigRecepB,OvSam);
        if Ploting
            plot(ff,20*log10(abs(fftshift(fft(SigRecepB)./length(SigRecepB)))));
            set(gcf,'units','normalized','outerposition',[0 0 1 1]);
        end
    case 'AM'
        [BitFiltEle,~] = FiltroGaussiano(ff,OBw/2,0,1);           %Creating filter for selection of the received signal
        BitFiltEle = fftshift(BitFiltEle);
        for jj=1:NumFra
            ModFilta(:,jj) = BitFiltEle;
        end
        SigRecepB   = SigRecepA.*cos(-2*pi*Ofc*tta);
        if Ploting
            figure;
            hold all;
            plot(ff,20*log10(abs(fftshift(fft(SigRecepB)./length(SigRecepB)))));
            plot(ff,20*log10(abs(fftshift(BitFiltEle))));
        end
        if SelModFilt
            SigRecepB = ifft(fft(SigRecepB).*ModFilta);
        end
        SigRecep  = intdump(SigRecepB,OvSam);
        if Ploting
            plot(ff,20*log10(abs(fftshift(fft(SigRecepB)./length(SigRecepB)))));
            set(gcf,'units','normalized','outerposition',[0 0 1 1]);
        end
        a=1;
    case 'AMSSB'
        [BitFiltEle,~] = FiltroGaussiano(ff,OBw/2,0,1);           %Creating filter for selection of the received signal
        BitFiltEle = fftshift(BitFiltEle);
        for jj=1:NumFra
            ModFilta(:,jj) = BitFiltEle;
        end
        SigRecepB   = SigRecepA.*exp(-1j*2*pi*Ofc*tta);
        if Ploting
            figure;
            hold all;
            plot(ff,20*log10(abs(fftshift(fft(SigRecepB)./length(SigRecepB)))));
            plot(ff,20*log10(abs(fftshift(BitFiltEle))));
        end
        if SelModFilt
            SigRecepB = ifft(fft(SigRecepB).*ModFilta);
        end
        SigRecep  = intdump(SigRecepB,OvSam);              %Over sampling
        if Ploting
            plot(ff,20*log10(abs(fftshift(fft(SigRecepB)./length(SigRecepB)))));
            set(gcf,'units','normalized','outerposition',[0 0 1 1]);
        end
        a=1;
    otherwise
        SigRecepB = demod(SigRecepA,Ofc,Ofs,'pm');
        SigRecep  = intdump(SigRecepB,OvSam);              %Over sampling
end
%             SigRecep = reshape(SigRecep,NumFra,NFFT);                      %Reshaping the signal when multiple frames were transmited
SigRecep1 = fft(SigRecep,NFFT);                              %FFT demultiplexing process
if UsingHermitian
    SigRecep2 = SigRecep1(2+Zt-ZtC:Ns+1+Zt-ZtC,:);                               %Taking the actual data
else
%     SigRecep2 = SigRecep1(1+(NFFT-size(TxSymbH,1))/2:size(TxSymbH,1)+(NFFT-size(TxSymbH,1))/2,:);
    SigRecep2 = [SigRecep1(1:size(TxSymbH,1)/2,:);SigRecep1(end+1-size(TxSymbH,1)/2:end,:)];
end
SigRecep3 = SigRecep2;
%             SigRecep3 = SigRecep3(:).';                                    %Changing from parralel to serial
%             SigRecep3 = SigRecep2.';
if Ploting
    switch OfdMod                                                  %Modulation Tx signal for future comparison
        case 'qam'
            TxDataA = TxDataMat(1,:);
            TxDataA = reshape(TxDataA,NumFra,length(TxDataA)/NumFra);
            TxDataA = TxDataA.';
            for CarDmtM=1:length(DmtMve)
                M = DmtMve(CarDmtM);
                TxSigToPlot(CarDmtM,:)  = qammod(TxDataA(CarDmtM,:),M);            %Modulating information
            end
        otherwise
            TxDataA = TxDataMat(1,:);
            TxDataA = reshape(TxDataA,NumFra,length(TxDataA)/NumFra);
            TxDataA = TxDataA.';
            for CarDmtM=1:length(DmtMve)
                M = DmtMve(CarDmtM);
                TxSigToPlot(CarDmtM,:)  = dpskmod(TxDataA(CarDmtM,:),M);            %Modulating information
            end
            %                     SigRecepA = reshape(TxDataMat(ThisCarr,:),NumFra,NFFT);
            %                     TxSigToPlot = dpskmod(TxDataA,M);                  %Modulating information
            %                     TxSigToPlot = TxSigToPlot(:);
    end
    %         TxSigToPlot = dpskmod(TxDataMat(ThisCarr,2:end),M);
    if Ploting
        figure;
        txcolor = [0.2 0 1];
        rxcolor = [1 0.4 0];
        hold all;
        plot(TxSigToPlot(:),'*','LineWidth',2,'color',txcolor);
        plot(SigRecep3(:),'o','LineWidth',2,'color',rxcolor);
        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    end
end
SaveRxNotEq = [SaveRxNotEq SigRecep3(:)];

switch SelModTp
    case 'BP'
        switch OfdMod
            case 'qam'
                RxSigOfdmNoEq(CurrentTest,1,:) = SigRecep3(:).';
                clear SigRecep4;
                for CarDmtM=1:length(DmtMve)
                    M = DmtMve(CarDmtM);
                    SigRecep4(CarDmtM,:)  = qamdemod(SigRecep3(CarDmtM,:),M);            %Modulating information
                end
            otherwise
                clear SigRecep4;
                for CarDmtM=1:length(DmtMve)
                    M = DmtMve(CarDmtM);
                    SigRecep4(CarDmtM,:)  = dpskdemod(SigRecep3(CarDmtM,:),M);            %Modulating information
                end
        end
    case 'AMSSB'
        switch OfdMod
            case 'qam'
                RxSigOfdmNoEq(CurrentTest,1,:) = SigRecep3(:).';
                equa = ChanelEqualizer(1,:);
                equ = reshape(equa,length(equa)/NumFra,NumFra);
                if ~ChanAwgn
                    SigRecep3 = ((SigRecep3)./equ);%Fixing phase rotation and amplitude variation with the equalizer
                end
                clear SigRecep4;
                for CarDmtM=1:length(DmtMve)
                    M = DmtMve(CarDmtM);
                    SigRecep4(CarDmtM,:)  = qamdemod(SigRecep3(CarDmtM,:),M);            %Modulating information
                end
            otherwise
                clear SigRecep4;
                for CarDmtM=1:length(DmtMve)
                    M = DmtMve(CarDmtM);
                    SigRecep4(CarDmtM,:)  = dpskdemod(SigRecep3(CarDmtM,:),M);            %Modulating information
                end
        end
    otherwise
        switch OfdMod
            case 'qam'
                RxSigOfdmNoEq(CurrentTest,1,:) = SigRecep3(:).';
                equa = ChanelEqualizer(1,:);
                equ = reshape(equa,length(equa)/NumFra,NumFra);
                SigRecep3 = ((SigRecep3)./equ);%Fixing phase rotation and amplitude variation with the equalizer
                clear SigRecep4;
                for CarDmtM=1:length(DmtMve)
                    M = DmtMve(CarDmtM);
                    SigRecep4(CarDmtM,:)  = qamdemod(SigRecep3(CarDmtM,:),M);            %Modulating information
                end
            otherwise
                clear SigRecep4;
                for CarDmtM=1:length(DmtMve)
                    M = DmtMve(CarDmtM);
                    SigRecep4(CarDmtM,:)  = dpskdemod(SigRecep3(CarDmtM,:),M);            %Modulating information
                end
        end
end
EvmMatRec(1,:) = SigRecep3(:);%Taking just the middle samples as references
[EvmDB(CurrentTest,1), EvmPer(CurrentTest,1), EvmRms(CurrentTest,1) ] = ...
                                  EvmCalc( EvmMatRef(1,:),SigRecep3(:).' );
[EvmDBJ(CurrentTest,1),EvmPerJ(CurrentTest,1),EvmRmsJ(CurrentTest,1)] = ...
                              evm1(M,OfdMod,EvmMatRef(1,:),SigRecep3(:).');
%         SigRecep4 = dpskdemod(SigRecep3,M);
SigRecep4 = SigRecep4.';
SigRecep4 = SigRecep4(:).';
RxSigOfdm(CurrentTest,1,:) = SigRecep4;
SaveRxEq = [SaveRxEq SigRecep4];
if Ploting
    figure(111);
    txcolor = [0.2 0 1];
    rxcolor = [1 0.4 0];
    hold all;
    plot(TxSigToPlot(:),'*','LineWidth',2,'color',txcolor);
    plot(SigRecep3(:),'o','LineWidth',2,'color',rxcolor);
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    figure;
    hold all;
    plot(TxDataMat(1,:));
    plot(SigRecep4);
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
end
SymbErr = ~(TxDataMat(1,:)==SigRecep4);%Measuring system bit error ration
DmtMvep = repmat(DmtMve.',1,NumFra);
DmtKvep = log2(DmtMvep(:).');
%             BerOFDM(CurrentTest,ThisCarr) = sum(SymbErr)/length(SymbErr);
BerOFDM(CurrentTest,1) = sum(SymbErr.*DmtKvep)/sum(DmtKvep);
% BerOFDM
a=0;
% close all;