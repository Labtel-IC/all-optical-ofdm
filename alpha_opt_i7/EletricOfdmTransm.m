%%
            TxData = zeros(length(DmtMve),NumFra,CurTesExtSiz);
            TxSymb = zeros(length(DmtMve),NumFra,CurTesExtSiz);
            TxDataMat = zeros(1,NumFra*length(DmtMve)*CurTesExtSiz);
            %U1t = 0;
            UniDmtMve = unique(DmtMve);
            UniDmtMve = fliplr(UniDmtMve);
            for CarDmtM=1:length(UniDmtMve)
                M = UniDmtMve(CarDmtM);
                DaSiAu = sum(ismember(DmtMve,M));
                TxData = randi([0 M-1],DaSiAu,NumFra,CurTesExtSiz);%Generation random information
                TxDataMat(1,1+CurTesExtSiz*DaSiAu*NumFra*(CarDmtM-1):CurTesExtSiz*DaSiAu*NumFra*CarDmtM) = TxData(:);%Saving data for future comparison
                switch OfdMod                                              %Sellect which modulation was used
                    case 'qam'
                        TxSymb(1+DaSiAu*(CarDmtM-1):DaSiAu*CarDmtM,1:NumFra,1:CurTesExtSiz)  = qammod(TxData,M);                    %Modulating information by QAM
                    otherwise
                        TxSymb(1+DaSiAu*(CarDmtM-1):DaSiAu*CarDmtM,1:NumFra,1:CurTesExtSiz)  = dpskmod(TxData,M);                   %Modulating information DPSK as default
                end
            end
%             TxSymb      = TxSymbAux.';                                %Converting from serial to paralel
%             EvmMatRef(1,:) = TxSymb(:);
            %Construction the Hermitian Simetry
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
            TxSymbConj      = conj(flipud(TxSymb));                    %Singnal conjugate format at reverse order
            if UsingHermitian
                TxSymbH = [zeros(1,NumFra,CurTesExtSiz);TxSymb;zeros(1,NumFra,CurTesExtSiz);TxSymbConj];%Hermitian Symetri
                TxSymbC = zeros(NFFT,NumFra,CurTesExtSiz);
                TxSymbC(1+Zt-ZtC:Zt-ZtC+size(TxSymbH,1)/2,:,:) = TxSymbH(1:size(TxSymbH,1)/2,:,:);
                TxSymbC(end-Zt+ZtC-(size(TxSymbH,1)/2)+1:end-Zt+ZtC,:,:) = TxSymbH(1+size(TxSymbH,1)/2:end,:,:);
            else
                TxSymbH = TxSymb;
                %     TxSymbC = [zeros((NFFT-size(TxSymbH,1))/2,NumFra);TxSymbH;zeros((NFFT-size(TxSymbH,1))/2,NumFra)];
                TxSymbC = [TxSymbH(1:size(TxSymbH,1)/2,:,:);zeros((NFFT-size(TxSymbH,1)),NumFra,CurTesExtSiz);TxSymbH(end+1-size(TxSymbH,1)/2:end,:,:)];
            end
            %Here was implemented the modulation through FFT process, it basica mix-up
            %all infomation following Fourier transform.
            TxSymbMod  = ifft(TxSymbC,NFFT,1);
            TxSymbMod  = rectpulse(TxSymbMod,OvSam);
            TxSymbMod  = reshape(TxSymbMod,NFFT*OvSam,NumFra,CurTesExtSiz);
            %Sellecting whether the OFDM signal will be transmited on
            %Base band or not
            tt = linspace(0,1*Te,size(TxSymbMod,1));
            ff = time2freq(tt);
            tta = repmat(tt.',1,NumFra,CurTesExtSiz);
            switch SelModTp
                case 'BP'                                              %Base Band transmission
%                     TxSymbModA  = rectpulse(TxSymbModA,OvSam);         %Over sampling
                    %                         TxSymbMod = TxSymbModA;
%                     TxSymbMod   = TxSymbMod;
                    if SelecGaus==1
                        [ ModFilt ] = FiltroGaussiano(ff,OBw/2,0,1);
                    else
                        [ ModFilt ] = Filtro_Retangular(OBw,0,ff);
                    end
                    ModFilt = fftshift(ModFilt);
                    if  (Ploting)
                        figure;hold all;
                        plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1,1))./length(TxSymbMod(:,1,1))))));
                    end
                    ModFilta = repmat(ModFilt.',1,NumFra,CurTesExtSiz);
                    if SelModFilt==1
                        TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
                    end
                    if  (Ploting)
                        %                             for jj=1:NumFra
                        %                                 ff2(:,jj) = ff;
                        %                             end
                        plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1,1))./length(TxSymbMod(:,1,1))))));
                        plot(ff,20*log10(abs(fftshift(ModFilt))));
                        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                    end
                case 'AM'                                              %Out of base band
%                     TxSymbModA  = rectpulse(TxSymbModA,OvSam);         %Over sampling
                    TxSymbMod   = TxSymbMod.*cos(2*pi*Ofc*tta);
                    if SelecGaus==1
                        [ ModFilt ] = FiltroGaussiano(ff,OBw,Ofc,1);
                    else
                        [ ModFilt ] = Filtro_Retangular( OBw,Ofc,ff);
                    end
                    ModFilt = fftshift(ModFilt);
                    if  (Ploting)
                        figure;hold all;
                        plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1,1))./length(TxSymbMod(:,1,1))))));
                    end
                    ModFilta = repmat(ModFilt.',1,NumFra,CurTesExtSiz);
                    if SelModFilt==1
                        TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
                    end
                    if  (Ploting)
                        plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1,1))./length(TxSymbMod(:,1,1))))));
                        plot(ff,20*log10(abs(fftshift(ModFilt))));
                        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                    end
                case 'AMSSB'                                           %Out of base band
%                     TxSymbModA  = rectpulse(TxSymbModA,OvSam);         %Over sampling
                    TxSymbMod   = TxSymbMod.*exp(1j*2*pi*Ofc*tta);
                    if SelecGaus==1
                        [ ModFilt ] = FiltroGaussiano(ff,OBw,Ofc,1);
                    else
                        [ ModFilt ] = Filtro_Retangular( OBw,Ofc,ff);
                    end
                    ModFilt = fftshift(ModFilt);
                    if  (Ploting)
                        figure;hold all;
                        plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1,1))./length(TxSymbMod(:,1,1))))));
                    end
                    ModFilta = repmat(ModFilt.',1,NumFra,CurTesExtSiz);
                    if SelModFilt==1
                        TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
                    end
                    if  (Ploting)
                        plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1,1))./length(TxSymbMod(:,1,1))))));
                        plot(ff,20*log10(abs(fftshift(ModFilt))));
                        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                        a=a+1;
                    end
                otherwise
%                     TxSymbModA  = rectpulse(TxSymbModA,OvSam);         %Over sampling
                    [TxSymbMod,tt] = modulate(TxSymbMod,Ofc,Ofs,'pm');
                    if SelecGaus==1
                        [ ModFilt ] = FiltroGaussiano(ff,3*Ofc,0,1);
                    else
                        [ ModFilt ] = Filtro_Retangular( 6*Ofc,0,ff);
                    end
                    ModFilt = fftshift(ModFilt);
                    if  (Ploting)
                        figure;hold all;
                        plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1,1))./length(TxSymbMod(:,1,1))))));
                    end
                    ModFilta = repmat(ModFilt.',1,NumFra,CurTesExtSiz);
                    TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
                    if  (Ploting)
                        plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1,1))./length(TxSymbMod(:,1,1))))));
                        plot(ff,20*log10(abs(fftshift(ModFilt))));
                        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                    end
            end
            %
            %                 figure;plot(TxSymbMod(:,1));
            
            %AuxTxSymbMod(1,:)  = TxSymbMod(:);
            U1t  = [TxSymbMod(1+end-NPOFEX:end,:,:);TxSymbMod];   %Adding Ciclyc prefix
%             U1t  = rectpulse(TxSymbModA,NPPOF);                      %Over sampling
            
            TxSigA = reshape(U1t,NumFra*size(U1t,1),CurTesExtSiz);                                       %Serializing signal
            U1t = [TxSigA(1+end-NPSTUf:end,:);TxSigA];                         %Conforming vector sizes to the same length
            % CpPaFr = CpPaFr + 1;
            %Assigning the eletrical signal to one drive of the MZM - The aplitude of
            %the signal can be controlled by the variable DatGai, which can be
            %understood as an gain for the eletrical signal or an atenuation. The
            %second signal will be similar with the only difference a phase shift of pi
            NormFact = max(max(U1t));
            U1t = (U1t./NormFact);                     %Normalizing the sinal to mach with the MZM espect to receive.
            
            % if  (PlotingThisThis)
            %     figure;hold all;grid on;
            %     plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
            %     axis([-25e9 37.5e9 -200 0]);
            %     set(gcf,'units','normalized','outerposition',[0 0 1 1]);
            % end
            
            
            if ChanAwgn==1
                [~,SigPower] = MeasPower(U1t);
                SigPower2 = SigPower-30;%20*log10(SigPowerI);
                %     SNR = CarSNR + 10*log10(log2(max(DmtMve))) - 10*log10(1*Ts*OBw);
                SNR = CarSNR + 10*log10(log2(max(DmtMve))) - 10*log10(OvSam);
                %     SNR
                U1t = awgn(U1t,SNR,SigPower2);
                %     EoutAux = TxSig;
            else
%                 U1t = U1t;
            end