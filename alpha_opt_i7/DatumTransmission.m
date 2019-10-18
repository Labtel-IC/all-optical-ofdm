%c
%c                                                       ..'Ž`'..'Ž`..'Ž`..
%c       File: DatumTransmission File For the All Optical OFDM
%c(Main file to call each idividual parameters)
%c
%c     This main code is resposible to call and run the all the fuction to
%c related to this simulation, which is responsible to generate a random
%c strem of data, modulated it and transmite it accross a SSMF.
%c
%c
%c
%c                                           by P.Marciano LG
%c                                           29/10/2017
%c                                           Last UpDate
%c                                           09/05/2018
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
% The first step here is to configure all the parameters that will used
% hereafter by the transmiter. Hence, implementing functionalities such as
% generation of random data to fill the pay load, create the syncronism
% symble, modulating the incoming data acordingly, transmiting this
% information throughout a medium (fiber or B2B). All those functionalities
% will need variables to be first initialyzed in this following Input Data.
DatumTransmissionInputData;
%%               Creating Data for transmission
% Creating the information accordingly with the type of transmition

%It is important to mention that an stream of data will be formed for each
%individual carrier. For the last but not the least, the DownStream and the
%UpStream carriers will be interleaved within the same transmition pass
%band.For example, if the user choses carriers 1,3,5 and 7 as DownStream
%carriers, the UpStream will be formed by carriers 2,4,6, and 8. This
%manuver was suggested by Segatto to address the problem of ICI (Inter
%Carrier Interference). The transmited signal still an OFDM once the
%carriers from 1 to 8, for instance, still orthogonal to each other.
%Therefore, the for loop to generate the information to be transmited will
%have the passe of 2.

if SendingDowStr
    
%%  Selecting part of the UOFCS
%The generation fo the Ultra Optical Flat Comb Source is not perfect, which
%may drive the user to select a specific part of it for its actual use.
%As a result, this condition verify if the user wants to select a part of
%the UOFCS or use it as an whole.
%There are two ways, implemented, to split the incoming OFCS. One by
%filters e another by the optical FFT. But this step just need to be done
%once in the simulation.

%After a recent updated where the OIFFT was added to summup the modulated
%signal, it is important to keep track of the right order of the carriers.
%Because the OIFFT is very sensitive to the order of the income signals. As
%it has a periodic response, if the correct signal was not placed at the
%right port the response will destroy the OFDM signal.
    if ~exist('EoutTx','var')                                              %Verify is this step was previouly done
        if SelecSetUP                                                      %Vefify is the split process will be by filter...
            EoutTx=SelectEachCarrier(Eout,NumCarr,f,fin,FBWD,Order,fc);    %This function is resposible to split each carrier
            VetThisCarrTx = (RefCarr-1)+1:(RefCarr-1)+NumCarr;             %Keeping track of the right carrier sequence
        else                                                               %... or by the OFFT
%As the OFFT has a periodic response the OFC needs to be constrained other 
%whise carrier multiple carrier may interfir with other channels, This
%first selection was done with a broad pass-band filter. A higher OFFT
%order can also be used, although it may increase the computational time
%and my not be exactly feasible in the real world.
            if Selecting                                                    
                EoutT = ifft(fft(Eout).*SelecFilt);
%                 PrintInfo(Ploting*1,f,EoutT);
%                 a=0;
            else
                EoutT = Eout;
            end                                                                 
            [EoutTx,~,VetThisCarrTx]=OpticalFFTN(f,T,MaxNumStagTx,EoutT);  %This function is resposible to split each carrier
        end
%         PrintInfo(Ploting*51,EoutTx,f);                              %Printing for qualitative analizes.
%         axis([min(f) max(f) -400 0]);
%         a=0;
    end
    switch Modulation
        %The reason to implement the OFDM system was to work with ROF
        %architecture. In that area the signal will not be converted to
        %data before being transmited through optical fiber. The signal
        %must be the modulator of a optical carrier. At LTE system, the
        %signal received is OFDM electrical and it will also be the signal
        %format for the 5G technology. Therefore, it was important to study
        %how a OFDM signal would behaviour in the architecture here
        %proposed.
        case 'OFDM' 
            for kk=InitCarrDo:CarrPass:NumCarr                             %Generating different data for each carrier
                TxSigB = 0;
%                 CpPaFr = 1;
                %for CarrOffSet=-1*(NumFraPar-1)/2:1:(NumFraPar-1)/2
                    for CarDmtM=1:length(DmtMve)
                        M = DmtMve(CarDmtM);
                        TxData          = randi([0 M-1],NumFra,1);                %Generation random information
                        TxDataMat(kk,1+NumFra*(CarDmtM-1):NumFra*CarDmtM) = TxData(:);                               %Saving data for future comparison
                        switch OfdMod                                              %Sellect which modulation was used
                            case 'qam'
                                TxSymbAux(1:NumFra,CarDmtM)  = qammod(TxData,M);                    %Modulating information by QAM
                            otherwise
                                TxSymbAux(1:NumFra,CarDmtM)  = dpskmod(TxData,M);                   %Modulating information DPSK as default
                        end
                    end
                    TxSymb      = TxSymbAux.';                                %Converting from serial to paralel 
                    TxSigOfdm(CurrentTest,ObsCarrPos==kk,:) = TxSymb(:);
                    EvmMatRef(ObsCarrPos==kk,:) = TxSymb(:);                   %Taking just the middle samples as references
                    %Construction the Hermitian Simetry
                    TxSymbConj      = conj(flipud(TxSymb));                    %Singnal conjugate format at reverse order
                    
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
                    if UsingHermitian
                        TxSymbH = [zeros(1,NumFra);TxSymb;zeros(1,NumFra);TxSymbConj];%Hermitian Symetri
                        TxSymbC = zeros(NFFT,NumFra);
                        TxSymbC(1+Zt-ZtC:Zt-ZtC+size(TxSymbH,1)/2,:) = TxSymbH(1:size(TxSymbH,1)/2,:);
                        TxSymbC(end-Zt+ZtC-(size(TxSymbH,1)/2)+1:end-Zt+ZtC,:) = TxSymbH(1+size(TxSymbH,1)/2:end,:);
                    else
                        TxSymbH = TxSymb;
                    %     TxSymbC = [zeros((NFFT-size(TxSymbH,1))/2,NumFra);TxSymbH;zeros((NFFT-size(TxSymbH,1))/2,NumFra)];
                        TxSymbC = [TxSymbH(1:size(TxSymbH,1)/2,:);zeros((NFFT-size(TxSymbH,1)),NumFra);TxSymbH(end+1-size(TxSymbH,1)/2:end,:)];
                    end
    %Here was implemented the modulation through FFT process, it basica mix-up
    %all infomation following Fourier transform.
                    TxSymbModA       = ifft(TxSymbC,NFFT);
                    %Sellecting whether the OFDM signal will be transmited on
                    %Base band or not
                    switch SelModTp
                        case 'BP'                                              %Base Band transmission
                            TxSymbModA  = rectpulse(TxSymbModA,OvSam);         %Over sampling
    %                         TxSymbMod = TxSymbModA;
                            tt = linspace(0,1*Te,size(TxSymbModA,1));
                            for jj=1:NumFra
                                tta(:,jj) = tt;
                            end
                            TxSymbMod   = TxSymbModA;%.*exp(1j*2*pi*CarrOffSet*FramSpac*tta);
                            ff = time2freq(tt);
                            if SelecGaus
                                [ ModFilt ] = FiltroGaussiano(ff,OBw/2,0,1);%CarrOffSet*FramSpac,1);
                            else
                                [ ModFilt ] = Filtro_Retangular(OBw,0,ff);%CarrOffSet*FramSpac,ff);
                            end
                            ModFilt = fftshift(ModFilt);
                            if (kk==1)&&(Ploting)
                                figure;hold all;
                                plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                            end
                            for jj=1:NumFra
                                ModFilta(:,jj) = ModFilt;
                            end
                            if SelModFilt
                                TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
                            end
                            if (kk==1)&&(Ploting)
    %                             for jj=1:NumFra
    %                                 ff2(:,jj) = ff;
    %                             end
                                plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                                plot(ff,20*log10(abs(fftshift(ModFilt))));
                                set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                            end
                        case 'AM'                                              %Out of base band
                            TxSymbModA  = rectpulse(TxSymbModA,OvSam);         %Over sampling
                            tt = linspace(0,1*Te,size(TxSymbModA,1));
                            ff = time2freq(tt);
                            for jj=1:NumFra
                                tta(:,jj) = tt;
                            end
                            TxSymbMod   = TxSymbModA.*cos(2*pi*Ofc*tta);
                            if SelecGaus
                                [ ModFilt ] = FiltroGaussiano(ff,OBw,Ofc,1);
                            else
                                [ ModFilt ] = Filtro_Retangular( OBw,Ofc,ff);
                            end
                            ModFilt = fftshift(ModFilt);
                            if (kk==1)&&(Ploting)
                                figure;hold all;
                                plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                            end
                            for jj=1:NumFra
                                ModFilta(:,jj) = ModFilt;
                            end
                            if SelModFilt
                                TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
                            end
                            if (kk==1)&&(Ploting)
                                plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                                plot(ff,20*log10(abs(fftshift(ModFilt))));
                                set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                            end
                        case 'AMSSB'                                           %Out of base band
                            TxSymbModA  = rectpulse(TxSymbModA,OvSam);         %Over sampling
                            tt = linspace(0,1*Te,size(TxSymbModA,1));
                            ff = time2freq(tt);
                            for jj=1:NumFra
                                tta(:,jj) = tt;
                            end
                            TxSymbMod   = TxSymbModA.*exp(1j*2*pi*Ofc*tta);
                            if SelecGaus
                                [ ModFilt ] = FiltroGaussiano(ff,OBw,Ofc,1);
                            else
                                [ ModFilt ] = Filtro_Retangular( OBw,Ofc,ff);
                            end
                            ModFilt = fftshift(ModFilt);
                            if (kk==1)&&(Ploting)
                                figure;hold all;
                                plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                            end
                            for jj=1:NumFra
                                ModFilta(:,jj) = ModFilt;
                            end
                            if SelModFilt
                                TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
                            end
                            if (kk==1)&&(Ploting)
                                plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                                plot(ff,20*log10(abs(fftshift(ModFilt))));
                                set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                                a=1;
                            end
                        otherwise
                            TxSymbModA  = rectpulse(TxSymbModA,OvSam);         %Over sampling
                            [TxSymbMod,tt] = modulate(TxSymbModA,Ofc,Ofs,'pm');
                            ff = time2freq(tt);
                            if SelecGaus
                                [ ModFilt ] = FiltroGaussiano(ff,3*Ofc,0,1);
                            else
                                [ ModFilt ] = Filtro_Retangular( 6*Ofc,0,ff);
                            end
                            ModFilt = fftshift(ModFilt);
                            if (kk==1)&&(Ploting)
                                figure;hold all;
                                plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                            end
                            for jj=1:NumFra
                                ModFilta(:,jj) = ModFilt;
                            end
                            TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
                            if (kk==1)&&(Ploting)
                                plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                                plot(ff,20*log10(abs(fftshift(ModFilt))));
                                set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                            end
                    end
    %                 
    %                 figure;plot(TxSymbMod(:,1));

                    AuxTxSymbMod(1,:)  = TxSymbMod(:);
                    TxSymbModA  = [TxSymbMod(1+end-NPOFEX:end,:);TxSymbMod];   %Adding Ciclyc prefix
                    SigTx  = rectpulse(TxSymbModA,NPPOF);                      %Over sampling

                    TxSigA = SigTx(:).';                                       %Serializing signal
                    TxSigB = TxSigB + [TxSigA(1:NPSTUf) TxSigA];                         %Conforming vector sizes to the same length
                    TxSig = real(TxSigB).*exp(1j.*(2*pi*FeqC.*t)) + imag(TxSigB).*exp(1j.*(2*pi*FeqC.*t + pi/2));
                    %CpPaFr = CpPaFr + 1;
                %end
%Assigning the eletrical signal to one drive of the MZM - The aplitude of 
%the signal can be controlled by the variable DatGai, which can be 
%understood as an gain for the eletrical signal or an atenuation. The 
%second signal will be similar with the only difference a phase shift of pi
                NormFact = max((TxSig));
                TxSig = 0.95*(TxSig./NormFact);                     %Normalizing the sinal to mach with the MZM espect to receive.
                U.U1t = TxSig;
                U.U2t = exp(-1j*pi).*TxSig;
%As both signals will have the mostrly the same characteristics with the 
%only difference the phase shift of 180 degress. The MZM-I will be working 
%on the Push-Pull configuration. It is necessary to reduce the Chirp noise
%to zero.
                [EoutModAux,~]=Mach_Zehnder_Modulator_simplificado_2(freqGHz,...
                  EoutTx(VetThisCarrTx==(RefCarr-1+kk),:),U,MZ_Input_File);%Modulating individual carriers
                EoutModTem(ObsCarrPos(kk),1:length(EoutModAux)) = ...
                                                                EoutModAux;%Adding the current Modulated output to the final OutPut
                if (kk==1)&&(Ploting)
                    figure;hold all;grid on;
                    plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
                    axis([-25e9 37.5e9 -200 0]);
                    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                end
%This section adds the unmodulated upstream carrier on the optical OFDM
%signal
                if (~mod(CarrPass,2))&&(SendingUp)
                    clear EoutModAux;
                    EoutModAux = EoutTx(VetThisCarrTx==(RefCarr-1+kk...
                                                                +1),:);
                if AddingNoiseO                                            %Adding optical noise, if needed
                    EsN0 = 10*log10(snr);
                    [~,sigpower] = MeasPower(EoutModAux);
                    EoutModAux = awgn(EoutModAux,osnr,sigpower-30);
                    [~,sigpower] = MeasPower(EoutModAux);
                end
                    EoutModTem(ObsCarrPos(kk+1),1:length(EoutModAux)...
                                                            ) = EoutModAux;%Adding the current Unmodulated output to the final OutPut
                end
            end
%At this section individual and parralel carriers will be assemble to
%compose the optical OFDM signal. It can be done in two way, with a simple
%adder, which gives a small ICI clearance. Or by the optical IFFT, which
%gives a better ICI clearance.
            if size(EoutModTem,1)<NumCarr
                EoutModTem(end+1,:) = 0;
            end
            if IfftOrSum                                                   %Use the OIFFT
                if size(EoutModTem,1)>1                                    %Just if the signal has more than 1 dimension
                    [EoutMod,~,~] = OpticalIFFTN(f,T,MaxNumStagT,EoutModTem);
                else
                    EoutMod = EoutModTem;                                  %otherwise just let the signal pass
                end
            else                                                           %Use a simple adder which don't display a better result
                if size(EoutModTem,1)>1
                    EoutMod = sum(EoutModTem);
                else
                    EoutMod = EoutModTem;
                end
            end
        case 'DPSK'
            %%            Generate the data DPSK
            for kk=InitCarrDo:CarrPass:NumCarr                             %Generating different data for each carrier
                if CarrUsedDo(kk)
                    if SendingData                                         % Frist it is chosen to transmite a random information...
                        TxData = (randi(2,1,NbDPSK)-1);                    %Creating the data stream to be transmited
                        TxData(1:JusPos) = JusVal;                         %Adding the Justification simble at the begining of the stream to sincronize received data frame
                        TxData(end-(JusPos-1):end) = JusValEnd;            %Adding the Justification simble at the end of the stream to sincronize received data frame
                    else                                                   %... or just one high pulse for testing the channel
                        TxData        = zeros(1,NbDPSK);                   %Creating the base fo the data stream
                        TxData(end/2) = 1;                                 %Adding the code in the data stream to produce the highest level possible
                    end
                    TxDataMat(kk,:)   = TxData;                            %Storring the transmited information for latter evaluation
%Once that we have the data to be transmited it needs to be converted to 
%the eletrical signal that hereafter will modulate our carrier through the 
%MZM-I
                    SigTx           = DpskEncodEq( TxData);
                    SigTx(SigTx==1) =  1.9;
                    SigTx(SigTx==0) = -1.9;
                    SigTx           = rectpulse(SigTx,NPPB);
                    % Adding CP to the data
                    if AddCP
                        TxAux = reshape(SigTx,NPPB,NbDPSK);
                        if SetCpSampZer
                            TxAux = [zeros(NumAmosCP,size(TxAux,2));...
                                   TxAux;zeros(NumAmosCP,size(TxAux,2))];
                        else
                            TxAux = [TxAux(1:NumAmosCP,:);TxAux;TxAux(...
                                                 end-(NumAmosCP-1):end,:)];
                        end
                        TxAux = reshape(TxAux,1,(2*NumAmosCP+NPPB)*NbDPSK);
                        SigTx = [TxAux TxAux(end-(StuffSampels-1):end)];
                    end
                
%Because of reasons, the real world does not have components capable 
%instantly change the voltage level. Thus, to emulate this behaviour the 
%eletrical PAM signal pass through a gaussian filter that will remove the 
%frequencies of higher order, which will result in small slope in the level 
%variation
                
                    [BitFilt,~] = FiltroGaussiano(f,BWD,CenFeq,FiltOrd);   %Creating filter for conformation of the input information
                    BitFilt     = fftshift(BitFilt);                       %Doing a shift on the vector for matching the transmited data
                    %TxSig      = ifft(fft(SigTx).*BitFilt);               %Conforming the information and Creating the modulation signal
                    TxSig       = SigTx;                                   %Adding a gain to the eletrical signal

%                     if (kk==1)&&(Ploting)%Conforming the information and Creating the modulation signal
%                     figure;hold all;grid on;
%                     plot(t,TxSig)
%                     end
                    if AddingNoiseE
                        EsN0 = 10*log10(snr);
                        [~,sigpower] = MeasPower(TxSig);
                        TxSig = awgn(TxSig,snr,sigpower-30);
                        [~,sigpower] = MeasPower(TxSig);
                    end
%                     if (kk==1)&&(Ploting)
%                       plot(t,TxSig)
%                       PrintInfo(Ploting*2,TxSig(kk,1:end-...
%                       StuffSampels),t(2*NumAmosCP+NPPB),2*NumAmosCP+...
%                                                                    NPPB);
%                       a=4;
%                     end

%Assigning the eletrical signal to one drive of the MZM - The aplitude of 
%the signal can be controlled by the variable DatGai, which can be 
%understood as an gain for the eletrical signal or an atenuation. The 
%second signal will be similar with the only difference a phase shift of pi
                    U.U1t = TxSig;
                    U.U2t = exp(-1j*pi).*TxSig;
%As both signals will have the mostrly the same characteristics with the 
%only difference the phase shift of 180 degress. The MZM-I will be working 
%on the Push-Pull configuration. It is necessary to reduce the Chirp noise
%to zero.
                    [EoutModAux,~]=Mach_Zehnder_Modulator_simplificado_2(freqGHz...
                    ,EoutTx(VetThisCarrTx==(RefCarr-1+kk),:),U,...
                                                            MZ_Input_File);%Modulating individual carriers
                                                        
%                     if (kk==1)&&(Ploting)
%                         figure;hold all;grid on;
%                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
%                     end
                    if AddingNoiseO
                        EsN0 = 10*log10(snr);
                        [~,sigpower] = MeasPower(EoutModAux);
                        EoutModAux = awgn(EoutModAux,osnr,sigpower-30);
                        [~,sigpower] = MeasPower(EoutModAux);
                    end
%                     if (kk==1)&&(Ploting)
%                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
%                         PrintInfo(Ploting*4,abs(EoutModAux(kk,1:end-...
%                         StuffSampels)).^2./max(abs(EoutModAux(kk,1:end-...
%                         StuffSampels)).^2),t(2*NumAmosCP+NPPB),2*...
%                         NumAmosCP+NPPB,abs(EoutModAux(1:end-StuffSampels...
%                         )).^2./max(abs(EoutModAux(kk,1:end-StuffSampels)...
%                                                                     ).^2));
%                         a=4;
%                     end
                    EoutModTem(ObsCarrPos(kk),1:length(EoutModAux)) = ...
                                                                EoutModAux;%Adding the current Modulated output to the final OutPut
                    %% Taking the sampling the EVM meassurement
                    PhaDel          = 0;
                    TimDel          = T;
                    [ESync1,ESync2] = DelayInterf(f,TimDel,PhaDel,EoutModAux); %Coverting phase shift to amplitude variation
                    ESync1 = ESync1.*conj(ESync1);
                    ESync2 = ESync2.*conj(ESync2);
                    Esync  = ESync2-ESync1;
                    PosAuxEout      = (2*NumAmosCP+NPPB)/2:(2*NumAmosCP+...
                                   NPPB):(length(Esync)-StuffSampels);%Varriable respossible to take just the samples at the middle of the symbol
                    Esync           = ifft(fft(Esync).*EvmFilt);           %Removing higher signal generated by the receiving process
%                     Esync  = Esync-min(Esync);
                    Esync  = Esync-mean(Esync);
                    Esync  = Esync./max(abs(Esync));
                    IxAux  = Esync(PosAuxEout);  
                    
                    TxSymbAmos(CurrentTest,ObsCarrPos==kk,:) = IxAux;%RxSymbAmosDo = [];
                    a=0;
                    EvmMatRef(ObsCarrPos==kk,:) = IxAux;                   %Taking just the middle samples as references
%                     if (kk==1)&&(Ploting)
%                         figure;hold all;grid on;
%                         plot(IxAux,zeros(1,length(IxAux)),'rx');
%                         set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%                         a=0;
%                     end
                    %%
                    if (~mod(CarrPass,2))&&(SendingUp)
                        clear EoutModAux;
                        EoutModAux = EoutTx(VetThisCarrTx==(RefCarr-1+kk...
                                                                    +1),:);
                                                                
%                     if (kk==1)&&(Ploting)
%                         figure;hold all;grid on;
%                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
%                     end
                    if AddingNoiseO
                        EsN0 = 10*log10(snr);
                        [~,sigpower] = MeasPower(EoutModAux);
                        EoutModAux = awgn(EoutModAux,osnr,sigpower-30);
                        [~,sigpower] = MeasPower(EoutModAux);
                    end
%                     if (kk==1)&&(Ploting)
% %                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
%                         PrintInfo(Ploting*4,abs(EoutModAux(kk,1:end-...
%                         StuffSampels)).^2./max(abs(EoutModAux(kk,1:end-...
%                         StuffSampels)).^2),t(2*NumAmosCP+NPPB),2*...
%                         NumAmosCP+NPPB,abs(EoutModAux(1:end-StuffSampels...
%                         )).^2./max(abs(EoutModAux(kk,1:end-StuffSampels)...
%                                                                     ).^2));
%                         a=4;
%                     end
                        EoutModTem(ObsCarrPos(kk+1),1:length(EoutModAux)...
                                                            ) = EoutModAux;%Adding the current Modulated output to the final OutPut
                    end
                else
                    EoutModAux = EoutTx(VetThisCarrTx==(RefCarr-1+kk),:);
                    %EoutMod = EoutMod + EoutModAux + EoutTx(kk+1,:);      %Adding the current Modulated output to the final OutPut
%                     if (kk==1)&&(Ploting)
%                         figure;hold all;grid on;
%                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
%                     end
                    if AddingNoiseO
                        EsN0 = 10*log10(snr);
                        [~,sigpower] = MeasPower(EoutModAux);
                        EoutModAux = awgn(EoutModAux,osnr,sigpower-30);
                        [~,sigpower] = MeasPower(EoutModAux);
                    end
%                     if (kk==1)&&(Ploting)
%                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
%                         PrintInfo(Ploting*4,abs(EoutModAux(kk,1:end-...
%                         StuffSampels)).^2./max(abs(EoutModAux(kk,1:end-...
%                         StuffSampels)).^2),t(2*NumAmosCP+NPPB),2*...
%                         NumAmosCP+NPPB,abs(EoutModAux(1:end-StuffSampels...
%                         )).^2./max(abs(EoutModAux(kk,1:end-StuffSampels)...
%                                                                     ).^2));
%                         a=4;
%                     end
                    EoutModTem(ObsCarrPos(kk),1:length(EoutModAux)) = ...
                                                                EoutModAux;%Adding the current Modulated output to the final OutPut
                    if ~mod(CarrPass,2)
                        clear EoutModAux;
                        EoutModAux = EoutTx(VetThisCarrTx==(RefCarr-1+kk...
                                                                    +1),:);
%                     if (kk==1)&&(Ploting)
%                         figure;hold all;grid on;
%                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
%                     end
                    if AddingNoiseO
                        EsN0 = 10*log10(snr);
                        [~,sigpower] = MeasPower(EoutModAux);
                        EoutModAux = awgn(EoutModAux,osnr,sigpower-30);
                        [~,sigpower] = MeasPower(EoutModAux);
                    end
%                     if (kk==1)&&(Ploting)
%                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
%                         PrintInfo(Ploting*4,abs(EoutModAux(kk,1:end-...
%                         StuffSampels)).^2./max(abs(EoutModAux(kk,1:end-...
%                         StuffSampels)).^2),t(2*NumAmosCP+NPPB),2*...
%                         NumAmosCP+NPPB,abs(EoutModAux(1:end-StuffSampels...
%                         )).^2./max(abs(EoutModAux(kk,1:end-StuffSampels)...
%                                                                     ).^2));
%                         a=4;
%                     end
                        EoutModTem(ObsCarrPos(kk+1),1:length(EoutModAux)...
                                                            ) = EoutModAux;%Adding the current Modulated output to the final OutPut
                    end
                end
            end
            if size(EoutModTem,1)<NumCarr
                EoutModTem(end+1,:) = 0;
            end
            if IfftOrSum
%                 Ekaux=sum(EoutModTem);
%                 figure;hold all;
%                 plot(f,20*log10(abs(fftshift(fft(Ekaux)./length(Ekaux)))));
%                 set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                if size(EoutModTem,1)>1
                    [EoutMod,~,~] = OpticalIFFTN(f,T,MaxNumStagT,EoutModTem);
                else
                    EoutMod = EoutModTem;
                end
%                 plot(f,20*log10(abs(fftshift(fft(EoutMod)./length(...
%                                                               EoutMod)))));
%                 axis([0.9*fin 1.2*NumCarr*fc+fin -140 -20]);
%                 a=1;
            else
                if size(EoutModTem,1)>1
                    EoutMod = sum(EoutModTem);
                else
                    EoutMod = EoutModTem;
                end
            end
        case 'DQPSK'
            %%            Generate the data DQPSK
            EoutMod    = 0;
            EoutModTem = 0;
            for kk=InitCarrDo:CarrPass:NumCarr                             %Generating different data for each carrier
                if CarrUsedDo(kk)
                    if SendingData == 1                                    % Frist it is chosen to transmite a random information...
                        TxData = (randi(2,1,NbDQPSK)-1);                   %Creating the data stream to be transmited
                        TxData(1:JusPos) = JusVal;                         %Adding the Justification simble at the begining of the stream to sincronize received data frame
                        TxData(end-(JusPos-1):end) = JusValEnd;            %Adding the Justification simble at the end of the stream to sincronize received data frame
                        PreBits  = [0 0];                                  %Seting the start point for the DQPSK maping
                    elseif SendingData == 2
%For previous test it was needed a known bit sequence. This
%section creates it. The [I Q] output will be:
%[1 1 0 0 0 1 0 0 0 0 1 0 0 1 0 1 1 1 1 1 0 1 1 0 1 1 1 0 1
% 0 0 0]
%This sequency covers all possible combinations of I and Q
%which is:
%I = 1 0 0 0 0 1 0 0 1 1 0 1 1 1 1 0
%Q = 1 0 1 0 0 0 1 1 1 1 1 0 1 0 0 0
                        TxDataSeg = [0 0; %01
                                     0 0; %02
                                     0 1; %03
                                     1 0; %04
                                     1 1; %05
                                     1 0; %06
                                     0 0; %07
                                     1 1; %08
                                     0 1; %09
                                     1 1; %10
                                     1 0; %11
                                     0 0; %12
                                     1 0; %13
                                     0 1; %14
                                     1 1; %15
                                     0 1];%16
                        TxData   = [];
                        PreBits  = [0 0];
                        countseg = 1;
                        for jj=1:Nb
                            TxData   = [TxData TxDataSeg(countseg,:)];
                            countseg = countseg + 1;
                            if countseg > 16
                                countseg = 1;
                            end
                        end
                    else                                                   %... or just one high pulse for testing the channel
                        TxData = zeros(1,NbDQPSK);                         %Creating the base fo the data stream
                        TxData((end/2)-1:end/2) = [1 1];                   %Adding the code in the data stream to produce the highest level possible
                        PreBits  = [0 0];
                    end
                    TxDataPos       = linspace(1,NbDQPSK,NbDQPSK);         %Auxiliar variable to split the TxData
                    [DataI,DataQ]   = DqpskEncodEq(TxData,PreBits);        %Maping the income data to the DQPSK format
%Converting the I and Q components to the polirezed NRZ format
                    DataI(DataI==1) =  1.9;
                    DataI(DataI==0) = -1.9;
                    DataQ(DataQ==1) =  1.9;
                    DataQ(DataQ==0) = -1.9;
                    TxOdd  = TxData(logical(mod(TxDataPos,2)));            %Spliting the information of odd positions
                    TxEven = TxData(~(mod(TxDataPos,2)));                  %Spliting the information of even positions
                    TxDataMat(kk,:)         = TxOdd;                       %Storring the transmited information for latter evaluation
                    TxDataMat(kk+NumCarr,:) = TxEven;                      %Storring the transmited information for latter evaluation
%Once that we have the data to be transmited it needs to be converted to 
%the eletrical signal that hereafter will modulate our carrier through the 
%MZM-I
                    SigI = rectpulse(DataI,NPPB);
                    SigQ = rectpulse(DataQ,NPPB);

                    % Adding CP to the data
                    if AddCP
                        TxAux1 = reshape(SigI,NPPB,NbDQPSK/2);
                        if SetCpSampZer
                            TxAux1 = [zeros(NumAmosCP,size(TxAux1,2));...
                                   TxAux1;zeros(NumAmosCP,size(TxAux1,2))];
                        else
                            TxAux1 = [TxAux1(1:NumAmosCP,:);TxAux1;...
                                          TxAux1(end-(NumAmosCP-1):end,:)];
                        end
                        TxAux1 = reshape(TxAux1,1,(2*NumAmosCP+NPPB)*...
                                                                NbDQPSK/2);
                        SigI = [TxAux1 TxAux1(end-(StuffSampels-1):end)];

                        TxAux2 = reshape(SigQ,NPPB,NbDQPSK/2);
                        if SetCpSampZer
                            TxAux2 = [zeros(NumAmosCP,size(TxAux2,2));...
                                   TxAux2;zeros(NumAmosCP,size(TxAux2,2))];
                        else
                            TxAux2 = [TxAux2(1:NumAmosCP,:);TxAux2;...
                                          TxAux2(end-(NumAmosCP-1):end,:)];
                        end
                        TxAux2 = reshape(TxAux2,1,(2*NumAmosCP+NPPB)*...
                                                                NbDQPSK/2);
                        SigQ = [TxAux2 TxAux2(end-(StuffSampels-1):end)];
                    end
                
%Because of reasons, the real world does not have components capable 
%instantly change the voltage level. Thus, to emulate this behaviour the 
%eletrical PAM signal pass through a gaussian filter that will remove the 
%frequencies of higher order, which will result in small slope in the level 
%variation
                    [BitFilt,~] = FiltroGaussiano(f,BWD,CenFeq,FiltOrd);   %Creating filter for conformation of the input information
                    BitFilt = fftshift(BitFilt);                           %Doing a shift on the vector for matching the transmited data
                    TxSigI = SigI;                                         %Conforming the information and Creating the modulation signal
                    TxSigQ = SigQ;                                         %Conforming the information and Creating the modulation signal
                    %TxSigI = ifft(fft(SigI).*BitFilt);                    %Conforming the information and Creating the modulation signal
                    %TxSigQ = ifft(fft(SigQ).*BitFilt);                    %Conforming the information and Creating the modulation signal
                    
%                     if (kk==1)&&(Ploting)
%                       figure;hold all;grid on;
%                       plot(t,TxSig1,t,TxSig2)
%                     end
                    if AddingNoiseE
                        EsN0 = 10*log10(snr);
                        [~,sigpower] = MeasPower(TxSigI);
                        [~,sigpower2] = MeasPower(TxSigQ);
                        TxSigI = awgn(TxSigI,snr,sigpower-30);
                        TxSigQ = awgn(TxSigQ,snr,sigpower2-30);
                        [~,sigpower] = MeasPower(TxSigI);
                        [~,sigpower2] = MeasPower(TxSigQ);
                    end
%                     if (kk==1)&&(Ploting)
%                       plot(t,TxSig1,t,TxSig2)
%                       PrintInfo(Ploting*3,TxSigI(1:end-StuffSampels),2*t(2*NumAmosCP+NPPB),2*NumAmosCP+NPPB,TxSigQ(1:end-StuffSampels));
%                       a=4;
%                     end
                    [EoutModAux] = IqMod(EoutTx(VetThisCarrTx==(RefCarr-...
                                            1+kk),:),TxSigI,TxSigQ,Vpi,V0);
%                     if (kk==1)&&(Ploting)
%                         figure;hold all;grid on;
%                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
%                     end
                    if AddingNoiseO
                        EsN0 = 10*log10(snr);
                        [~,sigpower] = MeasPower(EoutModAux);
                        EoutModAux = awgn(EoutModAux,osnr,sigpower-30);
                        [~,sigpower] = MeasPower(EoutModAux);
                    end
%                     if (kk==1)&&(Ploting)
% %                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
%                         PrintInfo(Ploting*4,abs(EoutModAux(kk,1:end-...
%                         StuffSampels)).^2./max(abs(EoutModAux(kk,1:end-...
%                         StuffSampels)).^2),t(2*NumAmosCP+NPPB),2*...
%                         NumAmosCP+NPPB,abs(EoutModAux(1:end-StuffSampels...
%                         )).^2./max(abs(EoutModAux(kk,1:end-StuffSampels)...
%                                                                     ).^2));
%                         a=4;
%                     end
                    EoutModTem(ObsCarrPos(kk),1:length(EoutModAux)) = ...
                                                                EoutModAux;%Adding the current Modulated output to the final OutPut
                    %% Taking the sampling the EVM meassurement
                    PhaDel = 1*pi/4;
                    TimDel = T;
                    [EoutA,EoutB] = DelayInterf(f,TimDel,PhaDel,EoutModAux);%Coverting phase shift to amplitude variation
                    PhaDel = -1*pi/4;
                    TimDel = T;
                    [EoutC,EoutD] = DelayInterf(f,TimDel,PhaDel,EoutModAux);%Coverting phase shift to amplitude variation
                    EoutA = EoutA.*conj(EoutA);
                    EoutB = EoutB.*conj(EoutB);
                    EoutC = EoutC.*conj(EoutC);
                    EoutD = EoutD.*conj(EoutD);                   %Shifting the filter for matching the received signal
                    EoutA = ifft(fft(EoutA).*EvmFilt);
                    EoutB = ifft(fft(EoutB).*EvmFilt);
                    EoutC = ifft(fft(EoutC).*EvmFilt);
                    EoutD = ifft(fft(EoutD).*EvmFilt);
                    EoutI = (EoutB - EoutA);
                    EoutQ = (EoutD - EoutC);
                    EoutI = EoutI-mean(EoutI);
                    EoutQ = EoutQ-mean(EoutQ);
                    EoutI = EoutI./max(abs(EoutI));                                %Normalizing the signal
                    EoutQ = EoutQ./max(abs(EoutQ));                                %Normalizing the signal
                    PosAuxEout      = (2*NumAmosCP+NPPB)/2:(2*NumAmosCP+...
                                   NPPB):(length(EoutModAux)-StuffSampels);%Varriable respossible to take just the samples at the middle of the symbol
                    IxAux = EoutI(PosAuxEout) + 1j.*EoutQ(PosAuxEout);  
                    
                    TxSymbAmos(CurrentTest,ObsCarrPos==kk,:) = IxAux;%RxSymbAmosDo = [];
                    
                    EvmMatRef(ObsCarrPos==kk,:) = IxAux;       %Taking just the middle samples as references
%                     if (kk==1)&&(Ploting)
%                         figure;hold all;grid on;
%                         plot(real(IxAux),imag(IxAux),'rx');
%                         set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%                         a=0;
%                     end
                    %%
                    if (~mod(CarrPass,2))&&(SendingUp)
                        clear EoutModAux;
                        EoutModAux = EoutTx(VetThisCarrTx==(RefCarr-1+kk+1),:);
%                     if (kk==1)&&(Ploting)
%                         figure;hold all;grid on;
%                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
%                     end
                    if AddingNoiseO
                        EsN0 = 10*log10(snr);
                        [~,sigpower] = MeasPower(EoutModAux);
                        EoutModAux = awgn(EoutModAux,osnr,sigpower-30);
                        [~,sigpower] = MeasPower(EoutModAux);
                    end
%                     if (kk==1)&&(Ploting)
%                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
%                         PrintInfo(Ploting*4,abs(EoutModAux(kk,1:end-...
%                         StuffSampels)).^2./max(abs(EoutModAux(kk,1:end-...
%                         StuffSampels)).^2),2*t(2*NumAmosCP+NPPB),2*...
%                         NumAmosCP+NPPB,abs(EoutModAux(1:end-StuffSampels...
%                         )).^2./max(abs(EoutModAux(kk,1:end-StuffSampels)...
%                                                                     ).^2));
%                         a=4;
%                     end
                        EoutModTem(ObsCarrPos(kk+1),1:length(EoutModAux)...
                                                            ) = EoutModAux;%Adding the current Modulated output to the final OutPut
                    end
                else
                    EoutModAux = EoutTx(VetThisCarrTx==(RefCarr-1+kk),:);
                    
%                     if (kk==1)&&(Ploting)
%                         figure;hold all;grid on;
%                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
%                     end
                    if AddingNoiseO
                        EsN0 = 10*log10(snr);
                        [~,sigpower] = MeasPower(EoutModAux);
                        EoutModAux = awgn(EoutModAux,osnr,sigpower-30);
                        [~,sigpower] = MeasPower(EoutModAux);
                    end
%                     if (kk==1)&&(Ploting)
% %                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
%                         PrintInfo(Ploting*4,abs(EoutModAux(kk,1:end-...
%                         StuffSampels)).^2./max(abs(EoutModAux(kk,1:end-...
%                         StuffSampels)).^2),t(2*NumAmosCP+NPPB),2*...
%                         NumAmosCP+NPPB,abs(EoutModAux(1:end-StuffSampels...
%                         )).^2./max(abs(EoutModAux(kk,1:end-StuffSampels)...
%                                                                     ).^2));
%                         a=4;
%                     end
                    EoutModTem(ObsCarrPos(kk),1:length(EoutModAux)) = ...
                                                                EoutModAux;%Adding the current Modulated output to the final OutPut
                    if ~mod(CarrPass,2)
                        clear EoutModAux;
                        EoutModAux = EoutTx(VetThisCarrTx==(RefCarr-1+kk...
                                                                    +1),:);
%                     if (kk==1)&&(Ploting)
%                         figure;hold all;grid on;
%                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
%                     end
                    if AddingNoiseO
                        EsN0 = 10*log10(snr);
                        [~,sigpower] = MeasPower(EoutModAux);
                        EoutModAux = awgn(EoutModAux,osnr,sigpower-30);
                        [~,sigpower] = MeasPower(EoutModAux);
                    end
%                     if (kk==1)&&(Ploting)
%                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
%                         PrintInfo(Ploting*4,abs(EoutModAux(kk,1:end-...
%                         StuffSampels)).^2./max(abs(EoutModAux(kk,1:end-...
%                         StuffSampels)).^2),t(2*NumAmosCP+NPPB),2*...
%                         NumAmosCP+NPPB,abs(EoutModAux(1:end-StuffSampels...
%                         )).^2./max(abs(EoutModAux(kk,1:end-StuffSampels)...
%                                                                     ).^2));
%                         a=4;
%                     end
                        EoutModTem(ObsCarrPos(kk+1),1:length(EoutModAux)...
                                                            ) = EoutModAux;%Adding the current Modulated output to the final OutPut
                    end
                end
            end
            if size(EoutModTem,1)<NumCarr
                EoutModTem(end+1,:) = 0;
            end
            if IfftOrSum
%                 Ekaux=sum(EoutModTem);
%                 figure;hold all;
%                 plot(f,20*log10(abs(fftshift(fft(Ekaux)./length(Ekaux)))));
%                 set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                if size(EoutModTem,1)>1
                    [EoutMod,~,~] = OpticalIFFTN(f,T,MaxNumStagT,EoutModTem);
                else
                    EoutMod = EoutModTem;
                end
%                 plot(f,20*log10(abs(fftshift(fft(EoutMod)./length(...
%                                                               EoutMod)))));
%                 axis([0.9*fin 1.2*NumCarr*fc+fin -140 -20]);
%                 a=1;
            else
                if size(EoutModTem,1)>1
                    EoutMod = sum(EoutModTem);
                else
                    EoutMod = EoutModTem;
                end
            end
        case '4PAM'
            %%        Generate the data 4PAM
            EoutMod    = 0;
            EoutModTem = 0;
            for kk=InitCarrDo:CarrPass:NumCarr                             %Generating different data for each carrier
                if CarrUsedDo(kk)                                          % Frist it is chosen to transmite a random information...
                    if SendingData
                        TxData = (randi(2,1,Nb4Pam)-1);                    %Creating the data stream to be transmited
                        TxData(1:JusPos) = JusVal;                         %Adding the Justification simble at the begining of the stream to sincronize received data frame
                        TxData(end-(JusPos-1):end) = JusValEnd;            %Adding the Justification simble at the end of the stream to sincronize received data frame
                    else                                                   %... or just one high pulse for testing the channel
                        TxData = zeros(1,Nb4Pam);                          %Creating the base fo the data stream
                        TxDataPos = linspace(1,Nb4Pam,Nb4Pam);             %Creating a vector to auxilliate the addresing process
                        TxData(~mod(TxDataPos,2))=1;                       %Adding the code in the data stream to produce the lowest level possible
                        TxData(end/2 - 1) = 1;                             %Adding the code in the data stream to produce the highest level possible
                    end
                    
                    TxDataMat(kk,:) = TxData;                              %Storring the transmited information for latter evaluation
%Once that we have the data to be transmited it needs to be
%converted to the eletrical signal that hereafter will modulate
%our carrier through the MZM-I
                    if Which4PAM                                           %Selecting if the PAM will be done on Electrical or Optical domain
                        %  Good for others kms
                        [Phi1,Phi2] = Maping4PamIq(TxData,Vmin,Vmax,...
                                         ModSchem,FiberLength,SetCpSampZer);%Generating the eletrical signal for the optical PAM4
                    else
                        [Phi1,Phi2] = Maping4Pam(TxData,VPI,Polirized,...
                                                               MaxAmp4PAM);%Generating the eletrical signal for the electrical PAM4
                    end
%The signal generated are not yet with the same number of samples as the 
%OFCS loaded. These nexte lines do oversampling of the electrical signal.
                    TxSig1 = rectpulse(Phi1,NPPB);
                    TxSig2 = rectpulse(Phi2,NPPB);
%Thus, if it would be required to add cycle prefix the number of samples 
%per symbol needs to change as well as some adjustments needs to be done 
%for the new signal match in size with the size of the vector time. This 
%problem just exist on simulation, at practice the main point is the 
%syncronism of the signals.
                    % Adding CP to the data
                    if AddCP
                        TxAux1 = reshape(TxSig1,NPPB,Nb4Pam/2);
                        if SetCpSampZer
                            TxAux1 = [zeros(NumAmosCP,size(TxAux1,2));...
                                   TxAux1;zeros(NumAmosCP,size(TxAux1,2))];
                        else
                            TxAux1 = [TxAux1(1:NumAmosCP,:);TxAux1;...
                                          TxAux1(end-(NumAmosCP-1):end,:)];
                        end
                        TxAux1 = reshape(TxAux1,1,(2*NumAmosCP+NPPB)*...
                                                                 Nb4Pam/2);
                        TxSig1 = [TxAux1 TxAux1(end-(StuffSampels-1):end)];
                        TxAux2 = reshape(TxSig2,NPPB,Nb4Pam/2);
                        if SetCpSampZer
                            TxAux2 = [zeros(NumAmosCP,size(TxAux2,2));...
                                   TxAux2;zeros(NumAmosCP,size(TxAux2,2))];
                        else
                            TxAux2 = [TxAux2(1:NumAmosCP,:);TxAux2;...
                                          TxAux2(end-(NumAmosCP-1):end,:)];
                        end
                        TxAux2 = reshape(TxAux2,1,(2*NumAmosCP+NPPB)*...
                                                                 Nb4Pam/2);
                        TxSig2 = [TxAux2 TxAux2(end-(StuffSampels-1):end)];
                    end
                    
%                     if (kk==1)&&(Ploting)
%                       figure;hold all;grid on;
%                       plot(t,TxSig1,t,TxSig2)
%                     end
                    if AddingNoiseE
                        EsN0 = 10*log10(snr);
                        [~,sigpower] = MeasPower(TxSig1);
                        [~,sigpower2] = MeasPower(TxSig2);
                        TxSig1 = awgn(TxSig1,snr,sigpower-30);
                        TxSig2 = awgn(TxSig2,snr,sigpower2-30);
                        [~,sigpower] = MeasPower(TxSig1);
                        [~,sigpower2] = MeasPower(TxSig2);
                    end
%                     if (kk==1)&&(Ploting)
%                       plot(t,TxSig1,t,TxSig2)
%                       PrintInfo(Ploting*4,TxSig1(1:end-StuffSampels),2*t(2*NumAmosCP+NPPB),2*NumAmosCP+NPPB,TxSig2(1:end-StuffSampels));
%                     end
%Because of reasons, the real world does not have components capable 
%instantly change the voltage level. Thus, to emulate this behaviour the 
%eletrical PAM signal pass through a gaussian filter that will remove the 
%frequencies of higher order, which will result in small slope in the level 
%variation                    
                    [BitFilt,~] = FiltroGaussiano(f,BWD,CenFeq,FiltOrd);   %Creating filter for conformation of the input information
                    BitFilt     = fftshift(BitFilt);                       %Doing a shift on the vector for matching the transmited data
                    %TxSig1     = ifft(fft(TxSig1).*BitFilt);              %Conforming the information and Creating the modulation signal
                    %TxSig2     = ifft(fft(TxSig2).*BitFilt);              %Conforming the information and Creating the modulation signal
%                     if (kk==1)&&(Ploting)
%                         PrintInfo(Ploting*4,TxSig1(kk,1:end-...
%                         StuffSampels),t(2*NumAmosCP+NPPB),2*NumAmosCP+...
%                                           NPPB,TxSig2(1:end-StuffSampels));
%                         a=4;
%                     end

%If the user choses to modulate all carriers at the same time there is no 
%need of this script to generate data for each individual carrier. 
%Therefore, this For loop can be halted
                    U.U1t = TxSig2;                                        %Assigning the electrical signal to one drive of the MZM
                    U.U2t = TxSig1;                                        %Assigning the electrical signal to another drive of the MZM
                    if Which4PAM
                        if ModSchem
                            [EoutModAux] = ...
                            Mach_Zehnder_Modulator_simplificado_2(freqGHz,EoutTx...
                            (VetThisCarrTx==(RefCarr-1+kk),:),U,...
                                                            MZ_Input_File);
                        else
                            [EoutModAux] = ...
                            Mach_Zehnder_Modulator_simplificado_2(freqGHz,EoutTx...
                            (VetThisCarrTx==(RefCarr-1+kk),:),U,...
                                                            MZ_Input_File);
                        end
                    else
                        [EoutModAux] = ...
                        Mach_Zehnder_Modulator_simplificado_2(freqGHz,EoutTx(...
                         VetThisCarrTx==(RefCarr-1+kk),:),U,MZ_Input_File);
                    end
                    
                    if AddingNoiseO
                        EsN0 = 10*log10(snr);
                        [~,sigpower] = MeasPower(EoutModAux);
                        EoutModAux = awgn(EoutModAux,osnr,sigpower-30);
                        [~,sigpower] = MeasPower(EoutModAux);
                    end
%                     if (kk==1)&&(Ploting)
% %                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
%                         PrintInfo(Ploting*4,abs(EoutModAux(kk,1:end-...
%                         StuffSampels)).^2./max(abs(EoutModAux(kk,1:end-...
%                         StuffSampels)).^2),2*t(2*NumAmosCP+NPPB),2*...
%                         NumAmosCP+NPPB,abs(EoutModAux(1:end-StuffSampels...
%                         )).^2./max(abs(EoutModAux(kk,1:end-StuffSampels)...
%                                                                     ).^2));
%                         a=4;
%                     end
                    EoutModTem(ObsCarrPos(kk),1:length(EoutModAux)) = ...
                                                                EoutModAux;%Adding the current Modulated output to the final OutPut
                    %% Taking the sampling the EVM meassurement
                    PosAuxEout = (2*NumAmosCP+NPPB)/2:(2*NumAmosCP+NPPB)...
                                        :(length(EoutModAux)-StuffSampels);%Varriable respossible to take just the samples at the middle of the symbol
                    Ix         = EoutModAux.*conj(EoutModAux);             %Recovering the signal that will be transmited
                    Ix         = ifft(fft(Ix).*EvmFilt);                   %Removing higher signal generated by the receiving process
                    Ix         = Ix - min(Ix);
                    Ix         = Ix./max(abs(Ix));                         %Normalizing the reference
                    IxAux      = Ix(PosAuxEout);
                    a=0;  
                    
                    TxSymbAmos(CurrentTest,ObsCarrPos==kk,:) = IxAux;%RxSymbAmosDo = [];
                    EvmMatRef(ObsCarrPos==kk,:) = IxAux;                   %Taking just the middle samples as references
%                     if (kk==1)&&(Ploting)
%                         figure;hold all;grid on;
%                         plot(t,Ix);plot(t(PosAuxEout),IxAux,'x');
%                         set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%                         figure;hold all;grid on;
%                         plot(IxAux,zeros(1,length(IxAux)),'rx');
%                         set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%                         a=0;
%                     end
                    %%
                    if (~mod(CarrPass,2))&&(SendingUp)
                        clear EoutModAux;
                        EoutModAux = EoutTx(VetThisCarrTx==(RefCarr-1+kk...
                                                                    +1),:);
%                     if (kk==1)&&(Ploting)
%                         figure;hold all;grid on;
%                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
%                     end
                    if AddingNoiseO
                        EsN0 = 10*log10(snr);
                        [~,sigpower] = MeasPower(EoutModAux);
                        EoutModAux = awgn(EoutModAux,osnr,sigpower-30);
                        [~,sigpower] = MeasPower(EoutModAux);
                    end
%                     if (kk==1)&&(Ploting)
%                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
%                         PrintInfo(Ploting*4,abs(EoutModAux(kk,1:end-...
%                         StuffSampels)).^2./max(abs(EoutModAux(kk,1:end-...
%                         StuffSampels)).^2),2*t(2*NumAmosCP+NPPB),2*...
%                         NumAmosCP+NPPB,abs(EoutModAux(1:end-StuffSampels...
%                         )).^2./max(abs(EoutModAux(kk,1:end-StuffSampels)...
%                                                                     ).^2));
%                         a=4;
%                     end
                        EoutModTem(ObsCarrPos(kk+1),1:length(EoutModAux)...
                                                            ) = EoutModAux;%Adding the current Modulated output to the final OutPut
                    end
%                     if (kk==1)&&(Ploting)
%                         figure;hold all;grid on;
%                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
%                     end
                    if AddingNoiseO
                        EsN0 = 10*log10(snr);
                        [~,sigpower] = MeasPower(EoutModAux);
                        EoutModAux = awgn(EoutModAux,osnr,sigpower-30);
                        [~,sigpower] = MeasPower(EoutModAux);
                    end
%                     if (kk==1)&&(Ploting)
%                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
%                         PrintInfo(Ploting*4,abs(EoutModAux(kk,1:end-...
%                         StuffSampels)).^2./max(abs(EoutModAux(kk,1:end-...
%                         StuffSampels)).^2),2*t(2*NumAmosCP+NPPB),2*...
%                         NumAmosCP+NPPB,abs(EoutModAux(1:end-StuffSampels...
%                         )).^2./max(abs(EoutModAux(kk,1:end-StuffSampels)...
%                                                                     ).^2));
%                         a=4;
%                     end
                else
                    EoutModAux = EoutTx(VetThisCarrTx==(RefCarr-1+kk),:);
%                     if (kk==1)&&(Ploting)
%                         figure;hold all;grid on;
%                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
%                     end
                    if AddingNoiseO
                        EsN0 = 10*log10(snr);
                        [~,sigpower] = MeasPower(EoutModAux);
                        EoutModAux = awgn(EoutModAux,osnr,sigpower-30);
                        [~,sigpower] = MeasPower(EoutModAux);
                    end
%                     if (kk==1)&&(Ploting)
% %                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
%                         PrintInfo(Ploting*4,abs(EoutModAux(kk,1:end-...
%                         StuffSampels)).^2./max(abs(EoutModAux(kk,1:end-...
%                         StuffSampels)).^2),t(2*NumAmosCP+NPPB),2*...
%                         NumAmosCP+NPPB,abs(EoutModAux(1:end-StuffSampels...
%                         )).^2./max(abs(EoutModAux(kk,1:end-StuffSampels)...
%                                                                     ).^2));
%                         a=4;
%                     end
                    EoutModTem(ObsCarrPos(kk),1:length(EoutModAux)) = ...
                                                                EoutModAux;%Adding the current Modulated output to the final OutPut
                    if ~mod(CarrPass,2)
                        clear EoutModAux;
                        EoutModAux = EoutTx(VetThisCarrTx==(RefCarr-1+kk...
                                                                    +1),:);
%                     if (kk==1)&&(Ploting)
%                         figure;hold all;grid on;
%                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
%                     end
                    if AddingNoiseO
                        EsN0 = 10*log10(snr);
                        [~,sigpower] = MeasPower(EoutModAux);
                        EoutModAux = awgn(EoutModAux,osnr,sigpower-30);
                        [~,sigpower] = MeasPower(EoutModAux);
                    end
%                     if (kk==1)&&(Ploting)
%                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
%                         PrintInfo(Ploting*4,abs(EoutModAux(kk,1:end-...
%                         StuffSampels)).^2./max(abs(EoutModAux(kk,1:end-...
%                         StuffSampels)).^2),t(2*NumAmosCP+NPPB),2*...
%                         NumAmosCP+NPPB,abs(EoutModAux(1:end-StuffSampels...
%                         )).^2./max(abs(EoutModAux(kk,1:end-StuffSampels)...
%                                                                     ).^2));
%                         a=4;
%                     end
                        EoutModTem(ObsCarrPos(kk+1),1:length(EoutModAux)...
                                                            ) = EoutModAux;%Adding the current Modulated output to the final OutPut
                    end
                end
            end
            if size(EoutModTem,1)<NumCarr
                EoutModTem(end+1,:) = 0;
            end
            if IfftOrSum
%                 Ekaux=sum(EoutModTem);
%                 figure;hold all;
%                 plot(f,20*log10(abs(fftshift(fft(Ekaux)./length(Ekaux)))));
%                 set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                if size(EoutModTem,1)>1
                    [EoutMod,~,~] = OpticalIFFTN(f,T,MaxNumStagT,EoutModTem);
                else
                    EoutMod = EoutModTem;
                end
%                 plot(f,20*log10(abs(fftshift(fft(EoutMod)./length(...
%                                                               EoutMod)))));
%                 axis([0.9*fin 1.2*NumCarr*fc+fin -140 -20]);
%                 a=1;
            else
                if size(EoutModTem,1)>1
                    EoutMod = sum(EoutModTem);
                else
                    EoutMod = EoutModTem;
                end
            end
        otherwise
            %%        Generate the data OOK
            EoutMod    = 0;
            EoutModTem = 0;
            for kk=InitCarrDo:CarrPass:NumCarr                             %Generating different data for each carrier
                if CarrUsedDo(kk)                                          % Frist it is chosen to transmite just one high pulse for
                    if AdjusData                                           %testing the channel...
                        TxData        = zeros(1,Nb - NumBitDesc);          %creating vector of zerros for correcting the delay caused by the transmission
                        TxData(end/2) = 1;                                 %Adding an strategical pulse for measuring the result.
                        TxDataMat(kk,:) = TxData;                          %Storring the transmited information for latter evaluation
                        if NRZPolarity
                            TxData(TxData==0) = NrzMin;
                            TxData(TxData==1) = NrzMax;
                        end
                    else                                                   %... or just a random information
                        
                        TxData = (randi(2,1,Nb - NumBitDesc)-1);           %Creating Random Information that will be loaded in each individual subcarrier
                        TxData(1:JusLen)           = JusVal;               %Making the First 4 bits equal to zero
                        TxData(end-(JusLen-1):end) = JusVal;               %Making the Last 4 bits equal to zero
                        TxDataMat(kk,:) = TxData;                          %Storring the transmited information for latter evaluation
                        if NRZPolarity
                            TxData(TxData==0) = NrzMin;
                            TxData(TxData==1) = NrzMax;
                        end
                    end

%The signal generated are not yet with the same number of samples as the 
%OFCS loaded. These nexte lines do oversampling
                    TxDataRes = rectpulse(TxData,NPPB);                    %Changing the length of the Data acordingly with the time vector
%Thus, if it would be required to add cycle prefix the number of samples 
%per symbol needs to change as well as some adjustments needs to be done 
%for the new signal match in size with the size of the vector time. This 
%problem just exist on simulation, at practice the main point is the 
%syncronism of the signals.
                    % Adding CP to the data
                    if AddCP
                        TxAux = reshape(TxDataRes,NPPB,Nb-NumBitDesc);
                        if SetCpSampZer
                            TxAux = [zeros(NumAmosCP,size(TxAux,2));...
                                     TxAux;zeros(NumAmosCP,size(TxAux,2))];
                        else
                            TxAux = [TxAux(1:NumAmosCP,:);TxAux;TxAux(...
                                                 end-(NumAmosCP-1):end,:)];
                        end
                        TxAux = reshape(TxAux,1,(2*NumAmosCP+NPPB)*(Nb-...
                                                              NumBitDesc));
                        TxDataRes = [TxAux TxAux(end-(StuffSampels-1):...
                                                                     end)];
                    end
                    
                    
%Because of reasons, the real world does not have components capable 
%instantly change the voltage level. Thus, to emulate this behaviour the 
%eletrical PAM signal pass through a gaussian filter that will remove the 
%frequencies of higher order, which will result in small slope in the level 
                                                                 %variation
                    
                    [BitFilt,~] = FiltroGaussiano(f,BWD,CenFeq,FiltOrd);   %Creating filter for conformation of the input information
                    BitFilt = fftshift(BitFilt);                           %Doing a shift on the vector for matching the transmited data
%                   TxSig = ifft(fft(TxDataRes).*BitFilt);                 %Conforming the information and Creating the modulation signal
                    TxSig = TxDataRes;     
%                     if (kk==1)&&(Ploting)%Conforming the information and Creating the modulation signal
%                     figure;hold all;grid on;
%                     plot(t,TxSig)
%                     end
                    if AddingNoiseE
                        EsN0 = 10*log10(snr);
                        [~,sigpower] = MeasPower(TxSig);
                        TxSig = awgn(TxSig,snr,sigpower-30);
                        [~,sigpower] = MeasPower(TxSig);
                    end
%                     if (kk==1)&&(Ploting)
%                     plot(t,TxSig)
%                     PrintInfo(Ploting*4,TxSig(1:end-StuffSampels),2*t(2*NumAmosCP+NPPB),2*NumAmosCP+NPPB,TxSig(1:end-StuffSampels));
%                     end
%                   TxSigMat(kk,:) = TxSig;                                %Storring the transmited information for latter evaluation
%                   if (kk==1)&&(Ploting)
%                       PrintInfo(Ploting*5,t,TxSigMat(kk,:),TxDataRes);
%                   end
%If the user choses to modulate all carriers at the same time there is no 
%need of this script to generate data for each individual carrier. 
%Therefore, this For loop can be halted


%Assigning the eletrical signal to one drive of the MZM - The aplitude of 
%the signal can be controlled by the variable DatGai, which can be 
%understood as an gain for the eletrical signal or an atenuation. The 
%second signal will be similar with the only difference a phase shift of pi
                    U.U1t = TxSig;
                    U.U2t = exp(-1j*pi).*TxSig;
%As both signals will have the mostrly the same characteristics with the 
%only difference the phase shift of 180 degress. The MZM-I will be working 
%on the Push-Pull configuration. It is necessary to reduce the Chirp noise
%to zero.
                    [EoutModAux,~]=Mach_Zehnder_Modulator_simplificado_2(freqGHz...
                    ,EoutTx(VetThisCarrTx==(RefCarr-1+kk),:),U,...
                                                            MZ_Input_File);%Modulating individual carriers
%                     if (kk==1)&&(Ploting)
%                         figure;hold all;grid on;
%                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
%                     end
                    if AddingNoiseO
                        EsN0 = 10*log10(snr);
                        [~,sigpower] = MeasPower(EoutModAux);
                        EoutModAux = awgn(EoutModAux,osnr,sigpower-30);
                        [~,sigpower] = MeasPower(EoutModAux);
                    end
%                     if (kk==1)&&(Ploting)
% %                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
%                         PrintInfo(Ploting*4,abs(EoutModAux(kk,1:end-...
%                         StuffSampels)).^2./max(abs(EoutModAux(kk,1:end-...
%                         StuffSampels)).^2),t(2*NumAmosCP+NPPB),2*...
%                         NumAmosCP+NPPB,abs(EoutModAux(1:end-StuffSampels...
%                         )).^2./max(abs(EoutModAux(kk,1:end-StuffSampels)...
%                                                                     ).^2));
%                         a=4;
%                     end
                    
                    EoutModTem(ObsCarrPos(kk),1:length(EoutModAux)) = ...
                                                                EoutModAux;%Adding the current Modulated output to the final OutPut
                    %% Taking the sampling the EVM meassurement
                    PosAuxEout = (2*NumAmosCP+NPPB)/2:(2*NumAmosCP+NPPB)...
                                        :(length(EoutModAux)-StuffSampels);%Varriable respossible to take just the samples at the middle of the symbol
                    Ix         = EoutModAux.*conj(EoutModAux);             %Recovering the signal that will be transmited
                    Ix         = ifft(fft(Ix).*EvmFilt);                   %Removing higher signal generated by the receiving process
                    Ix         = Ix - min(Ix);
                    Ix         = Ix./max(abs(Ix));                         %Normalizing the reference
                    IxAux      = Ix(PosAuxEout);                           %Normalizing the reference
                    a=0;  
                    
                    TxSymbAmos(CurrentTest,ObsCarrPos==kk,:) = IxAux;%RxSymbAmosDo = [];
                    EvmMatRef(ObsCarrPos==kk,:) = IxAux;                   %Taking just the middle samples as references
%                     if (kk==1)&&(Ploting)
%                         figure;hold all;grid on;
%                         plot(IxAux,zeros(1,length(IxAux)),'rx');
%                         set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%                         a=0;
%                     end
                     %%
%                     if kk~=1
                    if (~mod(CarrPass,2))&&(SendingUp)
                        clear EoutModAux;
                        EoutModAux = EoutTx(VetThisCarrTx==(RefCarr-1+kk+1),:);
%                         if (kk==1)&&(Ploting)
%                         figure;hold all;grid on;
%                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
%                         end
                        if AddingNoiseO
                            EsN0 = 10*log10(snr);
                            [~,sigpower] = MeasPower(EoutModAux);
                            EoutModAux = awgn(EoutModAux,osnr,sigpower-30);
                            [~,sigpower] = MeasPower(EoutModAux);
                        end
%                         if (kk==1)&&(Ploting)
%                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
%                         PrintInfo(Ploting*4,abs(EoutModAux(kk,1:end-...
%                         StuffSampels)).^2./max(abs(EoutModAux(kk,1:end-...
%                         StuffSampels)).^2),t(2*NumAmosCP+NPPB),2*...
%                         NumAmosCP+NPPB,abs(EoutModAux(1:end-StuffSampels...
%                         )).^2./max(abs(EoutModAux(kk,1:end-StuffSampels)...
%                                                                     ).^2));
%                         end
                        EoutModTem(ObsCarrPos(kk+1),1:length(EoutModAux)...
                                                            ) = EoutModAux;%Adding the current Modulated output to the final OutPut
                    end
%                     a=1;
%                     else
%                     if ~mod(CarrPass,2)
%                         EoutModTem(ObsCarrPos(kk+1),1:length(EoutModAux)) = 0;%Adding the current Modulated output to the final OutPut
%                     end
%                     end
                else
                    EoutModAux = EoutTx(VetThisCarrTx==(RefCarr-1+kk),:);
%                     if (kk==1)&&(Ploting)
%                         figure;hold all;grid on;
%                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
%                     end
                    if AddingNoiseO
                        EsN0 = 10*log10(snr);
                        [~,sigpower] = MeasPower(EoutModAux);
                        EoutModAux = awgn(EoutModAux,osnr,sigpower-30);
                        [~,sigpower] = MeasPower(EoutModAux);
                    end
%                     if (kk==1)&&(Ploting)
%                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
%                         PrintInfo(Ploting*4,abs(EoutModAux(kk,1:end-...
%                         StuffSampels)).^2./max(abs(EoutModAux(kk,1:end-...
%                         StuffSampels)).^2),t(2*NumAmosCP+NPPB),2*...
%                         NumAmosCP+NPPB,abs(EoutModAux(1:end-StuffSampels...
%                         )).^2./max(abs(EoutModAux(kk,1:end-StuffSampels)...
%                                                                     ).^2));
%                         a=4;
%                     end
                    
                    EoutModTem(ObsCarrPos(kk),1:length(EoutModAux)) = ...
                                                                EoutModAux;%Adding the current Modulated output to the final OutPut
                    if ~mod(CarrPass,2)
                        clear EoutModAux;
                        EoutModAux = EoutTx(VetThisCarrTx==(RefCarr-1+kk+1),:);
%                         if (kk==1)&&(Ploting)
%                         figure;hold all;grid on;
%                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
%                         end
                        if AddingNoiseO
                            EsN0 = 10*log10(snr);
                            [~,sigpower] = MeasPower(EoutModAux);
                            EoutModAux = awgn(EoutModAux,osnr,sigpower-30);
                            [~,sigpower] = MeasPower(EoutModAux);
                        end
%                         if (kk==1)&&(Ploting)
%                         plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
%                         PrintInfo(Ploting*4,abs(EoutModAux(kk,1:end-...
%                         StuffSampels)).^2./max(abs(EoutModAux(kk,1:end-...
%                         StuffSampels)).^2),t(2*NumAmosCP+NPPB),2*...
%                         NumAmosCP+NPPB,abs(EoutModAux(1:end-StuffSampels...
%                         )).^2./max(abs(EoutModAux(kk,1:end-StuffSampels)...
%                                                                     ).^2));
%                         end
                        EoutModTem(ObsCarrPos(kk+1),1:length(EoutModAux)...
                                                            ) = EoutModAux;%Adding the current Modulated output to the final OutPut
                    end
%                     EoutModTem(ObsCarrPos(kk),1:length(EoutModAux)) = 0;%Adding the current Modulated output to the final OutPut
%                     if ~mod(CarrPass,2)
%                         EoutModTem(ObsCarrPos(kk+1),1:length(EoutModAux)) = 0;%Adding the current Modulated output to the final OutPut
%                     end
                end
            end
            if size(EoutModTem,1)<NumCarr
                EoutModTem(end+1,:) = 0;
            end
            if IfftOrSum
%                 Ekaux=sum(EoutModTem);
%                 figure;hold all;
%                 plot(f,20*log10(abs(fftshift(fft(Ekaux)./length(Ekaux)))));
%                 set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                if size(EoutModTem,1)>1
                    [EoutMod,~,~] = OpticalIFFTN(f,T,MaxNumStagT,EoutModTem);
                else
                    EoutMod = EoutModTem;
                end
%                 plot(f,20*log10(abs(fftshift(fft(EoutMod)./length(...
%                                                               EoutMod)))));
%                 axis([0.9*fin 1.2*NumCarr*fc+fin -140 -20]);
%                 a=1;
            else
                if size(EoutModTem,1)>1
                    EoutMod = sum(EoutModTem);
                else
                    EoutMod = EoutModTem;
                end
            end
    end
else
%There are two ways, implemented, to split the incoming OFCS. One by
%filters e another by the optical FFT. But this step just need to be done
%once in the simulation.
    if ~exist('EoutTx','var')                                              %Verify is this step was previouly done
        if SelecSetUP                                                      %Vefify is the split process will be by filter...
            EoutTx=SelectEachCarrier(EoutT,NumCarr,f,fin,FBWD,Order,fc);   %This function is resposible to split each carrier
            VetThisCarrTx = (RefCarr-1)+1:(RefCarr-1)+NumCarr;
        else                                                               %... or by the OFFT
            [EoutTx,~,VetThisCarrTx]=OpticalFFTN(f,T,MaxNumStagTx,EoutT);  %This function is resposible to split each carrier
        end
%         PrintInfo(Ploting*51,EoutTx,f);                                %Printing for qualitative analizes.
    end
    EoutMod = 0;
    if IfftOrSum
        [EoutMod,~,~] = OpticalIFFTN(f,T,MaxNumStagT,EoutTx);
    else
        EoutMod = sum(EoutTx);
    end
end
%%   Transmission of the OFDM Symble through a channel
% Having data stored and ready to be sent to end user. At the stage this
% script is responsible to chose the medium where this signal will travel.
% It may be withing an optical fiber or Back-toBack transmission.
% [~,PdBm] = MeasPower(EoutMod)
% [~,PdBm] = MeasPower(EoutMod)
switch Medium
    case 'B2B'
        EoutRec = EoutMod;
    case 'Fiber'
        [EoutRec,~,~]=Fibra_Monomodo1(t,EoutMod,lambdac,T,FiberLength,0,f);
%         PrintInfo(Ploting*11,EoutRec.*conj(EoutRec)./max(abs(EoutRec...
%                                                  .*conj(EoutRec))),T,NPPB);
    otherwise
        EoutRec = EoutMod;
end

EoutRec = EoutRec/SplitRatio;

if AddingNoiseF
    EsN0 = 10*log10(snr);
    [~,sigpower] = MeasPower(EoutRec);
    EoutRec = awgn(EoutRec,osnrf,sigpower-30);
    [~,sigpower] = MeasPower(EoutRec);
end

clear TxSigMat TxSig TxSigI TxSigQ EoutAux1 VetThisCarr EoutMod Ix...
    EoutModAux Esync Esync1 Esync2 EoutA EoutB EoutC EoutD EoutModTem
% clearvars -except fc Rb NPPB fsample T Ta NumberOf_T FinalTime Nb NbDPSK...
% TotalSamples t f Eout EoutMod EoutRec TxDataMat NumCarr NbDQPSK ValsLev...
% Modulation Polirized MaxAmp4PAM NumbOfTest UsedModula TestBundle Nb4Pam ...
% ThisModula CurrentTest CurrentModula ThisTest ThisModula OfcName BerOOK ...
% Ber4PAM BerDQPSK BerDPSK FiberDelay Medium CurrentMedium OSNRPC AddCP V0...
% OcsToTest CurrentOCS PulseResp NumAmosCP NumBitDesc StuffSampels JusPos ...
% contcp Ber4PAM PAM4Type Vpi MZ_Input_File ModAll JusValEnd NRZPolarity ...
% VetElecPower AberLev fin FFTSplit PastTime PowRef SnrRef CarSNR ValsLev2...
% AberLevI ValsLevI AberLevQ ValsLevQ VetElecPowerQ VetElecPowerI SigAten ...
%                                      lambdac Atenuation VetOptiPower JusVal