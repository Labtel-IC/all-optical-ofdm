close all;clear;clc;
%%              Gerando o sinal a propagar pela fibra
Nb = 2^12;
Rb = 1e9;
Tb = 1/Rb;
NPPB = 2^10;
NumberOf_T = Nb;
tf = Nb*Tb;
TotalSamples = NPPB*Nb;
t = linspace(0,tf,TotalSamples);
f = time2freq(t);
spc = Rb;
df = f(2)-f(1);
aux1 = round(spc/df);
fc = aux1*df;
T=1/fc;
Gvar = 1;
FiberLength = 100;
PlurySignals = 0;
lambda      = 1550e-9;
% To select the type of signal use set SigSelec to:
%  'Gaus'  --->  For gausian signal
%  'Tone'  --->  For a complex signal with one tone (frequency)
%  'Ofdm'  --->  For a OFDM signal

SigSelec = 'Ofdm';

switch SigSelec
    case 'Gaus'
        kk = 1;
        sinal=puls_gau(t,t(end/2),Gvar*T,1);
        if PlurySignals
            for kk=1:4
                sinal = sinal + puls_gau(t,t(end/2),Gvar*T*kk,1);
            end
        end
        T = Gvar*T*kk;
    case 'Tone'
        sinal =  exp(1j*2*pi*fc*t);
        if PlurySignals
            for kk=1:4
                sinal = sinal + exp(1j*2*pi*fc*t*kk);
            end
        end
    case 'Ofdm'
        sinal =  exp(1j*2*pi*fc*t);
        kk=1;
        SelecGaus = 0;
        OfdMod    = 'qam';                                                 %Chossing the type of modulation for the electrical subcarriers
        SelModTp  = 'BP';                                               %Selecting with OFDM will be in Base Band "BP", Amplitude Modulation 
                                                                           %with double side band "AM" or with single side band "AMSSB" when a 
                                                                           %different option is choosen this script will implement phase modulation.
        TapN   = 1;                                                        %Number of equalizer used
        ExtSam = 1;                                                        %Set to the number of frames that will be used.
        CpLe   = 0.9;                                                        %Set the percentage of electrical carriers that will be used
        if CpLe <= 0.5                                                     
            error('CpLe can not be bellow 50%. More than half of available carriers must be transmited')
        end
        BW     = fc;                                                       %Signal available bandwidth
        NFFT   = 2^10;                                                     %Size of the electical FFT
        if NFFT>Nb
            NFFT = Nb;
        end
        DatSiz = NFFT;                                                     %Set to the total number of electrical carriers to be used
        NpEl   = 2^0;
        Te     = T*NpEl; 
        M      = 4;                                                        %Modulation Level
        Tu     = NFFT*Te;                                                   %Time that will actually be used to transmit information
        Tg     = 0.0*Tu;                                                   %Time for the guard band
        if ((NFFT*NpEl)==Nb)&&(Tg~=0)
            error('NFFT must be smaller than Nb to use some ciclic prefix');
        end                                                    %Number of over samples for out of base band modulations
        Ts     = Tu + Tg;                                                  %Time of the OFDM symbol
        g      = Tu/Ts;                                                    %Percentage of band guard
        Dtf    = 1/Tu;                                                     %Chanal signaling ratio
        N      = NFFT/2 - 1;                                               %Number of carriers found
        Ns     = ceil((CpLe*(DatSiz))/2) - 1;                              %Number of carriers used
        DifNsN = ceil((N - Ns)/1);                                         %Number of carriers unused
        k      = log2(M);                                                  %Number of bits per symbol
        OBw    = NFFT/Tu;                                                  %OFDM bandwidth
        switch SelModTp
            case 'BP'
                Ofc    = 2*OBw;                                            %Central frequency of out of base band transmission
                Ofs    = 4*OBw;                                            %Electrical sampling frequency
                OvSam  = Tu*OBw*2^2;%Tu*(Ofc*2^1);  
                OvSamS = Ts*OBw*2^2;%Ts*(Ofc*2^1);  
            otherwise
                Ofc    = 2*OBw;                                            %Central frequency of out of base band transmission
                Ofs    = 4*Ofc;                                            %Electrical sampling frequency
                OvSam  = Tu*(Ofs*2^1);  
                OvSamS = Ts*(Ofs*2^1);  
        end
        nRb    = (OBw/(N+2))*(Ns*k);                                       %New Rb found
        %As not all electrical carries may be used, the user have the
        %choise to select where the OFDM signal will be centralized within
        %all carrier available. Zt set were the OFDM signal will start. If
        %Zt was set to zero, the OFDM signal will start at possition 1.
        %Thus, all unsued carriers will placed after it. If Zt was set to
        %10 the OFDM signal will start at electrical carrier 11. This is
        %important when using Hermitiam symetri. Zt will determinate how
        %further away the OFDM signal will be from the electrical carrier.
        if DifNsN > 51
            Zt  = 8;%floor(DifNsN/2);                                      
        else 
            Zt  = floor(DifNsN/2);
        end
        NumFra = floor((Nb*NPPB)/(OvSamS*NFFT));                                  %Number of frames
        NumFrU = floor((Nb*NPPB)/(OvSam*NFFT));                                  %Number of frames
        NuAmOf = (OvSam*NumFra*NFFT);                                      %Total number of samples of the electrical OFDM signal
        NPOFEX = (NuAmOf/NumFra)*(Tg/Tu);                                  %Extra samples per OFDM frame. It is the number of CP
        NuAmTo = NuAmOf*(1+Tg/Tu);                                         %Total number of samples of the electrical OFDM signal added CP
        NPPOF = floor((Nb*NPPB)/NuAmTo);                                   %Over sampling for optical modulation transmission
        NPSTUf = Nb*NPPB-NPPOF*NuAmTo;                                     %Number of extra sample to make all vector to have the same length
        NumFrE = NumFrU - NumFra;                                          %Number of frames
        MZ_Input_File = 2;                                                 %Inputdata file to the MZM
%         NuSaTs = Nb*NPPB;
%         NuSaTu = ceil(g*NPPB);
%         NuSaTg = NuSaTs - NuSaTu*Nb;
        
        Testando_Fibra_Eq_Find_Equalizer;
        
        TxData          = randi([0 M-1],NumFra,Ns);                        %Generation random information
        TxDataPara      = TxData.';                                        %Converting from serial to paralel 
        switch OfdMod                                                      %Sellect which modulation was used
            case 'qam'
                TxSymb  = qammod(TxDataPara,M);                            %Modulating information by QAM
            otherwise
                TxSymb  = dpskmod(TxDataPara,M);                           %Modulating information DPSK as default
        end
        %Construction the Hermitian Simetry
        TxSymbConj      = conj(flipud(TxSymb));                            %Singnal conjugate format at reverse order
        TxSymbH = [zeros(1,NumFra);TxSymb;zeros(1,NumFra);TxSymbConj];     %Hermitian Symetri
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
        TxSymbC1 = [TxSymb;TxSymb(1:DifNsN,:)];
        TxSymbConjC      = conj(flipud(TxSymbC1));
        TxSymbC = zeros(NFFT,NumFra);
        TxSymbC(1+Zt:Zt+size(TxSymbH,1)/2,:) = TxSymbH(1:size(TxSymbH,1)/2,:);
        TxSymbC(end-Zt-(size(TxSymbH,1)/2)+1:end-Zt,:) = TxSymbH(1+size(TxSymbH,1)/2:end,:);
%Here was implemented the modulation through FFT process, it basica mix-up
%all infomation following Fourier transform.
        TxSymbModA       = ifft(TxSymbC,NFFT);
        %Sellecting whether the OFDM signal will be transmited on
        %Base band or not
        switch SelModTp
            case 'BP'                                                      %Base Band transmission
                TxSymbModA  = rectpulse(TxSymbModA,OvSam);                 %Over sampling
                TxSymbMod = TxSymbModA;
                tt = linspace(0,Ts,size(TxSymbModA,1));
%                 TxSymbMod = TxSymbMod.*(exp(1j*2*pi*Ofc*(tt.')));
                ff = time2freq(tt);
                if SelecGaus
                    [ ModFilt ] = FiltroGaussiano(ff,(OBw*(1-Tg/Tu))/2/NumFra,Ofc,1);
                else
                    [ ModFilt ] = Filtro_Retangular((OBw*(1-Tg/Tu))/NumFra,Ofc,ff);
                end
                ModFilt = fftshift(ModFilt);
                if kk==1
                    figure;hold all;
                    plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                end
                for jj=1:NumFra
                    ModFilta(:,jj) = ModFilt;
                end
                TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
                if kk==1
                    plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                    plot(ff,20*log10(abs(fftshift(ModFilt))));
                    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                end
                a=0;
            case 'AM'                                                      %Out of base band
                TxSymbModA  = rectpulse(TxSymbModA,OvSam);                 %Over sampling
                [TxSymbMod,tt] = modulate(TxSymbModA,Ofc,Ofs,'am');
                ff = time2freq(tt);
                if SelecGaus
                    [ ModFilt ] = FiltroGaussiano(ff,3*Ofc,0,1);
                else
                    [ ModFilt ] = Filtro_Retangular( OBw,Ofc,ff);
                end
                ModFilt = fftshift(ModFilt);
                if kk==1
                    figure;hold all;
                    plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                end
                for jj=1:NumFra
                    ModFilta(:,jj) = ModFilt;
                end
                TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
                if kk==1
                plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                plot(ff,20*log10(abs(fftshift(ModFilt))));
                set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                end
            case 'AMSSB'                                                   %Out of base band
                TxSymbModA  = rectpulse(TxSymbModA,OvSam);                 %Over sampling
                [TxSymbMod,tt] = modulate(TxSymbModA,Ofc,Ofs,'amssb');
                ff = time2freq(tt);
                if SelecGaus
                    [ ModFilt ] = FiltroGaussiano(ff,2*Ofc,0,1);
                else
                    [ ModFilt ] = Filtro_Retangular( 4*Ofc,0,ff);
                end
                ModFilt = fftshift(ModFilt);
                if kk==1
                    figure;hold all;
                    plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                end
                for jj=1:NumFra
                    ModFilta(:,jj) = ModFilt;
                end
                TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
                if kk==1
                    plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                    plot(ff,20*log10(abs(fftshift(ModFilt))));
                    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                end
            otherwise
                TxSymbModA  = rectpulse(TxSymbModA,OvSam);                 %Over sampling
                [TxSymbMod,tt] = modulate(TxSymbModA,Ofc,Ofs,'pm');
                ff = time2freq(tt);
                if SelecGaus
                    [ ModFilt ] = FiltroGaussiano(ff,3*Ofc,0,1);
                else
                    [ ModFilt ] = Filtro_Retangular( 6*Ofc,0,ff);
                end
                ModFilt = fftshift(ModFilt);
                if kk==1
                    figure;hold all;
                    plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                end
                for jj=1:NumFra
                    ModFilta(:,jj) = ModFilt;
                end
                TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
                if kk==1
                    plot(ff,20*log10(abs(fftshift(fft(TxSymbMod(:,1))./length(TxSymbMod(:,1))))));
                    plot(ff,20*log10(abs(fftshift(ModFilt))));
                    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                end
        end
        AuxTxSymbMod(1,:)  = TxSymbMod(:);
        TxSymbModA  = [TxSymbMod(1+end-NPOFEX:end,:);TxSymbMod];           %Adding Ciclyc prefix
        SigTx  = rectpulse(TxSymbModA,NPPOF);                              %Over sampling
        ffa  = rectpulse(ff,NPPOF);                              %Over sampling
%Actualy transmiting the information Downstream
        TxSigA = SigTx(:).';                                               %Serializing signal
        TxSig = [TxSigA(1:NPSTUf) TxSigA];                                 %Conforming vector sizes to the same length
        ffa = [ff(1:NPSTUf) ffa];                                 %Conforming vector sizes to the same length
        T = Te;
    otherwise
end

figure;
hold on;
grid on;
plot(f,10*log10(abs((fftshift(fft(sinal./length(sinal)))))));

%%                   Transmission of the Signal

switch SigSelec
    case 'Ofdm'
        TxSig = 0.95*(TxSig./max(abs(TxSig)));                                     %Normalizing the sinal to mach with the MZM espect to receive.
        U.U1t = TxSig;
        U.U2t = exp(-1j*pi).*TxSig;
        %As both signals will have the mostrly the same char
        [sinalMod,~]=Mach_Zehnder_Modulator_simplificado(t,sinal,U,MZ_Input_File);  %Modulating individual carriers
    otherwise
        sinalMod = sinal;
end

[Eout,LD,LNL] = Fibra_Monomodo1(t,sinalMod,lambda,t(end),FiberLength,0,f);

% Eout = Eout./max(abs(Eout));
figure;
plot(f,20*log10(abs(fftshift(fft(sinalMod)./length(sinalMod)))));
% axis([-10e7 10e7 -200 0]);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);

figure;
hold all;
plot(t,sinalMod);
plot(t,abs(Eout)./max(abs(Eout)));
set(gcf,'units','normalized','outerposition',[0 0 1 1]);

figure;
hold all;
plot(abs(sinalMod./max(abs(sinalMod))),atan(imag(sinalMod)./(real(sinalMod)))*180/pi);
plot(abs(Eout./max(abs(Eout))),atan(imag(Eout)./(real(Eout)))*180/pi);
xlabel('Magnitude');
ylabel('Phase');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
figure;
hold all
plot(f,atan(imag(sinalMod)./(real(sinalMod)))*180/pi);
plot(f,atan(imag(Eout)./(real(Eout)))*180/pi);
xlabel('Frequency');
ylabel('Phase');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
a=1;
%%                  Receiving the signal
switch SigSelec
    case 'Ofdm'
        EoutAux = Eout;
        %The current incoming signal is them converted from the optical
        %domain to the eletrical domain with the help of an photo
        %detector.
        Ix = EoutAux.*conj(EoutAux);
        [BitFilt,~] = FiltroGaussiano(f,fc,0,1);                           %Creating filter for selection of the received signal
        BitFilt = fftshift(BitFilt);                                       %Shifting the filter for matching the received signal
%         Ix = ifft(fft(Ix).*BitFilt);                                       %Filtering the signal for better detection
        Ix = Ix - min(Ix);                                                 %Removing any off-set that may exist
        Ix = Ix./max(abs(Ix));                                             %Normalizing the eletrical signal (amplifying what is needed)

        figure;
        hold all;
        grid on;
        plot(f,20*log10(abs(fftshift(fft(Ix)./length(Ix)))));
        axis([-25e9 37.5e9 -200 0]);
        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
        %This section is missing the syncronisation process

        Ix = Ix(1+NPSTUf:end);
        Ix = reshape(Ix,length(Ix)/NumFra,NumFra);
        Ix = intdump(Ix,NPPOF);                                            %Downsampling the income signal
        SigRecepA = Ix(1+NPOFEX:end,:);

        switch SelModTp
            case 'BP'
                SigRecep  = intdump(SigRecepA,OvSam);                      %Over sampling
            case 'AM'
                SigRecepA = demod(SigRecepA,Ofc,Ofs,'am');
                SigRecep  = intdump(SigRecepA,OvSam);                      %Over sampling
            case 'AMSSB'
                SigRecepA = demod(SigRecepA,Ofc,Ofs,'amssb');
                SigRecep  = intdump(SigRecepA,OvSam);                      %Over sampling
            otherwise
                SigRecepA = demod(SigRecepA,Ofc,Ofs,'pm');
                SigRecep  = intdump(SigRecepA,OvSam);                      %Over sampling
        end
        
        figure(1);hold all;
        plot(ff,20*log10(abs(fftshift(fft(SigRecepA(:,1))./length(SigRecepA(:,1))))));
        SigRecep1 = fft(SigRecep,NFFT);                                    %FFT demultiplexing process
        SigRecep2 = SigRecep1(2+Zt:Ns+1+Zt,:);                             %Taking the actual data                               
        SigRecep3 = SigRecep2;
        if 1
            switch OfdMod                                                  %Modulation Tx signal for future comparison
                case 'qam'
                    TxDataA = TxData;
                    TxDataA = reshape(TxDataA,NumFra,length(TxDataA)/NumFra);
                    TxDataA = TxDataA.';
                    TxSigToPlot  = qammod(TxDataA,M);                      %Modulating information
                otherwise
                    TxDataA = TxData;
                    TxDataA = reshape(TxDataA,NumFra,length(TxDataA)/NumFra);
                    TxDataA = TxDataA.';
                    TxSigToPlot = dpskmod(TxDataA,M);                      %Modulating information
            end
        %         TxSigToPlot = dpskmod(TxDataMat(ThisCarr,2:end),M);
        if 1
            figure;
            hold all;
            plot(TxSigToPlot(:),'*','LineWidth',2);
            plot(SigRecep3(:),'o','LineWidth',2);
            set(gcf,'units','normalized','outerposition',[0 0 1 1]);
        end
        end

        switch OfdMod
            case 'qam'
                equ = reshape(equa,length(equa)/NumFra,NumFra);
                SigRecep3 = ((SigRecep3)./equ);%Fixing phase rotation and amplitude variation with the equalizer
                SigRecep4  = qamdemod(SigRecep3,M);                    %Modulating information
            otherwise
                SigRecep4 = dpskdemod(SigRecep3,M);                        %Modulating information
        end
        SigRecep4 = SigRecep4.';
        SigRecep4 = SigRecep4(:).';
        if 1
            figure;
            hold all;
            plot(TxSigToPlot(:),'*','LineWidth',2);
            plot(SigRecep3(:),'o','LineWidth',2);
            set(gcf,'units','normalized','outerposition',[0 0 1 1]);
            figure;
            hold all;
            plot(TxData);
            plot(SigRecep4);
            set(gcf,'units','normalized','outerposition',[0 0 1 1]);
        end
        SymbErr = ~(TxData==SigRecep4);                                    %Measuring system bit error ration
        BerOFDM = sum(SymbErr)/length(SymbErr);
        a=0;
        close all;
    otherwise
end












