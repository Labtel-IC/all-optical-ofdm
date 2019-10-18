
TxData          = randi([0 M-1],NumFra,Ns);                        %Generation random information
TxDataParaEq      = TxData.';                                        %Converting from serial to paralel
switch OfdMod                                                      %Sellect which modulation was used
    case 'qam'
        TxSymb  = qammod(TxDataParaEq,M);                            %Modulating information by QAM
    otherwise
        TxSymb  = dpskmod(TxDataParaEq,M);                           %Modulating information DPSK as default
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
        tt = linspace(0,NumFra*Ts,size(TxSymbModA,1));
        ff = time2freq(tt);
        if SelecGaus
            [ ModFilt ] = FiltroGaussiano(ff,(OBw*(1-Tg/Tu))/2/NumFra,0,1);
        else
            [ ModFilt ] = Filtro_Retangular((OBw*(1-Tg/Tu))/NumFra,0,ff);
        end
        ModFilt = fftshift(ModFilt);
        for jj=1:NumFra
            ModFilta(:,jj) = ModFilt;
        end
        TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
    case 'AM'                                                      %Out of base band
        TxSymbModA  = rectpulse(TxSymbModA,OvSam);                 %Over sampling
        [TxSymbMod,tt] = modulate(TxSymbModA,Ofc,Ofs,'am');
        ff = time2freq(tt);
        if SelecGaus
            [ ModFilt ] = FiltroGaussiano(ff,3*Ofc,0,1);
        else
            [ ModFilt ] = Filtro_Retangular( 6*Ofc,0,ff);
        end
        ModFilt = fftshift(ModFilt);
        for jj=1:NumFra
            ModFilta(:,jj) = ModFilt;
        end
        TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
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
        for jj=1:NumFra
            ModFilta(:,jj) = ModFilt;
        end
        TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
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
        for jj=1:NumFra
            ModFilta(:,jj) = ModFilt;
        end
        TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
end
AuxTxSymbMod(1,:)  = TxSymbMod(:);
TxSymbModA  = [TxSymbMod(1+end-NPOFEX:end,:);TxSymbMod];           %Adding Ciclyc prefix
SigTx  = rectpulse(TxSymbModA,NPPOF);                              %Over sampling
%Actualy transmiting the information Downstream
TxSigA = SigTx(:).';                                               %Serializing signal
TxSig = [TxSigA(1:NPSTUf) TxSigA];                                 %Conforming vector sizes to the same length
T = Te;
%%                   Transmission of the Signal

TxSig = 0.95*(TxSig./max(abs(TxSig)));                                     %Normalizing the sinal to mach with the MZM espect to receive.
U.U1t = TxSig;
U.U2t = exp(-1j*pi).*TxSig;
%As both signals will have the mostrly the same char
[sinalMod,~]=Mach_Zehnder_Modulator_simplificado(t,sinal,U,MZ_Input_File);  %Modulating individual carriers

[Eout,LD,LNL] = Fibra_Monomodo1(t,sinalMod,lambda,T,FiberLength,0,f);

%%                  Receiving the signal
EoutAux = Eout;
%The current incoming signal is them converted from the optical
%domain to the eletrical domain with the help of an photo
%detector.
Ix = EoutAux.*conj(EoutAux);
[BitFilt,~] = FiltroGaussiano(f,fc,0,1);                           %Creating filter for selection of the received signal
BitFilt = fftshift(BitFilt);                                       %Shifting the filter for matching the received signal
% Ix = ifft(fft(Ix).*BitFilt);                                       %Filtering the signal for better detection
Ix = Ix - min(Ix);                                                 %Removing any off-set that may exist
Ix = Ix./max(abs(Ix));                                             %Normalizing the eletrical signal (amplifying what is needed)

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
SigRecep1 = fft(SigRecep,NFFT);                                    %FFT demultiplexing process
SigRecep2 = SigRecep1(2+Zt:Ns+1+Zt,:);                             %Taking the actual data
SigRecep3 = SigRecep2;
equa        = (SigRecep3)./(TxSymb);












