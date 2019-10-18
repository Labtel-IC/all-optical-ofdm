function [ ChanelEqualizer ] = FindChanelElectricalResponse(ControlVars,FiberDelay,equ1,equ2)
%%                  FindChanelResponse
%This function is resposible to find the equalizer for each carrier
if nargin<2
    FiberDelay = [];
end
OfdMod        = ControlVars.OfdMod;
M             = ControlVars.M;
NFFT          = ControlVars.NFFT;
Ns            = ControlVars.Ns;
NumFra        = ControlVars.NumFra;
MZ_Input_File = ControlVars.MZ_Input_File;
%     SelecFilt     = ControlVars.SelecFilt;
Medium        = ControlVars.DifNsN;
DifNsN        = ControlVars.DifNsN;
fin           = ControlVars.fin;
FBWD          = ControlVars.FBWD;
Order         = ControlVars.Order;
FibLen        = ControlVars.FibLen;
Ts            = ControlVars.Ts;
Tu            = ControlVars.Tu;
OBw           = ControlVars.OBw;
Zt            = ControlVars.Zt;
OvSam         = ControlVars.OvSam;
Ofc           = ControlVars.Ofc;
Ofs           = ControlVars.Ofs;
SelModTp      = ControlVars.SelModTp;
NPOFEX        = ControlVars.NPOFEX;
NPPOF         = ControlVars.NPPOF;
NPSTUf        = ControlVars.NPSTUf;
SelecGaus     = ControlVars.SelecGaus;
Te            = ControlVars.Te;
ZtC           = ControlVars.ZtC;
DmtMve        = ControlVars.DmtMve;
NumFraPar     = ControlVars.NumFraPar;
freqGHz       = ControlVars.freqGHz;
VetSnr        = ControlVars.VetSnr;
UsingHermitian= ControlVars.UsingHermitian;
SNR           = ControlVars.SNR;
SelModFilt    = ControlVars.SelModFilt;
OfcName       = ControlVars.OfcName;
InitCarrDo    = ControlVars.InitCarrDo;
InitCarrUp    = ControlVars.InitCarrUp;
NumCarr       = ControlVars.NumCarr;
SelecSetUP    = ControlVars.SelecSetUP;
RefCarr       = ControlVars.RefCarr;
CarrPass      = ControlVars.CarrPass;
ObsCarrPos    = ControlVars.ObsCarrPos;
SendingUp     = ControlVars.SendingUp;
IfftOrSum     = ControlVars.IfftOrSum;
FFTSplit      = ControlVars.FFTSplit;
ReceptorNoise = ControlVars.ReceptorNoise;
ChanAwgn      = ControlVars.ChanAwgn;
UsingGpu      = ControlVars.UsingGpu;
fgpu          = ControlVars.fgpu;
CurTesSiz     = ControlVars.CurTesSiz;
ModFilta      = ControlVars.ModFilta;
tt            = ControlVars.tt;
ff            = ControlVars.ff;
tta           = ControlVars.tta;

%     VetSnr = 1e9;
%% Loading data
%     MainInputData;
%FiberLength = FibLen;
%load([OfcName '.mat']);
%f             = repmat(f.',1,CurTesSiz);
%t             = repmat(t.',1,CurTesSiz);
%if ~exist('freqGHz','var')
	%tps        = t/1E-12;
	%freqTHz    = time2freq_lamb_2(tps);
	%freqGHz    = freqTHz*1e-3;			                                  % Frequencia em GHz
	%freqGHz    = -fftshift(freqGHz);
	%freqGHz(1) = freqGHz(2);
%end
%if UsingGpu == 1
	%fgpu = gpuArray(f);
%else
	%fgpu = 0;
%end

%     VetSnr = 1e9;
Selecting = 1;
%% Loading data
% MainInputData;

%% Transmiting DownStream
TxSig = 0;
UniDmtMve = unique(DmtMve);
UniDmtMve = fliplr(UniDmtMve);
DaSiAn = 0;
DaSiPo = 0;
for CarDmtM=1:length(UniDmtMve)
    M = UniDmtMve(CarDmtM);
    DaSiAu = sum(ismember(DmtMve,M));
    TxData          = randi([0 M-1],NumFra*CurTesSiz,DaSiAu);                %Generation random information
    TxDataMat(1,1+DaSiAn:DaSiAn + CurTesSiz*DaSiAu*NumFra) = TxData(:);%Saving data for future comparison
    DaSiAn = DaSiAn + CurTesSiz*DaSiAu*NumFra;
	switch OfdMod                                              %Sellect which modulation was used
		case 'qam'
			TxSymbAux(1:NumFra*CurTesSiz,1+DaSiPo:DaSiPo + DaSiAu)  = qammod(TxData,M);                    %Modulating information by QAM
		otherwise
			TxSymbAux(1:NumFra*CurTesSiz,1+DaSiPo:DaSiPo + DaSiAu)  = dpskmod(TxData,M);                   %Modulating information DPSK as default
    end
    DaSiPo = DaSiPo + DaSiAu;
end
%for CarDmtM=1:length(DmtMve)
    %M = DmtMve(CarDmtM);
    %TxData          = randi([0 M-1],NumFra,1);                %Generation random information
    %TxDataMat(1,1+NumFra*(CarDmtM-1):NumFra*CarDmtM) = TxData(:);                               %Saving data for future comparison
    %switch OfdMod                                              %Sellect which modulation was used
        %case 'qam'
            %TxSymbAux(1:NumFra,CarDmtM)  = qammod(TxData,M);                    %Modulating information by QAM
        %otherwise
            %TxSymbAux(1:NumFra,CarDmtM)  = dpskmod(TxData,M);                   %Modulating information DPSK as default
    %end
%end
TxSymb      = TxSymbAux.';                                %Converting from serial to paralel
%EvmMatRef(1,:) = TxSymb(:);
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
    TxSymbH = [zeros(1,NumFra*CurTesSiz);TxSymb;zeros(1,NumFra*CurTesSiz);TxSymbConj];%Hermitian Symetri
    TxSymbC = zeros(NFFT,NumFra*CurTesSiz);
    TxSymbC(1+Zt-ZtC:Zt-ZtC+size(TxSymbH,1)/2,:) = TxSymbH(1:size(TxSymbH,1)/2,:);
    TxSymbC(end-Zt+ZtC-(size(TxSymbH,1)/2)+1:end-Zt+ZtC,:) = TxSymbH(1+size(TxSymbH,1)/2:end,:);
else
    TxSymbH = TxSymb;
%     TxSymbC = [zeros((NFFT-size(TxSymbH,1))/2,NumFra*CurTesSiz);TxSymbH;zeros((NFFT-size(TxSymbH,1))/2,NumFra*CurTesSiz)];
    TxSymbC = [TxSymbH(1:size(TxSymbH,1)/2,:);zeros((NFFT-size(TxSymbH,1)),NumFra*CurTesSiz);TxSymbH(end+1-size(TxSymbH,1)/2:end,:)];
end
%Here was implemented the modulation through FFT process, it basica mix-up
%all infomation following Fourier transform.
TxSymbModA       = ifft(TxSymbC,NFFT);
%Sellecting whether the OFDM signal will be transmited on
%Base band or not
switch SelModTp
	case 'BP'                                              %Base Band transmission
		TxSymbModA  = rectpulse(TxSymbModA,OvSam);         %Over sampling
		TxSymbMod = TxSymbModA;
		if SelModFilt
			TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
		end
	case 'AM'                                              %Out of base band
		TxSymbModA  = rectpulse(TxSymbModA,OvSam);         %Over sampling
		TxSymbMod   = TxSymbModA.*cos(2*pi*Ofc*tta);
		if SelModFilt
			TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
		end
	case 'AMSSB'                                           %Out of base band
		TxSymbModA  = rectpulse(TxSymbModA,OvSam);         %Over sampling
		TxSymbMod   = TxSymbModA.*exp(1j*2*pi*Ofc*tta);
		if SelModFilt
			TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
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
		for jj=1:NumFra
			ModFilta(:,jj) = ModFilt;
		end
		TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
end

TxSymbModA  = [TxSymbMod(1+end-NPOFEX:end,:);TxSymbMod];
U1t  = rectpulse(TxSymbModA,NPPOF);%Over sampling
TxSigA = reshape(U1t,OvSam*NumFra*NFFT,CurTesSiz);
U1t = [TxSigA(1+end-NPSTUf:end,:);TxSigA];                         %Conforming vector sizes to the same length
NormFact = max(U1t)./0.95;
NormFact = repmat(NormFact,size(U1t,1),1);
U1t = (U1t./NormFact);


if ChanAwgn
    [~,SigPower] = MeasPower(U1t);
    SigPower2 = SigPower-30;%20*log10(SigPowerI);
%     SNR = CarSNR + 10*log10(log2(max(DmtMve))) - 10*log10(1*Tu*OBw);
    EoutAux = awgn(U1t,SNR,SigPower2);
else
    EoutAux = U1t;
end
%%  Receiving DownStream
% EoutAux = EoutRec;
%         sn=((1/2)*(((1/2)^MaxNumStag)-1))/((1/2)-1);
%         DiffPos = round((sn/(2/2^IfftOrSum))*T/Ta);
%         EoutAux = [EoutAux(DiffPos+1:end) EoutAux(1:DiffPos)];
%The current incoming signal is them converted from the optical
%domain to the eletrical domain with the help of an photo
%detector.
Ix = EoutAux.*NormFact;
% NormFact = max(Ix);
% NormFact = repmat(NormFact,size(Ix,1),1);
% Ix = Ix./NormFact;
if ReceptorNoise                                               %Verify whether the noise is to be added or not
    [~,SigPower] = MeasPower(Ix);
    SigPower2 = SigPower-30;%20*log10(SigPowerI);
    Ix = awgn(Ix,SNR,SigPower2);
end
Ix = Ix(1+NPSTUf:end,:);
Ix = reshape(Ix,OvSam*NFFT,NumFra*CurTesSiz);
% Ix = intdump(Ix,NPPOF);                            %Downsampling the income signal
SigRecepA = Ix(1+NPOFEX:end,:);
switch SelModTp
	case 'BP'
		SigRecep  = intdump(SigRecepA,OvSam);
	case 'AM'
		[BitFiltEle,~] = FiltroGaussiano(ff,OBw/2,0,1);           %Creating filter for selection of the received signal
		BitFiltEle = fftshift(BitFiltEle);
		ModFilta = repmat(BitFiltEle.',1,NumFra*CurTesSiz);
		SigRecepA   = SigRecepA.*cos(-2*pi*Ofc*tta);
		if SelModFilt
			SigRecepA = ifft(fft(SigRecepA).*ModFilta);
		end
		SigRecep  = intdump(SigRecepA,OvSam);              %Over sampling
	case 'AMSSB'
		[BitFiltEle,~] = FiltroGaussiano(ff,OBw/2,0,1);           %Creating filter for selection of the received signal
		BitFiltEle = fftshift(BitFiltEle);
		ModFilta = repmat(BitFiltEle.',1,NumFra*CurTesSiz);
		SigRecepA   = SigRecepA.*exp(-1j*2*pi*Ofc*tta);
		if SelModFilt
			SigRecepA = ifft(fft(SigRecepA).*ModFilta);
		end
		SigRecep  = intdump(SigRecepA,OvSam);              %Over sampling
	otherwise
		SigRecepA = demod(SigRecepA,Ofc,Ofs,'pm');
		SigRecep  = intdump(SigRecepA,OvSam);              %Over sampling
end
SigRecep1 = fft(SigRecep,NFFT);                              %FFT demultiplexing process
if UsingHermitian
    SigRecep2 = SigRecep1(2+Zt-ZtC:Ns+1+Zt-ZtC,:);                               %Taking the actual data
else
%     SigRecep2 = SigRecep1(1+(NFFT-size(TxSymbH,1))/2:size(TxSymbH,1)+(NFFT-size(TxSymbH,1))/2,:);
    SigRecep2 = [SigRecep1(1:size(TxSymbH,1)/2,:);SigRecep1(end+1-size(TxSymbH,1)/2:end,:)];
end
SigRecep3 = SigRecep2;
%         SigRecep2 = SigRecep2(:).';
switch OfdMod
    case 'qam'
		TxDataA = TxDataMat(1,:);
		TxDataA = reshape(TxDataA,NumFra*CurTesSiz,length(TxDataA)/(NumFra*CurTesSiz));
		TxDataA = TxDataA.';
		UniDmtMve = unique(DmtMve);
		UniDmtMve = fliplr(UniDmtMve);
        DaSiAn =  0;
        DaSiPo = 0;
        for CarDmtM=1:length(UniDmtMve)
            M = UniDmtMve(CarDmtM);
            DaSiAu = sum(ismember(DmtMve,M));
            TxSigToPlot(1+DaSiPo:DaSiPo + DaSiAu,:) = qammod(TxDataA(1+DaSiPo:DaSiPo + DaSiAu,:),M);%Saving data for future comparison
            DaSiAn = DaSiAn + CurTesSiz*DaSiAu*NumFra;
            DaSiPo = DaSiPo + DaSiAu;
            %TxSigToPlot(1+DaSiAu*(CarDmtM-1):DaSiAu*CarDmtM,:)  = qammod(TxDataA(1+DaSiAu*(CarDmtM-1):DaSiAu*CarDmtM,:),M);
        end
		%for CarDmtM=1:length(UniDmtMve)
			%M = UniDmtMve(CarDmtM);
			%DaSiAu = sum(ismember(DmtMve,M));
			%TxSigToPlot(1+DaSiAu*(CarDmtM-1):DaSiAu*CarDmtM,:)  = qammod(TxDataA(1+DaSiAu*(CarDmtM-1):DaSiAu*CarDmtM,:),M);
		%end
    otherwise
        TxDataA = TxDataMat(1,:);
        TxDataA = reshape(TxDataA,NumFra,length(TxDataA)/NumFra);
        TxDataA = TxDataA.';
        %                 TxSig = reshape(TxDataMat(ThisCarr,:),NumFra,NFFT);
        TxSigToPlot  = dpskmod(TxDataMat(1,:),M);                           %Modulating information
        %                 TxSigToPlot  = TxSigToPlot(:);
end
equ        = (SigRecep3)./(TxSigToPlot);                               % Coeficientes do Equalizador

figure;
hold all;
plot(TxSigToPlot,'b*','LineWidth',2);
plot(SigRecep2,'o','LineWidth',2,'color',[1 0.4 0]);
set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
% figure;
% hold all;
% plot(TxSigToPlot(:),'b','LineWidth',2);
% plot(SigRecep2(:),'LineWidth',2,'color',[1 0.4 0]);
% set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);

ChanelEqualizer(1,:)        = equ(:);

clear EoutModTem EoutModAux EoutMod;
% %% Transmiting UpStream
% TxSig = 0;
% for CarDmtM=1:length(DmtMve)
%     M = DmtMve(CarDmtM);
%     TxData          = randi([0 M-1],NumFra,1);                %Generation random information
%     TxDataMat(2,1+NumFra*(CarDmtM-1):NumFra*CarDmtM) = TxData(:);                               %Saving data for future comparison
%     switch OfdMod                                              %Sellect which modulation was used
%         case 'qam'
%             TxSymbAux(1:NumFra,CarDmtM)  = qammod(TxData,M);                    %Modulating information by QAM
%         otherwise
%             TxSymbAux(1:NumFra,CarDmtM)  = dpskmod(TxData,M);                   %Modulating information DPSK as default
%     end
% end
% TxSymb      = TxSymbAux.';                                %Converting from serial to paralel
% EvmMatRef(2,:) = TxSymb(:);
% %Construction the Hermitian Simetry
% %The Hermitian simetry was composed of the signal its conjugate format and
% %zeros to fill empty spaces. The signal is formed as show below:
% %                    signal     midle  signal conjugated
% %                       |        |       |
% %                       V        V       V
% %            TxH = [0 a + jb 0 0 0 0 0 a - jb 0];
% %
% %I thought in many ways to form this signal, replacing the centred zeros by
% %copy of the signal bay just geting its information and adding
% %respectively. Also, by creating the signal with redundancy at the exact
% %length needed to make the hermitian singnal. But all tests end up with
% %similar results, no improvement was noticed.
% TxSymbConj      = conj(flipud(TxSymb));                    %Singnal conjugate format at reverse order
% if UsingHermitian
%     TxSymbH = [zeros(1,NumFra);TxSymb;zeros(1,NumFra);TxSymbConj];%Hermitian Symetri
%     TxSymbC = zeros(NFFT,NumFra);
%     TxSymbC(1+Zt-ZtC:Zt-ZtC+size(TxSymbH,1)/2,:) = TxSymbH(1:size(TxSymbH,1)/2,:);
%     TxSymbC(end-Zt+ZtC-(size(TxSymbH,1)/2)+1:end-Zt+ZtC,:) = TxSymbH(1+size(TxSymbH,1)/2:end,:);
% else
%     TxSymbH = TxSymb;
% %     TxSymbC = [zeros((NFFT-size(TxSymbH,1))/2,NumFra);TxSymbH;zeros((NFFT-size(TxSymbH,1))/2,NumFra)];
%     TxSymbC = [TxSymbH(1:size(TxSymbH,1)/2,:);zeros((NFFT-size(TxSymbH,1)),NumFra);TxSymbH(end+1-size(TxSymbH,1)/2:end,:)];
% end
% %Here was implemented the modulation through FFT process, it basica mix-up
% %all infomation following Fourier transform.
% TxSymbModA       = ifft(TxSymbC,NFFT);
% %Sellecting whether the OFDM signal will be transmited on
% %Base band or not
% switch SelModTp
%     case 'BP'                                              %Base Band transmission
%         TxSymbModA  = rectpulse(TxSymbModA,OvSam);         %Over sampling
%         %                         TxSymbMod = TxSymbModA;
%         tt = linspace(0,1*Te,size(TxSymbModA,1));
%         for jj=1:NumFra
%             tta(:,jj) = tt;
%         end
%         TxSymbMod   = TxSymbModA;
%         ff = time2freq(tt);
%         if SelecGaus
%             [ ModFilt ] = FiltroGaussiano(ff,OBw/2,0,1);
%         else
%             [ ModFilt ] = Filtro_Retangular(OBw,0,ff);
%         end
%         ModFilt = fftshift(ModFilt);
%         for jj=1:NumFra
%             ModFilta(:,jj) = ModFilt;
%         end
%         if SelModFilt
%             TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
%         end
%     case 'AM'                                              %Out of base band
%         TxSymbModA  = rectpulse(TxSymbModA,OvSam);         %Over sampling
%         tt = linspace(0,1*Te,size(TxSymbModA,1));
%         ff = time2freq(tt);
%         for jj=1:NumFra
%             tta(:,jj) = tt;
%         end
%         TxSymbMod   = TxSymbModA.*cos(2*pi*Ofc*tta);
%         if SelecGaus
%             [ ModFilt ] = FiltroGaussiano(ff,OBw,Ofc,1);
%         else
%             [ ModFilt ] = Filtro_Retangular( OBw,Ofc,ff);
%         end
%         ModFilt = fftshift(ModFilt);
%         for jj=1:NumFra
%             ModFilta(:,jj) = ModFilt;
%         end
%         if SelModFilt
%             TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
%         end
%     case 'AMSSB'                                           %Out of base band
%         TxSymbModA  = rectpulse(TxSymbModA,OvSam);         %Over sampling
%         tt = linspace(0,1*Te,size(TxSymbModA,1));
%         ff = time2freq(tt);
%         for jj=1:NumFra
%             tta(:,jj) = tt;
%         end
%         TxSymbMod   = TxSymbModA.*exp(1j*2*pi*Ofc*tta);
%         if SelecGaus
%             [ ModFilt ] = FiltroGaussiano(ff,OBw,Ofc,1);
%         else
%             [ ModFilt ] = Filtro_Retangular( OBw,Ofc,ff);
%         end
%         ModFilt = fftshift(ModFilt);
%         for jj=1:NumFra
%             ModFilta(:,jj) = ModFilt;
%         end
%         if SelModFilt
%             TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
%         end
%     otherwise
%         TxSymbModA  = rectpulse(TxSymbModA,OvSam);         %Over sampling
%         [TxSymbMod,tt] = modulate(TxSymbModA,Ofc,Ofs,'pm');
%         ff = time2freq(tt);
%         if SelecGaus
%             [ ModFilt ] = FiltroGaussiano(ff,3*Ofc,0,1);
%         else
%             [ ModFilt ] = Filtro_Retangular( 6*Ofc,0,ff);
%         end
%         ModFilt = fftshift(ModFilt);
%         for jj=1:NumFra
%             ModFilta(:,jj) = ModFilt;
%         end
%         TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
% end
% %
% %                 figure;plot(TxSymbMod(:,1));
% 
% TxSymbModA  = [TxSymbMod(1+end-NPOFEX:end,:);TxSymbMod];   %Adding Ciclyc prefix
% SigTx  = rectpulse(TxSymbModA,NPPOF);                      %Over sampling
% 
% TxSigA = SigTx(:).';                                       %Serializing signal
% TxSig = TxSig + [TxSigA(1:NPSTUf) TxSigA];                         %Conforming vector sizes to the same length
% % CpPaFr = CpPaFr + 1;
% %Assigning the eletrical signal to one drive of the MZM - The aplitude of
% %the signal can be controlled by the variable DatGai, which can be
% %understood as an gain for the eletrical signal or an atenuation. The
% %second signal will be similar with the only difference a phase shift of pi
% NormFact = max(TxSig);
% TxSig = (TxSig./NormFact);                     %Normalizing the sinal to mach with the MZM espect to receive.
% 
% % if  (PlotingThisThis)
% %     figure;hold all;grid on;
% %     plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
% %     axis([-25e9 37.5e9 -200 0]);
% %     set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
% % end
% 
% 
% if ChanAwgn
%     [~,SigPower] = MeasPower(TxSig);
%     SigPower2 = SigPower-30;%20*log10(SigPowerI);
%     SNR = CarSNR + 10*log10(log2(max(DmtMve))) - 10*log10(1*Tu*OBw);
%     EoutAux = awgn(TxSig,SNR,SigPower2);
% else
%     EoutAux = TxSig;
% end
% 
% %% Receving UpStream
% % EoutAux = EoutRec;
% %         sn=((1/2)*(((1/2)^MaxNumStag)-1))/((1/2)-1);
% %         DiffPos = round((sn/(2/2^IfftOrSum))*T/Ta);
% %         EoutAux = [EoutAux(DiffPos+1:end) EoutAux(1:DiffPos)];
% %The current incoming signal is them converted from the optical
% %domain to the eletrical domain with the help of an photo
% %detector.
% Ix = EoutAux.*NormFact;
% if ReceptorNoise                                               %Verify whether the noise is to be added or not
%     [~,SigPower] = MeasPower(Ix);
%     SigPower2 = SigPower-30;%20*log10(SigPowerI);
%     Ix = awgn(Ix,SNR,SigPower2);
% end
% Ix = Ix(1+NPSTUf:end);
% Ix = reshape(Ix,length(Ix)/NumFra,NumFra);
% Ix = intdump(Ix,NPPOF);                            %Downsampling the income signal
% SigRecepA = Ix(1+NPOFEX:end,:);
% %         IxAux = Ix(NuSaTg+1:end);
% switch SelModTp
%     case 'BP'
%         [BitFiltEle,~] = FiltroGaussiano(ff,OBw/2,0,1);        %Creating filter for selection of the received signal
%         %                     [ BitFiltEle ] = Filtro_Retangular(OBw,0,ff);
%         BitFiltEle = fftshift(BitFiltEle);
%         for jj=1:NumFra
%             ModFilta(:,jj) = BitFiltEle;
%         end
%         SigRecepB   = SigRecepA;
%         if SelModFilt
%             SigRecepB = ifft(fft(SigRecepB).*ModFilta);
%         end
%         SigRecep  = intdump(SigRecepB,OvSam);
%     case 'AM'
%         [BitFiltEle,~] = FiltroGaussiano(ff,OBw/2,0,1);           %Creating filter for selection of the received signal
%         BitFiltEle = fftshift(BitFiltEle);
%         for jj=1:NumFra
%             ModFilta(:,jj) = BitFiltEle;
%         end
%         SigRecepB   = SigRecepA.*cos(2*pi*Ofc*tta);
%         if SelModFilt
%             SigRecepB = ifft(fft(SigRecepB).*ModFilta);
%         end
%         SigRecep  = intdump(SigRecepB,OvSam);
%     case 'AMSSB'
%         [BitFiltEle,~] = FiltroGaussiano(ff,OBw/2,0,1);           %Creating filter for selection of the received signal
%         BitFiltEle = fftshift(BitFiltEle);
%         for jj=1:NumFra
%             ModFilta(:,jj) = BitFiltEle;
%         end
%         SigRecepB   = SigRecepA.*exp(-1j*2*pi*Ofc*tta);
%         if SelModFilt
%             SigRecepB = ifft(fft(SigRecepB).*ModFilta);
%         end
%         SigRecep  = intdump(SigRecepB,OvSam);              %Over sampling
%     otherwise
%         SigRecepA = demod(SigRecepA,Ofc,Ofs,'pm');
%         SigRecep  = intdump(SigRecepA,OvSam);              %Over sampling
% end
% SigRecep1 = fft(SigRecep,NFFT);                              %FFT demultiplexing process
% if UsingHermitian
%     SigRecep2 = SigRecep1(2+Zt-ZtC:Ns+1+Zt-ZtC,:);                               %Taking the actual data
% else
% %     SigRecep2 = SigRecep1(1+(NFFT-size(TxSymbH,1))/2:size(TxSymbH,1)+(NFFT-size(TxSymbH,1))/2,:);
%     SigRecep2 = [SigRecep1(1:size(TxSymbH,1)/2,:);SigRecep1(end+1-size(TxSymbH,1)/2:end,:)];
% end
% SigRecep3 = SigRecep2;
% %         SigRecep2 = SigRecep2(:).';
% switch OfdMod
%     case 'qam'
%         TxDataA = TxDataMat(2,:);
%         TxDataA = reshape(TxDataA,NumFra,length(TxDataA)/NumFra);
%         TxDataA = TxDataA.';
%         for CarDmtM=1:length(DmtMve)
%             M = DmtMve(CarDmtM);
%             TxSigToPlot(CarDmtM,:)  = qammod(TxDataA(CarDmtM,:),M);            %Modulating information
%         end
%     otherwise
%         TxDataA = TxDataMat(2,:);
%         TxDataA = reshape(TxDataA,NumFra,length(TxDataA)/NumFra);
%         TxDataA = TxDataA.';
%         %                 TxSig = reshape(TxDataMat(ThisCarr,:),NumFra,NFFT);
%         TxSigToPlot  = dpskmod(TxDataMat(2,:),M);                           %Modulating information
%         %                 TxSigToPlot  = TxSigToPlot(:);
% end
% equ        = (SigRecep3)./(TxSigToPlot);                               % Coeficientes do Equalizador
% 
%         figure;
%         hold all;
%         plot(TxSigToPlot,'*','LineWidth',2);
%         plot(SigRecep2,'o','LineWidth',2);
%         set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
% 
% ChanelEqualizer(2,:)        = equ(:);
% 
% clear EoutModTem EoutModAux EoutMod;

end

