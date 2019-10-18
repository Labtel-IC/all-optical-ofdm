function [ ChanelEqualizer ] = FindChanelResponse(ControlVars,FiberDelay,equ1,equ2)
%%                  FindChanelResponse
    %This function is resposible to find the equalizer for each carrier
    if nargin<2
        FiberDelay = [];
    end
    L             = ControlVars.L;
    U0            = ControlVars.U0;
    U_pi1         = ControlVars.U_pi1;
    U_pi2         = ControlVars.U_pi2;
    nopt          = ControlVars.nopt;
    nel           = ControlVars.nel;
    C             = ControlVars.C;
    alfa0         = ControlVars.alfa0;
    OfdMod        = ControlVars.OfdMod;
    M             = ControlVars.M;
    NFFT          = ControlVars.NFFT;
    Ns            = ControlVars.Ns;
    NumFra        = ControlVars.NumFra;
    MZ_Input_File = ControlVars.MZ_Input_File;
    SelecFilt     = ControlVars.SelecFilt;
    Medium        = ControlVars.DifNsN;
    DifNsN        = ControlVars.DifNsN;
    fin           = ControlVars.fin;
    FBWD          = ControlVars.FBWD;
    Order         = ControlVars.Order;
    FibLen        = ControlVars.FibLen;
    Ts            = ControlVars.Ts;
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
	ObsCarrUsed   = ControlVars.ObsCarrUsed;
	ReceptorNoise = ControlVars.ReceptorNoise;
	ChanAwgn      = ControlVars.ChanAwgn;
	UsingGpu      = ControlVars.UsingGpu;
	fgpu          = ControlVars.fgpu;
	CurTesSiz     = ControlVars.CurTesSiz;
	ModFilta      = ControlVars.ModFilta;
	tt            = ControlVars.tt;
	ff            = ControlVars.ff;
	tta           = ControlVars.tta;
	SplitRatio    = ControlVars.SplitRatio;
	UseIQ         = ControlVars.UseIQ;
	FcIQ          = ControlVars.FcIQ;
% 	FeqC          = ControlVars.FeqC;
    
%     VetSnr = 1e9;
    Selecting = 1;
%% Loading data
%     MainInputData;
    FiberLength = FibLen;
    load([OfcName '.mat']);
	f             = repmat(f.',1,CurTesSiz);
	t             = repmat(t.',1,CurTesSiz);
	if ~exist('freqGHz','var')
		freqGHz = zeros(size(t,1),size(t,2));
		for kk=1:CurTesSiz
			tps           = t(:,kk)/1E-12;
			freqTHzA      = time2freq_lamb_2(tps);
			freqGHzA      = freqTHzA*1e-3;			                                  % Frequencia em GHz
			freqGHzA      = -fftshift(freqGHzA);
			freqGHzA(1)   = freqGHzA(2);
			freqGHz(:,kk) = freqGHzA;
		end
	end
    InitCarr     = InitCarrDo;
    MaxNumStagTx = nextpow2(NumCarr+1);                                    %Next number, power of 2, by the given number of carriers
    % MaxNumCarrTx = 2^MaxNumStagTx;                                       %With the number of stages the actual nubmer of carriers can be calculated
    MaxNumStagT  = nextpow2(NumCarr);                                      %Next number, power of 2, by the given number of carriers
    MaxNumStag = nextpow2(NumCarr);                                        %Next number, power of 2, by the given number of carriers
%     MaxNumCarr = 2^MaxNumStag;                                           %With the number of stages the actual nubmer of carriers can be calculated
    if SelecSetUP                                                          %Vefify is the split process will be by filter...
        EoutTx=SelectEachCarrier(Eout,NumCarr,f(:,1).',fin,FBWD,Order,fc);        %This function is resposible to split each carrier
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
		
		if UsingGpu==1
			EoutGpu = gpuArray(EoutT);
			[EoutAux1Gpu,~,VetThisCarrGpu]=OpticalFFTG(fgpu(1:Nb*NPPB),gpuArray(T),gpuArray(MaxNumStagTx),EoutGpu);
			EoutTx = gather(EoutAux1Gpu);
			VetThisCarrTx = gather(VetThisCarrGpu);
			clear EoutGpu EoutAux1Gpu VetThisCarrGpu;
		else
			[EoutTx,~,VetThisCarrTx]=OpticalFFTN(f(1:Nb*NPPB),T,MaxNumStagTx,EoutT);  %This function is resposible to split each carrier
		end
        %[EoutTx,~,VetThisCarrTx]=OpticalFFTN(f(1:Nb*NPPB),T,MaxNumStagTx,EoutT);      %This function is resposible to split each carrier
    end
%% Transmiting DownStream
    TxDataMat = zeros(NumCarr,CurTesSiz*NumFra*length(DmtMve));
    for kk=InitCarrDo:CarrPass:NumCarr %Generating different data for each carrier
%         CpPaFr = 1;
		UniDmtMve = unique(DmtMve);
		UniDmtMve = fliplr(UniDmtMve);
        DaSiAn    = 0;
        DaSiPo    = 0;
		for CarDmtM=1:length(UniDmtMve)
			M = UniDmtMve(CarDmtM);
			DaSiAu = sum(ismember(DmtMve,M));
			TxData = randi([0 M-1],NumFra*CurTesSiz,DaSiAu);%Generation random information
			TxDataMat(kk,1+DaSiAn:DaSiAn + CurTesSiz*DaSiAu*NumFra) = TxData(:);%Saving data for future comparison
            DaSiAn = DaSiAn + CurTesSiz*DaSiAu*NumFra;
			switch OfdMod%Sellect which modulation was used
				case 'qam'
					TxSymbAux(1:NumFra*CurTesSiz,1+DaSiPo:DaSiPo + DaSiAu)  = qammod(TxData,M);%Modulating information by QAM
				otherwise
					TxSymbAux(1:NumFra*CurTesSiz,1+DaSiPo:DaSiPo + DaSiAu)  = dpskmod(TxData,M);%Modulating information DPSK as default
            end
            DaSiPo = DaSiPo + DaSiAu;
		end
        TxSymb = TxSymbAux.';
        
        TxSymbConj = conj(flipud(TxSymb));
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
        TxSymbModA  = ifft(TxSymbC,NFFT);
        
        TxSymbModA  = rectpulse(TxSymbModA,OvSam);         %Over sampling

        if UseIQ
            RealTx = real(TxSymbModA);
            ImagTx = imag(TxSymbModA);
            TxSymbModA = RealTx.*cos(2*pi*FcIQ.*tta) + ImagTx.*sin(2*pi*FcIQ.*tta);
        end
 
        switch SelModTp
            case 'BP'                                              %Base Band transmission
%                 TxSymbModA  = rectpulse(TxSymbModA,OvSam);         %Over sampling
                TxSymbMod = TxSymbModA;
                if SelModFilt
                    TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
                end
            case 'AM'                                              %Out of base band
%                 TxSymbModA  = rectpulse(TxSymbModA,OvSam);         %Over sampling
                TxSymbMod   = TxSymbModA.*cos(2*pi*Ofc*tta);
                if SelModFilt
                    TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
                end
            case 'AMSSB'                                           %Out of base band
%                 TxSymbModA  = rectpulse(TxSymbModA,OvSam);         %Over sampling
                TxSymbMod   = TxSymbModA.*exp(1j*2*pi*Ofc*tta);
                if SelModFilt
                    TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
                end
            otherwise
%                 TxSymbModA  = rectpulse(TxSymbModA,OvSam);         %Over sampling
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
 
        %AuxTxSymbMod(1,:)  = TxSymbMod(:);
        TxSymbModA  = [TxSymbMod(1+end-NPOFEX:end,:);TxSymbMod];
        U1t  = rectpulse(TxSymbModA,NPPOF);%Over sampling
		TxSigA = reshape(U1t,OvSam*NumFra*NFFT,CurTesSiz);
		U1t = [TxSigA(1+end-NPSTUf:end,:);TxSigA];                         %Conforming vector sizes to the same length
%         U1t = real(U1t.*exp(1j.*(2*pi*FeqC.*t))) + imag(U1t.*exp(1j.*(2*pi*FeqC.*t + pi/2)));
%Assigning the eletrical signal to one drive of the MZM - The aplitude of 
%the signal can be controlled by the variable DatGai, which can be 
%understood as an gain for the eletrical signal or an atenuation. The 
%second signal will be similar with the only difference a phase shift of pi
        ConNor = max(U1t)./0.95;
        ConNor = repmat(ConNor,size(U1t,1),1);
        U1t = (U1t./ConNor);
        U.U1t = U1t;
        U.U2t = exp(-1j*pi).*U1t;
%As both signals will have the mostrly the same characteristics with the 
%only difference the phase shift of 180 degress. The MZM-I will be working 
%on the Push-Pull configuration. It is necessary to reduce the Chirp noise
%to zero.
		EoutTxAux = repmat(EoutTx(VetThisCarrTx==(RefCarr-1+kk),:).',1,CurTesSiz);
		[EoutModAux,~]=MZM(freqGHz,EoutTxAux,U,L,U0,U_pi1,U_pi2,nopt,nel,C,alfa0);%Modulating individual carriers
        EoutModTem(1:length(EoutModAux),:,ObsCarrPos(kk)) = EoutModAux;      %Adding the current Modulated output to the final OutPut
%             if kk==1
%                 figure;hold all;grid on;
%                 plot(f,20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
%                 set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%             end
        if (~mod(CarrPass,2))&&(SendingUp)
            clear EoutModAux;
			EoutModAux = repmat(EoutTx(VetThisCarrTx==(RefCarr-1+kk+1),:).',1,CurTesSiz);
            EoutModTem(:,:,ObsCarrPos(kk+1)) = EoutModAux;%Adding the current Modulated output to the final OutPut
        end
    end

    if IfftOrSum
        if size(EoutModTem,3)>1
			if UsingGpu==1
                StagGpu = gpuArray(MaxNumStagT);
                Tgpu = gpuArray(T);
                EoutMod = zeros(size(EoutModTem,1),size(EoutModTem,2));
                for jj=1:size(EoutModTem,2)
                    EoutModTemGpu = gpuArray(EoutModTem(:,jj,:));
                    [EoutAux1Gpu,~,~]=OpticalIFFTGT(fgpu,Tgpu,StagGpu,EoutModTemGpu);
                    EoutMod(:,jj) = gather(EoutAux1Gpu);
                    clear EoutModTemGpu EoutAux1Gpu;
                end
				clear StagGpu Tgpu;
			else
				[EoutMod,~,~] = OpticalIFFTNT(f,T,MaxNumStagT,EoutModTem);
			end
            %[EoutMod,~,~] = OpticalIFFTN(f,T,MaxNumStagT,EoutModTem);
        else
            EoutMod = EoutModTem;
        end
    else
        if size(EoutModTem,3)>1
            EoutMod = sum(EoutModTem,3);
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
    
    EoutRec = EoutRec/SplitRatio;
%%  Receiving DownStream
    
    if FFTSplit
		if UsingGpu==1
            Tgpu = gpuArray(T);
            StagGpu = gpuArray(MaxNumStag);
            EoutAux1 = zeros(size(EoutRec,1),size(EoutRec,2),NumCarr);
            for jj=1:size(EoutRec,2)
                EoutGpu = gpuArray(EoutRec(:,jj));
                [EoutAux1Gpu,~,VetThisCarrGpu]=OpticalFFTGT(fgpu,Tgpu,StagGpu,EoutGpu);
                EoutAux1(:,jj,:) = gather(EoutAux1Gpu);
                VetThisCarr = gather(VetThisCarrGpu);
                clear EoutGpu EoutAux1Gpu VetThisCarrGpu;
            end
			clear Tgpu StagGpu;
		else
			[EoutAux1,~,VetThisCarr]=OpticalFFTNT(f,T,MaxNumStag,EoutRec);
		end
        %[EoutAux1,~,VetThisCarr]=OpticalFFTN(f,T,MaxNumStag,EoutRec);
    else
        VetThisCarr = ObsCarrPos;
        EoutAux1    = SelectEachCarrier(EoutRec,NumCarr,f,fin,1.0*fc,11,fc);
    end
%     figure;hold all;
% for kk=1:size(EoutAux1,3)
%     plot(f(:,1),20*log10(abs(fftshift(fft(EoutAux1(:,1,VetThisCarr==kk)./length(EoutAux1(:,1,VetThisCarr==kk)))))));
% end
% a=1;
    for ThisCarr=InitCarr:CarrPass:NumCarr
        EoutAux = EoutAux1(:,:,VetThisCarr==ObsCarrUsed(ThisCarr));
        switch Medium
            case 'Fiber'
                EoutAux = ifft(fft(EoutAux).*exp(1j*2*pi*f*(...
                    FiberDelay(ThisCarr)*Ta)));
            otherwise
        end
%         sn=((1/2)*(((1/2)^MaxNumStag)-1))/((1/2)-1);
%         DiffPos = round((sn/(2/2^IfftOrSum))*T/Ta);
%         EoutAux = [EoutAux(DiffPos+1:end) EoutAux(1:DiffPos)];  
        %The current incoming signal is them converted from the optical
        %domain to the eletrical domain with the help of an photo
        %detector.
        Ix =EoutAux.*conj(EoutAux);
%         Ix = (Ix.*exp(-1j.*(2*pi*FeqC.*t))) + 1j.*(Ix.*exp(-1j.*(2*pi*FeqC.*t + pi/2)));
        [BitFilt,~] = FiltroGaussiano(f,fc/2,0,1);                           %Creating filter for selection of the received signal
        BitFilt = fftshift(BitFilt);                                       %Shifting the filter for matching the received signal
        Ix = ifft(fft(Ix).*BitFilt);                                     %Filtering the signal for better detection
        Ix = Ix - min(min(Ix));                                                 %Removing any off-set that may exist
        Ix = Ix.*ConNor;  
        NorAux = max(Ix);
        NorAux = repmat(NorAux,size(Ix,1),1);
        Ix = (Ix./NorAux);                                           %Normalizing the eletrical signal (amplifying what is needed)
        
        if ReceptorNoise                                               %Verify whether the noise is to be added or not
            [~,SigPower] = MeasPower(Ix);
            SigPower2 = SigPower-30;%20*log10(SigPowerI);
%             SNR2 = CarSNR + 10*log10(1) - 10*log10(OvSam) + 10*log10(10^0.36);
%             if TimeSys
%                 SNR = CarSNR + 10*log10(mean(DmtMve)) - 10*log10(1*Te*fc);
%             else
%                 SNR = CarSNR + 10*log10(mean(DmtMve)) - 10*log10(1*Ts*fc);
%             end
            Ix = awgn(Ix,SNR,SigPower2);
        end
        
        %Ix = reshape(Ix,Nb*NPPB,CurTesSiz);
        Ix = Ix(1+NPSTUf:end,:);
        Ix = reshape(Ix,OvSam*NFFT,NumFra*CurTesSiz);
        Ix = intdump(Ix,NPPOF);                            %Downsampling the income signal
        SigRecepA = Ix(1+NPOFEX:end,:);
        if UseIQ
            SigReal = SigRecepA.*cos(2*pi*FcIQ.*tta);
            SigImag = SigRecepA.*sin(2*pi*FcIQ.*tta);
            
            if SelModFilt
                SigReal = ifft(fft(SigReal).*ModFilta);
                SigImag = ifft(fft(SigImag).*ModFilta);
            end
            
            SigRecepA = SigReal + 1j.*SigImag;
        end
%         IxAux = Ix(NuSaTg+1:end);
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
%                         x_ofdm = x_ofdm_ssb;
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
%         SigRecep = reshape(SigRecep,NumFra,NFFT);
        SigRecep1 = fft(SigRecep,NFFT);
%         SigRecep2 = SigRecep1(2+Zt-ZtC:Ns+Zt-ZtC+1,:);
        if UsingHermitian
            SigRecep2 = SigRecep1(2+Zt-ZtC:Ns+1+Zt-ZtC,:);                               %Taking the actual data
        else
        %     SigRecep2 = SigRecep1(1+(NFFT-size(TxSymbH,1))/2:size(TxSymbH,1)+(NFFT-size(TxSymbH,1))/2,:);
            SigRecep2 = [SigRecep1(1:size(TxSymbH,1)/2,:);SigRecep1(end+1-size(TxSymbH,1)/2:end,:)];
        end
%         SigRecep2 = SigRecep2(:).';
        switch OfdMod
            case 'qam'
                TxDataA = TxDataMat(ThisCarr,:);
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
                TxDataA = TxDataMat(ThisCarr,:);
                TxDataA = reshape(TxDataA,NumFra*CurTesSiz,length(TxDataA)/NumFra*CurTesSiz);
                TxDataA = TxDataA.';
                TxSigToPlot  = dpskmod(TxDataMat(ThisCarr,:),M);                           %Modulating information
        end
        equ        = (SigRecep2)./(TxSigToPlot);                               % Coeficientes do Equalizador
        
        figure;
        hold all;
        plot(TxSigToPlot,'*','LineWidth',2);
        plot(SigRecep2,'o','LineWidth',2);
        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
          
        ChanelEqualizer(ObsCarrUsed(ThisCarr),:)        = equ(:);

    end
    clear EoutModTem EoutModAux EoutMod;
    %% Transmiting UpStream
    for kk=InitCarrUp:CarrPass:NumCarr                                     %Generating different data for each carrier
%         CpPaFr = 1;
		UniDmtMve = unique(DmtMve);
		UniDmtMve = fliplr(UniDmtMve);
        DaSiAn    = 0;
        DaSiPo    = 0;
		for CarDmtM=1:length(UniDmtMve)
			M = UniDmtMve(CarDmtM);
			DaSiAu = sum(ismember(DmtMve,M));
			TxData = randi([0 M-1],NumFra*CurTesSiz,DaSiAu);%Generation random information
			TxDataMat(kk,1+DaSiAn:DaSiAn + CurTesSiz*DaSiAu*NumFra) = TxData(:);%Saving data for future comparison
            DaSiAn = DaSiAn + CurTesSiz*DaSiAu*NumFra;
			switch OfdMod                                              %Sellect which modulation was used
				case 'qam'
					TxSymbAux(1:NumFra*CurTesSiz,1+DaSiPo:DaSiPo + DaSiAu)  = qammod(TxData,M);                    %Modulating information by QAM
				otherwise
					TxSymbAux(1:NumFra*CurTesSiz,1+DaSiPo:DaSiPo + DaSiAu)  = dpskmod(TxData,M);                   %Modulating information DPSK as default
			end
            DaSiPo = DaSiPo + DaSiAu;
		end
        TxSymb      = TxSymbAux.';
        
        
        TxSymbConj      = conj(flipud(TxSymb));
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
        TxSymbModA       = ifft(TxSymbC,NFFT);
        
        TxSymbModA  = rectpulse(TxSymbModA,OvSam);         %Over sampling

        if UseIQ
            RealTx = real(TxSymbModA);
            ImagTx = imag(TxSymbModA);
            TxSymbModA = RealTx.*cos(2*pi*FcIQ.*tta) + ImagTx.*sin(2*pi*FcIQ.*tta);
        end
 
        switch SelModTp
            case 'BP'                                              %Base Band transmission
%                 TxSymbModA  = rectpulse(TxSymbModA,OvSam);         %Over sampling
                TxSymbMod = TxSymbModA;
                if SelModFilt
                    TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
                end
            case 'AM'                                              %Out of base band
%                 TxSymbModA  = rectpulse(TxSymbModA,OvSam);         %Over sampling
                TxSymbMod   = TxSymbModA.*cos(2*pi*Ofc*tta);
                if SelModFilt
                    TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
                end
            case 'AMSSB'                                           %Out of base band
%                 TxSymbModA  = rectpulse(TxSymbModA,OvSam);         %Over sampling
                TxSymbMod   = TxSymbModA.*exp(1j*2*pi*Ofc*tta);
                if SelModFilt
                    TxSymbMod = ifft(fft(TxSymbMod).*ModFilta);
                end
            otherwise
%                 TxSymbModA  = rectpulse(TxSymbModA,OvSam);         %Over sampling
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
 
        %AuxTxSymbMod(1,:)  = TxSymbMod(:);
        TxSymbModA  = [TxSymbMod(1+end-NPOFEX:end,:);TxSymbMod];
        U1t  = rectpulse(TxSymbModA,NPPOF);%Over sampling
		TxSigA = reshape(U1t,OvSam*NumFra*NFFT,CurTesSiz);
		U1t = [TxSigA(1+end-NPSTUf:end,:);TxSigA];                         %Conforming vector sizes to the same length
%         U1t = real(U1t.*exp(1j.*(2*pi*FeqC.*t))) + imag(U1t.*exp(1j.*(2*pi*FeqC.*t + pi/2)));
%Assigning the eletrical signal to one drive of the MZM - The aplitude of 
%the signal can be controlled by the variable DatGai, which can be 
%understood as an gain for the eletrical signal or an atenuation. The 
%second signal will be similar with the only difference a phase shift of pi
        ConNor = max(U1t)./0.95;
        ConNor = repmat(ConNor,size(U1t,1),1);
        U1t = (U1t./ConNor);
        U.U1t = U1t;
        U.U2t = exp(-1j*pi).*U1t;
%As both signals will have the mostrly the same characteristics with the 
%only difference the phase shift of 180 degress. The MZM-I will be working 
%on the Push-Pull configuration. It is necessary to reduce the Chirp noise
%to zero.
		EoutTxAux = EoutAux1(:,:,VetThisCarr==(RefCarr-1+kk));
		[EoutModAux,~]=MZM(freqGHz,EoutTxAux,U,L,U0,U_pi1,U_pi2,nopt,nel,C,alfa0);%Modulating individual carriers
        EoutModTem(1:length(EoutModAux),:,ObsCarrPos(kk)) = EoutModAux;      %Adding the current Modulated output to the final OutPut
%             if kk==1
%                 figure;hold all;grid on;
%                 plot(f(:,1),20*log10(abs(fftshift(fft(EoutModAux)./length(EoutModAux)))));
%                 set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%                 figure;hold all;grid on;
%                 plot(f(:,1),20*log10(abs(fftshift(fft(EoutAux1(:,:,VetThisCarr==(RefCarr-1+kk)))./length(EoutAux1(:,:,VetThisCarr==(RefCarr-1+kk)))))));
%                 set(gcf,'units','normalized','outerposition',[0 0 1 1]);
% %             end
% a=1;
    end
    if IfftOrSum
        if size(EoutModTem,3)>1
			if UsingGpu==1
                Tgpu = gpuArray(T);
                StagGpu = gpuArray(MaxNumStagT);
                EoutMod = zeros(size(EoutRec,1),size(EoutRec,2));
                for jj=1:size(EoutRec,2)
				EoutModTemGpu = gpuArray(EoutModTem(:,jj,:));
				[EoutAux1Gpu,~,~]=OpticalIFFTGT(fgpu,Tgpu,StagGpu,EoutModTemGpu);
				EoutMod(:,jj) = gather(EoutAux1Gpu);
				clear EoutModTemGpu EoutAux1Gpu;
                end
				clear Tgpu StagGpu;
			else
				[EoutMod,~,~] = OpticalIFFTNT(f,T,MaxNumStagT,EoutModTem);
			end
            %[EoutMod,~,~] = OpticalIFFTN(f,T,MaxNumStagT,EoutModTem);
        else
            EoutMod = EoutModTem;
        end
    else
        if size(EoutModTem,3)>1
            EoutMod = sum(EoutModTem,3);
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
    if FFTSplit
		if UsingGpu==1
            Tgpu = gpuArray(T);
            StagGpu = gpuArray(MaxNumStag);
            EoutAux1 = zeros(size(EoutRec,1),size(EoutRec,2),NumCarr);
            for jj=1:size(EoutRec,2)
                EoutGpu = gpuArray(EoutRec(:,jj));
                [EoutAux1Gpu,~,VetThisCarrGpu]=OpticalFFTGT(fgpu,Tgpu,StagGpu,EoutGpu);
                EoutAux1(:,jj,:) = gather(EoutAux1Gpu);
                VetThisCarr = gather(VetThisCarrGpu);
                clear EoutGpu EoutAux1Gpu;
            end
			clear StagGpu Tgpu;
		else
			[EoutAux1,~,VetThisCarr]=OpticalFFTNT(f,T,MaxNumStag,EoutRec);
		end
        %[EoutAux1,~,VetThisCarr]=OpticalFFTN(f,T,MaxNumStag,EoutRec);
    else
        VetThisCarr = ObsCarrPos;
        EoutAux1    = SelectEachCarrier(EoutRec,NumCarr,f,fin,1.0*fc,11,fc);
    end
    InitCarr = InitCarrUp;
    for ThisCarr=InitCarr:CarrPass:NumCarr
        EoutAux = EoutAux1(:,:,VetThisCarr==ObsCarrUsed(ThisCarr));
        switch Medium
            case 'Fiber'
                EoutAux = ifft(fft(EoutAux).*exp(1j*2*pi*f*(...
                    FiberDelay(ThisCarr)*Ta)));
            otherwise
        end
%         sn=((1/2)*(((1/2)^MaxNumStag)-1))/((1/2)-1);
%         DiffPos = round((sn/(2/2^IfftOrSum))*T/Ta);
%         EoutAux = [EoutAux(DiffPos+1:end) EoutAux(1:DiffPos)];  
        %The current incoming signal is them converted from the optical
        %domain to the eletrical domain with the help of an photo
        %detector.
        Ix =EoutAux.*conj(EoutAux);
%         Ix = (Ix.*exp(-1j.*(2*pi*FeqC.*t))) + 1j.*(Ix.*exp(-1j.*(2*pi*FeqC.*t + pi/2)));
        [BitFilt,~] = FiltroGaussiano(f,fc/2,0,1);                           %Creating filter for selection of the received signal
        BitFilt = fftshift(BitFilt);                                       %Shifting the filter for matching the received signal
        Ix = ifft(fft(Ix).*BitFilt);                                     %Filtering the signal for better detection
        Ix = Ix - min(min(Ix));                                                 %Removing any off-set that may exist
        Ix = Ix.*ConNor;    
        NorAux = max(Ix);
        NorAux = repmat(NorAux,size(Ix,1),1);
        Ix = (Ix./NorAux);                                         %Normalizing the eletrical signal (amplifying what is needed)
        
        if ReceptorNoise                                               %Verify whether the noise is to be added or not
            [~,SigPower] = MeasPower(Ix);
            SigPower2 = SigPower-30;%20*log10(SigPowerI);
            Ix = awgn(Ix,SNR,SigPower2);
        end
        
        %Ix = reshape(Ix,Nb*NPPB,CurTesSiz);
        Ix = Ix(1+NPSTUf:end,:);
        Ix = reshape(Ix,OvSam*NFFT,NumFra*CurTesSiz);
        Ix = intdump(Ix,NPPOF);                            %Downsampling the income signal
        SigRecepA = Ix(1+NPOFEX:end,:);
        
        if UseIQ
            SigReal = SigRecepA.*cos(2*pi*FcIQ.*tta);
            SigImag = SigRecepA.*sin(2*pi*FcIQ.*tta);
            
            if SelModFilt
                SigReal = ifft(fft(SigReal).*ModFilta);
                SigImag = ifft(fft(SigImag).*ModFilta);
            end
            
            SigRecepA = SigReal + 1j.*SigImag;
        end
        
%         IxAux = Ix(NuSaTg+1:end);
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
%                         x_ofdm = x_ofdm_ssb;
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
%         SigRecep = reshape(SigRecep,NumFra,NFFT);
        SigRecep1 = fft(SigRecep,NFFT);
%         SigRecep2 = SigRecep1(2+Zt-ZtC:Ns+Zt-ZtC+1,:);
        if UsingHermitian
            SigRecep2 = SigRecep1(2+Zt-ZtC:Ns+1+Zt-ZtC,:);                               %Taking the actual data
        else
        %     SigRecep2 = SigRecep1(1+(NFFT-size(TxSymbH,1))/2:size(TxSymbH,1)+(NFFT-size(TxSymbH,1))/2,:);
            SigRecep2 = [SigRecep1(1:size(TxSymbH,1)/2,:);SigRecep1(end+1-size(TxSymbH,1)/2:end,:)];
        end
%         SigRecep2 = SigRecep2(:).';
        switch OfdMod
            case 'qam'
                TxDataA = TxDataMat(ThisCarr,:);
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
                TxDataA = TxDataMat(ThisCarr,:);
                TxDataA = reshape(TxDataA,NumFra*CurTesSiz,length(TxDataA)/NumFra*CurTesSiz);
                TxDataA = TxDataA.';
%                 TxSig = reshape(TxDataMat(ThisCarr,:),NumFra,NFFT);
                TxSigToPlot  = dpskmod(TxDataMat(ThisCarr,:),M);                           %Modulating information
%                 TxSigToPlot  = TxSigToPlot(:);
        end
        equ        = (SigRecep2)./(TxSigToPlot);                               % Coeficientes do Equalizador
        
        figure;
        hold all;
        plot(TxSigToPlot,'*','LineWidth',2);
        plot(SigRecep2,'o','LineWidth',2);
        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
          
        ChanelEqualizer(ObsCarrUsed(ThisCarr),:)        = equ(:);
    end
end

