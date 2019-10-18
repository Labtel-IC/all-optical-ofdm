function [Delay,PowRef] = PARProgationDelay(NumOfCarrier,temp,lambda,To,L,...
                              regime,f,Ploting,FiberInpuData,FiberFunction)
    %%
    if nargin < 7
        error('Not Enough Input Argument');
    elseif nargin<8
        Ploting = 0;
        FiberInpuData = 'Fibra_Monomodo_Input_Data';
        FiberFunction = 'Fibra_Monomodo1';
    elseif nargin<9
        FiberInpuData = 'Fibra_Monomodo_Input_Data';
        FiberFunction = 'Fibra_Monomodo1';
    elseif nargin<10
        FiberFunction = 'Fibra_Monomodo1';
    end
    %%
    BWD      = 1/To;
    CenFeq   = 0;
    FiltOrd  = 1;
    [Filt,~] = FiltroGaussiano(f,BWD,CenFeq,FiltOrd);                          %Creating filter for selection of the received signal
    Filt     = fftshift(Filt);                                                 %Shifting the filter for matching the received signal 

    NPPB = round(To/(temp(2)-temp(1)));
    Nb = round(temp(end)/To);

    VetPul = zeros(1,Nb);
    VetPul(round(end/2)) = 1;

    TxPul = rectpulse(VetPul,NPPB);


    VetDif = length(TxPul)-length(temp);

    if VetDif > 0 
        TxPul = [TxPul zeros(1,VetDif)];
    elseif VetDif < 0
        TxPul = TxPul(1:end-VetDif);
    end
    Delay  = zeros(1,NumOfCarrier);
    PowRef = zeros(1,NumOfCarrier);
    parfor kk=1:NumOfCarrier
        TxPul2 = ifft(fft(TxPul).*Filt);
        PiloPuls = TxPul2.*exp(1j*2*pi*(kk/To).*temp);
        %%
        fh = str2func(FiberFunction);

        [EoutRec,~,~]=fh(temp,PiloPuls,lambda,To,L,regime,f,FiberInpuData);

        Ix =EoutRec.*conj(EoutRec);
        Ix = ifft(fft(Ix).*Filt);
        Ix = Ix - min(Ix);                                                         %Removing any off-set that may exist 
        Ix = Ix./max(abs(Ix));                                         %Normalizing the eletrical signal (amplifying what is needed)
        %% Ploting the result for qualitative analizes
%         if Ploting
%             figure;hold all;
%             plot(linspace(1,length(temp),length(temp)),abs(PiloPuls)./max(abs(PiloPuls)),linspace(1,length(temp),length(temp)),Ix);
%             title(['Comparison between transmited and received ' ...
%                           'signal'],'FontSize',16,'FontWeight','bold');
%             xlabel('Time [s]','FontSize',14);%,'FontWeight','bold');
%             ylabel('Amplitude [u]','FontSize',14);%,'FontWeight','bold');
%             legend({'Tx','Rx'},'FontSize',12,'Location','best',...
%                                       'FontWeight','bold','Box','off');
%             grid on;
%         end
        %%
        [~,SyncedPos] = findpeaks(abs(PiloPuls),'SortStr','descend','NPeaks',1);
        [~,PosToSync] = findpeaks(Ix,'SortStr','descend','NPeaks',1);

        Delay(kk)  = PosToSync - SyncedPos;
        PowRef(kk) = sum(abs(Ix).^2)./length(Ix);
        
%         a=1;
    end
end