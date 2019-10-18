close all;clear;clc;
%%
UseHerSym = true;
UsePasBan = true;
useSiSiBa = true;
UseExp = false;
UseIq = false;
Bw = 10e6;
NFFT = 512;
if UseHerSym
    Ns = (NFFT/2) - 1;
else
    Ns = NFFT-4;
end
% Ns = 64;
Fu = Bw/512;
Tu = 1/Fu;
NumFra = 1;
if UsePasBan
%     Fc = 5e6;
    Fc = Bw;
else
    Fc = Bw/2;
end
Fa = Fc + 7*Bw;
Ta = 1/Fa;
tf = NumFra*Tu;
tn = linspace(0,tf,NFFT*NumFra);
fn = time2freq(tn);
[ FiltroN ] = Filtro_Retangular( Bw,Fc,fn);
FiltroN = fftshift(FiltroN);
TotalSampPerSymb = Tu/Ta;
% TotalSampPerSymb = NFFT;
OveSamp = ceil(TotalSampPerSymb/NFFT);
TotalSampPerSymb = OveSamp*NFFT;
log2(TotalSampPerSymb)
TotalSimuSamp = TotalSampPerSymb*NumFra;
log2(TotalSimuSamp)
t = linspace(0,tf,TotalSimuSamp);
f = time2freq(t);
[ Filtro ] = Filtro_Retangular( Bw,0,f);
Filtro = fftshift(Filtro);
[ FiltroRx ] = Filtro_Retangular( Bw,0,f);
FiltroRx = fftshift(FiltroRx);
t_u = linspace(0,Tu,TotalSampPerSymb);
f_u = time2freq(t_u);
M = 4;
k = log2(M);
Mult = 1;
%%

Data = randi([0 M-1],Ns,NumFra);
% Data = zeros(Ns,NumFra);
% Data(end/4) = M-1;

DataMod = qammod(Data,M);

tx = zeros(NFFT,NumFra);
if UseHerSym
    DataConj = conj(DataMod);
    DataFlip = flipud(DataConj);
    tx(2:1+size(DataMod,1),:) = DataMod;
    tx(2+end/2:1+(end/2) + size(DataFlip,1),:) = DataFlip;
else
    tx(2:1+size(DataMod,1)/2,:) = DataMod(1:end/2,:);
    tx(2+end/2:1+(end/2) + size(DataMod,1)/2,:) = DataMod(1+end/2:end,:);
%     tx(2:1+size(DataMod,1),:) = DataMod;
end
Tx = ifft(tx,NFFT);
if OveSamp > 1
Tx_Over = rectpulse(Tx,OveSamp);
else
    Tx_Over = Tx;
end
Tx_A = reshape(Tx_Over,1,OveSamp*NFFT*NumFra);
Tx_B = real(Tx_A);
Tx_C = imag(Tx_A);
if UsePasBan
    if useSiSiBa
        if UseExp
            Tx_Sig = Tx_B.*exp(1j.*(2*pi*Fc.*t)) + Tx_C.*exp(1j.*(2*pi*Fc.*t + pi/2));
        else
            Tx_Sig = Tx_A.*exp(1j.*(2*pi*Fc.*t));
        end
    else
        if UseExp
            Tx_Sig = Tx_B.*cos(2*pi*Fc*t) + Tx_C.*sin(2*pi*Fc*t);
        else
            Tx_Sig = Tx_A.*cos(2*pi*Fc*t);
        end
    end
else
    if UseIq
        Tx_Sig = Tx_B.*cos(2*pi*Fc.*t) + Tx_C.*sin(2*pi*Fc.*t);
        Tx_Sig = ifft(fft(Tx_Sig).*Filtro);
    else
        Tx_Sig = Tx_A;
    end
end
figure;
hold all;
TxPlot = fft(Tx);
plot(fn,10*log10(abs(fftshift(TxPlot(:)./length(TxPlot(:))))));
plot(fn,10*log10(abs(fftshift(FiltroN))));
set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
figure
hold all;

plot(f,10*log10(abs(fftshift(fft(Tx_Sig)./length(Tx_Sig)))));

plot(f_u,10*log10(abs(fftshift(fft(Tx_Sig(1+TotalSampPerSymb*(Mult-1):Mult*TotalSampPerSymb)./length(Tx_Sig(1+TotalSampPerSymb*(Mult-1):Mult*TotalSampPerSymb)))))))
plot(f,10*log10(abs(fftshift(Filtro))));

set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
a=1;
figure
hold all;
plot(t,abs(Tx_Sig))
set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
a=1;

Tx_Sig = Tx_Sig.*exp(1j.*(2*pi*Fc*t));
figure;
hold all;
RxPlot = fft(Tx_Sig);
plot(f,10*log10(abs(fftshift(RxPlot(:)./length(RxPlot(:))))));
set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
a=1;
%%
Rx_Sig = Tx_Sig;

%%
Rx_Sig = Rx_Sig.*conj(Rx_Sig);
figure;
hold all;
RxPlot = fft(Rx_Sig);
plot(f,10*log10(abs(fftshift(RxPlot(:)./length(RxPlot(:))))));
set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
a=1;
if UsePasBan
    if useSiSiBa
        if UseExp
            Rx_Real = Rx_Sig;
            Rx_Imag = Rx_Sig;
            Rx_B = Rx_Sig.*exp(-1j.*(2*pi*Fc.*t));
            Rx_C = Rx_Sig.*exp(-1j.*(2*pi*Fc.*t + pi/2));
            Rx_A = Rx_B + 1j.*Rx_C;
        else
            Rx_Real = Rx_Sig;
            Rx_Imag = Rx_Sig;
            Rx_A = Rx_Sig;%.*exp(-1j.*(2*pi*Fc.*t));
        end
    else
        if UseExp
            Rx_Real = Rx_Sig;
            Rx_Imag = Rx_Sig;
            Rx_B = Rx_Sig.*cos(2*pi*Fc*t);
            Rx_C = Rx_Sig.*sin(2*pi*Fc*t);
            Rx_A = Rx_B + 1j.*Rx_C;
        else
            Rx_Real = Rx_Sig;
            Rx_Imag = Rx_Sig;
            Rx_A = Rx_Sig.*cos(2*pi*Fc*t);
        end
    end
else
    if UseIq
        Rx_Real = Rx_Sig.*cos(2*pi*Fc*t);
%         Rx_Real = ifft(fft(Rx_Real).*FiltroRx);
        Rx_Imag = Rx_Sig.*sin(2*pi*Fc*t);
%         Rx_Imag = ifft(fft(Rx_Imag).*FiltroRx);
        Rx_A = Rx_Real + 1j.*Rx_Imag;
%         Rx_A = Rx_A.*exp(-1j.*(2*pi*Fc*t+pi));
%         Rx_A = fftshift(Rx_A);

    else
        Rx_Real = Rx_Sig;
        Rx_Imag = Rx_Sig;
        Rx_A = Rx_Sig;
    end
end
figure;
hold all;
RxPlot = fft(Rx_Real);
plot(f,10*log10(abs(fftshift(RxPlot(:)./length(RxPlot(:))))));
RxPlot = fft(Rx_Imag);
plot(f,10*log10(abs(fftshift(RxPlot(:)./length(RxPlot(:))))));
plot(f,10*log10(abs(fftshift(FiltroRx))));
% set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
% a=1;
% figure;
% hold all;
RxPlot = fft(Rx_A);
plot(f,10*log10(abs(fftshift(RxPlot(:)./length(RxPlot(:))))));
plot(f,10*log10(abs(fftshift(FiltroRx))));
set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
a=1;
Rx_Over = reshape(Rx_A,OveSamp*NFFT,NumFra);

if OveSamp > 1
Rx = intdump(Rx_Over,OveSamp);
else
    Rx = Rx_Over;
end
figure;
hold all;
RxPlot = fft(Rx);
plot(fn,10*log10(abs(fftshift(RxPlot(:)./length(RxPlot(:))))));
plot(fn,10*log10(abs(fftshift(FiltroN))));
set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
a=1;

rx = fft(Rx,NFFT);
if UseHerSym
    rx_mod = rx(2:1+Ns,:);
else
    rx_mod = [rx(2:1+Ns/2,:);rx(2+end/2:1+end/2+Ns/2,:)];
end
rx_data = qamdemod(rx_mod,M);

figure;
hold all;
posaux = zeros(1,size(rx_mod,1)*size(rx_mod,2));
plot(real(rx_mod(:)),imag(rx_mod(:)),'o','color',[1 0.4 0]);
plot(real(DataMod(:)),imag(DataMod(:)),'*','color',[0 0 0.6]);
set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);

figure;
hold all;
plot(Data(:));
plot(rx_data(:));
set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
a=1;

CompSymb = ~(Data(:)==rx_data(:));
SymbErr = sum(CompSymb);
BitErr = SymbErr*k;
BER = BitErr/(size(Data,1)*size(Data,2)*k);

StrToPrint = sprintf('BER = %0.6f',BER);

display(StrToPrint);

a=1;












