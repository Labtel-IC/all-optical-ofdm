%%
FeqC = 0;
SelModFilt= 0;
SelecGaus = 1;
OfdMod    = 'qam';                                                %Chossing the type of modulation for the electrical subcarriers
SelModTp  = 'AMSSB';                                                  %Selecting with OFDM will be in Base Band "BP", Amplitude Modulation
%with double side band "AM" or with single side band "AMSSB" when a
%different option is choosen this script will implement phase modulation.
UsingHermitian = 0;
% if ~UseOfdmElect
%     UsingHermitian = 1;
% end
ForcOvSam = 0;
UsingParCar= 0;
NuDmSuDi  = 1;
DmtPas    = 2;
DmtMinM   = 3;
TapN      = 1;                                                     %Number of equalizer used
ExtSam    = 1;                                                     %Set to the number of frames that will be used.
NpEl      = 2^9;                                                  %Number of Electrical carriers
NFFT      = NpEl/1;                                                %Size of the electical FFT
OveT      = 2^9;                                                   %Variable to set the eletrical carrier period to be OveT time the symbol period
M         = 4;                                                     %Modulation Level
% Not all electrical carriers may be used. For those unused
%carriers, zeros are placed in their locations. As a result, there
%are free spaces to place the used carriers as the Zt variable
%control how many zeros will be inserted on the FFT frame to offset
%the used carriers from the beginning. The value of Zt is
%approximately half of the number of the free carrier. Thus,
%ControlZt variable controls the percentage that Zt will influence
%the FFT offset. As a result, more flexibility was possible for
%placing the used carriers on the available spectrum.
ControlZt = 1;
CpLe   = 1;                                                      %Set the percentage of electrical carriers that will be used
if CpLe <= 0.5
    error(['CpLe can not be bellow 50%. More than half of ' ...
        'available carriers must be transmited']);
end
BW     = fc;                                                       %Signal available bandwidth
%         NumFra = floor(BW*t(end));                                       %Number of frames
%         if NFFT>Nb
%             NFFT = Nb;
%         end
DatSiz = NFFT;                                                     %Set to the total number of electrical carriers to be used
Te     = OveT*T;
Tu     = Te;                                                       %Time that will actually be used to transmit information
Tg     = 0.0*Tu;                                                     %Time for the guard band
if ((NFFT*NpEl)==Nb)&&(Tg~=0)
    error('NFFT must be smaller than Nb to use some ciclic prefix');
end
Ts     = Tu + Tg;                                                  %Time of the OFDM symbol
g      = Tu/Ts;                                                    %Percentage of band guard
Dtf    = 1/Tu;                                                     %Chanal signaling ratio
if UsingHermitian
    N      = NFFT/2 - 1;                                           %Number of carriers found
    Ns     = ceil((CpLe*(DatSiz))/2) - 1;                          %Number of carriers used
else
    N      = NFFT;                                                 %Number of carriers found
    Ns     = ceil((CpLe*(DatSiz)));                                %Number of carriers used
    if mod(Ns,2)
        Ns = Ns - 1;
    end
end
DifNsN = ceil((N - Ns)/1);                                         %Number of carriers unused

k      = log2(M);                                                  %Number of bits per symbol
MZ_Input_File = 2;                                                 %Inputdata file to the MZM
NuSaTs = Nb*NPPB;
NuSaTu = ceil(g*NPPB);
NuSaTg = NuSaTs - NuSaTu*Nb;
OBw    = NFFT/Tu;                                                  %OFDM bandwidth
Ofc    = 1.0*OBw;                                                  %Central frequency of out of base band transmission
Ofs    = 2*Ofc;                                                    %Electrical sampling frequency
OvSam = round(Te/t(2)/NFFT);                                       %Number of over samples for out of base band modulations
if UsingHermitian
    NuPart = floor(Ns/NuDmSuDi);
    NuSuSe = floor(Ns/NuPart);
    DmtMve = ones(1,Ns);
    for ElCaM=NuSuSe:-1:2
        DmtMve(1+(NuSuSe-ElCaM)*NuPart:(NuSuSe-ElCaM+1)*NuPart)= 2^(ElCaM);%2^(ElCaM*DmtPas-1);
    end
else
    NuPart = floor(Ns/NuDmSuDi);
    NuSuSe = floor(Ns/NuPart);
    DmtMve = ones(1,Ns/2);
    DmtMve(:) = 2^DmtMinM;
    for ElCaM=NuSuSe:-1:2
        DmtMve(1+(NuSuSe-ElCaM)*round(NuPart/2):(NuSuSe-ElCaM+1)*round(NuPart/2))= 2^(ElCaM);%2^(ElCaM*DmtPas-1);
    end
    DmtMve = [DmtMve fliplr(DmtMve)];
%     figure;hold on;plot(1:Ns,DmtMve)
end
switch OfdMod
    case 'dpsk'
        DmtMve(:) = 2^DmtMinM;
    otherwise
        if NuDmSuDi==1
            DmtMve(:) = 2^DmtMinM;
        else
            DmtMve(1+(NuSuSe-ElCaM+1)*NuPart:end)= 2^DmtMinM;
        end
end

if UsingHermitian
    nRb    = (OBw/(2*(N+2)))*(Ns*mean(log2(DmtMve)));                                       %New Rb found
else
    nRb    = (OBw/(N))*(Ns*mean(log2(DmtMve)));                                       %New Rb found
end

nRbAgg = nRb*NumCarr/2;                                            %New Aggregate Symbol ratio found
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
    if DifNsN < 1
        Zt = 0;
    else
        Zt  = floor(DifNsN/2);
    end
end
ZtC = floor(Zt*ControlZt);
switch SelModTp
    case 'AM'
        NumFraPar = 1;
        FramSpac = 1;
        ForcOvSam = 0;
    case 'AMSSB'
        NumFraPar = 1;
        FramSpac = 1;
        ForcOvSam = 0;
    otherwise
        OBw    = NFFT/Tu;                                          %OFDM bandwidth
        Ofc    = 0*OBw;                                            %Central frequency of out of base band transmission
        Ofs    = 2*Ofc;                                            %Electrical sampling frequency
        nRb    = (OBw/(N+2))*(Ns*k);                               %New Rb found
        FramSpac = 2*OBw;
        if UsingParCar
            NumFraPar = floor(Rb/((OBw+FramSpac)*1));
            %                 NumFraPar = 3;
            if (NumFraPar>1)&&(~mod(NumFraPar,2))
                NumFraPar = NumFraPar - 1;
            end
            if NumFraPar>0
                %                     OvSam  = round(Te*2*OBw*2^((NumFraPar-1)/2));                   %Number of over samples for out of base band modulations
            else
                %                     OvSam  = round(Te*2*OBw);
                NumFraPar = 1;
            end
        else
            NumFraPar = 1;
        end
        OvSam = round(Te/t(2)/NFFT);
end
if ForcOvSam
    OvSam = ForcOvSam;
end
NumFra    = floor((Nb*NPPB)/((OvSam*NFFT)+(OvSam*NFFT*Tg/Tu)));
NumFraTem = NumFra/NumFraPar;               %Number of frames
%         if NumFra>1
%             if mod(NumFra,2)
%                 NumFra = NumFra - 1;
%             end
%         end
NumFrU = floor((Nb*NPPB)/(OvSam*NFFT));                            %Number of frames
NuAmOf = (OvSam*NumFra*NFFT);                                      %Total number of samples of the electrical OFDM signal
NPOFEX = (NuAmOf/NumFra)*(Tg/Tu);                                  %Extra samples per OFDM frame. It is the number of CP
NuAmTo = NuAmOf*(1+Tg/Tu);                                         %Total number of samples of the electrical OFDM signal added CP
NPPOF = floor((Nb*NPPB)/NuAmTo);                                   %Over sampling for optical modulation transmission
if (NuAmOf<=0)||(NPPOF<=0)
    error('OFDM paramters outof range. Please double check them.');
end
if NumFra/OBw > Nb*T
    error(['The simulation is exiding the maximum allowed time '...
        '.Please refrain the amount of serial frames used.']);
end
NPSTUf = Nb*NPPB-NPPOF*NuAmTo;                                     %Number of extra sample to make all vector to have the same length
NumFrE = NumFrU - NumFra;                                          %Number of frames

% SNR = CarSNR + 10*log10(log2(max(DmtMve))) - 10*log10(1*Tu*OBw);
SNR = CarSNR + 10*log10(log2(max(DmtMve))) - 10*log10(OvSam);

%         NPOFEX = floor(NumFrE*OvSam*NFFT/NumFra);
%         NPPOF  = floor((2^19)/((NPOFEX*NumFra)+(OvSam*NumFra*NFFT)));
%         NPSTUf = 2^19 - NPPOF*NumFra*OvSam*NFFT;
%Variables for the Selection of carriers - Wave Divided.
fin           = RefCarr*fc;                                        %The frequency of the first carrier of the optical COMB
FBWD          = 1.0*fc;                                            %The filter BandWidth to split each individual carrier
Order         = 5;                                                 %The order of the filter used to split each carrier
a=1;