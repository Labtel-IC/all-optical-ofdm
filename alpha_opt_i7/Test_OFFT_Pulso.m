close all;clear;clc;
thisone = 2;
if thisone==1    
    fc           = 12.5e9;%2^30;%center frequency of Electical signal ^10 kilo ^20 mega ^30 giga
    Rb           = 12.5e9;
%     f_max        = (2^0)*;                                                      %max frequency that matlab can reproduce
    NPPB         = 2^5;
    fsample      = NPPB*Rb;                                                   %Sampling frequency
    T            = 1/fc;                                                       %period of the filter
    Ta           = 1/fsample;                                                       %period of the filter
    Rb           = fc;                                                     %bit rate
    Tb           = 1/Rb;                                                       %time of one bit NRZ
    NPPB_This    = 2^0;                                            %Number of Point Per Bit
    NumberOf_T   = 2^13;                                                        %number of periods on the simulation 2^18 is the m�ximum value for numberEyes = 2
    FinalTime    = NumberOf_T*T;                                               %Final time of simulation
    Nb           = floor(FinalTime/Tb);                                               %number of bits
    TotalSamples = NPPB_This*NumberOf_T*(T/Ta);                                                   %Total Number of Samples
%     NPPB         = NPPB_This*(T/Ta);

    t            = linspace(0,FinalTime,round(TotalSamples));                         %Time vector for simulation
    f            = time2freq(t); 
    df = f(2)-f(1);
    aux1 = round(fc/df);
    spc = aux1*df
elseif thisone==2
    Nb = 2^13;
    Rb = 12.5e9;
    Tb = 1/Rb;
    NPPB = 2^5;
    NumberOf_T = Nb;
    tf = Nb*Tb;
    TotalSamples = NPPB*Nb;
    t = linspace(0,tf,TotalSamples);
    f = time2freq(t);
    spc = 12.5e9;
    df = f(2)-f(1);
    aux1 = round(spc/df);
    fc = aux1*df;
    T=1/fc;
else
    % Basic specification for the simulation
    fc   = 12.5e9;   % RF frequency [Hz]
    % fc=f_RF;
    n_sc   = 30;    % Number of required subcarriers
    spare  = 10;    % Frequency spare [%]
    N      = 2^18;  % Number of points in t and (of coarse) in f

    %% Finding the apropriate frequency and sampling time
    f_max = (1 + spare/100)*n_sc*(4*fc);                                      %Maximum frequency that this code 
                                                                               %will represent, according to Nyquist
    df    = 2*f_max/(N - 1);                                                  %Interval between samples frequencies
    dt2   = ((N - 1)/N)/(2*f_max);                                           %Time interval between samples
    t_max = (N - 1)*dt2;                                                      %End time for the time vector
    %
    t    = 0:dt2:(N - 1)*dt2;                                               %Time vector
    f    = time2freq(t);  
    dff = f(2)-f(1);
    aux1 = round(fc/dff);
    spc = aux1*df
    T=1/spc;
    NPPB = length((t<=T)~=0);
end
sinal=puls_gau(t,t(1),1e-14,1);
% sinal = zeros(1,length(t));
% sinal(1) = 1;
% sinal = exp(1j*2*pi*fc*t);

figure;
hold on;
grid on;
plot(t,sinal);
% axis([0 5e-10 -0.1 1.1]);
figure;
hold on;
grid on;
plot(f,10*log10(abs((fftshift(fft(sinal./length(sinal)))))));
SymPer=T;
% SymPer=NPPB;
MaxPort = 8;
MaxNumStag=ceil(log2(MaxPort));
Count=ones(MaxNumStag,1);
[EoutAux1,~,VetMap]=OpticalFFTN(t,SymPer,MaxNumStag,sinal);
VetMap
% [ EoutAux1,~,~ ] = OpticalFFT(t,SymPer,MaxNumStag,Count,sinal);
% [ EoutAux1,~,~ ] = OpticalFFTExp(t,SymPer,MaxNumStag,Count,sinal);

figure;
hold all;
% plot(f,db(abs((fftshift(fft(EoutT))))/length(EoutT)));
for kk=1:size(EoutAux1,1)
    subplot(MaxPort,1,kk);
    hold on;
    grid on;
    if mod(kk,2)
        plot(f,20.*log10(abs((fftshift(fft(EoutAux1(VetMap==kk,:)/length(EoutAux1(VetMap==kk,:))))))),'r','LineWidth',2);
    else
        plot(f,20.*log10(abs((fftshift(fft(EoutAux1(VetMap==kk,:)/length(EoutAux1(VetMap==kk,:))))))),'b','LineWidth',2);
    end
    if thisone==1
%         axis([-2e11 -0.25e11 -150 -108]);
    elseif thisone==2
%         axis([-2e11 -0.25e11 -150 -108]);
    else
%         axis([-1.65e12 -1.062e12 -140 -80]);
    end
    a=1;
end
figure;
hold all;
conting = 1;
for kk=1:size(EoutAux1,1)
    subplot(MaxPort/2,1,conting);
    hold on;
    grid on;
    if mod(kk,2)
        plot(f,20.*log10(abs((fftshift(fft(EoutAux1(VetMap==kk,:)/length(EoutAux1(VetMap==kk,:))))))),'r','LineWidth',2);
        conting = conting + 1;
        if conting>MaxPort/2
            break;
        end
    end
    if thisone==1
%         axis([-2e11 -0.25e11 -150 -108]);
    elseif thisone==2
%         axis([-2e11 -0.25e11 -150 -108]);
    else
%         axis([-1.65e12 -1.062e12 -140 -80]);
    end
    a=1;
axis([0 1e11 -200 -100]);
end
figure;
hold all;
conting = 1;
for kk=1:size(EoutAux1,1)
    subplot(MaxPort/2,1,conting);
    hold on;
    grid on;
    if ~mod(kk,2)
        plot(f,20.*log10(abs((fftshift(fft(EoutAux1(VetMap==kk,:)/length(EoutAux1(VetMap==kk,:))))))),'b','LineWidth',2);
        conting = conting + 1;
        if conting>MaxPort/2
            break;
        end
    end
    if thisone==1
%         axis([-2e11 -0.25e11 -150 -108]);
    elseif thisone==2
%         axis([-2e11 -0.25e11 -150 -108]);
    else
%         axis([-1.65e12 -1.062e12 -140 -80]);
    end
    a=1;
axis([1.25e10 1.125e11 -200 -100]);
end

    figure;
    hold on;
    grid on;
for kk=1:size(EoutAux1,1)
    plot(f,20.*log10(abs((fftshift(fft(EoutAux1(kk,:)/length(EoutAux1(kk,:))))))));
    if thisone==1
%         axis([-2e11 -0.25e11 -150 -108]);
    elseif thisone==2
%         axis([-2e11 -0.25e11 -150 -108]);
    else
%         axis([-1.65e12 -1.062e12 -140 -80]);
    end
    a=1;
end
% axis([-0.1e11 (size(EoutAux1,1))*1.1*fc -250 10]);

a=1;