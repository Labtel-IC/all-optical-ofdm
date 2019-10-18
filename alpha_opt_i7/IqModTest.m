function [Eout1,Eout2] = IqModTest (t,Ein,Isig,Qsig,U_pi2,Vbias)

%This test will implemente the IQ modulator as theoretically as possible.
%The CW signal will pass through an optical coupler that will divided the
%signal in two. Each part will then pass through a MZM. At the end of each
%MZM the signal will be combined with an optical couples as well. One of
%thos signals will be out phased by 90 degrees. Finaly anoter optical
%coupler will combine those two signal to form the  output.

%%
% CoupConst = 3*pi/4;    %Couples Constant of the DI
% ExpCoef   = exp(-1j*1);%Exponential component from Ramaswami equation
% MutCoef   = 1;         %Coeficient to control polarity
% E0 = 0;
E11 = Ein./2;
E12 = Ein./2;
% E11 = ExpCoef.*(cos(CoupConst).*Ein + MutCoef*1j*sin(CoupConst).*E0);        
% E12 = ExpCoef.*(cos(CoupConst).*E0 + MutCoef*1j*sin(CoupConst).*Ein);
% 
% E11 = ExpCoef.*(cos(CoupConst).*E11 + MutCoef*1j*sin(CoupConst).*E0);        
% E12 = ExpCoef.*(cos(CoupConst).*E0 + MutCoef*1j*sin(CoupConst).*E12);
%%
% 
% L          = 10;
U0         = Vbias;
% eta1       = 89;% 89;
% eta2       =-12.7;%-12.7;
% %
% nopt       =  2.17;
% nel        =  2.60;
% alfa_ins   =  5.1;
% phi_0      =  0.0;
% alfa0      =  0.55;
% 
% C     = (eta1+eta2)/(eta1-eta2);       	%  Parametro de chirp
% alfa0 = 10^(alfa0/20); 			% [1/cm.GHz^0.5]

U1t        = Isig;
U2t        = exp(-1j*pi).*Isig;
U_pi       = U_pi2;

%
% tps        = t/1E-12;
% ccmns      = 30;			     	% Velocidade da luz [cm/ns]
% freqTHz    = time2freq_lamb(tps);
% freqGHz    = freqTHz*1e-3;			% Frequencia em GHz
% freqGHz    = -fftshift(freqGHz);
% freqGHz(1) = freqGHz(2);
% n          = length(U1t); 			% Tamanho de U1t
% %
% 
% %%%%%%%  FUNCAO DE TRANSFEWRENCIA ELETRICA DO MODULADOR  %%%%%%%%%%%
% %
% alfaf  = 0.5*alfa0*(abs(freqGHz).^0.5);
% gamaf  = 2*pi*abs(freqGHz).*(nel - nopt)/ccmns;
% atn    = alfaf + 1j*gamaf;
% H      = (1./(atn*L)).*(1 - exp(-atn*L));
% 
% U1f   = fft(U1t);
% U1t   = ifft(U1f.*H);
% U2f   = fft(U2t);
% U2t   = ifft(U2f.*H);
EoutI = E11.*cos((pi/2).*(U1t - U2t - U0)/U_pi).*exp(-1j*(pi/2).*((U1t + U2t)/U_pi));
% EoutI = E11.*cos((pi/2).*(U1t + U0)/U_pi);
% EoutI = exp(1j*pi/2).*EoutI;
%%

% L          = 10;
U0         = Vbias;
% eta1       = 89;% 89;
% eta2       =-12.7;%-12.7;
% %
% nopt       =  2.17;
% nel        =  2.60;
% alfa_ins   =  5.1;
% phi_0      =  0.0;
% alfa0      =  0.55;
% 
% C     = (eta1+eta2)/(eta1-eta2);       	%  Parametro de chirp
% alfa0 = 10^(alfa0/20); 			% [1/cm.GHz^0.5]

U1t        = Qsig;
U2t        = exp(-1j*pi).*Qsig;
U_pi       = U_pi2;

%
% tps        = t/1E-12;
% ccmns      = 30;			     	% Velocidade da luz [cm/ns]
% freqTHz    = time2freq_lamb(tps);
% freqGHz    = freqTHz*1e-3;			% Frequencia em GHz
% freqGHz    = -fftshift(freqGHz);
% freqGHz(1) = freqGHz(2);
% n          = length(U1t); 			% Tamanho de U1t
% %

%%%%%%%  FUNCAO DE TRANSFEWRENCIA ELETRICA DO MODULADOR  %%%%%%%%%%%
%
% alfaf  = 0.5*alfa0*(abs(freqGHz).^0.5);
% gamaf  = 2*pi*abs(freqGHz).*(nel - nopt)/ccmns;
% atn    = alfaf + 1j*gamaf;
% H      = (1./(atn*L)).*(1 - exp(-atn*L));
% 
% U1f   = fft(U1t);
% U1t   = ifft(U1f.*H);
% U2f   = fft(U2t);
% U2t   = ifft(U2f.*H);
EoutQ = E12.*cos((pi/2).*(U1t - U2t - U0)/U_pi).*exp(-1j*(pi/2).*((U1t + U2t)/U_pi));
% EoutQ = E12.*cos((pi/2).*(U1t + U0)/U_pi);
%%    
% [EoutI,~] = Mach_Zehnder_Modulator_simplificado(t,Ein,Isig,Fid);
% [EoutQ1,~] = Mach_Zehnder_Modulator_simplificado(t,Ein,Qsig,Fid);
EoutQ = exp(1j*1*pi/2).*EoutQ;

% Eout1 = ExpCoef.*(cos(CoupConst).*EoutI + MutCoef*1j*sin(CoupConst).*EoutQ);
% Eout2 = ExpCoef.*(cos(CoupConst).*EoutQ + MutCoef*1j*sin(CoupConst).*EoutI);

Eout = (EoutI + EoutQ).*exp(1j*pi);
% Eout = (EoutI + EoutQ);
figure;plot(abs(Eout));hold all;grid on;plot(real(Eout));plot(imag(Eout));axis([0 1023 -0.8 0.8]);
Eout1 = Eout;
Eout2 = Eout;

figure;plot(t,abs(Eout),t,real(Eout),t,imag(Eout));
axis([0 t(16*64) 1.1*min(real(Eout)) 1.1*max(real(Eout))]);
    
figure;
hold all;
grid on;
% plot(t,abs(Eout1)./max(abs(Eout1)),t,abs(EoutI)./max(abs(EoutI)),t,abs(EoutQ)./max(abs(EoutQ)),...
%      t,real(Eout1)./max(real(Eout1)),t,real(EoutI)./max(real(EoutI)),t,real(EoutQ)./max(real(EoutQ)),...
%      t,imag(Eout1)./max(imag(Eout1)),t,imag(EoutI)./max(imag(EoutI)),t,imag(EoutQ)./max(imag(EoutQ)));
% plot(t,abs(Eout2)./max(abs(Eout2)),t,abs(EoutI)./max(abs(EoutI)),t,abs(EoutQ)./max(abs(EoutQ)),...
%      t,real(Eout2)./max(real(Eout2)),t,real(EoutI)./max(real(EoutI)),t,real(EoutQ)./max(real(EoutQ)),...
%      t,imag(Eout2)./max(imag(Eout2)),t,imag(EoutI)./max(imag(EoutI)),t,imag(EoutQ)./max(imag(EoutQ)));
plot(t,abs(Eout)./max(abs(Eout)),t,abs(EoutI)./max(abs(EoutI)),t,abs(EoutQ)./max(abs(EoutQ)),...
     t,real(Eout)./max(real(Eout)),t,real(EoutI)./max(real(EoutI)),t,real(EoutQ)./max(real(EoutQ)),...
     t,imag(Eout)./max(imag(Eout)),t,imag(EoutI)./max(imag(EoutI)),t,imag(EoutQ)./max(imag(EoutQ)));
a=1;
