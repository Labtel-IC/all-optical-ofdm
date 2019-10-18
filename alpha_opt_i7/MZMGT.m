function  [Eout,H] = MZMGT(freqGHz,Ein,U1t,U2t,L,U0,U_pi1,U_pi2,nopt,nel,C,alfa0)
% L          = 10;
% U0         = 1.900;
% U_pi1      = 3.8;
% U_pi2      = 3.8;
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
%c
%c function [Eout,H] = Mach_Zehnder_Modulator(t,Ein,U1t,U2t,MZ_Input_File);
%c
%c This script is responsible to reproduce the output field of a MZM
%c modulator of one or two arms depending on the inputs.
%c
%c
%c                                           Updated by P.Marciano LG
%c                                           18/09/2017
%c                                           pablorafael.mcx@gmail.com
%c
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                    G L O B A L  V A R I A B L E S                      %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                   L I S T  O F  V A R I A B L E S                      %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c
%c	Modulador Mach-Zehnder
%c
%c  Entrada:
%c 	t  			  : Vetor de tempo                                      [s]
%c  Ein			  : Campo Eletrico de Entrada                         [V/m]  
%c                  (dominio do tempo)
%c  U1t 		  : Tensao de entrada no braco 1 do modulador           [V]  
%c                  (dominio do tempo)
%c  U2t 		  : Tensao de entrada no braco 2 do modulador           [V] 
%c                  (dominio do tempo)
%c  MZ_Input_File :Arquivo de entrada contendo parametros do modulador
%c
%c  Saida
%c  Eout          : Campo Eletrico optico modulado no tempo           [V/m]
%c  H             : Funcao de tranferencia eletrica do modulador
%c
%c
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                             S T R U C T S                              %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c
%c   Usa: time2freq_lamb
%%
% if (nargin <= 3),
%   MZ_Input_File  ='MZ_Input_Data';
%   load([MZ_Input_File '.m']);
% else
%     [L,U0,U_pi1,U_pi2,~,~,nopt,nel,~,~,C,alfa0] = Import_Data_MZM (...
%                                                             MZ_Input_File);%Import the data from files 
% end
if size(U2t,1)>1
  electrodes = 2;  % Two arms modulator
  U_pi       = U_pi2;
else
  electrodes = 1;  % One arm modulator
  U_pi       = U_pi1;
end
%
% tps        = t/1E-12;
ccmns      = 30;			     	% Velocidade da luz [cm/ns]
% freqTHz    = time2freq_lamb(tps);
% freqGHz    = freqTHz*1e-3;			% Frequencia em GHz
% freqGHz    = -fftshift(freqGHz);
% freqGHz(1) = freqGHz(2);
n          = length(U1t); 			% Tamanho de U1t
%

%%%%%%%  FUNCAO DE TRANSFEWRENCIA ELETRICA DO MODULADOR  %%%%%%%%%%%
%
alfaf  = 0.5*alfa0*(abs(freqGHz).^0.5);
gamaf  = 2*pi*abs(freqGHz).*(nel - nopt)/ccmns;
atn    = alfaf + 1j*gamaf;
H      = (1./(atn*L)).*(1 - exp(-atn*L));
%
if electrodes == 1
  U1f   = fft(U1t);
  U1t   = real(ifft(U1f.*H));
  exp1  = exp(1j*C*(pi/2).*(U1t/U_pi));
  if C~=0
      Eout  = Ein.*cos((pi/2).*(U1t - U0)/U_pi).*exp1;
  else
      Eout  = Ein.*cos((pi/2).*(U1t - U0)/U_pi);
  end
else
  U1f   = fft(U1t);
  U1t   = ifft(U1f.*H);
  U2f   = fft(U2t);
  U2t   = ifft(U2f.*H);
  %The final part of the following equation is the Chirp parameter and it
  %should not be removed, because it will be the key parameter to create
  %the optical carriers at the same level.
  Eout = Ein.*cos((pi/2).*(U1t - U2t - U0)/U_pi).*exp(-1j*(pi/2).*((U1t + U2t)/U_pi));
end
