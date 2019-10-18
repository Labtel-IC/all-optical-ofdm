function  [Eout,H] = Mach_Zehnder_Modulator_DP(t,Ein,U,MZ_Input_File,Local)
%c
%c function [Eout,H] = Mach_Zehnder_Modulator_DP(t,Ein,U,MZ_Input_File,...
%c                                                                   Local)
%c
%c This script is responsible to reproduce the output field of a DP-MZM
%c modulator of one or two arms depending on the inputs.
%c
%c
%c                                           Updated by P.Marciano LG
%c                                           14/10/2017
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
%c  MZ_Input_File : Arquivo de entrada contendo parametros do modulador
%c  Local         : Caminho para o diretório que contém o arquivo de
%c                  configuração do DP-MZM
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
if (nargin < 4),
  MZ_Input_File  ='MZ_Input_Data_DP';
  load([MZ_Input_File '.m']);
elseif nargin < 5
    Local = [pwd '\'];
    [L,U10,U20,U_pi1,U_pi2,~,~,V1pi0,V2pi0,Vphi0,nopt,nel,~,phi_0,C, ...   %Import the data from files 
                         alfa0] = Import_Data_MZM_DP (MZ_Input_File,Local);
else
    [L,U10,U20,U_pi1,U_pi2,~,~,V1pi0,V2pi0,Vphi0,nopt,nel,~,phi_0,C, ...   %Import the data from files 
                         alfa0] = Import_Data_MZM_DP (MZ_Input_File,Local);
end
if isstruct(U),
  electrodes = 2;  % Two arms modulator
  U1t        = U.U1t;
  U2t        = U.U2t;
  U1_pi       = U_pi1;
  U2_pi       = U_pi2;
else
  electrodes = 1;  % One arm modulator
  U1t        = U;
  U2t        = zeros(1,length(U));  
  U1_pi       = U_pi1;
  U2_pi       = U_pi2;
end
%
tps        = t/1E-12;
ccmns      = 30;			     	% Velocidade da luz [cm/ns]
freqTHz    = time2freq_lamb(tps);
freqGHz    = freqTHz*1e-3;			% Frequencia em GHz
freqGHz    = -fftshift(freqGHz);
freqGHz(1) = freqGHz(2);
n          = length(U1t); 			% Tamanho de U1t
%

%%%%%%%  FUNCAO DE TRANSFEWRENCIA ELETRICA DO MODULADOR  %%%%%%%%%%%
%
alfaf  = 0.5*alfa0*(abs(freqGHz).^0.5);
gamaf  = 2*pi*abs(freqGHz).*(nel - nopt)/ccmns;
atn    = alfaf + 1j*gamaf;
H      = (1./(atn*L)).*(1 - exp(-atn*L));
%
% if (electrodes == 1),
U1f   = fft(U1t);
U1t   = real(ifft(U1f.*H));
U2f   = fft(U2t);
U2t   = ifft(U2f.*H);
exp1  = exp(1j*C*(pi/2).*(U1t/U1_pi));
exp2  = exp(1j*C*(pi/2).*(U2t/U2_pi));
Eout1  = (Ein./2).*cos((pi/2).*(U1t - U10  - V1pi0)/U1_pi).*exp1;
Eout2  = (Ein./2).*cos(((pi/2).*(U2t - U20 - V2pi0)/U2_pi)).*exp2;
% figure;plot(t,real(Eout2),t,imag(Eout2));hold all;
if phi_0~=0
    Eout2  = Eout2.*exp(-1j*((pi/2).*((phi_0 + Vphi0)/1.95)));
end
%figure;hold all;f_RF=1e9;n_sc=30;spare=10;N=2^19;f_max=(1+spare/100)*n_sc*(4*f_RF);df=2*f_max/(N-1);dt=((N-1)/N)/(2*f_max);t_max=(N-1)*dt;t=0:dt:(N-1)*dt;f=time2freq(t);grid on;plot(t,abs((Eout1-min(Eout1))./max((Eout1-min(Eout1)))));plot(t,((real(Eout2)+1)-min(real(Eout2)+1))./max(((real(Eout2)+1)-min(real(Eout2)+1))));
Eout = Eout1 + Eout2;


