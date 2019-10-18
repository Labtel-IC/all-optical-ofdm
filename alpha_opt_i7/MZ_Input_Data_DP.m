%c
%c                                                       ..'Ž`'..'Ž`..'Ž`..                                                   
%c
%c  	Arquivo de entrada usado pela funcao Mach_Zehnder_Modulator. Contem
%c informacoes sobre as caracteristicas fisicas do  modulador Mach-Zehnder.
%c
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                    G L O B A L  V A R I A B L E S                      %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c   
%c
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                   L I S T  O F  V A R I A B L E S                      %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c
%c	L	      : Comprimento do dispositivo [cm]
%c  U10       : Tensao de polarizacao do braço superior [V]
%c  U20       : Tensao de polarizacao do braço inferior [V]
%c  eletrodes : Numero de eletrodos (1,2). Definido pelo 
%c              numero de variaveis de entrada na funcao
%c              Mach_Zehnder_Modulator
%c  U_pi1     : Tensao de chaveamento para o eletrodo 1 [V]
%c  U_pi2     : Tensao de chaveamento para o eletrodo 2 [V]
%c  eta1      : Sensibilidade no caminho 1  [1/V.m]
%c  eta2      : Sensibilidade no caminho 2  [1/V.m]
%c  V1pi0	  : Constante da Tensao de chaveamento para o eletrodo 1 [V]
%c  V2pi0	  : Constante da Tensao de chaveamento para o eletrodo 2 [V]
%c  Vphi0	  : Constante da Tensao de chaveamento para controle de fase
%c                                                                      [V]
%c  nopt	  : Indice de refracao optico
%c  nel	      : Indice de refracao eletrico
%c  alfa_ins  : Perda por insercao [dB]
%c  phi_0 	  : Constante de fase entre os dois caminhos
%c  alfa0	  : Perda condutiva [dB/cm.GHz^0.5]
%c
%c                     by M.Segatto, S. Cani, B. Jesus e A. Togneri
%c                                          04/09/2002
%c
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                             S T R U C T S                              %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c   
%

L          = 10;
U10        = 1.550;
U20        = 1.625;
U_pi1      = 3.1;
U_pi2      = 3.25;
eta1       = 89;% 89;
eta2       =-12.7;%-12.7;
V1pi0      = 0.350;
V2pi0      = 1.200;
Vphi0      = 0.100;
%
nopt       =  2.17;
nel        =  2.60;
alfa_ins   =  5.1;
phi_0      =  1.85;
alfa0      =  0.55;

C     = (eta1+eta2)/(eta1-eta2);       	%  Parametro de chirp
alfa0 = 10^(alfa0/20); 			% [1/cm.GHz^0.5]


% 
% L          =  10;
% U0         =  0;
% U_pi1      =  0;
% U_pi2      =  0;
% eta1       =  89;
% eta2       = -12.7;
% %
% nopt       =  1;
% nel        =  2;
% alfa_ins   =  5.1;
% phi_0      =  0.0;
% alfa0      =  0.55;
% %
% C          =  (eta1+eta2)/(eta1-eta2);       	%  Parametro de chirp
% alfa0      =  10^(alfa0/20); 			% [1/cm.GHz^0.5]
