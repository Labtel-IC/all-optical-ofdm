function y = EDFA_SP(x,t,GAIN_dB,NF,PWsat,BW,c,lambdac)
%==========================================================================
% function y = EDFA_NF(x,t,GAIN_dB,noised,PWsat,BW,c,lambdac);
%
% Amplifica o sinal optico e gera ruido de acordo com a equacao da ASE
%   
%    ENTRADA:
%        x        : campo de entrada
%        t        : vetor dos tempos
%        BW       : largura de banda (frequencia de amostragem)
%        c        : velocidade da luz
%        lambdac  :comprimento de onda da portadora
%        GAIN_dB  : Ganho do amplificador [dB]
%        NF       : Figura de ruido do amplificador
%    SAIDA:
%
%        y        : saida do amplificador      
%
%                                             by Diogo Coelho
%                                            08/10/2015
%
%                                             Last update
%                                            03/02/2016
%==========================================================================

% global po 

[~,power_in]     = measure_power(x);
h                = 6.62607004e-34;   % constante de planck

%%%%%%%%%%%%%%%% saturacao %%%%%%%%%%%%%%

if power_in >= PWsat 

    gainLin       = 1;
    noise      = 0;
                
else

    gainLin = 10^(GAIN_dB/20);
    
    nf      = 10^(NF/10);
    
    nase    = (0.5).*nf.*gainLin.*h.*(c)/(lambdac.*1e-9);
    
    Nase    = 10*log10(nase);
    
    noise   = noise_gen(Nase,t,BW); % Geração do ruído (gaussiano) em I e em Q 

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y    = gainLin.*x + noise;


