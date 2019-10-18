function [E,Ei,Eq]=qi_MZM(E, sinal_i, sinal_q,V_bias,V_pi,t)
% QI Mach-Zehnder modulator

%==========================================================================
% function E=qi_MZM(E, sinal_i, sinal_q,V_pi,V_bias,ext_ratio)
%
% Modulador MZM IQ
%   
%    ENTRADA:
%        Ein            : Campo de entrada
%        sinali         : sinal eletrico em fase de entrada
%        sinalq         : sinal eletrico em quadratura de entrada
%        V_pi           : Vpi do MZM
%        V_bias         : tensao de polarizacao
%        ext_ratio      : razao de extincao
%        
%    SAIDA:
%
%        E              : Campo de saida      
%
%                                             by Diogo Coelho
%                                            09/10/2015
%==========================================================================


Ei = sqrt(2).*mzm1arm(E,sinal_i,V_pi,V_bias);
Eq = sqrt(2).*mzm1arm(E,sinal_q,V_pi,V_bias).*exp(1i*pi/2);

E = Ei + Eq;  
 
if nargin==6
    figure;plot(t,abs(E),t,real(E),t,imag(E));
    axis([0 16*1/(12.5e9) 1.1*min(real(E)) 1.1*max(real(E))]);
end