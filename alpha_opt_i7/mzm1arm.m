function Eout = mzm1arm(Ein,sinal1,V_pi,V_bias)

%==========================================================================
% function Eout = mzm1arm(Ein,sinal1,V_pi,V_bias)
%
% Modula o campo eletrico como um modulador MZM de um braço
%   
%    ENTRADA:
%        Ein            : Campo de entrada
%        sinal1         : sinal eletrico de entrada
%        V_pi           : Vpi do MZM
%        V_bias         : tensao de polarizacao
%        
%    SAIDA:
%
%        Eout            : Campo de saida      
%
%                                             by Diogo Coelho
%                                            09/10/2015
%==========================================================================
    
    digits(64);
    format longEng

   
    V1 = sinal1 + V_bias;     
    phi = (V1).*pi/(2*V_pi);
    
    Eout = Ein .*cos(phi);
    
end