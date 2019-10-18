function [power,power_dBm] = measure_power(E)

%==========================================================================
% function [power,power_dBm] = measure_power(E);
%
% Calcula a potencia do campo eletrico em Watts e dBm.
%
%
%       
%                                             by Diogo Coelho
%                                            08/10/2015
%==========================================================================

power_dBm = 30+10*log10(mean(abs(E.*conj(E))));

power = 10.^(power_dBm/10);