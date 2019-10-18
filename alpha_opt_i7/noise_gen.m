function [out] = noise_gen(noise_d,t,BWO)


%==========================================================================
% function [out] = noise_gen(in,noise_d,t,BWO)
%
% Gera ruido gaussiano ao longo de todo o espectro.
%
%
%       
%                                             by Diogo Coelho
%                                            30/11/2015
%==========================================================================
warning off;

noise_dens       = (10 .^ (noise_d./ 10));


% Geração do ruído (gaussiano) em I e em Q 

randn('state', sum(100*clock));                                        % Garante a aleatoriedade do ruído guassiano
noise_fase = sqrt(noise_dens .* BWO) .* randn(1,length(t));
% noise_quad = sqrt(noise_dens .* BWO) .* randn(1,length(t));
noise_quad = 0;

% Adição do ruído ao sinal
    
out = noise_fase + 1i.*noise_quad;

