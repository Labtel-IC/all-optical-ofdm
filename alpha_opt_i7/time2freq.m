function freq = time2freq(T)
%
%           function freq = time2freq(T);
%
%   Gera o vetor de freq��ncia para um dado vetor tempo para ser usado junto com a FFT.
%
%   ENTRADA:
%       T            : Vetor tempo
%
%   SA�DA:
%       freq         : Vetor freq��ncia

dT    = mean(diff(T));
nT    = length(T);
nT2   = nT/2;
df    = inv(dT*nT);
%
freq   = df*(-nT2:nT2 - 1);