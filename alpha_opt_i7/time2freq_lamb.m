function [freq,lamb] = time2freq_lamb(T,lc)
%c
%c           function [freq,lamb] = time2freq_lamb(T,lc);
%c
%c   Gera os vetores de frequencia e comprimento de onda para um dado vetor
%c   tempo.
%c
%c                                        by   :           M.Segatto
%c                                        Date :                2004
%c                                        email: segatto@ele.ufes.br
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                    G L O B A L  V A R I A B L E S                      %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c   
%c
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                   L I S T  O F  V A R I A B L E S                      %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c
%c   Entrada:
%c   T            : Vetor tempo                                        [ps]
%c   lc           : Comprimento de onda central                        [nm]
%c
%c   Saida:
%c   freq         : Vetor frequencia                                  [THz]
%c   lamb         : Vetor comprimento de onda                          [nm]
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                             S T R U C T S                              %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c   
lamb = 0;
cnmps = 2.998E5;                            % Velocidade da luz [nm/ps]

dT    = mean(diff(T));
nT    = length(T);
df    = inv(dT*nT);                        	 % [THz]
switch nargin
  case 1
    fc = 0;
    for n = 1:nT,
      freq(n) = fc - (n - 1 - nT/2)*df;	   	% [THz] vetor de freq
    end
  case 2
    fc    = cnmps/lc;	                        % [THz]
    for n = 1:nT,
      freq(n) = fc - (n - 1 - nT/2)*df;	   	% [THz] vetor de freq
    end
    lamb  = (cnmps./freq);                      % Gera vetor de lambdas [nm]
%
 otherwise error('Numero de parametros de entrada nao confere...');
end
