%                                                                         %
function [EVM_dB, EVM_porc,EVMrms] = evm1(M,ThisModula,simb_tx,simb_rx)              %
%                                                                         %
% ======================================================================= %
%       --------  EVM - Error Vector Magnitude (RMS)  --------            %
%  Calcula a M�trica EVMrms em [dB] e [%] de diagramas de constela��o     % 
%  de sistemas de trasnmiss�o digital multiportadoras (OFDM, DMT, etc).   %
%                                                                         %
%  Entradas:                                                              %
%     M       - �ndice de Modula��o da Constela��o (M = 16 para 16QAM)    %
%     simb_tx - Vetor de S�mbolos da Constela��o Transmitidos             %
%     simb_rx - Vetor de S�mbolos da Constela��o Recebidos                %
%                                                                         %
%  Sa�das:                                                                %
%     EVM_dB   - EVM em dB                                                %
%     EVM_porc - EVM em porcentagem                                       %
%     EVMrms   - EVM em rms                                               %
%                                                                         %
% by Jair Silva                                                           %
% Mar�o 2009                                                              %
%                                                                         %
%                                                                         %
% UpDate by Pablo Marciano  - Adicionado sa�da EVMrms, corre��o do nome   %
% da fun��o                                                               %
% Junho 2018                                                              %
% ======================================================================= %

% Gera e Normaliza uma constela��o
switch ThisModula
    case 'qam'
        constel = qammod([0:M-1],M);           % Gera a constela��o
    case 'dpsk'
        constel = dpskmod([0:M-1],M);           % Gera a constela��o
    case 'pam'
        constel = pammod([0:M-1],M);           % Gera a constela��o
    otherwise
        constel = pammod([0:M-1],M);           % Gera a constela��o
end

potmed  = mean(abs(constel).^2);       % Potencia m�dia dos s�mbolos �nicos
escala  = modnorm(constel,'avpow',1);  % fator de escala
constel = escala.*constel;             % Normaliza a constela��o ideal.

% Normaliza os s�mbolos trasnmitidos
escala  = modnorm(simb_tx,'avpow',1);  % fator de escala
simb_tx = escala.*simb_tx;             % Normaliza os simbolos unicos.

% Normaliza os s�mbolos recebidos
escala  = modnorm(simb_rx,'avpow',1);  % fator de escala
simb_rx = escala.*simb_rx;             % Normaliza os simbolos unicos.

% Calcula a EVM
num     = mean((abs(simb_tx - simb_rx)).^2);
den     = mean(abs(constel).^2);
EVMrms  = sqrt(num/den);

% Determina os valores em dB e porcentagem (%)
EVM_dB   = 20*log10(EVMrms);           % em [dB]
% EVM_porc = (EVMrms*100)/potmed;      % em [%]
EVM_porc = (EVMrms*100);        % em [%]

%
% -------------------------------------------------------------------------
%