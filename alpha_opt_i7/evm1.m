%                                                                         %
function [EVM_dB, EVM_porc,EVMrms] = evm1(M,ThisModula,simb_tx,simb_rx)              %
%                                                                         %
% ======================================================================= %
%       --------  EVM - Error Vector Magnitude (RMS)  --------            %
%  Calcula a Métrica EVMrms em [dB] e [%] de diagramas de constelação     % 
%  de sistemas de trasnmissão digital multiportadoras (OFDM, DMT, etc).   %
%                                                                         %
%  Entradas:                                                              %
%     M       - Índice de Modulação da Constelação (M = 16 para 16QAM)    %
%     simb_tx - Vetor de Símbolos da Constelação Transmitidos             %
%     simb_rx - Vetor de Símbolos da Constelação Recebidos                %
%                                                                         %
%  Saídas:                                                                %
%     EVM_dB   - EVM em dB                                                %
%     EVM_porc - EVM em porcentagem                                       %
%     EVMrms   - EVM em rms                                               %
%                                                                         %
% by Jair Silva                                                           %
% Março 2009                                                              %
%                                                                         %
%                                                                         %
% UpDate by Pablo Marciano  - Adicionado saída EVMrms, correção do nome   %
% da função                                                               %
% Junho 2018                                                              %
% ======================================================================= %

% Gera e Normaliza uma constelação
switch ThisModula
    case 'qam'
        constel = qammod([0:M-1],M);           % Gera a constelação
    case 'dpsk'
        constel = dpskmod([0:M-1],M);           % Gera a constelação
    case 'pam'
        constel = pammod([0:M-1],M);           % Gera a constelação
    otherwise
        constel = pammod([0:M-1],M);           % Gera a constelação
end

potmed  = mean(abs(constel).^2);       % Potencia média dos símbolos únicos
escala  = modnorm(constel,'avpow',1);  % fator de escala
constel = escala.*constel;             % Normaliza a constelação ideal.

% Normaliza os símbolos trasnmitidos
escala  = modnorm(simb_tx,'avpow',1);  % fator de escala
simb_tx = escala.*simb_tx;             % Normaliza os simbolos unicos.

% Normaliza os símbolos recebidos
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