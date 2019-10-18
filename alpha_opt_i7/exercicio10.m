
% =========================================================================
%                                                                         %
%   *************  C�digo Transceptor OFDM para Wireless  ************    %
%                                                                         %
%            ----- Objetivo: Gerar Curvas de BER x SNR  -----             %  
%                                                                         %
%         Para a Disserta��o de Mestrado do Leonardo Ribas Castor         %
%                                                                         %
% by Jair Silva                                                           %
% Para Comunica��o Digital da P�s-gradua��o                               % 
%    Novembro/2013                                                        %
% =========================================================================

% Inicializa��o
clc, clear all, close all, 

disp ('                                                                 ');
disp ('=================================================================');
disp ('    Sistema de Medi��o de BER versus SNR para Wireless           ');
disp ('=================================================================');
disp ('                                                                 ');

% ......... Parametros de entrada (Algum Padr�o OFDM ) ....................
pe.Bw    = 10e6;                 % Largura de banda do sistema
pe.Nfft  = 1024;                 % Pontos da IFFT/FFT
pe.Ns    = 200;                  % Qtd Subportadoras de dados 
pe.Ng    = 32;                   % Qtd de amostras do Prefixo C�clico (IG)
pe.IG    = round(pe.Ng/pe.Nfft); % Fator do intervalo de guarda IG

% ..... Parametros Calculados para a implementa��o .......................
pe.fs    = pe.Bw/4;             % Taxa de Amostragem do AD/DA
pe.N     = (pe.Nfft - 2)/2;     % Qtd Subportadoras considerando Hermitiana
pe.NsH   = 2*pe.Ns + 2;         % Qtd de subportadoras ap�s Hermitiana 
pe.Nz    = pe.Nfft - pe.NsH;    % Qtd subportadoras zeradas (zero padding)
pe.M     = 4;                   % �ndice de modula��o --> Escolhido
pe.n     = log2(pe.M);          % Qtd de bits por subportadora
pe.frame = 10;                  % Qtd de s�mbolos OFDM no frame - Escolhido

% ..... Parametros Calculados para Informes ..............................
pe.DF    = pe.Bw/pe.Nfft;             % Espa�amento em freq. entre subport.
pe.Bwe   = pe.Ns*pe.DF;               % Largura de banda efetiva
pe.Rb    = (pe.DF*pe.Ns*pe.n)/(1+pe.IG);  % Taxa de transmiss�o REAL 
pe.Tu    = 1/pe.DF;                   % Dura��o �til do sinal OFDM
pe.Tg    = pe.Ng*(1/pe.fs);         % Dura��o do intervalo de guarda (CP)
pe.Ts    = pe.Tu + pe.Tg;           % Dura��o total sinal OFDM ( = t(end))
pe.Perda = 10*log10(1-pe.Tg/pe.Ts); % Perda gerada pela inser��o do IG.

% ..... Parametros do Loop de Simula��o ..................................
vsnr  = [0:2:16];        % Vetor SNR em dB
max_erros = 100;         % Computa BER depois de no maximo 100 erros
max_loop  = 1000;        % Numero maximo de itera��es do loop interno


% .................. Transmissor OFDM  ...................................

% Gera��o dos dados
Dadostx = randint(pe.Ns,pe.frame,pe.M);  

% Modula��o QAM (Mapeamento)
DadosModtx = qammod(Dadostx,pe.M);

% Simetria Hermitiana --> antes da centraliza��o do espectro
Hermit_tx = [zeros(1,pe.frame); DadosModtx; zeros(1,pe.frame); conj(flipud(DadosModtx))];  

% Zero padding --> Centraliza o espectro
pe.TH   = length(Hermit_tx);       
info_tx = zeros(pe.Nfft,pe.frame);   
info_tx(1:pe.TH/2,:)         = Hermit_tx(1:end/2,:);
info_tx(end-(pe.TH/2-1):end,:) = Hermit_tx(end/2+1:end,:);

% "Modula��o" via IDFT (ifft)
ofdm_symbol_tx = ifft(info_tx,pe.Nfft);

% Intervalo de Guarda
ofdm_symbol_ig = [ofdm_symbol_tx(end+1-pe.Ng:end,:); ofdm_symbol_tx];

% Concatena todos os s�mbolos OFDM
ofdm_signal  = reshape(ofdm_symbol_ig,1,pe.frame*length(ofdm_symbol_ig));

% Normaliza o sinal
pe.fat_norm = (max(ofdm_signal));
ofdm_signal = ofdm_signal./pe.fat_norm;

% Calcula a pot�ncia do sinal OFDM gerado
potw   = sum(ofdm_signal.*conj(ofdm_signal))/length(ofdm_signal);
potdBm = 10*log10(potw/1e-3); 
pe.pot = potdBm;

% Calcula o PAPR do sinal OFDM gerado
num  = max((abs(ofdm_signal)).^2);
den  = mean((abs(ofdm_signal)).^2);
PAPR = 10*log10(num/den);
pe.PAPR = PAPR;

% .........  Gera vetores tempo e frequencia .............................
t  = linspace(0,pe.frame*pe.Ts,length(ofdm_signal));
f  = linspace(0,pe.Bw,length(ofdm_signal));


% .... Gera o perfil de pot�ncia do canal (taps no dom�nio do tempo) ......
h = [(randn+j*randn) (randn+j*randn)/2 (randn+j*randn)/4 (randn+j*randn)/4 ];


% ================  Inicia o Loop BER x SNR ===============================
figure;
v_ber = [];                 % inicializa��o
ii = 0;
for snr = vsnr              % loop no vetor de SNR
    vv_ber = [];            % inicializa��o
    loop = 1;               % inicializa o loop interno
    n_erros = 0;            % reseta o contador de erros
    while loop > 0          % loop interno

        %  ################# Canal de Comunica��o #########################

        % Convolu��o com o canal
        ofdm_signal_noise = conv(ofdm_signal, h);

        % Retira as amostras excedentes devido � convolu��o
        ofdm_signal_noise = ofdm_signal_noise(1:length(ofdm_signal));

        % Insere Ru�do AWGN
        ofdm_signal_noise = awgn(ofdm_signal_noise,snr,'measured');

        % #################################################################



        % ..................... Receptor OFDM  ............................

        % Retira a normaliza��o
        ofdm_signal_rx = ofdm_signal_noise.*pe.fat_norm;

        % Convers�o serie - Paralelo (Matriz de sinais OFDM)
        ofdm_signal_rx = reshape(ofdm_signal_rx,length(ofdm_signal_rx)/pe.frame,pe.frame);

        % Retita o Intervalo de Guarda (IG)
        ofdm_symbol_rx = ofdm_signal_rx(pe.Ng+1:end,:);

        % Demodula��o via DFT (fft)
        info_rx = fft(ofdm_symbol_rx,pe.Nfft);

        % Reconhecimento do Canal
        HH = info_rx(:,1:3)./info_tx(:,1:3);  % Coeficientes do Equalizador
        H = (mean(HH.')).';
        [lin,col] = size(info_rx);
        HH = repmat(H,1,col);

        % Equaliza��o via one tap equalizer (Zero forcing)
        info_rx = info_rx./HH;
        % info_rx = info_rx;

        % Desfaz o zero padding
        Hermit_rx = [info_rx(1:pe.Ns+1,:); info_rx(end-pe.Ns:end,:)];

        % Remove a simetria Hermitiana
        DadosModrx = Hermit_rx(2:pe.Ns+1,:);

        % Detec��o
        Dadosrx = qamdemod(DadosModrx,pe.M);

        % Calcula a BER
        [n_erros,BER] = biterr(Dadosrx(:,4:end),Dadostx(:,4:end),pe.n);
        vv_ber = [vv_ber BER];

        % ********** Analisa os Crit�rios de Parada ***********************
        %
        if (loop >= max_loop)|(n_erros >= max_erros)
            n_loop = loop;
            loop = 0;
        else
            loop = loop + 1;
        end;
    end;
    ber = mean(vv_ber);
    v_ber = [v_ber ber];
    if ber == 0
        ii = ii + 1;
    end
    semilogy(snr,ber,'bo');
    axis([min(vsnr) max(vsnr)  1e-6 1]);
    drawnow;
    hold on;
    grid on;
    
    % ********* Mostra Resultados Provis�rios na Tela *********************
    disp(sprintf('SNR: %4.1f dB BER: %5.1e Bitstream: %7d N_erros: %7d',...
    snr,ber,n_loop*pe.frame,n_erros));
end
    

%........................... Plota as Curvas ..............................
% Simulada
clf
v_ber = v_ber(1:end-ii);
vsnr = vsnr(1:end-ii);
snr_min = min(vsnr);
snr_max = max(vsnr);
vsnr_fit = snr_min:0.25:snr_max;  % interpola os valores de SNR 
ber_fit  = berfit(vsnr,v_ber,vsnr_fit);  % ajusta a curva
% figure;
semilogy(vsnr,v_ber,'bo',vsnr_fit,ber_fit,'b-');
grid on;
title('BER x SNR')
xlabel('SNR');
ylabel('BER');
%
% Te�rica (Emp�rica)
hold on;
for j=1:length(vsnr)
    EbNo = vsnr(j)-10*log10(log2(pe.M));
    ber_awgn_teorico(j) = berawgn(EbNo, 'qam', pe.M);
    ber_rayl_teorico(j) = berfading(EbNo, 'qam', pe.M, 1);
end
semilogy(vsnr,ber_awgn_teorico,'k-*')
hold on, semilogy(vsnr,ber_rayl_teorico,'r*')
grid on;
legend('Simulado','Simulado','Emp�rico',3,'location','NorthEastOutside');


% ......................... Comandos de Sa�da ............................
disp ('                                                                 ');
fprintf(' Taxa de amotragem [Sps]              : %0.7g \n',pe.fs);
fprintf(' Largura de Banda Total [Hz]          : %0.7g \n',pe.Bw);
fprintf(' Taxa de Transfer�ncia [bps]          : %0.7g \n',pe.Rb);
fprintf(' Largura de Banda Efetiva [Hz]        : %0.7g \n',pe.Bwe);
fprintf(' Subportadoras de Informa��o          : %0.7g \n',pe.Ns);
fprintf(' Espa�amento entre Subportadoras [Hz] : %0.7g \n',pe.DF);
fprintf(' N�vel de Modula��o por Subportadora  : %0.7g QAM \n',pe.M);
fprintf(' Quantidade de Pontos da IFFT/FFT     : %0.7g \n',pe.Nfft);
fprintf(' Dura��o Total do S�mbolo OFDM [s]    : %0.7g \n',pe.Ts);
fprintf(' Dura��o �til do S�mbolo OFDM [s]     : %0.7g \n',pe.Tu);
fprintf(' Dura��o do Intervalo de Guarda [s]   : %0.7g \n',pe.Tg);
fprintf(' Perda SNR por inser��o do IG [dB]    : %0.7g \n',pe.Perda);
fprintf(' Potencia Sinal OFDM [dBm]            : %0.7g \n',pe.pot);
fprintf(' PAPR Raz�o Pot. M�x. e Pot M�d. [dB] : %0.7g \n',pe.PAPR);
fprintf(' Taxe de Erro de Bits BER             : %0.7g \n',BER);
disp ('                                                                 ');
disp ('=================================================================');
%