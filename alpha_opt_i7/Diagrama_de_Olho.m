function [Abertura,Tolho,olho]=Diagrama_de_Olho(T,pulse,...
    n_bits,offset,Bit_Rate,Sample_Rate,Sequence_Length);
%
%    function [Abertura,Tolho,olho] = Diagrama_de_Olho(T,pulse,...
%                      n_bits,offset,Bit_Rate,Sample_Rate,Sequence_Length);
%
%       Calcula o diagrama e a abertura de olho para uma dada 
%  sequencia de bits. Se o numero de saidas for nulo, plota o
%  diagrama e o histograma vertical da sequencia.
%
%   ENTRADA:
%    T                  : Vetor tempo [s]
%    pulse              : Vetor contendo o sinal [u.a]
%    n_bits             : Numero de bits a serem mostrados no diagrama
%   !! Importante       => n_bits deve ser divisivel por Sequence_Length
%    offset             : Offset [s]
%    Bit_Rate           : Taxa de Transmissao [bit/s]
%    Sample_Rate        : Numero de pontos por bit
%    Sequence_Length    : Numero de bits da sequencia
%
%   SAIDA:
%    Tolho              : Vetor tempo para o diagrama de olho [s]
%    olho               : Matriz contendo o diagrama de olho
%    Abertura           : Abertura do olho
%
%                 by M.Segatto
%                   18/02/2002
%                 Atualizado por
%                 Arnaldo P. Togneri 06/09/03   <- Inclusao de 'Abertura'
%                 Atualizado por
%                 Carlos E. Castellani 31/01/07
%                                       Cálculo da abertura pelo histograma  
%
if (offset < 0),
  offset = -1*offset;
end
%
n_bits    = round(n_bits);
TB        = inv(Bit_Rate);                 % [s]
noffset   = round((offset/TB)*Sample_Rate);
N         = Sample_Rate*Sequence_Length;
%
if (noffset == 0),
  noffset = 1;
end
%
%  Verifica se eh possivel fazer o offset. Caso nao, zera o offset
%
if noffset > N
  msg  = 'Valor muito alto para o offset. offset foi zerado!';
  warning(msg);
  noffset = 1;
end
pulse1(1:noffset+1) = pulse(N-noffset:N);
%
for ii = noffset + 2:N,
  pulse1(ii) = pulse(ii - noffset - 1);
end
%
%  Verifica se eh possivel fazer o reshape
%
if rem(Sequence_Length,n_bits) ~= 0
  msg  = 'O numero de bits do diagrama nao eh possivel! ';
  msg  = [msg 'O Diagrama serah calculado com 1 bit.'];
  warning(msg);
  n_bits = 1;
end
%
N1     = Sample_Rate*n_bits;
N2     = Sequence_Length/n_bits;
Tolho  = T(1:N1);
olho   = reshape(pulse1,N1,N2);
%
% mvalue = mean(pulse);          % valor medio flutuante
%
% Calculo da abertura do diagrama de olho
%%
for ii = 1:N1
  media(ii) = mean(olho(ii,:));
  nhi       = find(olho(ii,:) > media(ii));
  nlow      = find(olho(ii,:) < media(ii));
  hi(ii)    = mean(olho(ii,nhi));
  low(ii)   =  mean(olho(ii,nlow));
%
  Lhi(ii)   = mean(hi)  - std(olho(ii,nhi));     
  Llow(ii)  = mean(low) + std(olho(ii,nlow));   
  Ab(ii)    = Lhi(ii) - Llow(ii);
end
Abertura   = mean(Ab);
%
if (nargout == 0),
  maxolho = max(max(pulse));
  minolho = min(min(pulse));
  x       = linspace(minolho,maxolho,50);
  y       = hist(pulse,x);
  subplot('position',[0.10 0.09 0.15 0.85]),barh(x,y)
  grid on,ylabel('Amplitude [u.a]')
  subplot('position',[0.33 0.09 0.60 0.85])
  plot(Tolho/1E-12,olho,'r','LineWidth',1)
  grid on,xlabel('Tempo [ps]'),...
  ylabel(sprintf('Abertura do diagrama de olho:  %2.5f [u.a]',Abertura))
end
% 
%    Essa eh a maneiro de fazer o calculo rapidamente, porem com resultados
% inferiores.
%
% nhi  = find(pulse > mvalue);
% nlow = find(pulse < mvalue);
% %
% hi   = pulse(nhi);
% low  = pulse(nlow);
% %
% Lhigh     = mean(hi)-std(hi);     % media dos maximos - Desvio padrao
% Llow      = mean(low)+std(low);   % media dos minimos + Desvio padrao
% Abertura  = Lhigh - Llow;