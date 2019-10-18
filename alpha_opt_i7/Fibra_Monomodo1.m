function [uf,LD,LNL] = Fibra_Monomodo1(temp,u,lambda,To,L,regime,frequency,Fibra_Monomodo_Input_File)
%%
% function [uf,LD,LNL] = Fibra_Monomodo(t,u,lambda,To,L,regime,Fibra_Monomodo_Input_File);
%
%     Calcula a resposta de uma fibra optica monomodo nao linear
%  atraves do metodo Split-Step Fourier.
%
%  ENTRADA:
%    t                         : Vetor tempo [s]
%    u                         : Campo eletrico a ser aplicado na fibra (envelope)
%    lambda                    : Comprimento de onda do sinal optico  [m]
%    To                        : Largura do pulso TDM [s]
%    L                         : Comprimento da fibra [km]
%    regime                    : 1 -> roda o split-step  /  0-> usa a funcao
%                                de transf. linear da fibra
%    Fibra_Monomodo_Input_File : Arquivo com os parametros da fibra
%
%  SAIDA:
%    uf                      : Campo eletrico na sa?da da fibra (dominio do tempo)
%    LD                      : Comprimento de dispersao  [km]
%    LNL                     : Comprimento de nao linearidade  [km]
%
%       Ref.: G. Agrawal, "Nonlinear Fiber Optics,"
%             chapter 2, pp. 44-49, 1989.
%                                        by M.Segatto e A. Togneri
%                                          08/12/2004
%                Modificado do original nlse2.m escrito por
%            M. Freitas, M. Ribeiro, e M. Segatto
%                         DEL/CT/UFES - 1996/1997
%
%
if (nargin <= 7),
  Fibra_Monomodo_Input_File = 'Fibra_Monomodo_Input_Data';
end
eval(Fibra_Monomodo_Input_File);

%
%    Calculos intermediarios
%
L           = L*1e3;      % agora em [m]
To          = To/1e-12;   % agora em [ps]
n           = length(temp);
% window      = 64;
a           = n - 1;
%   Eixo de frequencia
%
c           = 3e8;
freq        = c./(lambda);
omegaf      = 2*pi*freq;
d_omegaf_r  = omegaf-omegaf(lambdar);
%
%   Constante de atenuacao
%

 alfa        = (alfadB/(10*log10(exp(1))))*1e-3;    % Transformaçao dB/km para Np/m

%
%   Dispersao
%

Dref        = S.*(lambda(lambdar)-1550e-9) + D1550;
D           = S.*(lambda-lambda(lambdar))+Dref ;      % [s/(m-m)]
%
beta2       = -((lambda.^2).*D/(2*pi*c));             % [s^2/m]
% beta2       = - 20*1e-27;      % [s^2/m]

beta3       = (D.*lambda.^3)./((2*pi*c)^2);           % [s^3/m]
d           = beta2.*d_omegaf_r ;                     % Walk-off [m]
%
%   Nao linearidade
%
gama        = n2*omegaf/(c*Aeff);

% % Nao linearidade tirada do Pulse evolution (19/02/2010)
% gama = 3;  		  % Nonlinear coefficient em 1/(W-Km)
% gama = gama*1E-3; % Agora em 1/(W-m)

%
fmax        = 0.5*n/temp(n);
d_freq      = 2*fmax/(n-1);
% frequency   = (-fmax:d_freq:fmax);
% frequency   = freq;
omega       = 2*pi*frequency;
omega       = fftshift(omega);
%
  LD   = To^2/abs(beta2*1E27);                 % [km]
%   LNL  = (1./(gama.*max(abs(u))^2))/1000;          % [km]
  
  Pw_ent_fibra = sum(abs(u).^2)/length(u);
  LNL  = (1./(gama.*Pw_ent_fibra))/1000;
%
switch regime
%
case 0,
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     F U N C A O  D E  T R A N S F .  L I N E A R      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
  beta21     = beta2*1E27;         % [ps^2/km]
  omega1     = (omega)*1E-12;
%
  U          = fft(u);
%
% MODIFICADO em 07/2009 no IST
  % Para Potência
%   UF1         = U.*exp((1j/2)*beta21*(L/1000)*omega1.^2)*exp(-(alfa)*L);
  
  % Para Campo
%   UF1         = U.*exp((1j/2)*beta21*(L/1000)*omega1.^2)*exp(-(alfa/2)*L);
  
  % Pedido do Cartaxo
  UF1         = U.*exp(1j*((D*lambda^2)/(4*pi*c))*(((2*pi)^2)*(frequency.^2)*L))*exp(-(alfa/2)*L);
%
  uf          = ifft(UF1);
%
case 1,
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         S P L I T - S T E P  F O U R I E R            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
  nc        = length(d);
  omega2    = omega.^2;
  omega3    = omega2.*omega;
  v6jomega  = 6*j*omega;
  v3jomega2 = 3*j*omega2;
%
  Az        = u;
%
  for ii=1:nc,
    Dw(ii,1:n)   = (1/6)*((-3*alfa)-(v6jomega*d(ii))+(v3jomega2*beta2(ii))-(j*beta3(ii)*omega3));
  end
%
%   inicio do Split-Step
%
  Z     = 0;
  soma  = sum(2*abs(Az).^2,1);
  for ii = 1:nc,
    v1      = j*gama(ii);
    v2      = abs(Az(ii,:));
    N2(ii,:) = v1*(soma-(v2).^2); % Operador nao linear
  end
%
  Azw   = fft(Az);
  while (Z < L),
    h1      = 3e-3/(2*mean(abs(N2),2));
    h       = min(h1);
   
    LZ      = L - Z;
    if (LZ < h) && (LZ ~= 0),
      h = LZ;
    end
    Z       = Z + h;
    h2      = h/2;
    exph2Dw = exp(h2*Dw);
    soma1   = sum(2*abs(Az).^2,1);
%
    Azw     = exph2Dw.*Azw;    % Dispersao no intervalo h/2
    Az1d    = ifft(Azw);
%
    soma2   = 0;
    for i = 1:nc,
      N1(i,:)  = j*gama(i).*( soma1- abs(Az(i,:)).^2 );
      Azl(i,:) = exp(h*N1(i,:)).*Az1d(i,:);
    end
    Azwl    = fft(Azl);
    Azwl    = exph2Dw.*Azwl;
    Azl     = ifft(Azwl);
    soma2   = sum(2*abs(Azl).^2, 1);
    for i = 1:nc,
      N2(i,:)  = j*gama(i)*(soma2 - abs(Azl(i,:)).^2);
      Azl(i,:) = exp( h*(N1(i,:)+N2(i,:))/2 ).*Az1d(i,:);
    end
    soma2   = zeros;
    Azwl    = fft(Azl);
    Azwl    = exph2Dw.*Azwl;
    Azl     = ifft(Azwl);
    soma2  = sum(2*abs(Azl).^2,1);
    for i=1:nc,
      N2(i,:)  = j*gama(i)*(soma1-(abs(Azl(i,:))).^2);
      Azl(i,:) = exp(h*(N1(i,:)+N2(i,:))/2).*Az1d(i,:);
    end
%
    Azw = fft(Azl,[],2);
    Azw = exph2Dw.*Azw;
    Az  = ifft(Azw,[],2);
  end;
uf = Az;
%
otherwise,
  error('Opcao invalida. favor verificar a variavel "regime"')
end

