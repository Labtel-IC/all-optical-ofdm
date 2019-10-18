function [I,Q] = DQPSK_encoder(u,v,Io,Qo)

%==========================================================================
% function [Iout,Qout] = DQPSK_encoder(in,quad)
%
% Codifica dados de entrada para um sistema DQPSK óptico
%   
%    ENTRADA:
%        u : dados em fase de entrada
%        v : dados em quadratura de entrada
%        
%    SAIDA:
%
%        I: dados em fase diferencial
%        Q: dados em quadratura diferencial
%                                             by Diogo Coelho
%                                            17/02/2017
%                                            
%   
%==========================================================================


% baseado no artigo:

% "Implementation of Differential Precoder for High-Speed Optical DQPSK Modulation"


aux1 = length(u);
aux2 = length(v);

if aux1~=aux2

    error('dados de entrada devem ter o mesmo tamanho')
    
end

% Io = 0;   % estado inicial para fase
% Qo = 0;   % estado inicial para quadratura


II = zeros(1,aux1);
QQ = zeros(1,aux1);


for ii = 1:1:aux1
   
    
%    a = and(not(q),xor(i,v(ii)));
%    b = and(q,not(xor(i,u(ii))));
%    
%    
%    I(ii) = or(a,b);
%    Q(ii) = not(xor(xor(q,u(ii)),v(ii)));


   a = and(and(not(u(ii)),not(v(ii))),not(Io));
   b = and(and(not(u(ii)),v(ii)),Qo);
   c = and(and(u(ii),v(ii)),Io);
   d = and(and(u(ii),not(v(ii))),not(Qo));
    
   e = and(and(not(u(ii)),not(v(ii))),not(Qo));
   f = and(and(not(u(ii)),v(ii)),not(Io));
   g = and(and(u(ii),v(ii)),Qo);
   h = and(and(u(ii),not(v(ii))),Io);
   
   
   
   II(ii) = or(or(a,b),or(c,d));
   QQ(ii) = or(or(e,f),or(g,h));
   
     
   Io = II(ii);
   Qo = QQ(ii);
    
end

I = zeros(1,aux1);
Q = zeros(1,aux1);

for iii = 1:aux1
    
   if II(iii) == 0
    I(iii) = -1;
   else
    I(iii) = 1;   
   end 
   
   if QQ(iii) == 0
    Q(iii) = -1;
   else
    Q(iii) = 1; 
   end 
end

clear l
