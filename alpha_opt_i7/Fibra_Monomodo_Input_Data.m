%
%       Dados de entrada para a funcao Fibra_Monomodo.
%  Contem os parametros fisicos da fibra e o tipo de
%  simulacao a ser realizada.
%
%
%     alfadB  : Atenuacao da fibra  [dB/km]
%     D1550   : Dispersao da fibra em pwd = 1550 nm  [ps/km.nm].
%                 O valor de beta2 eh calculado a partir deste valor.
%     S       : Inclinacao da dispersao (slope)  [ps/(km.nm)^2]
%     n2      : Indice de refracao nao linear    [m^2/W]
%     Aeff    : Area efetiva da fibra [m^2]
%
%                                        by M.Segatto e A. Togneri
%                                          08/12/2004
%
%    Atenuacao
%
alfadB	    = 0.2;       % [dB/km]
% alfadB	    = alfadB/2;       % [dB/km]
%
%    Dispersao
%
D1550 	    = 16;        % [ps/(nm.km)]
% D1550	    = D1550/2;       % [dB/km]
S	        = 0.09;      % [ps/(nm.km)^2]
lambdar     = 1;
%
%   Nao linearidade
%
n2          = 3.2E-20;
% n2          = 2.5e-20;
% Aeff        = 70e-12;    % [m^2]
Aeff        = 80e-12;    % [m^2]
%
%   Calculos intermediarios
%
D1550       = D1550*1e-6;               % [s/(m-m)]    (p/ lambda 1550e-9)
S           = S*1e3;                    % [s/(m^2-m)]

