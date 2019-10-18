function [ Filtro ] = Filtro_Retangular( FBW,fc,f)
%c
%c function [ Filtro ] = Filtro_Retangular( FBW,fc,f);
%c
%c
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                   L I S T  O F  V A R I A B L E S                      %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c
%c	Filtro de janela retangular, essa função considera que todas os
%c paramentros passados como entrada estarão na mesma unidade e grandeza.
%c
%c  Entrada:
%c 	FBW    : Largura de banda do filtro [xHz]
%c  fc     : Frequência central do filtro [xHz]  
%c                  (dominio do tempo)
%c  f 	   : Vetor frequência [xHz]  
%c
%c  Saida
%c  Filtro : Filtro retangular de saida
%c
%c
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                             S T R U C T S                              %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c
%%
Filtro = zeros(1,length(f));
% Filtro = linspace(-99999,-99999,length(f));
Filtro(((fc - FBW/2)<f)&(f<(fc + FBW/2))) = 1;

end

