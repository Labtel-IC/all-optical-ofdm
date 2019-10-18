function [L,U10,U20,U_pi1,U_pi2,eta1,eta2,V1pi0,V2pi0,Vphi0,nopt,nel,...
             alfa_ins,phi_0,C,alfa0] = Import_Data_MZM_DP (FileIndex,Local)
%% Import data from text file.
%c function [L,U10,U20,U_pi1,U_pi2,eta1,eta2,V1pi0,V2pi0,Vphi0,nopt,nel,...
%c           alfa_ins,phi_0,C,alfa0] = Import_Data_MZM_DP (FileIndex,Local)
%c Script for importing data from the following text file:
%c
%c    [pwd '\input_files\AA_MZ_Input_Data_t' num2str(FileIndex) '.m'];
%c
%c To extend the code to different selected data or a different text file,
%c generate a function instead of a script.
%c
%c Auto-generated by MATLAB on 2017/10/04 11:05:03
%c
%c
%c
%c                                           Updated by P.Marciano LG
%c                                           18/10/2017
%c                                           pablorafael.mcx@gmail.com
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                    G L O B A L  V A R I A B L E S                      %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                   L I S T  O F  V A R I A B L E S                      %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c
%c    Import_Data_MZM
%c
%c  Input:
%c 	FileIndex : Indentifier of the file to be load into the workspace   [-]
%c  Local     : Current path for the location of the files              [-]
%c
%c  Output:
%c	L		  : Comprimento do dispositivo                             [cm]
%c	U10       : Tensao de polarizacao do MZM superior                   [V]
%c	U20       : Tensao de polarizacao do MZM inferior                   [V]
%c  eletrodes : Numero de eletrodos (1,2). Definido pelo numero de 
%c              variaveis de entrada na funcao Mach_Zehnder_Modulator     
%c  U_pi1     : Tensao de chaveamento para o eletrodo 1                 [V]
%c  U_pi2     : Tensao de chaveamento para o eletrodo 2                 [V]
%c  eta1      : Sensibilidade no caminho 1                          [1/V.m]
%c  eta2      : Sensibilidade no caminho 2                          [1/V.m]
%c  V1pi0	  : Constante da Tensao de chaveamento para o eletrodo 1    [V]
%c  V2pi0	  : Constante da Tensao de chaveamento para o eletrodo 2    [V]
%c  Vphi0	  : Constante da Tensao de chaveamento para controle de fase[V]
%c  nopt	  : Indice de refracao optico                               [-]
%c  nel	      : Indice de refracao eletrico                             [-]
%c  alfa_ins  : Perda por insercao                                     [dB]
%c  phi_0 	  : Constante de fase entre os dois caminhos                [-]
%c  C         : Parametro de chirp                                      [?]
%c  alfa0     : Perda condutiva                             [dB/cm.GHz^0.5]
%c
%c
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                             S T R U C T S                              %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c
%%
%% Initialize variables.
if nargin<2                                                                %Check if there is a specific laction to look for the files
    filename = [pwd '\input_files\AA_MZ_Input_Data_t' num2str(FileIndex)...%If not assume that there is a input_files folder
                                                                     '.m'];
else
    filename = [Local 'AA_MZ_Input_Data_t' num2str(FileIndex) '.m'];       %If yes, open the files on the given path
end
delimiter = ' ';

%% Format string for each line of text:
%   column1: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, ...
                     'MultipleDelimsAsOne', true,  'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
VarName1 = dataArray{:, 1};
                 
L          = VarName1(1);
U10        = VarName1(2);
U20        = VarName1(3);
U_pi1      = VarName1(4);
U_pi2      = VarName1(5);
eta1       = VarName1(6);
eta2       = VarName1(7);
V1pi0      = VarName1(8);
V2pi0      = VarName1(9);
Vphi0      = VarName1(10);
nopt       = VarName1(11);
nel        = VarName1(12);
alfa_ins   = VarName1(13);
phi_0      = VarName1(14);
alfa0      = VarName1(15);

C     = (eta1+eta2)/(eta1-eta2); %  Parametro de chirp
alfa0 = 10^(alfa0/20); 			 % [1/cm.GHz^0.5]


%% Clear temporary variables
clearvars filename delimiter formatSpec fileID dataArray ans;