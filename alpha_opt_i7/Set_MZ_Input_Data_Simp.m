function [ FileName ] = Set_MZ_Input_Data( FileIdex,L,U0,U_pi1,U_pi2,...
                             eta1,eta2,nopt,nel,alfa_ins,phi_0,alfa0,Local)
%% %c                 Generate Input Files for MZM
%c
%c function [ FileName ] = Set_MZ_Input_Data( FileIdex,L,U0,U_pi1,U_pi2,...
%c                           eta1,eta2,nopt,nel,alfa_ins,phi_0,alfa0,Local)
%c
%c                                           by P.Marciano LG
%c                                           5/10/2017
%c                                           pablorafael.mcx@gmail.com
%c 
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                   L I S T  O F  V A R I A B L E S                      %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c
%c	Generator of Input File for the Modulador Mach-Zehnder
%c
%c  Input:
%c 	FileIdex    : Number of the file in the sequence given by variable S
%c	L		    : Comprimento do dispositivo [cm]
%c	U0          : Tensao de polarizacao [V]
%c  eletrodes	: Numero de eletrodos (1,2). Definido pelo numero de 
%c                variaveis de entrada na funcao Mach_Zehnder_Modulator     
%c  U_pi1       : Tensao de chaveamento para 1 eletrodo  [V]
%c  U_pi2       : Tensao de chaveamento para 2 eletrodos [V]
%c  eta1        : Sensibilidade no caminho 1  [1/V.m]
%c  eta2        : Sensibilidade no caminho 2  [1/V.m]
%c  nopt	    : Indice de refracao optico
%c  nel	        : Indice de refracao eletrico
%c  alfa_ins	: Perda por insercao [dB]
%c  phi_0 		: Constante de fase entre os dois caminhos
%c  alfa0		: Perda condutiva [dB/cm.GHz^0.5]
%c  Local       : Variav�l onde os arquivos ser�o salvos
%c
%c  Output:
%c  FileName    : Nome do arquivo gerado
%c
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                             S T R U C T S                              %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c
%%
if nargin==12 %Verify if the storage folder was passed to this function
    Local = [pwd '\'];%if not the files will be saved on the current folder
elseif nargin<12%Verify if there is all parameters needed
    error('Not enough input arguments for Set_MZ_Input_Data');
end

%Location and name of the input file to be created
Name = [Local 'AA_MZ_Input_Data_t' num2str(FileIdex) '.m'];
fid = fopen( Name, 'wt' );
ContTimeOut = toc;
while(fid<3)
    if ((toc-ContTimeOut)/60)>5
        break;
    end
    fid = fopen( Name, 'wt' );
end
%opening the file or creating it
%writing in the file the parameters passed to this function
fprintf( fid, '%6.2f\n',L);
fprintf( fid, '%6.3f\n',U0);
fprintf( fid, '%6.2f\n',U_pi1);
fprintf( fid, '%6.2f\n',U_pi2);
fprintf( fid, '%6.2f\n',eta1);
fprintf( fid, '%6.2f\n',eta2);
fprintf( fid, '%6.2f\n',nopt);
fprintf( fid, '%6.2f\n',nel);
fprintf( fid, '%6.2f\n',alfa_ins);
fprintf( fid, '%6.2f\n',phi_0);
fprintf( fid, '%6.2f\n',alfa0);

fclose(fid);%closing the file

FileName = ['MZ_Input_Data_t' num2str(FileIdex)];%Passing out the file
end