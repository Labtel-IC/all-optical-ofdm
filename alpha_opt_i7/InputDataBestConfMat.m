%c
%c                                                       ..'Ž`'..'Ž`..'Ž`..                                                   
%c       File: InputDataBestConfMat 
%c
%c       This code is used to set up the input parameters that will be used
%c in the whole simualation.
%c
%c
%c
%c
%c
%c                                           by P.Marciano LG
%c                                           13/10/2017
%c                                           pablorafael.mcx@gmail.com
%c 
%c
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                    G L O B A L  V A R I A B L E S                      %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                   L I S T  O F  V A R I A B L E S                      %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c
%c
%c
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                             S T R U C T S                              %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c   Here it will be added the functions that this code will call
%c   
%c
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Start of the Program                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      OneArmLoopFc10G_   22-Oct-2017
load([pwd '\save_files\OneArmLoopFc10G_M10' '24-Oct-2017' '.mat']);           %Load the result of the past simulation for evaluation

AvaRecFil = 0;%Set this variable to 1 if the file to be evaluate is from the selection with the rectangular filter
LinSta = 1;

if exist('SFName','var')
    if AvaRecFil
        LocalSaved = [pwd '\save_files\' 'MatResRectFil' SFName];          %Local where the results will be saved
    else
        LocalSaved = [pwd '\save_files\' 'MatRes' SFName];                 %Local where the results will be saved
    end
else
    if AvaRecFil
        LocalSaved = [pwd '\save_files\' 'MatrizResultadosRectFil'];       %Local where the results will be saved
    else
        LocalSaved = [pwd '\save_files\' 'MatrizResultados'];              %Local where the results will be saved
    end
end