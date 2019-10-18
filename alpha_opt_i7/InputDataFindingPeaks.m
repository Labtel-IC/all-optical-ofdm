%c
%c                                                       ..'Ž`'..'Ž`..'Ž`..                                                   
%c       File: InputDataFindingPeaks (Detection of only the interest 
%c point for evaluation)
%c
%c       This code is used to ...
%c 
%c
%c
%c
%c
%c                                           by P.Marciano LG
%c                                           26/09/2017
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

%%   Loads the selected parameter to be evaluated
load([pwd '\save_files\OneArmLoopFc10G_M10' '24-Oct-2017' '.mat']);                %Load the result of the past simulation for evaluation
Ploting = 1;                                                               %Flag to indicate if the scrip must plot the results or not
% load([pwd '\save_files\TwoArmNoLoopFc1G_' '13-Oct-2017' '.mat']);          %Load the result of the past simulation for evaluation
%%   Creates tools for evaluation
Best_Dif_Amp    = 2;                                                      %Variation of peaks power amplitude of the final COMB desired in [db]
Dif_Amp         = Best_Dif_Amp + 1;                                        %For an initial analizes it is selected an higher variation of peaks
                                                                           %power amplitude of the final COMB desired in [db]
SubCarrPerGroup = 4;                                                       %Number of subcarriers to be considered inside an group 
                                                                           %for planiciti evaluation
MinObsGroup     = 2;                                                       %Minimum number of carriers group observed for further classification.
                                                                           %Plus 1 because the center carrier also counts.
if exist('SFName','var')
    LocalSaved = [pwd '\save_files\' 'Ava' SFName];                        %Local where the results will be saved
else
    LocalSaved = [pwd '\save_files\' 'Avaliando'];                         %Local where the results will be saved
end

%%    Filter for better evaluation of the planicity
FBW            = SubCarrPerGroup*f_RF;                                     %BandWidth of the Retangular Filter
fc_RetFil      = 0;                                                        %Central frequency for the retangular filter
[ Filtro_Sqr ] = Filtro_Retangular( FBW,fc_RetFil,f);                      %Creating the retangular filter


