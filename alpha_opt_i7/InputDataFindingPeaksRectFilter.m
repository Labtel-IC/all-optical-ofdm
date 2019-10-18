%c
%c                                                       ..'Ž`'..'Ž`..'Ž`..                                                   
%c       File: InputDataFindingPeaksRectFilter (Detection of only the interest 
%c point for evaluation)
%c
%c       This code is used to create the different types of possible
%c combinations with the variables of interest and apply then to a MZM with
%c one single arm and evaluate the number of peaks created and their
%c flatness.
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
load([pwd '\save_files\DpmzmNoRing' '15-Oct-2017' '.mat']);                %Load the result of the past simulation for evaluation
Ploting = 0;
% load([pwd '\save_files\TwoArmNoLoopFc1G_' '13-Oct-2017' '.mat']);          %Load the result of the past simulation for evaluation
%%   Creates tools for evaluation
MinCentDiff    = 20;
MaxSideDiff    = 10;
Best_Dif_Amp   = 40;                                                        %Variation of peaks power amplitude of the final COMB desired in [db]
Num_Subcari    = 1;                                                        %Minimun number of subcariers available at the COMB
Obs_Port       = 2 + 1;                                          %Number of subcariers observed for measuring the planicity of the comb.
                                                                           %Plus 1 because the center carrier also counts.
if exist('SFName','var')
    LocalSaved = [pwd '\save_files\' 'AvaRecFil' SFName];                  %Local where the results will be saved
else
    LocalSaved = [pwd '\save_files\' 'AvaliandoRecFil'];                   %Local where the results will be saved
end

FBW            = Obs_Port*f_RF;                                            %BandWidth of the Retangular Filter
fc_RetFil      = 0;                                                        %Central frequency for the retangular filter
[ Filtro_Sqr ] = Filtro_Retangular( FBW,fc_RetFil,f);                      %Creating the retangular filter