%c
%c                                                       ..'Ž`'..'Ž`..'Ž`..                                                   
%c       File: Make_MZ_Input_Files_Simp
%c(Creating the input file necessary for the MZM)
%c
%c     This main code is resposible to call and load the right parameters
%c to be saved in a file for latter be used by the MZM
%c
%c
%c                                           by P.Marciano LG
%c                                           18/08/2017
%c                                           pablorafael.mcx@gmail.com
%c 
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                    G L O B A L  V A R I A B L E S                      %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c   
%c
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                   L I S T  O F  V A R I A B L E S                      %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c   
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                             S T R U C T S                              %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c   Here it will be added the functions that this code will call
%c   
%c
%c   Set_MZ_Input_Data. - Function that will create and load the files with
%c the given in data.
%c
%c
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Start of the Program                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('L','var')
    MZ_Input_Data_DP;                                                      %Loading the basic data to be saved
end
S=1;                                                                       %Variabel that give the name for the input file
for Inc0=1:phi_0_steps 
    for Inc1=1:V1bias_steps                                                %Secundary loop for the possible combinations
        for Inc2=1:V2bias_steps                                            %Secundary loop for the possible combinations
            phi_0 = phi_0_vet(Inc0);
            U10               = V1bias(Inc1);                              %Variation of the Vbias
            U20               = V2bias(Inc2);                              %Variation of the Vbias
            [MZ_Input_Sdata] = Set_MZ_Input_Data_DP(S,L,U10,U20,U_pi1,...  %Functionresponsible to create the input files
            U_pi2,eta1,eta2,V1pi0,V2pi0,Vphi0,nopt,nel,alfa_ins,phi_0,...
                                                           alfa0,Local_Dp);
            S                = S+1;
        end
    end
end
% clear L U0 U_pi1 U_pi2 eta1 eta2 nopt nel alfa_ins phi_0 alfa0 C alfa0