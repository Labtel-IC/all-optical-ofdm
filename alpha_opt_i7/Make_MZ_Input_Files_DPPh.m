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
%c   Set_MZ_Input_Data_DP. - Function that will create and load the files
%c with the given in data.
%c
%c
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Start of the Program                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('L','var')
    MZ_Input_Data_DP; %Loading the basic data to be saved
end
% U_pi1      = 0.1;
% U_pi2      = 0.1;
S=1; %Variabel that give the name for the input file
for Inc0=1:phi_0_steps 
    phi_0 = phi_0_vet(Inc0);
    for Inc=1:phase_steps                                                  %Main loop for all possible combinations of the input
        U10              = 0;                                              %Variation of the Vbias
        U20              = 0;                                              %Variation of the Vbias
        [MZ_Input_Sdata] = Set_MZ_Input_Data_DP(S,L,U10,U20,U_pi1,U_pi2,...%Functionresponsible to create the input files
        eta1,eta2,V1pi0,V2pi0,Vphi0,nopt,nel,alfa_ins,phi_0,alfa0,LocalPh);
        S                = S+1;
    end
end