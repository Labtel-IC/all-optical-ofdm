%c
%c                                                       ..'Ž`'..'Ž`..'Ž`..                                                   
%c       File: Make_MZ_Input_Files 
%c(Creating the input file necessary for the MZM)
%c
%c     This main code is resposible to call and load the right parameters
%c to be saved in a file for latter be used by the MZM
%c
%c
%c                                           by P.Marciano LG
%c                                           18/08/2018
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
MZ_Input_Data; %Loading the basic data to be saved

for Inc=1:phase_steps%Main loop for all possible combinations of the input
    for Inc2=1:Vbias_steps%Secundary loop for the possible combinations
        U0               = Vbias(Inc2);%Variation of the Vbias
        [MZ_Input_Sdata] = Set_MZ_Input_Data(S,L,U0,U_pi1,U_pi2,eta1,...   %Functionresponsible to create the input files
                                 eta2,nopt,nel,alfa_ins,phi_0,alfa0,Local);
        S                = S+1;
    end
end
% clear L U0 U_pi1 U_pi2 eta1 eta2 nopt nel alfa_ins phi_0 alfa0 C alfa0