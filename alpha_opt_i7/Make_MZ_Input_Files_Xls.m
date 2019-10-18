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
% MZ_Input_Data; %Loading the basic data to be saved
% phase_steps = 1000
% Vbias_steps = 1000
L          = 1;
U_pi1      = 4.8;
U_pi2      = 4.8;
eta1       = 0.5;
eta2       = -0.5;
nopt       =  0;
nel        =  0;
alfa_ins   =  0;
phi_0      =  0;
alfa0      =  0;
S=1; %Variabel that give the name for the input file
for Inc=1:phase_steps%Main loop for all possible combinations of the input
    for Inc2=1:Vbias_steps%Secundary loop for the possible combinations
        U0               = Vbias(Inc2);%Variation of the Vbias
        aux = S - 1;
        done = 1;
        counting = [];
        while done
            counting = [mod(aux,26) counting];
            aux2 = (aux - mod(aux,26))/26;
            if aux2 < 1
                done =0;
            end
            aux = aux2-1;
        end
        Index = char(counting+65);
        Write_MZ_Input_Data(Index,L,U0,U_pi1,U_pi2,eta1,eta2,nopt,nel,...
                                               alfa_ins,phi_0,alfa0,Local);
%         sprintf('For %i is %s',S,Index);
        S = S + 1;
        (S-1)*100/(phase_steps*Vbias_steps)
    end
end
% clear L U0 U_pi1 U_pi2 eta1 eta2 nopt nel alfa_ins phi_0 alfa0 C alfa0