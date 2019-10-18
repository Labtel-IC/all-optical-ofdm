%c
%c                                                       ..'Ž`'..'Ž`..'Ž`..                                                   
%c       File:  BestConfMat
%c(Main File to find the COMB with the selected characteristcs)
%c
%c     This main code is resposible to call and run the all the fuction to
%c run this simulation. Here it is possible to change any configuration
%c that was previouly stated on the Input data file of this simulation.
%c
%c      Eout = Ein*cos( (phi_1-phi_2)/2 )*e^[j*( (phi_1+phi_2)/2 )]
%c                              
%c
%c          phi_1 = (pi/V_pi)*V*sin(2*pi*f*t)
%c
%c          phi_1 = (pi/V_pi)*V*sin(2*pi*f*t)
%c
%c                                           by P.Marciano LG
%c                                           08/09/2017
%c                                           pablorafael.mcx@gmail.com
%c 
%c     References:
%c@article{kim2002chirp,
%c  title={Chirp characteristics of dual-drive. Mach-Zehnder modulator with a finite DC extinction ratio},
%c  author={Kim, Hoon and Gnauck, Alan H},
%c  journal={IEEE Photonics Technology Letters},
%c  volume={14},
%c  number={3},
%c  pages={298--300},
%c  year={2002},
%c  publisher={IEEE}
%c}
%c
%c@phdthesis{togneri2005analise,
%c  title={An{\'a}lise de Sistemas de Multiplexa{\c{c}}{\~a}o por Subportadora-SCM},
%c  author={Togneri, Arnaldo Paterline},
%c  year={2005},
%c  school={UNIVERSIDADE FEDERAL DO ESP{\'I}RITO SANTO}
%c}
%c
%c@article{oliveiralarge,
%c  title={Large Signal Analysis of Mach-Zehnder Modulator Intensity Response in a Linear Dispersive Fiber},
%c  author={Oliveira, JMB and Salgado, HM and Rodrigues, MRD}
%c}
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
%c   RFS_DSB_INPUT_DATA. - Keeps and create all the Variables and the
%c                         Vast majority of signal that will be used
%c                         .Does not have input parameters
%c
%c
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Start of the Program                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear ; clc; close all;
InputDataBestConfMat;
load([pwd '\save_files\DpmzmNoRing' '15-Oct-2017' '.mat']);          %Load the result of the past simulation for evaluation
% load([pwd '\save_files\TwoArmNoLoopFc1G_' '13-Oct-2017' '.mat']);          %Load the result of the past simulation for evaluation
% load([pwd '\save_files\' 'AvaliandoRectFilt.mat']);                      %Load the result of the past simulation for evaluation
%%
% Local = [pwd '\input_files\'];%variable that shows where the input files for
                                                     %the mzm will be saved
display('Making Files...');
Make_MZ_Input_Files_DP;    %Create the files for the MZM input data
display('Files Ready');
drawnow;            

    S=1;                     %variable to control the different Simulations
                             %Setings for comparation. For example if it
                             %want to simulate just the sistem response
                             %with Vbias = 2.8 [vol] and phase = 90 [deg]
                             %the value of S will be 1 through the
                             %simulation. But if you want to compare
                             %different Setings such as Vbias = 2.8 and
                             %Vbias = 3.8 with phase = 90. Then S will
                             %first assume the value 1 for the first seting
                             %then S will have the value 2 for the secont
                             %seting. The variables will be a matix with N
                             %lines as the number of comparisons that you
                             %want to do. For instance in this last case
                             %the Eout will have two lines the first will
                             %store the simulation result of Eout for Vbias
                             %=2.8 phase 90 and a second line for the
                             %result of Eout for Vbias = 3.8 and phase = 90
%% Storing the Results in a matrix for further evaluation
if exist('SFName','var')
    if AvaRecFil
        load([pwd '\save_files\' 'AvaRecFil' SFName '.mat']);              %Load the results of the previous simulation
    else
        load([pwd '\save_files\' 'Ava' SFName '.mat']);                    %Load the results of the previous simulation
    end
else
    if AvaRecFil
        load([pwd '\save_files\' 'AvaliandoRecFil.mat']);                  %Load the results of the previous simulation
    else
        load([pwd '\save_files\' 'Avaliando.mat']);                        %Load the results of the previous simulation
    end
end
MatPos = 1;
for kk = melhorespontos
    BesResMat(1,MatPos)  = To_Eval{kk}.Ome/(2*pi);
    BesResMat(2,MatPos)  = To_Eval{kk}.Pha;
    BesResMat(3,MatPos)  = To_Eval{kk}.VA2;
    BesResMat(4,MatPos)  = To_Eval{kk}.VA1;
    BesResMat(5,MatPos)  = To_Eval{kk}.V1b;
    BesResMat(6,MatPos)  = To_Eval{kk}.V2b;
    BesResMat(7,MatPos)  = To_Eval{kk}.VPh;
    BesResMat(8,MatPos)  = To_Eval{kk}.Gai;
    BesResMat(9,MatPos)  = length(To_Eval{kk}.Pea);
    BesResMat(10,MatPos) = BestCentDiff(MatPos);
    BesResMat(11,MatPos) = BestSideDiff(MatPos);
%     BesResMat(8,MatPos) = BestPlanicity(MatPos);
    MatPos = MatPos + 1;
end
MatPos = MatPos - 1;
%Location and name of the input file to be created

%% main logic to adress data on excel files
%Those next following lines cost me 4 days to develop. Here is the main
%idea of how to get any number and convert it to Excel alpha numeric collum
%indexation. The idea is to change the base of the given number. Let us
%assum that the number given to be converted is 10 based which means it
%varies from 0 to 9 (ten numbers in total). Now the Brazilian alphabet,
%which is based on the Roman alphabet added the k Y W letters. It is
%composed of 26 different characters. Therefore it can be used as an number
%26 besed which means it varies from A to Z (26 numbers in total). The
%point here is to conver a number 10 based in another number 26 based. The
%main idea is to divide the given number, for example 824, by 26 while the
%quotient be greater than 1. The new number 26 based will be composed by
%the rest of the division, for instance:
%               824/26 = 31 with 18 of rest
%The quotient without the rest is again divided by 26:
%                31/26 = 1 with 5 of rest
%The quotient without the rest is again divided by 26:
%                 1/26 = 0 with 1 of rest
%Finally, the new number 26 based is then 1 5 18. But in 10 based number
%the number 1 can be writen as 01 or yet 00001. We for convention don't
%represent the number zero at left. But in the 26 based number (represented
%by letter) the 0 is A. Therefore if we wanted to represent 1 in 26 based
%number we needed to AAAAAAAAAAAAAAAAAAA...B. It would be needed infinit
%numbers As before the number 1 10 based with is number B in 26 based. For
%convinience one unity will be remove from each part of the new number for
%better representation. Thus, the new number 26 based will be 0 4 17.
%Therefore, when converted to the 26 based notation (Alphabetic) it will be
%AER. If you need to check if the convertion is correct you just need to
%solve the equation:
%   ((nth10based + 1)*26^(nth-1)+...(1th10based + 1)*26^(1-1))
%For instance let us take the number AER for exemple:
%A is 0, E is 4 and R is 17 thus:
% (0+1)*26^2  + (4+1)*26^1  +  (17+1)*26^0 = 676 + 130 + 18 = 824 !!!
%As we can see the convertion was succefull.
OffSet = 1;                                                                %Offset from the first positon( starting from B or C or D... and so on)
aux = MatPos - (1-OffSet);                                                 %Variable that will recive the given number for convertion
done = 1;                                                                  %Flag to show that the Convertion is done
counting = [];                                                             %Vector to store the result from convertion
while done                                                                 %Main loop
    counting = [mod(aux,26) counting];                                     %For each interaction the rest from the division is sotred
    aux2 = (aux - mod(aux,26))/26;                                         %The rest is removed from the number and then divided for the next interaction
    if aux2 < 1                                                            %Check if the quotient is less than 1
        done =0;                                                           %If it is the convertion get to an end
    end                                                                    %If not the process continue
    aux = aux2-1;                                                          %The next value is removed one unity for better representation (as previously explained)
end
Index = char(counting+65);                                                 %Converting the new number to alphabetic pattern
if AvaRecFil                                                               %Selecting the name of the file to be saved
    Name = [pwd '\save_files\' 'RecFilBestConf' SFName '.xlsx'];
else
    Name = [pwd '\save_files\' 'BestConf' SFName '.xlsx'];
end
sheet = 1;                                                                 %Selecting the sheet in the Excel file that will receive the data
xlRange = [char((OffSet)+65) num2str(LinSta) ':' Index num2str((LinSta- ...
                                                    1)+size(BesResMat,1))];                              %Saving the range in which this data will be saved
xlswrite(Name,BesResMat,sheet,xlRange);                                    %Writing the data in the Excel file
%%
a=1;
