%c
%c                                                       ..'Ž`'..'Ž`..'Ž`..                                                   
%c       File: CharacAnaliticMzmDp
%c(Main File to Run the scrip to characterize the MZM_DP function)
%c
%c     This main script is used to verify the analitical equation for the
%c Dual Parallel Mach Zehnder Modulator and the results from the laboratory
%c with the same device.
%c
%c
%c                                           by P.Marciano LG
%c                                           14/10/2017
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Start of the Program                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear ; clc; close all;
DPMZM_Expe_Resul;       %Load the results from laboratory for comparison. 
                        %OBS.: this line of code must be called before the 
                        %inicializatino of any variable because in it are 
                        %place the comands clc, close anc clear from matlab

%% First Instance: Responsible to generate the basic signal and constants
% This is the first point of the program whre the file that loads the main
% configuration is loaded. This file is used to make the major definitions 
% needed for the simulation. Therefore, whatever modification needed will
% be done inside of it. No other part of this code must be changed. If you 
% do, please save this file with a different name.
InputDataCharcAnaliticMzmDp;

if ~exist([LocSav SFName '.mat'],'file')
    %% Second Instance: Calling main Function to Execute the Simulation
    display('Making Files...');
    Make_MZ_Input_Files_DPV1;    %Create the files for the MZM input data
    Make_MZ_Input_Files_DPV2;    %Create the files for the MZM input data
    Make_MZ_Input_Files_DPPh;    %Create the files for the MZM input data
    display('Files Ready');
    drawnow;                %It needs to be done because the Matlab do not
                            %update work space automaticaly

        S=1;                %variable to control the different Simulations
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

%% Loop For the Variation on Arm 1
%Fristly, it will vary the voltage on the Bias voltage for the Intensity
%MZM located at the up arm of the DP-MZM. Also, the others voltages (phase
%valtage and Bias voltage on the second arm) will be at zero.

    for kk=1:V1bias_steps%and a Vbais
    %% Loading variables for testing arm 1 and others in zerro
        EleSig2 = zeros(1,length(t));                                      %RF signal at the second arm
        EleSig1 = zeros(1,length(t));                                      %RF signal at the first arm
        EleSig.U1t =EleSig1;                                               %There will not have any signal going to the MZM model
        EleSig.U2t =EleSig2;                                               %There will not have any signal going to the MZM model
        MZ_Input_Sdata = (kk);                                             %The Index for the file to be loaded by the MZM
    %% Mach Zhender Modulator
       [EoutA1,~]=Mach_Zehnder_Modulator_DP(t,CW,EleSig,MZ_Input_Sdata,...
                                                                  LocalV1);%Script to implement the DP-MZM
    %% End of the current Interaction
       power_A1(kk) = mean(abs(EoutA1).^2);                                %The avarege power of the output signal
       clc;display(kk*100/(V1bias_steps),'%');                                 %the simulation still runing.
    end
%% Loop For the Variation on Arm 2
%Secondly, it will vary the voltage on the Bias voltage for the Intensity
%MZM located at the Lower arm of the DP-MZM. Also, the others voltages (
%phase  valtage and Bias voltage on the first arm) will be at zero.
for kk=1:V2bias_steps%and a Vbais
    %% Loading variables for testing arm 1 and others in zerro
        EleSig2 = zeros(1,length(t));                                      %RF signal at the second arm
        EleSig1 = zeros(1,length(t));                                      %RF signal at the first arm
        EleSig.U1t =EleSig1;                                               %There will not have any signal going to the MZM model
        EleSig.U2t =EleSig2;                                               %There will not have any signal going to the MZM model
        MZ_Input_Sdata = (kk);                                             %The Index for the file to be loaded by the MZM
    %% Mach Zhender Modulator
       [EoutA2,~]=Mach_Zehnder_Modulator_DP(t,CW,EleSig,MZ_Input_Sdata,...
                                                                  LocalV2);%Script to implement the DP-MZM
    %% End of the current Interaction
       power_A2(kk) = mean(abs(EoutA2).^2);                                %The avarege power of the output signal
       clc;display(kk*100/(V2bias_steps),'%');                                 %the simulation still runing.
    end
%% Loop For the Variation on Phase
%Finally, it will vary the voltage on the Phase controler among each arm of
%the DP-MZM (between Up and Lower arm) . Also, the others voltages (Bias
%valtage on the first and second arm) will be at zero.
    for kk=1:phi_0_steps
    %% Loading variables for testing Phase and others in zerro
        EleSig2 = zeros(1,length(t));                                      %RF signal at the second arm
        EleSig1 = zeros(1,length(t));                                      %RF signal at the first arm
        EleSig.U1t =EleSig1;                                               %There will not have any signal going to the MZM model
        EleSig.U2t =EleSig2;                                               %There will not have any signal going to the MZM model
        MZ_Input_Sdata = (kk);                                             %The Index for the file to be loaded by the MZM
    %% Mach Zhender Modulator
       [EoutPh,~]=Mach_Zehnder_Modulator_DP(t,CW,EleSig,MZ_Input_Sdata,...
                                                                   LocalPh);%Script to implement the DP-MZM
    %% End of the current Interaction
       power_Ph(kk) = mean(abs(EoutPh).^2);                                %The avarege power of the output signal
       clc;display(kk*100/(phi_0_steps),'%');                                  %the simulation still runing.
    end
    %% Acquiering Results from Laboratory
    PoutArm1AveNor = (mat_dp_mzm(:,PotRis1_mW)'+mat_dp_mzm(:,PotFal1_mW...
                                                                    )')./2;%Output power result from variation on Vbias 1
    PoutArm1AveNor = PoutArm1AveNor./max(PoutArm1AveNor);

    PoutArm2AveNor = (mat_dp_mzm(:,PotRis2_mW)'+mat_dp_mzm(:,PotFal2_mW...
                                                                    )')./2;%Output power result from variation on Vbias 2
    PoutArm2AveNor = PoutArm2AveNor./max(PoutArm2AveNor);

    PoutFaseAveNor = (mat_dp_mzm(:,PotRisFas_mW)' + mat_dp_mzm(:,...
                                                        PotFalFas_mW)')./2;%Output power result from variation on Fase
    PoutFaseAveNor = PoutFaseAveNor./max(PoutFaseAveNor);
    Vbias_Exp = mat_dp_mzm(:,Vbias_Index)';
    %% End of Script
    save([LocSav SFName]);                                                 %When a simulation is finished all the data it is saved.
    else
    load([LocSav SFName '.mat']);
end
%% Ploting results for Comparison
figure;
hold on;
grid on;
plot(V1bias,power_A1./max(power_A1),'-','color',[0 0.6 0],'LineWidth',2);
plot(Vbias_Exp,PoutArm1AveNor,'*','color',[0 0.6 0],'LineWidth',2);
plot(V2bias,power_A2./max(power_A2),'-','color',[0 0.4 0.8],'LineWidth',2);
plot(Vbias_Exp,PoutArm2AveNor,'o','color',[0 0.4 0.8],'LineWidth',2);
plot(V2bias,power_Ph./max(power_Ph),'-','color',[1 0.8 0],'LineWidth',2);
plot(Vbias_Exp,PoutFaseAveNor,'d','color',[1 0.8 0],'LineWidth',2);
title('Potência de Saída Vs Tensão de Polarização','FontSize',20);
xlabel('Tensão de Polarização [V]','FontSize',16);
ylabel('Potência de Saída [W]','FontSize',16);
legend({'Vbias1-Sim','Vbias1-Exp','Vbias2-Sim','Vbias2-Exp','Vfase-Sim',...
'Vfase-Exp'},'FontSize',14,'Location','best','FontWeight','bold','Box',...
                                                                    'off');

