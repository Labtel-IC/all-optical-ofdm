%c
%c                                                       ..'Ž`'..'Ž`..'Ž`..                                                   
%c       File: Caracte_MZM_Model
%c(Main File to find the Characteristical curve of the MZM model)
%c
%c     Each MZM device has its own characteristics. At this code we aim to
%c understand the bahaviour of the analitical model that we have and to
%c evaluate if it matches with the experiment developed in laboratory for
%c further comparison with an practical experiment.
%c
%c Where the analitical expression is givenn by:
%c      Eout = Ein*cos( (phi_1-phi_2)/2 )*e^[j*( (phi_1+phi_2)/2 )]
%c                              
%c
%c          phi_1 = (pi/V_pi)*V*sin(2*pi*f*t)
%c
%c          phi_1 = (pi/V_pi)*V*sin(2*pi*f*t)
%c
%c in addition the theoritical equation is given by:
%c Eout_theo = sqrt(Po).*cos((pi/2).*(-Vbias+2*S(t)+AV1)./(Vpi1))
%c
%c                                           by P.Marciano LG
%c                                           29/09/2017
%c                                           pablorafael.mcx@gmail.com
%c 
%c     References:
%c@article{pereira2015impact,
%c  title={Impact of optical power in the guard-band reduction of an optimized DDO-OFDM system},
%c  author={Pereira, Esequiel da V and Helder, R de O and Nunes, Reginaldo B and Segatto, Marcelo EV and Silva, Jair AL},
%c  journal={Journal of Lightwave Technology},
%c  volume={33},
%c  number={23},
%c  pages={4717--4725},
%c  year={2015},
%c  publisher={IEEE}
%c}
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
%c   Eout_d     :Analitical output of the MZM                           [?]
%c   U_pi1      :Polarization characteristic voltage of the MZM (Vpi)   [V]
%c   Eout_theo  :Theoretical output of the MZM                          [?]
%c   power_theo :Power of the signal Eout_theo                          [W]
%c   power_Aved :Power of the signal Eout_d                             [W]
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                             S T R U C T S                              %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c   Here it will be added the functions that this code will call
%c   
%c
%c
%c
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Start of the Program                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear ; clc; close all; %Cleaning all variables and close all extra windows
                        %for a fresh start of the program
DPMZM_Expe_Resul;       %Load the results from laboratory for comparison. 
                        %OBS.: this line of code must be called before the 
                        %inicializatino of any variable because in it are 
                        %place the comands clc, close anc clear from matlab
SDMZM_Expe_Resul;
SDMZM10_Expe_Resul;
inputdata_Caracte_MZM_Model;

%% Criating input files for Simulation

display('Making Files...');
Make_MZ_Input_Files_Simp;  %Create the files for the MZM input data
display('Files Ready');
drawnow;                   %It needs to be done because the Matlab do not
                           %update work space automaticaly

%% Main Loop
if ~exist([pwd '\save_files\' 'Caract_MZM' date '.mat'],'file')
    for kk=1:Vbias_steps %for each Vbias it is measured the output power   
    %% main logic to adress data on excel files
    %     aux = kk - 1;
    %     done = 1;
    %     counting = [];
    %     while done
    %         counting = [mod(aux,26) counting];
    %         aux2 = (aux - mod(aux,26))/26;
    %         if aux2 < 1
    %             done =0;
    %         end
    %         aux = aux2-1;
    %     end
    %     Index = char(counting+65);
    %%
        [Eout_d,~] = Mach_Zehnder_Modulator_simplificado(t,CW,EleSig1,kk); %This is to call the right MZM-Input-Data for 
                                                                           %the simulation it's name vary with values of 
                                                                           %phase and Vbias
        [Eout_E,~] = Mach_Zehnder_Modulator_Ezequiel(t,CW,EleSig1,kk);     %This is to call the right MZM-Input-Data for 
                                                                           %the simulation it's name vary with values of 
                                                                           %phase and Vbias
        [Eout_M,~] = Mach_Zehnder_Modulator_Mach10(t,CW,EleSig1,kk);     %This is to call the right MZM-Input-Data for 
                                                                           %the simulation it's name vary with values of 
                                                                           %phase and Vbias
        [~,~,U_pi1,~,~,~,~,~,~,~,~,~] = Import_Data_MZM (kk);              %Loads just the polarization voltage for the 
                                                                           %current interaction. This value will be used 
                                                                           %on the theoretical equation
        Eout_theoE = sqrt(Po).*cos((pi/2).*(-Vbias(kk)+2*EleSig1+AV1)./...
                                                                  (4.7));%Electical output field of the MZM founded by 
                                                                           %the theoretical equation
        power_theoE(kk) = mean(abs(Eout_theoE).^2);                          %Theoretical_Average Power of optical signal
        
        Eout_theo = sqrt(Po).*cos((pi/2).*(-Vbias(kk)+2*EleSig1+0)./...
                                                                  (U_pi1));%Electical output field of the MZM founded by 
                                                                           %the theoretical equation
        power_theo(kk) = mean(abs(Eout_theo).^2);                          %Theoretical_Average Power of optical signal
        power_Aved(kk) = mean(abs(Eout_d).^2);                             %Average Power of optical signal
        power_AveE(kk) = mean(abs(Eout_E).^2);                             %Average Power of optical signal
        power_AveM(kk) = mean(abs(Eout_M).^2);                             %Average Power of optical signal
        U_pi1 = 3.9;
        Eout_theoDpPha = sqrt(Po).*cos((pi/2).*(-Vbias(kk)+2*EleSig1-0.14)./...
                                                                  (U_pi1));%Electical output field of the MZM founded by 
                                                                           %the theoretical equation
        U_pi1 = 4.3;
        Eout_theoDpVb1 = sqrt(Po).*cos((pi/2).*(-Vbias(kk)+2*EleSig1-0.4)./...
                                                                  (U_pi1));%Electical output field of the MZM founded by 
                                                                           %the theoretical equation
        U_pi1 = 4.75;
        Eout_theoDpVb2 = sqrt(Po).*cos((pi/2).*(-Vbias(kk)+2*EleSig1-1.2)./...
                                                                  (U_pi1));%Electical output field of the MZM founded by 
        U_pi1 = 7.4;
        Eout_theoM10 = sqrt(Po).*cos((pi/2).*(-Vbias(kk)+2*EleSig1+4.1)./...
                                                                  (U_pi1));%Electical output field of the MZM founded by 
                                                                           %the theoretical equation
%         power_theo(kk) = mean(abs(Eout_theo).^2);                          %Theoretical_Average Power of optical signal
        power_theoDpPha(kk) = mean(abs(Eout_theoDpPha).^2);                          %Theoretical_Average Power of optical signal
        power_theoDpVb1(kk) = mean(abs(Eout_theoDpVb1).^2);                          %Theoretical_Average Power of optical signal
        power_theoDpVb2(kk) = mean(abs(Eout_theoDpVb2).^2);                          %Theoretical_Average Power of optical signal
        power_theoM10(kk) = mean(abs(Eout_theoM10).^2);                          %Theoretical_Average Power of optical signal
        power_Aved(kk) = mean(abs(Eout_d).^2);                             %Average Power of optical signal
        clc;(kk)*100/Vbias_steps
    end
else
    load([pwd '\save_files\' 'Caract_MZM' date '.mat']);
end
%% Acquiering Results from Laboratory
PoutArmAveNor = (mat_sd_mzm(:,PotRis_mW)');                                %(mat_sd_mzm(:,PotRis_mW)'+mat_sd_mzm(:,PotFal_mW)')./2; %Output power result from variation on Vbias 1
PoutArmAveNor = PoutArmAveNor./max(PoutArmAveNor);

PoutArmAveNorM10 = (mat_mach10(:,M10PotRis1_mW)' + mat_mach10(:,...
M10PotFal1_mW)'+ mat_mach10(:,M10PotRis2_mW)' + mat_mach10(:,...
                                                       M10PotFal2_mW)')./4; %Output power result from variation on Vbias 1
PoutArmAveNorM10 = PoutArmAveNorM10./max(PoutArmAveNorM10);

PoutArm1AveNor = (mat_dp_mzm(:,PotRis1_mW)'+mat_dp_mzm(:,PotFal1_mW)')./2; %Output power result from variation on Vbias 1
PoutArm1AveNor = PoutArm1AveNor./max(PoutArm1AveNor);

PoutArm2AveNor = (mat_dp_mzm(:,PotRis2_mW)'+mat_dp_mzm(:,PotFal2_mW)')./2; %Output power result from variation on Vbias 2
PoutArm2AveNor = PoutArm2AveNor./max(PoutArm2AveNor);

PoutFaseAveNor = (mat_dp_mzm(:,PotRisFas_mW)' + mat_dp_mzm(:,...
                                                        PotFalFas_mW)')./2;%Output power result from variation on Fase
PoutFaseAveNor = PoutFaseAveNor./max(PoutFaseAveNor);
Vbias_Exp = mat_dp_mzm(:,Vbias_Index)';
%% Ploting results for evaluation

figure;
hold on;
%% Result from EZ MZM
plot(Vbias_Exp,PoutArmAveNor,'-*','color',[0 1 0.4],'LineWidth',2);%experiment
plot(Vbias,power_AveE,'-','color',[0 1 1],'LineWidth',2);%simulation
plot(Vbias,power_theoE,'-','color',[1 0 0],'LineWidth',2);%theoretical
%% Result from the model in study
plot(Vbias,power_Aved,'-o','color',[0 0 0],'LineWidth',2);%simulation
plot(Vbias,power_theo,'-','color',[1 0 0],'LineWidth',2);%theoretical
%% Result From the DP-MZM
%Phase
plot(Vbias_Exp,PoutFaseAveNor,'-d','color',[0 0.8 0.8],'LineWidth',2);%experiment
plot(Vbias,power_theoDpPha,'-','color',[1 0 0],'LineWidth',2);%theoretical
%Vbias 1
plot(Vbias_Exp,PoutArm1AveNor,'-.','color',[0 0.8 0],'LineWidth',2);%experiment
plot(Vbias,power_theoDpVb1,'-','color',[1 0 0],'LineWidth',2);%theoretical
%Vbias 2
plot(Vbias_Exp,PoutArm2AveNor,'--','color',[0.4 0 0.4],'LineWidth',2);%experiment
plot(Vbias,power_theoDpVb2,'-','color',[1 0 0],'LineWidth',2);%theoretical
%% Results from the MZM 10
plot(Vbias_Exp,PoutArmAveNorM10,'-*','color',[0 1 0.1],'LineWidth',2);%experiment
plot(Vbias,power_AveM,'-','color',[1 1 0],'LineWidth',2);%simulation
plot(Vbias,power_theoM10,'-','color',[1 0 0],'LineWidth',2);%theoretical

title('Potência de Saída Vs Tensão de Polarização','FontSize',20);
xlabel('Tensão de Polarização [V]','FontSize',16);
ylabel('Potência de Saída [W]','FontSize',16);
legend({'SD-MZM-EZ-Exp','SD-MZM-EZ-Sim','SD-MZM-The','MZM-Module-Sim',...
'MZM-Module-The','DP-MZM-Pha-Exp','DP-MZM-Pha-Theo','DP-MZM-Vb1-Exp',...
'DP-MZM-Vb1-The','DP-MZM-Vb2-Exp','DP-MZM-Vb2-The','MZM10-Exp',...
'MZM10-Sim','MZM10-The'},'FontSize',14,'Location','best',...
                                          'FontWeight','bold','Box','off');
grid on;
print(['-f' num2str(1)],['figura' num2str(1)],'-dpng','-r0');
a=1;
save([pwd '\save_files\' 'Caract_MZM' date]);                              %When a simulation is finished all the data is is saved.