%c
%c                                                       ..'�`'..'�`..'�`..                                                   
%c       File: GenerateOFC
%c(Main File to Run The Generator of Ultra Optical Flat Comb Source)
%c
%c     This main code is resposible to call and run the all the fuction to
%c related to this simulation. Here it is possible to change any 
%c configuration that was previouly stated on the Input data file of this 
%c simulation.
%c
%c      Eout = Ein*cos( (phi_1-phi_2)/2 )*e^[j*( (phi_1+phi_2)/2 )]
%c                              
%c
%c          phi_1 = (pi/V_pi)*V*sin(2*pi*f*t)
%c
%c          phi_1 = (pi/V_pi)*V*sin(2*pi*f*t)
%c
%c                                           by P.Marciano LG
%c                                           28/10/2017
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
%c  S     - This is a variable to store the overall of interaction of the 
%c          simulation.
%c
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                             S T R U C T S                              %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c   Here it will be added the functions that this code will call
%c   
%c
%c   MAIN_COMB_INPUT_DATA. - Keeps and create all the Variables and the
%c                         Vast majority of signal that will be used
%c                         .Does not have input parameters
%c
%c   Make_MZ_Input_Files_Simp - Creates the files that will be the input 
%c                        data to the Mach-Zehnder Modulator
%c                    .Does not have input parameters
%c
%c  Mach_Zehnder_Modulator_simplificado - MZM model for simulation
%c
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Start of the Program                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear ; clc; close all;

%% First Instance: Responsible to generate the basic signal and constants
% This is the first point of the program whre the file that loads the main
% configuration is loaded. This file is used to make the major definitions 
% needed for the simulation. Therefore, whatever modification needed will
% be done inside of it. No other part of this code must be changed. If you 
% do, please save this file with a different name.
GenerateOfcInputDataExp;
%% Second Instance: Calling main Function to Execute the Simulation
display('Making Files...');
Make_MZ_Input_Files_Simp;    %Create the files for the MZM input data
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
%% Third Instance: Generation of the COMB
% At this point this script will save the result of the DP-MZM output for
% each possible combination of variables previouly chosen.

%Main Loop
% if Ploting
%     figure;
%     hold all;
%     grid on;
%     plot(f,ML*log10(abs(fftshift(Filtro))));
% end

for Va1=1:Varm1_steps                                                      %For each possible value of Amplitude of the RF signal at the Uper arm
    for Va2=1:Varm2_steps                                                  %For each possible value of Amplitude of the RF signal at the Lower arm
        for Ph12=1:phase_steps                                             %For each Phase differency between the RF sgnals
            for Gai=1:Gain_steps                                           %For each gain given to the output optical signal
                for Inc0=1:Vbias_steps                                     %For each Phase differency between Uper and Lower arm
                                                                           %This main loop will run the simulation for every possible combination of those variables
                    %% Loading variables  
                    Gain        = Gain_vet(Gai);                           %The Gain That will be applied at the output of the MZM
                    phase_is    = phase(Ph12);
                    EleSig1     = Varm1(Va1)*sin(2*pi*fc*t);               %The first RF signal at the Uper arm.
                    EleSig2     = Varm2(Va2)*sin(2*pi*fc*t + phase_is);    %The second RF singal at the Lower arm.
                                                                           %will have a phase shift different of 
                                                                           %the first arm.
                    EleSig.U1t  = EleSig1;                                 %The MZM model require that the voltages
                    EleSig.U2t  = EleSig2;                                 %The MZM model require that the voltages
                    MzmInputInd = Inc0;                                    %The Index for the file to be loaded by the MZM
                %% Mach Zhender Modulator
%                    OFC_RIM_RFS;                                          %Script to implement the optical ring
%                    OFC_RIMExpIQ;                                             %Script to implement the optical ring
                   OFC_RIMExp;                                             %Script to implement the optical ring
                %% Finding the peaks
                    ponto_max =max(abs(fftshift(fft(Eout)./length(Eout))));%Before anything else it is necessary firt to 
                                                                           %know the maximum power value of the signal
                    limiar = ponto_max/Dif_Amp;                            %Finding the threshold
                    [peaks_aux,locs_aux] = findpeaks(abs(fftshift(...      %Now it is possible to find how many peaks the 
                         fft(Eout)./length(Eout))),'MinPeakHeight',limiar);%signal has within our first specification whith the Dif_Amp varialbe
                    peaks(S,1:length(peaks_aux)) = peaks_aux;              %After finding the number of peaks they are stored for latter evaluation.
                    locs(S,1:length(locs_aux)) = locs_aux;                 %Possitions of the frequency refent of the founded peaks
                    [~,U0,~,~,~,~,~,~,~,~,~,~] = Import_Data_MZM(Inc0);%Import the data from files
                    %% Saving main parameters of the simulation
                    To_Eval{S}.Ome = Rad_f;                                % Store the angular frequency of the signals generated
                    To_Eval{S}.Pha = phase_is;                             % Stores the phase difference between each RF signal
                    To_Eval{S}.VA1 = Varm1(Va1);                           % Stores the Amplitude of the RF signal on the first arm
                    To_Eval{S}.VA2 = Varm2(Va2);                           % Stores the Amplitude of the RF signal on the second arm
                    To_Eval{S}.Vbi = U0;                                   % Stores the bias voltage applied
                    To_Eval{S}.Pea = peaks_aux;                            % Stores the peaks founded on the simulation
                    To_Eval{S}.Loc = locs_aux;                             % Stored the location of the frequency respective to the peaks founded
                    To_Eval{S}.Fin = MzmInputInd;                          % Stored the Index number for the MZM configuration file
                    To_Eval{S}.Gai = Gain;                                 % Stored the Gain applied in the optical ring
                    %% End of the current interaction
                    S = S + 1;                                             %counter variable
                end
                clc;
                display((S-1)*100/(N_combination),'%');                    %the simulation still runing.
            end
        end
    end
end

% ThisFigure = 7;
% CriandoFiguras;
%% Saving Results and closing the simulation
save([LocSav SFName]);                                                     %When a simulation is finished all the data it is saved.

