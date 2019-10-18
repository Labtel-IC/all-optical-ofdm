%c
%c                                                       ..'�`'..'�`..'�`..                                                   
%c       File: MAIN_RFS_DSB 
%c(Main File to Run The Recirculating Frequency Shifting Double Sided Band)
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
%c                                           01/09/2016
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
%c  Sfile - This is a variable to store the name of a file that the data/or
%c          simulation will be Stored. It have by Default the value
%c          'Simulacao_1'
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
%c   Make_MZ_Input_Files - Creates the files that will be the input data
%c                         to the Mach-Zehnder Modulator
%c                    .Does not have input parameters
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
RFS_DSB_INPUT_DATA;
%% Second Instance: Calling main Function to Execute the Simulation
display('Making Files...');
Make_MZ_Input_Files;    %Create the files for the MZM input data
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
% At this point the program will make a comparison as setup earlier on the 
% the program. For every Phase set up it will simulate for all Combination
% ofVbias previouly set up. The number of simulations that it will run is
% given by (phase_steps*Vbias_steps)^2 because the voltage of each arm on 
% the MZM also will be changed.

%Main Loop
for Inc=1:(phase_steps)   %For EACH phase
    phase_is = phase(Inc);%Salved here to be used on the main Loop
    for Inc2=1:Vbias_steps%and a Vbais
        [MZ_Input_Sdata] = ['AA_MZ_Input_Data_t' ...
                                      num2str(Inc2 + Vbias_steps*(Inc-1))];%This is to call the right MZM-Input-Data for 
                                                                           %the simulation it's name vary with values of 
                                                                           %phase and Vbias
         for kk=1:Vbias_steps   %This loop is responsible for vary the 
                                %voltage at the second arm.
            for jj=1:Vbias_steps %While this loop on the First arm.
            % Electrical parameters
                EleSig2 = Vbias(kk)*sin(Rad_f*t + phase_is);               %The second arm will have a phase shift different of the first arm
                EleSig1 = Vbias(jj)*sin(Rad_f*t); %Set up of 1st arm
                EleSig.U1t=EleSig1;%The MZM model require that the voltages
                EleSig.U2t=EleSig2;%to be stores at an structure variable
            % Mach Zhender Modulator
                [Eout,H]=Mach_Zehnder_Modulator_simplificado(t,CW,...
                                                    EleSig,MZ_Input_Sdata);
            % Finding the peaks
                ponto_max_db = max(ML*log10(abs(fftshift(fft(Eout)./...    %Before anything else it is necessary firt to 
                                                          length(Eout)))));%know the maximum power value of the signal
                limiar_db = ponto_max_db - Dif_Amp;%Finding the threshold
                [peaks_aux,locs_aux] = findpeaks(ML*log10(abs(fftshift(... %Now it is possible to find how many peaks the 
                     fft(Eout)./length(Eout)))),'MinPeakHeight',limiar_db);%signal has within our first specification whith the Dif_Amp varialbe
                peaks(S,1:length(peaks_aux)) = peaks_aux;                  %After finding the number of peaks they are stored for latter evaluation.
                locs(S,1:length(locs_aux)) = locs_aux;                     %Possitions of the frequency refent of the founded peaks
                S = S + 1;%conter variable
            end
         end
        display(S*100/(phase_steps*Vbias_steps*phase_steps*Vbias_steps),...%Used to show that
                                                                      '%');%the simulation still runing.
    end
end
save('final');%When a simulation is finished all the data it is saved.
