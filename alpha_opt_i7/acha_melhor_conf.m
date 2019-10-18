%c
%c                                                       ..'Ž`'..'Ž`..'Ž`..                                                   
%c       File: achar_melhor_conf 
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
%c                                           08/0/2017
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

                   %% First Instance of the Program Responsible to generate the basic signal and constants
RFS_DSB_INPUT_DATA;          %This is the file that loads the main
                             %configuration of this program. Use this file
                             %to make the major definitions of this
                             %simulation. But if it is needed you can
                             %change some parameters (just those related 
                             %with 'on' and 'off')

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
%% Generation of the COMB
% At this point the program will make a comparison as setup earlier on the 
% the program. For every Phase set up it will simulate for all Combination
% ofVbias previouly set up. The number of simulations that it will run is
% given by (phase_steps*Vbias_steps)^2 because the voltage of each arm on 
% the MZM also will be changed.
melhor_phi = [];%Best phase delay founded
melhor_vbi = [];%Best polarization voltage founded
melhor_va1 = [];%Best voltage applyed to the first arm founded
melhor_va2 = [];%Best voltage applyed to the second arm founded
melhor_inp = [];%Best MZM input data founded
load('osmelhorespontos.mat');%Load the results of the previous simulation
% Finding the optmal point
    for Inc=1:(phase_steps)  %For EACH phase
        phase_is = phase(Inc);%Salved here to be used on the main Loop
        for Inc2=1:Vbias_steps%and a Vbais
            Vbias_is = Vbias(Inc2);%that will also be saved to be used on 
                                   %the main loop
             for kk=1:Vbias_steps%For each voltage apply at the second arm
                for jj=1:Vbias_steps%and at the first arm
                    if find(ismember(melhorespontos,S))                    %Check if it is a optimum point
                        melhor_phi = [melhor_phi phase_is];                %If "Yes" the variables are stored
                        melhor_vbi = [melhor_vbi Vbias_is];
                        melhor_va1 = [melhor_va1 Vbias(jj)];
                        melhor_va2 = [melhor_va2 Vbias(kk)];
                        melhor_inp = [melhor_inp Inc2 + Vbias_steps*(Inc-1)];
                    end                                                    %If "No" it does nothing.
                    S = S + 1;%Counter
                end
             end
            display(S*100/(phase_steps*Vbias_steps*phase_steps*...
                                                         Vbias_steps),'%');%Used to show that
                                                                           %the simulation still runing.
        end        
    end
 % Qualitative evaluation of the optimal aspects founded
    contador = 1;%Counter
    for cc=1:length(melhorespontos)%run each best realization and plot it
        EleSig2 = melhor_va2(cc).*sin(Rad_f*t + melhor_phi(cc));           %The second arm will have a phase shift different of the first arm
        EleSig1 = melhor_va1(cc).*sin(Rad_f*t);%Set up of 1st arm
        EleSig.U1t=EleSig1;%The MZM model require that the voltages
        EleSig.U2t=EleSig2;%to be stores at an structure variable
        [MZ_Input_Sdata] = ['AA_MZ_Input_Data_t' num2str(melhor_inp(cc))]; %This is to call the right MZM-Input-Data for 
                                                                           %the simulation it's name vary with values of 
                                                                           %phase and Vbias
        % Mach Zhender Modulator
        [Eout(contador,:),~]=Mach_Zehnder_Modulator_simplificado(t,CW,...
                                                    EleSig,MZ_Input_Sdata);
        % Ploting results for evaluation
        plot(f,ML*log10(abs(fftshift(fft(Eout(contador,:)./length(Eout(...
                                                         contador,:)))))));
        hold all;
        contador = contador + 1;
    end
    save('melhoresconfig'); %When a simulation is finished all the data is 
                            %is saved.


