%c
%c                                                       ..'Ž`'..'Ž`..'Ž`..                                                   
%c       File: main_two_arm_loop_mzm_DP
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
%c   Make_MZ_Input_Files - Creates the files that will be the input data
%c                         to the Mach-Zehnder Modulator
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
MAIN_COMB_INPUT_DATA_DP;
%% Second Instance: Calling main Function to Execute the Simulation
display('Making Files...');
Make_MZ_Input_Files_DP;    %Create the files for the MZM input data
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
if Ploting
    figure;
    hold all;
    grid on;
    plot(f,ML*log10(abs(fftshift(Filtro))));
end
To_Eval = [];
distcomp.feature( 'LocalUseMpiexec', false );
for Va1=1:Varm1_steps                                                       %For each possible value of Amplitude of the RF signal at the Uper arm
    for Va2=1:Varm2_steps                                                  %For each possible value of Amplitude of the RF signal at the Lower arm
        for Ph12=1:phase_steps                                             %For each Phase differency between the RF sgnals
            for Gai=1:Gain_steps                                           %For each gain given to the output optical signal
%                 To_Eval_aux = struct([]);
                %% Loading variables  
                Gain        =Gain_vet(Gai);                    %The Gain That will be applied at the output of the MZM
                phase_is    = phase(Ph12);
                EleSig1     = Varm1(Va1)*sin(Rad_f*t);         %The first RF signal at the Uper arm.
                EleSig2     = Varm2(Va2)*sin(Rad_f*t+phase_is);%The second RF singal at the Lower arm.
                                                               %will have a phase shift different of 
                                                               %the first arm.
                EleSig.U1t  = EleSig1;                         %The MZM model require that the voltages
                EleSig.U2t  = EleSig2;                         %The MZM model require that the voltages
%                 tic
                parfor Inc0=1:phi_0_steps*V2bias_steps*V1bias_steps           %For each Phase differency between Uper and Lower arm
                                                                           %For each possible value of the polarization voltage at Uper arm
                                                                           %For each possible value of the polarization voltage at Lower arm
                                                                           %This main loop will run the simulation for every possible combination of those variables
                    MzmInputInd = (Inc0);                          %The Index for the file to be loaded by the MZM
                    %% Mach Zhender Modulator
                    
                    E1=zeros(1,length(CW));
                    Eout=zeros(1,length(E1));
                    %%
                    for N=1:N_volta


                       %Choose to use a infinit Loop or to have a output
                       if InfinitLoop
                           E1 = Eout + CW;
                       else      
                           E1 = CuplerGain*(Eout + CW);
                       end

                        % Mach Zhender Modulator
                      [Eout,~]=Mach_Zehnder_Modulator_DP(t,E1,EleSig,MzmInputInd,Local);
                    %   plot(f,10*log10(abs(fftshift(fft(Eout)./length(Eout)))))
                       %Use the amplifyer to give a bust on the signal
                       if Amplify
                           Eout = Gain*Eout;
                       end

                       %Use the filter on the Loop signal to limit the Nmb of Carriers
                       if Filtering
                           Eout = ifft(fft(Eout).*Filtro);    
                       end

                       % add the noise to the Loop signal before re-enter into the MZM
                       if AddNoise
                          EsN0 = EbN0 + 3 + 10*log10(k) - 10*log10(nsamp);%Calculate the 
                                                                         %signal to noise ratio
                          Eout = awgn(Eout,EsN0);
                       end

                        if Ploting
                    %         figure;
                            plot(f,ML*log10(abs(fftshift(fft(Eout)./length(Eout)))));
                    %         ponto_max =max(abs(fftshift(fft(Eout)./length(Eout))));%Before anything else it is necessary firt to 
                    %                                                                %know the maximum power value of the signal
                    %         limiar = ponto_max/Dif_Amp;                            %Finding the threshold
                    %         [peaks_aux,locs_aux] = findpeaks(abs(...               %Now it is possible to find how many peaks the 
                    %         fftshift(fft(Eout)./length(Eout))),...                 %signal has within our first specification whith
                    %                                'MinPeakHeight',limiar);
                    %         plot(f(locs_aux),ML*log10(abs(fftshift(fft(peaks_aux)./length(peaks_aux)))));
                        end
                    end
                    %%
%                     OFC_RIM_DP;
                    if Ploting
                        plot(f,ML*log10(abs(fftshift(fft(Eout)./length(Eout)))));
                    end
                    %% Finding the peaks
                    ponto_max =max(abs(fftshift(fft(Eout)./...
                                                   length(Eout))));%Before anything else it is necessary firt to 
                                                                   %know the maximum power value of the signal
                    limiar = ponto_max/Dif_Amp;                    %Finding the threshold
                    [peaks_aux,locs_aux] = findpeaks(abs(...       %Now it is possible to find how many peaks the 
                    fftshift(fft(Eout)./length(Eout))),...         %signal has within our first specification whith
                                           'MinPeakHeight',limiar);% the Dif_Amp varialbe
%                     peaks(Inc0,1:length(peaks_aux)) = peaks_aux;      %After finding the number of peaks they are stored for latter evaluation.
%                     locs(Inc0,1:length(locs_aux)) = locs_aux;         %Possitions of the frequency refent of the founded peaks

                    [~,U10,U20,~,~,~,~,~,~,~,~,~,~,phi_0,~,~] = ...%[L,U10,U20,U_pi1,U_pi2,eta1,eta2,V1pi0,V2pi0,Vphi0,nopt,nel,alfa_ins,phi_0,C,alfa0]
                             Import_Data_MZM_DP(MzmInputInd,Local);%Import the data from files
                    %% Saving main parameters of the simulation
                    To_Eval_aux(Inc0)= struct('Ome',Rad_f,'Pha',phase_is...
                    ,'VA1',Varm1(Va1),'VA2',Varm2(Va2),'V1b',U10,'V2b',...
                    U20,'VPh',phi_0,'Gai',Gain,'Fin',MzmInputInd,'Pea',...
                                                 peaks_aux,'Loc',locs_aux);
%                     To_Eval_aux{Inc0}.Ome = Rad_f;                        % Store the angular frequency of the signals generated
%                     To_Eval_aux{Inc0}.Pha = phase_is;                     % Stores the phase difference between each RF signal
%                     To_Eval_aux{Inc0}.VA1 = Varm1(Va1);                   % Stores the Amplitude of the RF signal on the first arm
%                     To_Eval_aux{Inc0}.VA2 = Varm2(Va2);                   % Stores the Amplitude of the RF signal on the second arm
%                     To_Eval_aux{Inc0}.V1b = U10;                          % Stores the bias voltage applied at the Uper arm
%                     To_Eval_aux{Inc0}.V2b = U20;                          % Stores the bias voltage applied at the Lower arm
%                     To_Eval_aux{Inc0}.VPh = phi_0;                        % Stored the voltage applied at the phase controler between Uper and Lower arms
%                     To_Eval_aux{Inc0}.Gai = Gain;                         % Stored the Gain applied in the optical ring
%                     To_Eval_aux{Inc0}.Fin = MzmInputInd;                  % Stored the Index number for the MZM configuration file
%                     To_Eval_aux{Inc0}.Pea = peaks_aux;                    % Stores the peaks founded on the simulation
%                     To_Eval_aux{Inc0}.Loc = locs_aux;                     % Stored the location of the frequency respective to the peaks founded
                end
%                 toc
                To_Eval = [To_Eval To_Eval_aux];
                %% End of the current interaction
                S = phi_0_steps*V2bias_steps*V1bias_steps*(Gai + ...
                Gain_steps*(Ph12-1) + Gain_steps*phase_steps*(Va2-1) + ...
                               Gain_steps*phase_steps*Varm2_steps*(Va1-1));                                     %counter variable
            end
            display((S-1)*100/(N_combination),'%');                        %the simulation still runing.
        end
    end
end
%% Saving Results and closing the simulation
save([LocSav SFName]);                                                     %When a simulation is finished all the data it is saved.
