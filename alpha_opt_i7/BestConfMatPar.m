%c
%c                                                       ..'�`'..'�`..'�`..                                                   
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
% load([pwd '\save_files\' 'AvaliandoRectFilt.mat']);                      %Load the result of the past simulation for evaluation
%%
% Local = [pwd '\input_files\'];%variable that shows where the input files for
                                                     %the mzm will be saved
display('Making Files...');
Make_MZ_Input_Files_Simp;    %Create the files for the MZM input data
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
    BesResMat(1,MatPos) = To_Eval(kk).Ome/(2*pi);
    BesResMat(2,MatPos) = To_Eval(kk).Pha;
    BesResMat(3,MatPos) = To_Eval(kk).VA1;
    BesResMat(4,MatPos) = To_Eval(kk).VA2;
    BesResMat(5,MatPos) = To_Eval(kk).Vbi;
    BesResMat(6,MatPos) = To_Eval(kk).Fid;
    BesResMat(7,MatPos) = length(To_Eval(kk).Pea);
    BesResMat(8,MatPos) = To_Eval(kk).Gai;
    BesResMat(9,MatPos) = BestPlanicity(MatPos);
    MatPos = MatPos + 1;
end
MatPos = MatPos - 1;
Vet_Lab = [{'Frequ�ncia Central [Hz]'} {'Fase entre Eletrodos'} {'Tens�o Bra�o 1 [V]'} {'Tens�o Bra�o 2 [V]'} {'Tens�o de Polariza��o [V]'} {'Index dos Arquivos'} {'Picos Encontrados'} {'Ganho'}];
for kk=1:size(BesResMat,1)-1
    [Vet_aux,Ind_aux] = sort(BesResMat(kk,:));
    figure;grid on;
    if kk==2
        plot(Vet_aux.*(180/pi),BesResMat(end,Ind_aux),'-o');
    else
        plot(Vet_aux,BesResMat(end,Ind_aux),'-o');
    end
    title('Varia��o da Planicidade');
    xlabel(Vet_Lab(kk));
    ylabel('Planicidade');
    a=1;
end
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
xlRange = [char((OffSet)+65) '1:' Index num2str(size(BesResMat,1))];                              %Saving the range in which this data will be saved
xlswrite(Name,BesResMat,sheet,xlRange);                                    %Writing the data in the Excel file
xlswrite(Name,[{'FC '};{'Pha'};{'VA1'};{'VA2'};{'Vbi'};{'Fid'};{'NPe'};...
                                           {'Gai'};{'Pla'}],sheet,'A1:A9');                                    %Writing the data in the Excel file

toplot = 0;
while ~toplot
    clc;verif = input('Do you to finish the simulation? y or n\n','s');
    if upper(verif) ~=('Y')
        clc;PlotIt = input('Insert the number of the best result to be plot.\n','s');
        if (isnan(str2double(PlotIt)))||(str2double(PlotIt)>size(BesResMat,1))
            clc;[~]=input('Not a valid input. Press enter to continue.\n');
        else
%             EleSig1 = BesResMat(4,round(str2double(PlotIt)))*sin(Rad_f*t);
%             EleSig = EleSig1;
%             MZ_Input_Sdata = BesResMat(6,round(str2double(PlotIt)));
            Gain = BesResMat(8,str2double(PlotIt));
%%

            L          = 10;
            U0         = BesResMat(5,round(str2double(PlotIt)));
            U_pi1      = 3.8;
            U_pi2      = 3.8;
            eta1       = 89;% 89;
            eta2       =-12.7;%-12.7;
            %
            nopt       =  2.17;
            nel        =  2.60;
            alfa_ins   =  5.1;
            phi_0      =  0.0;
            alfa0      =  0.55;

            C     = (eta1+eta2)/(eta1-eta2);       	%  Parametro de chirp
            alfa0 = 10^(alfa0/20); 			% [1/cm.GHz^0.5]
%             if isstruct(U),
              electrodes = 2;  % Two arms modulator
              U1t        = BesResMat(3,round(str2double(PlotIt)))*sin(Rad_f*t);
              U2t        = BesResMat(4,round(str2double(PlotIt)))*sin(Rad_f*t+BesResMat(2,round(str2double(PlotIt))));
              U_pi       = U_pi2;
%             else
%               electrodes = 1;  % One arm modulator
%               U1t        = U;  
%               U_pi       = U_pi1;
%             end
            %
            tps        = t/1E-12;
            ccmns      = 30;			     	% Velocidade da luz [cm/ns]
            freqTHz    = time2freq_lamb(tps);
            freqGHz    = freqTHz*1e-3;			% Frequencia em GHz
            freqGHz    = -fftshift(freqGHz);
            freqGHz(1) = freqGHz(2);
            n          = length(U1t); 			% Tamanho de U1t
            %

            %%%%%%%  FUNCAO DE TRANSFEWRENCIA ELETRICA DO MODULADOR  %%%%%%%%%%%
            %
            alfaf  = 0.5*alfa0*(abs(freqGHz).^0.5);
            gamaf  = 2*pi*abs(freqGHz).*(nel - nopt)/ccmns;
            atn    = alfaf + 1j*gamaf;
            H      = (1./(atn*L)).*(1 - exp(-atn*L));
            %
            if (electrodes == 1),
              U1f   = fft(U1t);
              U1t   = real(ifft(U1f.*H));
              exp1  = exp(1j*C*(pi/2).*(U1t/U_pi));
              if C~=0
                  Eout  = CW.*cos((pi/2).*(U1t - U0)/U_pi).*exp1;
              else
                  Eout  = CW.*cos((pi/2).*(U1t - U0)/U_pi);
              end
            else
              U1f   = fft(U1t);
              U1t   = ifft(U1f.*H);
              U2f   = fft(U2t);
              U2t   = ifft(U2f.*H);
              if C ~= 0
                  Eout = CW.*cos((pi/2).*(U1t - U2t - U0)/U_pi).*exp(-1j*(pi/2).*((U1t + U2t)/U_pi));
              else
                  Eout = CW.*cos((pi/2).*(U1t - U2t - U0)/U_pi);
              end
            end
%%            
            plot(f,ML*log10(abs(fftshift(fft(Eout)./length(Eout)))))
            hold on
            plot(f(To_Eval(melhorespontos(str2double(PlotIt))).Loc),ML*...
              log10(To_Eval(melhorespontos(str2double(PlotIt))).Pea),'-o');
        end
    else
        toplot = 1;
    end
end
save(LocalSaved,'BesResMat');                         %Saving the results
%%
a=1;
