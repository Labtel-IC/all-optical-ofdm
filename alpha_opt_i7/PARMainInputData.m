%c
%c                                                       ..'Ž`'..'Ž`..'Ž`..                                                   
%c       File: MainInputDAta File With the Main Script Parameters 
%c(Main file to call each idividual parameters)
%c
%c     This main code is resposible to  load all necessary parameters of 
%c the MAIN script, just in this file changes should be done. Do not change
%c the MAIN file. 
%c                                           by P.Marciano LG
%c                                           28/10/2017
%c                                           last updated
%c                                           02/04/2018
%c                                           pablorafael.mcx@gmail.com
%c 
%c     References:
%c@article{hillerkuss2010simple,
%c  title={Simple all-optical FFT scheme enabling Tbit/s real-time signal processing},
%c  author={Hillerkuss, D and Winter, M and Teschke, M and Marculescu, A and Li, J and Sigurdsson, G and Worms, K and Ezra, S Ben and Narkiss, N and Freude, W and others},
%c  journal={Optics express},
%c  volume={18},
%c  number={9},
%c  pages={9324--9340},
%c  year={2010},
%c  publisher={Optical Society of America}
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
%c
%c
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

%%                    OFC Generator Parameters
ThisPath = pwd;                                                            %Path for the folter where the OFC is saved
UofcsLocal = '\save_files\';                                               %Files name for the location of the OFC
if exist('CurrentOCS','var')
    OfcName = [ThisPath UofcsLocal 'OCS_' num2str(OcsToTest(CurrentOCS))];
else
    OfcName = [ThisPath UofcsLocal 'OCS1_08-Mar-2018'];
end

%%                 Control Variable
SendingDowStr = 1;        %Set 1 to send data downstream 0 to not send.
SigAten       = 0;        %Set 1 to atenuate the OFCS signal. Set to zero 
                          %to turn off attenuation
SnrRef        = 0;        %Set 1 to select the own signal for the noise 
                          %reference at receptor. Set 0 for the pilote 
                          %carrier as noise reference at receptor.
FFTSplit      = 1;        %Set 1 to use the OFFT to split the income 
                          %signal on receptor. Set 0 to use filter to 
                          %split the income signal at receptor.
CarSNR        = 100;      %Set the snr to be applied at receptor.
PastTime      = 0;        %Storage for the time flow on simulation
Atenuation    = 10^(4.0); %Set how much the OFCS will be attenuated
NumbOfTest    = 1000;     %The number of repetition that one test will be 
                          %done
TestBundle    = 1;        %Set if the scrip will be run once or many times
CurrentMedium = 1;        %Sellection of channel 1 for Fiber other of B2B
PAM4Type      = [1 0];    %Sellection the type of 4PAM that will be tested
                          % 1 1  ----> 4PAM for IQ-MZM - Optical 4PAM
                          % 1 0  ----> 4PAM for DD-MZM - Optical 4PAM
                          % 0 x  ----> 4PAM for DD-MZM - Electrical 4PAM
%
DecLevDef3 = 0.6953;
DecLevDef2 = 0.4013;
DecLevDef1 = 0.1877;
%This will be the variable that will select the type of modulation used by
%this script. So far, this the options goes from 1 to 4, where:
%       1   -->    OOK
%       2   -->    4PAM
%       3   -->    DQPSK
%       4   -->    DPSK
% ThisModula = 1;                                            
UsedModula = [2];
UpStModula = [1];
%The following verification will determinate which type of test will be on 
%running, a test bundle or a single test. In the section the user may 
%change with variables will be loaded for testing.
if TestBundle                                                              %Test Bundle selected
    ThisTest   = NumbOfTest;                                               %The number of tests per type of modulation
    ThisModula = UsedModula;                                               %The modulations under test
else
    ThisTest   = 1;                                                        %The number of tests per type of modulation
    ThisModula = UsedModula(1);                                            %The modulations under test
end









