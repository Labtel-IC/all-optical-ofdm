%c
%c                                                       ..'Ž`'..'Ž`..'Ž`..                                                   
%c       File: SDMZM_Expe_Resul
%c(Loading the results from acquired in laboratory for further analizes)
%c
%c     This main code is resposible to load a matrix with the results that 
%c     was acquired.
%c
%c                                           by P.Marciano LG
%c                                           22/10/2017
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
%c Vbias_Index   :  %Polarization voltage
%c PotRis_mW    :  %Linear power measured when the first arm Bias Voltage 
%c                     %was changed in an rissing trend
%c PotRis_dBm   :  %Logarithm power measured when the first arm Bias 
%c                     %Voltage was changed in an rissing trend
%c PotFal_mW    :  %Linear power measured when the first arm Bias Voltage 
%c                     %was changed in an falling trend
%c PotFal_dBm   :  %Logarithm power measured when the first arm Bias 
%c                     %Voltage was changed in an falling trend
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
% clear, close all, clc;
%% Load the results from experiments
% Following matrix has its columns representing those values as noted:
%Vbias,	Potência Rising [mW],	Potência [dBm],	Potência Falling [mW],	
%Potência [dBm],	Potência Rising [mW],	Potência [dBm],	
%Potência Falling [mW],	Potência [dBm],	Potência Rising [mW],	
%Potência [dBm],	Potência Falling [mW] and	Potência [dBm] on this
%order
%% Those variables give the lable for each colun of the matrix with its 
%  respective meaning.
% It is important to mention that for all curves acquired the other inputs
% where granded for stability of the output signal.
M10Vbias_Index   = 1;  %Polarization voltage
M10PotRis1_mW    = 4;  %Linear power measured when the first arm Bias Voltage 
                    %was changed in an rissing trend
M10PotRis1_dBm   = 5;  %Logarithm power measured when the first arm Bias 
                    %Voltage was changed in an rissing trend
M10PotFal1_mW    = 2;  %Linear power measured when the first arm Bias Voltage 
                    %was changed in an falling trend
M10PotFal1_dBm   = 3;  %Logarithm power measured when the first arm Bias 
                    %Voltage was changed in an falling trend
M10PotRis2_mW    = 8;  %Linear power measured when the first arm Bias Voltage 
                    %was changed in an rissing trend
M10PotRis2_dBm   = 9;  %Logarithm power measured when the first arm Bias 
                    %Voltage was changed in an rissing trend
M10PotFal2_mW    = 6;  %Linear power measured when the first arm Bias Voltage 
                    %was changed in an falling trend
M10PotFal2_dBm   = 7;  %Logarithm power measured when the first arm Bias 
                    %Voltage was changed in an falling trend
%% Results from the Single Driven Mach-Zehnder Modulator
mat_mach10 = [...
7.00,0.59,-2.30,0.50,-3.01,0.53,-2.72,0.55,-2.56;
6.80,0.54,-2.70,0.45,-3.45,0.50,-3.02,0.51,-2.89;
6.60,0.49,-3.10,0.40,-4.01,0.44,-3.56,0.46,-3.33;
6.40,0.44,-3.54,0.35,-4.6,0.39,-4.12,0.43,-3.70;
6.20,0.39,-4.13,0.31,-5.14,0.32,-4.99,0.37,-4.33;
6.00,0.34,-4.74,0.27,-5.76,0.27,-5.72,0.33,-4.85;
5.80,0.29,-5.42,0.23,-6.41,0.22,-6.53,0.28,-5.53;
5.60,0.24,-6.19,0.19,-7.2,0.18,-7.41,0.23,-6.29;
5.40,0.19,-7.12,0.16,-7.9,0.14,-8.49,0.20,-7.04;
5.20,0.15,-8.10,0.13,-8.77,0.11,-9.66,0.16,-8.01;
5.00,0.12,-9.28,0.10,-9.8,0.08,-11.00,0.13,-8.99;
4.80,0.09,-10.67,0.08,-11.01,0.06,-12.56,0.10,-10.01;
4.60,0.06,-12.31,0.05,-12.6,0.04,-14.45,0.08,-11.20;
4.40,0.04,-14.39,0.04,-14.3,0.02,-16.79,0.06,-12.40;
4.20,0.02,-16.97,0.02,-16.2,0.01,-19.84,0.04,-13.80;
4.00,0.01,-20.39,0.01,-18.32,0.00,-23.25,0.03,-15.80;
3.80,0.00,-24.31,0.01,-21.06,0.00,-24.54,0.02,-18.06;
3.60,0.00,-23.70,0.01,-22.8,0.01,-21.65,0.01,-21.07;
3.40,0.01,-19.92,0.00,-24.4,0.01,-18.44,0.00,-24.01;
3.20,0.02,-16.67,0.00,-24.0,0.03,-15.90,0.00,-24.72;
3.00,0.04,-14.25,0.01,-21.02,0.04,-13.87,0.00,-23.48;
2.80,0.06,-12.38,0.01,-18.6,0.06,-12.26,0.01,-21.99;
2.60,0.08,-10.88,0.02,-17.02,0.08,-10.95,0.01,-18.25;
2.40,0.11,-9.66,0.02,-16.7,0.10,-9.84,0.03,-16.01;
2.20,0.14,-8.66,0.03,-14.84,0.13,-8.91,0.04,-14.38;
2.00,0.16,-7.84,0.05,-13.2,0.15,-8.14,0.05,-12.97;
1.80,0.19,-7.18,0.07,-11.8,0.18,-7.45,0.07,-11.75;
1.60,0.22,-6.54,0.08,-10.88,0.21,-6.86,0.09,-10.45;
1.40,0.25,-5.99,0.10,-9.9,0.23,-6.37,0.12,-9.39;
1.20,0.28,-5.49,0.13,-8.95,0.26,-5.82,0.14,-8.52;
1.00,0.31,-5.08,0.16,-8.06,0.29,-5.45,0.17,-7.68;
0.80,0.35,-4.60,0.18,-7.38,0.32,-4.99,0.20,-6.89;
0.60,0.38,-4.24,0.22,-6.56,0.35,-4.59,0.23,-6.34;
0.40,0.41,-3.89,0.25,-5.99,0.37,-4.28,0.27,-5.68;
0.20,0.44,-3.54,0.29,-5.43,0.40,-3.95,0.25,-6.07;
0.00,0.47,-3.25,0.30,-5.22,0.43,-3.68,0.28,-5.55;
-0.20,0.49,-3.12,0.34,-4.72,0.45,-3.48,0.31,-5.09;
-0.40,0.53,-2.77,0.37,-4.36,0.48,-3.22,0.35,-4.62;
-0.60,0.54,-2.7,0.41,-3.85,0.51,-2.90,0.38,-4.22;
-0.80,0.58,-2.37,0.44,-3.57,0.53,-2.77,0.42,-3.80;
-1.00,0.61,-2.16,0.48,-3.2,0.54,-2.69,0.43,-3.69;
-1.20,0.64,-1.94,0.50,-3.02,0.58,-2.38,0.46,-3.36;
-1.40,0.69,-1.6,0.52,-2.8,0.61,-2.12,0.48,-3.16;
-1.60,0.71,-1.47,0.55,-2.56,0.64,-1.92,0.52,-2.85;
-1.80,0.77,-1.16,0.59,-2.26,0.67,-1.73,0.54,-2.69;
-2.00,0.78,-1.08,0.63,-2.04,0.70,-1.56,0.59,-2.32;
-2.20,0.81,-0.92,0.65,-1.9,0.73,-1.39,0.61,-2.12;
-2.40,0.83,-0.8,0.71,-1.5,0.75,-1.23,0.64,-1.95;
-2.60,0.85,-0.7,0.72,-1.42,0.77,-1.14,0.68,-1.70;
-2.80,0.86,-0.65,0.76,-1.2,0.78,-1.09,0.71,-1.51;
-3.00,0.89,-0.5,0.82,-0.88,0.80,-0.95,0.74,-1.32;
-3.20,0.88,-0.55,0.84,-0.74,0.81,-0.90,0.76,-1.20;
-3.40,0.88,-0.54,0.87,-0.59,0.82,-0.88,0.78,-1.07;
-3.60,0.90,-0.48,0.86,-0.65,0.84,-0.75,0.80,-0.98;
-3.80,0.92,-0.36,0.87,-0.6,0.81,-0.89,0.80,-0.98;
-4.00,0.91,-0.39,0.90,-0.45,0.82,-0.86,0.83,-0.79;
-4.20,0.91,-0.42,0.88,-0.55,0.81,-0.94,0.80,-0.95;
-4.40,0.87,-0.6,0.91,-0.41,0.81,-0.89,0.84,-0.78;
-4.60,0.88,-0.56,0.90,-0.48,0.82,-0.85,0.82,-0.86;
-4.80,0.84,-0.78,0.85,-0.7,0.78,-1.07,0.80,-0.95;
-5.00,0.84,-0.77,0.84,-0.78,0.77,-1.14,0.78,-1.07;
-5.20,0.79,-1.02,0.81,-0.89,0.77,-1.13,0.76,-1.21;
-5.40,0.74,-1.29,0.79,-1.05,0.73,-1.37,0.71,-1.49;
-5.60,0.70,-1.55,0.73,-1.34,0.70,-1.58,0.68,-1.66;
-5.80,0.68,-1.67,0.70,-1.57,0.67,-1.76,0.63,-2.00;
-6.00,0.63,-1.98,0.65,-1.9,0.63,-2.04,0.61,-2.18;
-6.20,0.58,-2.39,0.58,-2.37,0.60,-2.25,0.56,-2.50;
-6.40,0.52,-2.8,0.53,-2.77,0.54,-2.68,0.50,-3.00;
-6.60,0.48,-3.2,0.47,-3.25,0.49,-3.10,0.45,-3.48;
-6.80,0.41,-3.83,0.41,-3.86,0.44,-3.59,0.40,-3.95;
-7.00,0.36,-4.4,0.36,-4.47,0.38,-4.20,0.36,-4.48];