%c
%c                                                       ..'Ž`'..'Ž`..'Ž`..                                                   
%c       File: SDMZM_Expe_Resul
%c(Loading the results from acquired in laboratory for further analizes)
%c
%c     This main code is resposible to load a matrix with the results that 
%c     was acquired.
%c
%c                                           by P.Marciano LG
%c                                           08/10/2017
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
Vbias_Index   = 1;  %Polarization voltage
PotRis_mW    = 2;  %Linear power measured when the first arm Bias Voltage 
                    %was changed in an rissing trend
PotRis_dBm   = 3;  %Logarithm power measured when the first arm Bias 
                    %Voltage was changed in an rissing trend
PotFal_mW    = 4;  %Linear power measured when the first arm Bias Voltage 
                    %was changed in an falling trend
PotFal_dBm   = 5;  %Logarithm power measured when the first arm Bias 
                    %Voltage was changed in an falling trend
%% Results from the Single Driven Mach-Zehnder Modulator
mat_sd_mzm = [...
-7.0	1.109	0.450	1.091	0.380;
-6.8	1.107	0.440	1.072	0.300;
-6.6	1.096	0.400	1.047	0.200;
-6.4	1.079	0.330	0.989	-0.050;
-6.2	1.057	0.240	0.966	-0.150;
-6.0	1.026	0.110	0.912	-0.400;
-5.8	0.989	-0.050	0.851	-0.700;
-5.6	0.946	-0.240	0.794	-1.000;
-5.4	0.895	-0.480	0.724	-1.400;
-5.2	0.841	-0.750	0.646	-1.900;
-5.0	0.782	-1.070	0.575	-2.400;
-4.8	0.716	-1.450	0.501	-3.000;
-4.6	0.649	-1.880	0.417	-3.800;
-4.4	0.577	-2.390	0.355	-4.500;
-4.2	0.507	-2.950	0.282	-5.500;
-4.0	0.438	-3.590	0.214	-6.700;
-3.8	0.370	-4.320	0.166	-7.800;
-3.6	0.304	-5.170	0.116	-9.340;
-3.4	0.241	-6.180	0.074	-11.300;
-3.2	0.182	-7.400	0.045	-13.500;
-3.0	0.132	-8.800	0.018	-17.500;
-2.8	0.089	-10.500	0.004	-23.500;
-2.6	0.054	-12.700	0.001	-29.900;
-2.4	0.028	-15.600	0.007	-21.300;
-2.2	0.009	-20.500	0.024	-16.200;
-2.0	0.001	-30.500	0.050	-13.050;
-1.8	0.004	-24.500	0.083	-10.800;
-1.6	0.016	-18.000	0.124	-9.050;
-1.4	0.037	-14.300	0.179	-7.470;
-1.2	0.066	-11.800	0.237	-6.250;
-1.0	0.105	-9.800	0.294	-5.310;
-0.8	0.148	-8.300	0.356	-4.490;
-0.6	0.200	-7.000	0.420	-3.770;
-0.4	0.263	-5.800	0.495	-3.050;
-0.2	0.331	-4.800	0.562	-2.500;
+0.0	0.398	-4.000	0.710	-1.490;
+0.2	0.495	-3.050	0.719	-1.430;
+0.4	0.562	-2.500	0.847	-0.720;
+0.6	0.631	-2.000	0.906	-0.430;
+0.8	0.708	-1.500	0.957	-0.190;
+1.0	0.776	-1.100	1.002	0.010;
+1.2	0.832	-0.800	1.038	0.160;
+1.4	0.891	-0.500	1.067	0.280;
+1.6	0.955	-0.200	1.086	0.360;
+1.8	1.000	0.000	1.096	0.400;
+2.0	1.047	0.200	1.099	0.410;
+2.2	1.072	0.300	1.091	0.380;
+2.4	1.096	0.400	1.072	0.300;
+2.6	1.104	0.430	1.050	0.210;
+2.8	1.109	0.450	1.012	0.050;
+3.0	1.102	0.420	0.971	-0.130;
+3.2	1.086	0.360	0.923	-0.350;
+3.4	1.059	0.250	0.869	-0.610;
+3.6	1.023	0.100	0.813	-0.900;
+3.8	0.977	-0.100	0.741	-1.300;
+4.0	0.923	-0.350	0.676	-1.700;
+4.2	0.863	-0.640	0.610	-2.150;
+4.4	0.785	-1.050	0.543	-2.650;
+4.6	0.716	-1.450	0.473	-3.250;
+4.8	0.643	-1.920	0.417	-3.800;
+5.0	0.569	-2.450	0.347	-4.600;
+5.2	0.494	-3.060	0.288	-5.400;
+5.4	0.417	-3.800	0.229	-6.400;
+5.6	0.347	-4.600	0.186	-7.300;
+5.8	0.282	-5.500	0.141	-8.500;
+6.0	0.214	-6.700	0.105	-9.800;
+6.2	0.162	-7.900	0.071	-11.500;
+6.4	0.098	-10.100	0.045	-13.500;
+6.6	0.058	-12.400	0.025	-16.000;
+6.8	0.030	-15.300	0.012	-19.200;
+7.0	0.010	-19.800	0.010	-20.000];