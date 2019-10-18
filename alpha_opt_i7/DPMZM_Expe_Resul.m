%c
%c                                                       ..'Ž`'..'Ž`..'Ž`..                                                   
%c       File: DPMZM_Expe_Resul
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
%c PotRis1_mW    :  %Linear power measured when the first arm Bias Voltage 
%c                     %was changed in an rissing trend
%c PotRis1_dBm   :  %Logarithm power measured when the first arm Bias 
%c                     %Voltage was changed in an rissing trend
%c PotFal1_mW    :  %Linear power measured when the first arm Bias Voltage 
%c                     %was changed in an falling trend
%c PotFal1_dBm   :  %Logarithm power measured when the first arm Bias 
%c                     %Voltage was changed in an falling trend
%c PotRis2_mW    :  %Linear power measured when the second arm Bias Voltage 
%c                     %was changed in an rissing trend
%c PotRis2_dBm   :  %Logarithm power measured when the second arm Bias 
%c                     %Voltage was changed in an rissing trend
%c PotFal2_mW    :  %Linear power measured when the second arm Bias Voltage 
%c                     %was changed in an falling trend
%c PotFal2_dBm   :  %Logarithm power measured when the second arm Bias 
%c                     %Voltage was changed in an falling trend
%c PotRisFas_mW  :  %Linear power measured when the fase voltage was 
%c                     %changed in a rissing trend
%c PotRisFas_dBm :  %Logarithm power measured when the fase voltage was 
%c                     %changed in a rissing trend
%c PotFalFas_mW  :  %Linear power measured when the fase voltage was 
%c                     %changed in a falling trend
%c PotFalFas_dBm :  ;%Logarithm power measured when the fase voltage was 
%c                     %changed in a falling trend
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
clear, close all, clc;
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
PotRis1_mW    = 2;  %Linear power measured when the first arm Bias Voltage 
                    %was changed in an rissing trend
PotRis1_dBm   = 3;  %Logarithm power measured when the first arm Bias 
                    %Voltage was changed in an rissing trend
PotFal1_mW    = 4;  %Linear power measured when the first arm Bias Voltage 
                    %was changed in an falling trend
PotFal1_dBm   = 5;  %Logarithm power measured when the first arm Bias 
                    %Voltage was changed in an falling trend
PotRis2_mW    = 6;  %Linear power measured when the second arm Bias Voltage 
                    %was changed in an rissing trend
PotRis2_dBm   = 7;  %Logarithm power measured when the second arm Bias 
                    %Voltage was changed in an rissing trend
PotFal2_mW    = 8;  %Linear power measured when the second arm Bias Voltage 
                    %was changed in an falling trend
PotFal2_dBm   = 9;  %Logarithm power measured when the second arm Bias 
                    %Voltage was changed in an falling trend
PotRisFas_mW  = 10; %Linear power measured when the fase voltage was 
                    %changed in a rissing trend
PotRisFas_dBm = 11; %Logarithm power measured when the fase voltage was 
                    %changed in a rissing trend
PotFalFas_mW  = 12; %Linear power measured when the fase voltage was 
                    %changed in a falling trend
PotFalFas_dBm = 13 ;%Logarithm power measured when the fase voltage was 
                    %changed in a falling trend
%% Results from the Dual Parallel Mach-Zehnder Modulator
mat_dp_mzm = [...
-7.0	0.011	-19.70	0.011	-19.62	0.004	-23.50	0.005	-23.08	...
0.585	-2.33	0.585	-2.33;
-6.8	0.011	-19.50	0.011	-19.54	0.005	-23.33	0.006	-22.34	...
0.540	-2.68	0.543	-2.65;
-6.6	0.011	-19.55	0.011	-19.62	0.005	-22.71	0.008	-21.05	...
0.495	-3.05	0.506	-2.96;
-6.4	0.010	-19.91	0.010	-20.18	0.007	-21.51	0.011	-19.40	...
0.457	-3.40	0.465	-3.33;
-6.2	0.009	-20.65	0.008	-21.18	0.010	-19.95	0.017	-17.66	...
0.403	-3.95	0.424	-3.73;
-6.0	0.007	-21.8	0.005	-22.89	0.015	-18.23	0.030	-15.3	...
0.343	-4.65	0.370	-4.32;
-5.8	0.004	-23.7	0.003	-25.56	0.022	-16.52	0.035	-14.5	...
0.294	-5.32	0.316	-5.00;
-5.6	0.002	-26.75	0.001	-30.33	0.032	-14.92	0.052	-12.83	...
0.239	-6.22	0.267	-5.73;
-5.4	0.001	-32.3	0.000	-39.21	0.045	-13.42	0.071	-11.5	...
0.189	-7.23	0.217	-6.63;
-5.2	0.000	-38.6	0.001	-29.77	0.062	-12.05	0.094	-10.28	...
0.141	-8.5	0.170	-7.70;
-5.0	0.001	-28.95	0.004	-23.78	0.083	-10.8	0.121	-9.17	...
0.101	-9.96	0.129	-8.89;
-4.8	0.005	-23.45	0.010	-19.95	0.108	-9.65	0.152	-8.18	...
0.066	-11.81	0.088	-10.55;
-4.6	0.010	-19.82	0.019	-17.12	0.136	-8.65	0.188	-7.25	...
0.040	-14.03	0.056	-12.48;
-4.4	0.019	-17.15	0.032	-14.98	0.170	-7.7	0.226	-6.45	...
0.019	-17.12	0.032	-14.98;
-4.2	0.032	-14.99	0.049	-13.11	0.207	-6.85	0.270	-5.68	...
0.008	-21.2	0.015	-18.15;
-4.0	0.048	-13.2	0.070	-11.52	0.248	-6.05	0.316	-5.00	...
0.005	-23.3	0.007	-21.68;
-3.8	0.068	-11.65	0.095	-10.22	0.288	-5.41	0.365	-4.38	...
0.010	-19.8	0.007	-21.55;
-3.6	0.093	-10.33	0.124	-9.05	0.339	-4.7	0.412	-3.85	...
0.025	-16.08	0.016	-18.05;
-3.4	0.120	-9.2	0.158	-8.02	0.385	-4.15	0.475	-3.23	...
0.047	-13.28	0.032	-14.93;
-3.2	0.152	-8.18	0.194	-7.13	0.432	-3.65	0.519	-2.85	...
0.076	-11.19	0.057	-12.45;
-3.0	0.188	-7.25	0.234	-6.31	0.481	-3.18	0.562	-2.5	...
0.112	-9.5	0.089	-10.53;
-2.8	0.226	-6.45	0.274	-5.62	0.535	-2.72	0.603	-2.2	...
0.153	-8.15	0.127	-8.95;
-2.6	0.267	-5.74	0.317	-4.99	0.570	-2.44	0.649	-1.88	...
0.200	-7.0	0.172	-7.65;
-2.4	0.308	-5.11	0.359	-4.45	0.624	-2.05	0.689	-1.62	...
0.246	-6.09	0.219	-6.60;
-2.2	0.352	-4.54	0.403	-3.95	0.665	-1.77	0.706	-1.51	...
0.296	-5.28	0.265	-5.76;
-2.0	0.393	-4.06	0.447	-3.5	0.689	-1.62	0.733	-1.35	...
0.347	-4.6	0.316	-5.00;
-1.8	0.435	-3.62	0.484	-3.15	0.733	-1.35	0.760	-1.19	...
0.398	-4.0	0.367	-4.35;
-1.6	0.474	-3.24	0.519	-2.85	0.755	-1.22	0.759	-1.2	...
0.441	-3.56	0.419	-3.78;
-1.4	0.511	-2.92	0.551	-2.59	0.774	-1.11	0.773	-1.12	...
0.494	-3.06	0.467	-3.31;
-1.2	0.541	-2.67	0.574	-2.41	0.767	-1.15	0.767	-1.15	...
0.520	-2.84	0.504	-2.98;
-1.0	0.566	-2.47	0.596	-2.25	0.774	-1.11	0.767	-1.15	...
0.569	-2.45	0.543	-2.65;
-0.8	0.593	-2.27	0.605	-2.18	0.767	-1.15	0.736	-1.33	...
0.589	-2.3	0.578	-2.38;
-0.6	0.603	-2.2	0.614	-2.12	0.750	-1.25	0.711	-1.48	...
0.627	-2.03	0.605	-2.18;
-0.4	0.617	-2.1	0.618	-2.09	0.733	-1.35	0.689	-1.62	...
0.638	-1.95	0.617	-2.10;
-0.2	0.612	-2.13	0.610	-2.15	0.700	-1.55	0.649	-1.88	...
0.644	-1.91	0.624	-2.05;
+0.0	0.608	-2.16	0.596	-2.25	0.668	-1.75	0.598	-2.23	...
0.631	-2.0	0.631	-2.00;
+0.2	0.603	-2.2	0.577	-2.39	0.624	-2.05	0.562	-2.5	...
0.624	-2.05	0.624	-2.05;
+0.4	0.578	-2.38	0.550	-2.6	0.596	-2.25	0.519	-2.85	...
0.610	-2.15	0.610	-2.15;
+0.6	0.556	-2.55	0.522	-2.82	0.537	-2.7	0.465	-3.33	...
0.596	-2.25	0.577	-2.39;
+0.8	0.531	-2.75	0.490	-3.1	0.495	-3.05	0.417	-3.8	...
0.556	-2.55	0.550	-2.60;
+1.0	0.495	-3.05	0.451	-3.46	0.444	-3.53	0.369	-4.33	...
0.528	-2.77	0.513	-2.90;
+1.2	0.459	-3.38	0.409	-3.88	0.396	-4.02	0.320	-4.95	...
0.484	-3.15	0.476	-3.22;
+1.4	0.417	-3.8	0.368	-4.34	0.347	-4.6	0.274	-5.63	...
0.437	-3.6	0.426	-3.71;
+1.6	0.373	-4.28	0.327	-4.85	0.302	-5.2	0.231	-6.36	...
0.394	-4.05	0.376	-4.25;
+1.8	0.331	-4.8	0.282	-5.5	0.254	-5.95	0.194	-7.12	...
0.343	-4.65	0.323	-4.91;
+2.0	0.290	-5.38	0.243	-6.15	0.214	-6.69	0.158	-8.00	...
0.288	-5.4	0.272	-5.65;
+2.2	0.248	-6.06	0.204	-6.91	0.176	-7.54	0.126	-9.00	...
0.238	-6.23	0.223	-6.51;
+2.4	0.207	-6.85	0.166	-7.8	0.143	-8.45	0.098	-10.07	...
0.191	-7.18	0.177	-7.53;
+2.6	0.170	-7.7	0.132	-8.79	0.112	-9.51	0.074	-11.28	...
0.147	-8.33	0.133	-8.75;
+2.8	0.136	-8.68	0.102	-9.9	0.086	-10.65	0.054	-12.65	...
0.106	-9.75	0.095	-10.22;
+3.0	0.105	-9.8	0.077	-11.13	0.060	-12.19	0.039	-14.05	...
0.071	-11.5	0.062	-12.05;
+3.2	0.078	-11.09	0.056	-12.55	0.043	-13.64	0.027	-15.68	...
0.043	-13.65	0.037	-14.35;
+3.4	0.056	-12.53	0.037	-14.26	0.030	-15.27	0.018	-17.5	...
0.022	-16.52	0.018	-17.45;
+3.6	0.038	-14.25	0.024	-16.22	0.020	-17.08	0.011	-19.55	...
0.010	-20.18	0.008	-21.15;
+3.8	0.024	-16.25	0.014	-18.62	0.012	-19.15	0.007	-21.81	...
0.005	-22.93	0.006	-22.45;
+4.0	0.013	-18.79	0.007	-21.62	0.007	-21.63	0.004	-24.2	...
0.009	-20.45	0.012	-19.25;
+4.2	0.006	-21.99	0.003	-25.5	0.004	-24.15	0.002	-26.52	...
0.021	-16.75	0.026	-15.80;
+4.4	0.002	-26.08	0.001	-29.57	0.002	-26.51	0.001	-28.25	...
0.041	-13.83	0.048	-13.15;
+4.6	0.001	-29.93	0.001	-29.28	0.001	-28.32	0.001	-29.03	...
0.069	-11.61	0.078	-11.09;
+4.8	0.001	-28.56	0.003	-25.95	0.001	-29.1	0.001	-29.15	...
0.104	-9.85	0.114	-9.44;
+5.0	0.003	-25.15	0.005	-23.33	0.001	-29.19	0.001	-29.13	...
0.144	-8.43	0.155	-8.09;
+5.2	0.005	-22.6	0.007	-21.55	0.001	-29.19	0.001	-29.35	...
0.190	-7.22	0.202	-6.95;
+5.4	0.008	-21.02	0.009	-20.33	0.001	-29.45	0.001	-29.77	...
0.239	-6.22	0.248	-6.05;
+5.6	0.010	-19.95	0.011	-19.58	0.001	-29.91	0.001	-30.05	...
0.292	-5.35	0.301	-5.21;
+5.8	0.012	-19.34	0.012	-19.15	0.001	-30.00	0.001	-29.33	...
0.343	-4.65	0.355	-4.50;
+6.0	0.012	-19.06	0.012	-19.08	0.001	-29.77	0.002	-27.33	...
0.391	-4.08	0.403	-3.95;
+6.2	0.012	-19.15	0.012	-19.27	0.002	-26.35	0.003	-24.82	...
0.442	-3.55	0.437	-3.60;
+6.4	0.011	-19.6	0.010	-19.85	0.004	-23.55	0.006	-22.2	...
0.485	-3.14	0.495	-3.05;
+6.6	0.009	-20.41	0.008	-20.77	0.008	-21.00	0.010	-19.84	...
0.533	-2.73	0.531	-2.75;
+6.8	0.007	-21.66	0.006	-22.11	0.014	-18.62	0.017	-17.76	...
0.569	-2.45	0.562	-2.50;
+7.0	0.004	-23.82	0.004	-24.09	0.022	-16.55	0.025	-16.00	...
0.600	-2.22	0.592	-2.28;];