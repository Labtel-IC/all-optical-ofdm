function [ Eout1,Eout2 ] = DelayInterfG( f,TimDel,PhaDel,E1,E0,CoupConst,...
                                                          ExpCoef,MutCoef )
%%                        Delay Interferometer
%c function[Eout1,Eout2] = DelayInterf( t,TimDel,PhaDel,E1,E0,CoupConst,...
%c                                                        ExpCoef,MutCoef )
%c This function is responsible for implementing an Delay Interferometer 
%c (DI). There many ways that this device can be implemented such as MZM,
%c passive components and couples with diffe light paths. This script will
%c inplement a DI ussing 3dB couples with 2-input and 2-outputs. The delay
%c given will be implement by light paths with different lengths and the
%c phase delay will be achieve with a phase controler. This will be the
%c principal and only device used to implement an All Optica FFT. Aiming to
%c separate subcarriers from an OFDM signal and convert from serial to
%c paralel.
%c
%c
%c                                           Created by P.Marciano LG
%c                                           02/11/2017
%c                                           pablorafael.mcx@gmail.com
%c
%c Refences:
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
%c
%c
%c
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                    G L O B A L  V A R I A B L E S                      %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                   L I S T  O F  V A R I A B L E S                      %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c
%c    DelayInterf
%c
%c  Input:
%c  t         : Time vector of the whole simulation [s]
%c  TimDel    : Time delay to be implemented by the DI
%c  PhaDel    : Phase delay to be implemented by the DI
%c  E1        : Principal input signal to be implemented
%c  E0        : Secundary input signal (usualy it is zero)
%c  CoupConst : Couples Constant of the DI
%c  ExpCoef   : Exponential component from Ramaswami equation
%c  MutCoef   : Coeficient to control polarity
%c  
%c  Output:
%c  Eout1     : Output1 from the DI interfometric response
%c  Eout2     : Output2 from the DI interfometric response
%c  
%c
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                             S T R U C T S                              %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c
%%
% At the very first part of this iscript, it is necessary to be sure that
% all variable that will be used are needed to be correctly initialized.
%% Delay Interferometers
% it was fundamental to use a teoretical model for the couple. Because just
% adding signals together or subtracting them those will not result as an 
% optical coupler. The model used at this script was taken from Optical
% Networks by Ramaswami. 

%Verifing the variables used on the optical coupler were correctly
%initialized.
if nargin<5
%     DelayInterfInputData;
    CoupConst = gpuArray(3*pi/4);    %Couples Constant of the DI
    ExpCoef   = gpuArray(exp(-1j*1));%Exponential component from Ramaswami equation
    MutCoef   = gpuArray(1);         %Coeficient to control polarity
    E0 = gpuArray(0);
elseif nargin<6
    CoupConst = gpuArray(3*pi/4);    %Couples Constant of the DI
    ExpCoef   = gpuArray(exp(-1j*1));%Exponential component from Ramaswami equation
    MutCoef   = gpuArray(1);         %Coeficient to control polarity
else
    error(['Input Arguments are not enougth. Please check this function'...
                                                                ' Call!']);
end
Pimag = gpuArray(1j);
Nimag = gpuArray(-1j);
G2    = gpuArray(2);
Gpi   = gpuArray(pi);
%% * Stage: 1 spliting the principal input signal
% Fristly it is needed to divide the signal in two given that by the length
% difference of the optical pahts and the phase delay btween them the
% signals will interfere with one another.
 % split
 % An 3dB coupler with 2 inputs and 2 outputs is used to slplit the input
 % signal.
 %The couples will output two signals that are a combination of E1 and E0.
 %As E0 = 0 the out put will be just the E1 signal slpit in two. The name 
 %notation E11 means the E1 signal at the arm 1 (uper arm) hence the 
 %notation name E12 means the E1 signal at the arm 2 (Lower arm).
E11 = ExpCoef.*(cos(CoupConst).*E1 + MutCoef*Pimag*sin(CoupConst).*E0);        
E12 = ExpCoef.*(cos(CoupConst).*E0 + MutCoef*Pimag*sin(CoupConst).*E1);       

%% * Stage: 2 Changing the signal at the lower arm
 %Dellay and change phase
  %Firstly it is taken in acont what means for the input signal a time 
  %delay by measuring the input time vector
% aux1 = E12(t<=TimDel);
% aux1 = E12(1:TimDel);
% f=time2freq(t);
  %Secondly the signal will be delayed by rotation its possition. That 
  %means a circular rith-shift will be performed. The number of possition 
  %to be deslocated will be given by the number of point that represent the 
  %delay based on the input time vector given as input. The time delay
  %time accuracy depends on the time vector accuracy.
% E12d = [E12(end-(length(aux1)-1):end) E12(1:end-(length(aux1)))];
E12d = ifft(fft(E12).*exp(Nimag*G2*Gpi*TimDel.*f));

  %In practice the phase delay will be perfomed by an phase controler. In
  %the simulation it will be done by multiplying the input signal by a
  %complex neperian exponential where the argument was received as an input
  %parameters.
  %The nema notation means that the E12 was Delayed thus E12d
E12d = E12d.*exp(Nimag*PhaDel);

%% * Stage: 3 Recombining both signals
 %Combining singnals
 %Similarly to te input step. An 3dB coupler with 2 inputs and 2 outputs is
 %used to combine the input signals from the arm1 and arm2. The couples 
 %will output two signals that are a combination of E12 and E12d.
 %Thus the Eout1 is the combination of E11 and E12d wich will be out put
 %at the port 1 of the couples. Hence, Eout2 is the combination of E11 and
 %E12d wich will be out put at the port 2.
Eout1 = ExpCoef.*(cos(CoupConst).*E11 + MutCoef*Pimag*sin(CoupConst).*E12d);
Eout2 = ExpCoef.*(cos(CoupConst).*E12d + MutCoef*Pimag*sin(CoupConst).*E11);

end

