%% Scrip to calculate the basic parameter to construct a power divider
close all,clc,clear;
f = 18e9;%Frequency operation in Hz
F = f/1e9;%Operation frequency in GHz
% W = 0.1;%Width of the cuper track mm
D = 1.55;%Thickness of the Board in mm
Er = 4.4;%Relative dieletrical value of the material
Z0 = 50;%The disered impedance of the wave guide in ohms
C = 3e8;%Speed of light m/s
K0 = (2*pi*f)*sqrt((pi*4e-7)*(8.854187817e-12));%Propagation constant for free space
A = (Z0/60)*sqrt((Er +1)/2) + ((Er-1)/(Er+2))*(0.23+(0.11/Er));
B = (377*pi)/(2*Z0*sqrt(Er));
W1 = (8*exp(A))/(exp(2*A)-2)
W2 = D*(2/pi)*(B-1-log(2*B-1)+((Er-1)/(2*Er))*(log(B-1)+0.39-(0.61/Er)))
Eeff1 = ((Er+1)/2) + (((Er-1)/2)*(1/(sqrt(1+((12*D)/W1)))))%Effective dielectric constant for the substrate
Eeff2 = ((Er+1)/2) + (((Er-1)/2)*(1/(sqrt(1+((12*D)/W2)))))%Effective dielectric constant for the substrate
vp1 = C/sqrt(Eeff1)%phase velocity
vp2 = C/sqrt(Eeff2)%phase velocity
Beta1 = K0*sqrt(Eeff1)%Propagation constant for the micro strip
Beta2 = K0*sqrt(Eeff2)%Propagation constant for the micro strip
LambdaM1 = 300/(F*sqrt(Eeff1));%The wave length in milimeters
LambdaM2 = 300/(F*sqrt(Eeff2));%The wave length in milimeters
LambdaM1_4 = LambdaM1/4
LambdaM2_4 = LambdaM2/4
W1/D
W2/D


f = 18e9;%Frequency operation in Hz
F = f/1e9;%Operation frequency in GHz
% W = 0.1;%Width of the cuper track mm
D = 1.55;%Thickness of the Board in mm
Er = 4.4;%Relative dieletrical value of the material
Z0 = sqrt(2)*50;%The disered impedance of the wave guide in ohms
C = 3e8;%Speed of light m/s
K0 = (2*pi*f)*sqrt((pi*4e-7)*(8.854187817e-12));%Propagation constant for free space
A = (Z0/60)*sqrt((Er +1)/2) + ((Er-1)/(Er+2))*(0.23+(0.11/Er));
B = (377*pi)/(2*Z0*sqrt(Er));
W1 = (8*exp(A))/(exp(2*A)-2)
W2 = D*(2/pi)*(B-1-log(2*B-1)+((Er-1)/(2*Er))*(log(B-1)+0.39-(0.61/Er)))
Eeff1 = ((Er+1)/2) + (((Er-1)/2)*(1/(sqrt(1+((12*D)/W1)))))%Effective dielectric constant for the substrate
Eeff2 = ((Er+1)/2) + (((Er-1)/2)*(1/(sqrt(1+((12*D)/W2)))))%Effective dielectric constant for the substrate
vp1 = C/sqrt(Eeff1)%phase velocity
vp2 = C/sqrt(Eeff2)%phase velocity
Beta1 = K0*sqrt(Eeff1)%Propagation constant for the micro strip
Beta2 = K0*sqrt(Eeff2)%Propagation constant for the micro strip
LambdaM1 = 300/(F*sqrt(Eeff1));%The wave length in milimeters
LambdaM2 = 300/(F*sqrt(Eeff2));%The wave length in milimeters
LambdaM1_4 = LambdaM1/4
LambdaM2_4 = LambdaM2/4
W1/D
W2/D
a=1;
% Err = 1;%Measurement error
% counter = 1;
% while ((W/D)<2) || Err
%     
%     if abs(Waux-W)<0.01
%         Err = 0;
%     end
%     if Waux>W
%         W = W + 0.01;
%     else
%         W = W - 0.01;
%     end
%     plot(counter,Waux,'o');
%     hold on;
%     counter = counter + 1;
%     drawnow;
%     a=1;
% end