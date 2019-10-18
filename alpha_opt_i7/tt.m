close all;clear;clc;
MZ_Input_File = 1;
fc = 12.5e9;
NPPB = 2^5;
NB   = 2^10;
amto = NB*NPPB;
t = linspace(0,NB/fc,amto);
EoutTx = ones(1,size(t,2));
U.U1t = zeros(1,size(t,2));
U.U2t = zeros(1,size(t,2));
Vbias_steps = 100;
Vbias = linspace(-3.8,3.8,Vbias_steps);
L        =  10.00;
% U0       =  2.390;
U_pi1    =   3.80;
U_pi2    =   3.80;
eta1     =  89.00;%1;%
eta2     = -12.70;%-1;%
nopt     =   2.17;
nel      =   2.60;
alfa_ins =   5.10;
phi_0    =   0.00;
alfa0    =   1.07;

p = pwd;                    %Geting the current path of the main program
Local = [p '\input_files\'];%variable that shows where the input files for
                                                     %the mzm will be saved
S=1; %Variabel that give the name for the input file
for Inc2=1:Vbias_steps%Secundary loop for the possible combinations
    U0               = Vbias(Inc2);%Variation of the Vbias
    [MZ_Input_Sdata] = Set_MZ_Input_Data_Simp(S,L,U0,U_pi1,U_pi2,eta1,...  %Functionresponsible to create the input files
                                 eta2,nopt,nel,alfa_ins,phi_0,alfa0,Local);
    S                = S+1;
end
figure;
hold all;
grid on;

for kk=1:Vbias_steps
    [EoutModAux] = Mach_Zehnder_Modulator_simplificado(t,EoutTx,U,kk);
    power_theo(kk) = mean(abs(EoutModAux).^2);  
    plot(Vbias(kk),power_theo(kk),'x');
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    drawnow;
%     pause(0.5);
end