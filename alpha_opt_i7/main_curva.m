%% Código para levantar a curva caracteristica do MZM de Vbias por PO
clear,clc,close all;

phase_steps = 1;
Vbias_steps = 100;
Vbias = linspace(-10,10,Vbias_steps);
% Make_MZ_Input_Files;

f_RF   = 1E6;   % RF frequency [Hz]
n_sc   = 10;    % Number of required subcarriers
spare  = 1;    % Frequency spare [%]
N      = 2^10;  % Number of points in t and (of coarse) in f
%
f_max  = (1 + spare/100)*n_sc*(4*f_RF);
df     = 2*f_max/(N - 1);
dt     = ((N - 1)/N)/(2*f_max);
t_max  = (N - 1)*dt;
%
t      = 0:dt:(N - 1)*dt;
f      = time2freq(t);
A = 0.5*exp(1j*(2*pi*f_RF*t));
B = 0.5 + 1j*0.5;
C = A + B;
aa=1;
EleSig2 = sin(2*pi*f_RF*t);
EleSig1 = zeros(1,length(t));

EleSig.U1t=EleSig1;
EleSig.U2t=EleSig2;

% CW = C;
CW = ones(1,length(t));
Eout = zeros(1,length(t));
H = zeros(1,length(t));
Powerdb = zeros(1,Vbias_steps);
PowerW = zeros(1,Vbias_steps);
PowerI = zeros(1,Vbias_steps);
PowerIdb = zeros(1,Vbias_steps);
%%
for kk=1:Vbias_steps
    [MZ_Input_Sdata] = ['AA_MZ_Input_Data_t' num2str(kk)];
    [Eout,H] = Mach_Zehnder_Modulator_simplificado(t,CW,EleSig,...
                                                           MZ_Input_Sdata);
    Powerdb(kk) = mean(20.*log10(sqrt(Eout.*conj(Eout))));
    PowerIdb(kk) = mean(20.*log10(sqrt(CW.*conj(CW))))-Powerdb(kk);
    PowerW(kk) = mean((sqrt(Eout.*conj(Eout))));
    PowerI(kk) = mean((sqrt(CW.*conj(CW))))-PowerW(kk);
end
Power = [PowerW(1:31) -1.*PowerW(32:69) PowerW(70:end)];
figure;
plot(Vbias,Power,'o','color',[0 0.6 1],'LineWidth',2);
hold on;
plot(Vbias,(PowerW),'color',[1 0.4 0],'LineWidth',2);
plot(Vbias,(PowerI),'color',[0 0 0],'LineWidth',2);
title('Potência de Saída Vs Tensão de Polarização','FontSize',20);
xlabel('Tensão de Polarização [V]','FontSize',16);
ylabel('Potência de Saída [W]','FontSize',16);
grid on;
% Power = [PowerW(1:31) -1.*PowerW(32:69) PowerW(70:end)];
% figure;
% plot(Vbias,(1+Power)./max((1+Power)),'o','color',[0 0.6 1],'LineWidth',2);
% hold on;
% plot(Vbias,(1+PowerW)./max((1+PowerW)),'color',[1 0.4 0],'LineWidth',2);
% title('Potência de Saída Vs Tensão de Polarização','FontSize',20);
% xlabel('Tensão de Polarização [V]','FontSize',16);
% ylabel('Potência de Saída [W]','FontSize',16);
% grid on;
print(['-f' num2str(1)],['figura' num2str(1)],'-dpng');
a=1;