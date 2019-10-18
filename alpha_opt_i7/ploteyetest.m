close all;clear;clc;
load('savingforplotingeye240.mat');
Nb = 1024;
TotalSamples = Nb*NPPB;
% t
PerUtil         = 0.122;                                               %Percentage for the ultil period of symbol
NumAmosCP       = ceil(PerUtil*NPPB);                                  %Half of the number of samples for the CP
NewTotalSamples = (2*NumAmosCP+NPPB)*Nb;                               %New total samples of the simulation
ExtraSamples    = NewTotalSamples - TotalSamples;                      %Extra samples created
NumBitDesc      = ceil(ExtraSamples/(2*NumAmosCP+NPPB));               %Number of bits excided
if mod(NumBitDesc,2)                                                   %Verify if the amount of bits to remove is odd
    NumBitDesc = NumBitDesc + 1;                                       %This simulation works just with even number of bit stream
end                                                                    %Thus, adding one more unit to be removed is needed
StuffSampels    = TotalSamples - (2*NumAmosCP+NPPB)*(Nb - NumBitDesc); %Number of samples for stuffing process
% Tcp             = t(2)*(2*NumAmosCP+NPPB);
    
AuxPlot = [];
figure;
hold all;
for kk=1:size(EyeToPlot,1)
    AuxPlot = [AuxPlot EyeToPlot(kk,1+2*(2*NumAmosCP+NPPB):end-StuffSampels-0)];
%     mod(size(AuxPlot,2),(2*NumAmosCP+NPPB))
%     plot(linspace(0,T*(size(EyeToPlot,2)/NPPB),size(EyeToPlot,2)),EyeToPlot(kk,:));
end
size(AuxPlot,2)
tmax = 3;
mod(size(AuxPlot,2),(tmax*(2*NumAmosCP+NPPB)))
Olho(AuxPlot,T,(2*NumAmosCP+NPPB),1,tmax);
axis([0.45e-10 2.05e-10 -1 1]);

clear;

load('savingforplotingeye10.mat');
Nb = 1024;
TotalSamples = Nb*NPPB;
% t
PerUtil         = 0.122;                                               %Percentage for the ultil period of symbol
NumAmosCP       = ceil(PerUtil*NPPB);                                  %Half of the number of samples for the CP
NewTotalSamples = (2*NumAmosCP+NPPB)*Nb;                               %New total samples of the simulation
ExtraSamples    = NewTotalSamples - TotalSamples;                      %Extra samples created
NumBitDesc      = ceil(ExtraSamples/(2*NumAmosCP+NPPB));               %Number of bits excided
if mod(NumBitDesc,2)                                                   %Verify if the amount of bits to remove is odd
    NumBitDesc = NumBitDesc + 1;                                       %This simulation works just with even number of bit stream
end                                                                    %Thus, adding one more unit to be removed is needed
StuffSampels    = TotalSamples - (2*NumAmosCP+NPPB)*(Nb - NumBitDesc); %Number of samples for stuffing process
% Tcp             = t(2)*(2*NumAmosCP+NPPB);
    
AuxPlot = [];
figure;
hold all;
for kk=1:size(EyeToPlot,1)
    AuxPlot = [AuxPlot EyeToPlot(kk,1+2*(2*NumAmosCP+NPPB):end-StuffSampels-0)];
%     mod(size(AuxPlot,2),(2*NumAmosCP+NPPB))
%     plot(linspace(0,T*(size(EyeToPlot,2)/NPPB),size(EyeToPlot,2)),EyeToPlot(kk,:));
end
size(AuxPlot,2)
tmax = 3;
mod(size(AuxPlot,2),(tmax*(2*NumAmosCP+NPPB)))
Olho(AuxPlot,T,(2*NumAmosCP+NPPB),1,tmax);
axis([0.45e-10 2.05e-10 -1 1]);

clear;

load('savingforplotingeye110.mat');
Nb = 1024;
TotalSamples = Nb*NPPB;
% t
PerUtil         = 0.122;                                               %Percentage for the ultil period of symbol
NumAmosCP       = ceil(PerUtil*NPPB);                                  %Half of the number of samples for the CP
NewTotalSamples = (2*NumAmosCP+NPPB)*Nb;                               %New total samples of the simulation
ExtraSamples    = NewTotalSamples - TotalSamples;                      %Extra samples created
NumBitDesc      = ceil(ExtraSamples/(2*NumAmosCP+NPPB));               %Number of bits excided
if mod(NumBitDesc,2)                                                   %Verify if the amount of bits to remove is odd
    NumBitDesc = NumBitDesc + 1;                                       %This simulation works just with even number of bit stream
end                                                                    %Thus, adding one more unit to be removed is needed
StuffSampels    = TotalSamples - (2*NumAmosCP+NPPB)*(Nb - NumBitDesc); %Number of samples for stuffing process
% Tcp             = t(2)*(2*NumAmosCP+NPPB);
    
AuxPlot = [];
figure;
hold all;
for kk=1:size(EyeToPlot,1)
    AuxPlot = [AuxPlot EyeToPlot(kk,1+2*(2*NumAmosCP+NPPB):end-StuffSampels-0)];
%     mod(size(AuxPlot,2),(2*NumAmosCP+NPPB))
%     plot(linspace(0,T*(size(EyeToPlot,2)/NPPB),size(EyeToPlot,2)),EyeToPlot(kk,:));
end
size(AuxPlot,2)
tmax = 3;
mod(size(AuxPlot,2),(tmax*(2*NumAmosCP+NPPB)))
Olho(AuxPlot,T,(2*NumAmosCP+NPPB),1,tmax);
axis([0.4e-10 2e-10 0 1]);

clear;

load('savingforplotingeye60.mat');
Nb = 1024;
TotalSamples = Nb*NPPB;
% t
PerUtil         = 0.122;                                               %Percentage for the ultil period of symbol
NumAmosCP       = ceil(PerUtil*NPPB);                                  %Half of the number of samples for the CP
NewTotalSamples = (2*NumAmosCP+NPPB)*Nb;                               %New total samples of the simulation
ExtraSamples    = NewTotalSamples - TotalSamples;                      %Extra samples created
NumBitDesc      = ceil(ExtraSamples/(2*NumAmosCP+NPPB));               %Number of bits excided
if mod(NumBitDesc,2)                                                   %Verify if the amount of bits to remove is odd
    NumBitDesc = NumBitDesc + 1;                                       %This simulation works just with even number of bit stream
end                                                                    %Thus, adding one more unit to be removed is needed
StuffSampels    = TotalSamples - (2*NumAmosCP+NPPB)*(Nb - NumBitDesc); %Number of samples for stuffing process
% Tcp             = t(2)*(2*NumAmosCP+NPPB);
    
AuxPlot = [];
figure;
hold all;
for kk=1:size(EyeToPlot,1)
    AuxPlot = [AuxPlot EyeToPlot(kk,1+2*(2*NumAmosCP+NPPB):end-StuffSampels-0)];
%     mod(size(AuxPlot,2),(2*NumAmosCP+NPPB))
%     plot(linspace(0,T*(size(EyeToPlot,2)/NPPB),size(EyeToPlot,2)),EyeToPlot(kk,:));
end
size(AuxPlot,2)
tmax = 3;
mod(size(AuxPlot,2),(tmax*(2*NumAmosCP+NPPB)))
Olho(AuxPlot,T,(2*NumAmosCP+NPPB),1,tmax);
axis([0.4e-10 2e-10 0 1]);
a=1;