%% Test Bundle for different OCS created
% This is an simple script where will be evaluate different OCS created
% with different input power and different OSNR Values
close all; clear; clc;tic;
OcsToTest = [0 10 20 30 40 50 60 70 80 90 100 108 116 124 132 140 ...
                      156 164 168 200 205 210 215 220 225 230 235 240 245];%Vector to store the index of the OCS to be tested
%%Main Loop of this script
for CurrentOCS=1:length(OcsToTest)
    clearvars -except CurrentOCS OcsToTest BerDPSK BerDQPSK Ber4PAM BerOOK
    MAIN;
    toc
    CurrentOCS*100/length(OcsToTest)
    tic
end

save(['BerVarForB2B100_' date]);
s=1;