clear all;close all;clc;

FilesToOpen = [{'GenerateOFC.m'},{'GenerateOfcInputData.m'},{'OFC_RIM.m'}];
for kk=1:size(FilesToOpen,2)
    open(char(FilesToOpen(kk)))
end