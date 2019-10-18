close all;clear;clc;
files_vet = [1 2 3 4 5 6 7 8 9];
VetBer = [];
VetErec = [];
for FileToLoad = files_vet
    load(['Teste' num2str(FileToLoad) '_4PAM_50km.mat']);
    VetBer  = [VetBer;Ber4PAM];
    for kk=1:size(ValsLev,3)
        for jj=1:size(ValsLev,1)
            auxcont = 1;
            for ii=1:2:size(ValsLev,2)
                MeanLevVal(jj,auxcont)=ValsLev(jj,ii,kk);
                auxcont = auxcont + 1;
            end
        end
    end
%     VetErec = [VetErec;EoutRec];
end

ValsPos = linspace(1,size(MeanLevVal,2),size(MeanLevVal,2));
MeanVetBer = mean(VetBer);
CarPos = 1:2:128;
LimMark = linspace(3e-3,3e-3,length(CarPos));
figure;
plot(CarPos,MeanVetBer(berpos),CarPos,LimMark);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
figure;
hold all;
for kk=1:size(MeanLevVal,1)
    plot(ValsPos,MeanLevVal(kk,:),[1 size(MeanLevVal,2)], [mean(MeanLevVal(kk,:)) mean(MeanLevVal(kk,:))]);
end
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
a=1;