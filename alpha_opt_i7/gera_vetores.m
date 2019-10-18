clear;close all;clc;
%% Gera vetores de forma automática
% Vetor de arestas
VetAre = [2 4 4 6 6 8 3 5 6 7 7 12 12 13 16 18 12 15 15 17];
% Número maximo de modens entre arestas
N = 6;
% Número aleatório de roteadores entre arestas
NumRot = randi(N+1,[1 10])-1
% Vetor que representa o grafo da rede
NewVet = [];
count = 1;
for kk=1:length(NumRot)
%     if mod(kk,2)
    NewVet = [NewVet [VetAre(count) ones(1,NumRot(kk)) VetAre(count+1)]];
    count = count + 2;
%     end
end
NewVet
a=1;