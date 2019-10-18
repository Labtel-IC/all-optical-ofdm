close all;clear;clc;
NormAmos = 3e6;
SeleAmos = 3;
SampSize = NormAmos+SeleAmos;
BoxSize = 2^10;
count=1;
tic
for kk=1:BoxSize-2
    for jj=kk+1:BoxSize-1
        for ll=jj+1:BoxSize
            p1=SeleAmos/(SampSize-(kk-1));
            p2=(SeleAmos-1)/(SampSize-(jj-1));
            p3=(SeleAmos-2)/(SampSize-(ll-1));
            p(count) = p1*p2*p3;
            count=count+1;
        end
    end
end
sum(p)
a=1;
