close all;clear;clc;
prebit = [0 0 0 1 1 1 1 0];
txdata = [0 0 0 1 1 1 1 0];% 0 1 1 0 1 1];
for kk=1:2:length(prebit)
    for jj=1:2:length(txdata)
        prebit(kk:kk+1)
        txdata(jj:jj+1)
        [I,Q] = DqpskEncodEq(txdata(jj:jj+1),prebit(kk:kk+1))
        a=1;
    end
end

for kk=1:2:length(prebit)
    for jj=1:2:length(txdata)
        prebit(kk:kk+1)
        txdata(jj:jj+1)
        [I,Q] = DqpskEncodNow(txdata(jj:jj+1),prebit(kk:kk+1))
        a=1;
    end
end

for kk=1:2:length(prebit)
    for jj=1:2:length(txdata)
        prebit(kk:kk+1)
        txdata(jj:jj+1)
        [I,Q] = DqpskEncod(txdata(jj:jj+1),prebit(kk:kk+1))
        a=1;
    end
end