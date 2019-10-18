function [ EoutI,EoutQ ] = DqpskEncod( TxData,PreBits )
    if nargin<2
        PreBits = [0 0];
    end
    EoutI = [];
    EoutQ = [];
    if mod(length(TxData),2)
        TxData = [TxData 0];
    end
    
    ResLen = linspace(1,length(TxData),length(TxData));
    Res2 = TxData(~mod(ResLen,2));
    Res1 = TxData(mod(ResLen,2)==1);
    Resized = (Res1.*2^1 + Res2.*2^0) + 1;
    
    for kk=1:length(Resized)
        LookTable = [num2str((PreBits(1)*2^1+PreBits(2)*2^0)+1) num2str(Resized(kk))];
        switch LookTable
            case '11'
                EoutI = [EoutI 0];
                EoutQ = [EoutQ 0];
                PreBits(1:2) = [ 0,0]; 
            case '12'
                EoutI = [EoutI 0];
                EoutQ = [EoutQ 1];
                PreBits(1:2) = [ 0,1];
            case '13'
                EoutI = [EoutI 1];
                EoutQ = [EoutQ 0];
                PreBits(1:2) = [ 1,0];
            case '14'
                EoutI = [EoutI 1];
                EoutQ = [EoutQ 1];
                PreBits(1:2) = [ 1,1];
            case '21'
                EoutI = [EoutI 0];
                EoutQ = [EoutQ 1];
                PreBits(1:2) = [ 0,1];
            case '22'
                EoutI = [EoutI 1];
                EoutQ = [EoutQ 1];
                PreBits(1:2) = [ 1,1];
            case '23'
                EoutI = [EoutI 0];
                EoutQ = [EoutQ 0];
                PreBits(1:2) = [ 0,0];
            case '24'
                EoutI = [EoutI 1];
                EoutQ = [EoutQ 0];
                PreBits(1:2) = [ 1,0];
            case '31'
                EoutI = [EoutI 1];
                EoutQ = [EoutQ 0];
                PreBits(1:2) = [ 1,0];
            case '32'
                EoutI = [EoutI 0];
                EoutQ = [EoutQ 0];
                PreBits(1:2) = [ 0,0];
            case '33'
                EoutI = [EoutI 1];
                EoutQ = [EoutQ 1];
                PreBits(1:2) = [ 1,1];
            case '34'
                EoutI = [EoutI 0];
                EoutQ = [EoutQ 1];
                PreBits(1:2) = [ 0,1];
            case '41'
                EoutI = [EoutI 1];
                EoutQ = [EoutQ 1];
                PreBits(1:2) = [ 1,1];
            case '42'
                EoutI = [EoutI 1];
                EoutQ = [EoutQ 0];
                PreBits(1:2) = [ 1,0];
            case '43'
                EoutI = [EoutI 0];
                EoutQ = [EoutQ 1];
                PreBits(1:2) = [ 0,1];
            case '44'
                EoutI = [EoutI 0];
                EoutQ = [EoutQ 0];
                PreBits(1:2) = [ 0,0];
            otherwise
                error('Could not find the combination');
        end
                
    end
end