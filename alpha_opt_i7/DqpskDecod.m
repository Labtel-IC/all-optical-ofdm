function [ Eout ] = DqpskDecod( RxRecI,RxRecQ,PreBits )
    if nargin<3
        PreBits = [0 0];
    end
    Eout = [];
    
    for kk=1:length(RxRecI)
        LookTable = [num2str((PreBits(1)*2^1+PreBits(2)*2^0)+1) num2str((RxRecI(kk)*2^1+RxRecQ(kk)*2^0)+1)];
        switch LookTable
            case '11'
                Eout = [Eout 0 0];
                PreBits(1:2) = [0,0]; 
            case '12'
                Eout = [Eout 0 1];
                PreBits(1:2) = [0,1]; 
            case '13'
                Eout = [Eout 1 0];
                PreBits(1:2) = [1,0]; 
            case '14'
                Eout = [Eout 1 1];
                PreBits(1:2) = [1,1]; 
            case '21'
                Eout = [Eout 1 0];
                PreBits(1:2) = [1,0]; 
            case '22'
                Eout = [Eout 0 0];
                PreBits(1:2) = [0,0];
            case '23'
                Eout = [Eout 1 1];
                PreBits(1:2) = [1,1]; 
            case '24'
                Eout = [Eout 0 1];
                PreBits(1:2) = [0,1]; 
            case '31'
                Eout = [Eout 0 1];
                PreBits(1:2) = [0,1]; 
            case '32'
                Eout = [Eout 1 1];
                PreBits(1:2) = [1,1]; 
            case '33'
                Eout = [Eout 0 0];
                PreBits(1:2) = [0,0]; 
            case '34'
                Eout = [Eout 1 0];
                PreBits(1:2) = [1,0]; 
            case '41'
                Eout = [Eout 1 1];
                PreBits(1:2) = [1,1]; 
            case '42'
                Eout = [Eout 1 0];
                PreBits(1:2) = [1,0]; 
            case '43'
                Eout = [Eout 0 1];
                PreBits(1:2) = [0,1]; 
            case '44'
                Eout = [Eout 0 0];
                PreBits(1:2) = [0,0]; 
            otherwise
                error('Could not find the combination');
        end
                
    end
end