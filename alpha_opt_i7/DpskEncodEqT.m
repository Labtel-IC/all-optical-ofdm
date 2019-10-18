function [ Eout ] = DpskEncodEqT( TxData,PreBits )
    if nargin<2
        PreBits = zeros(1,size(TxData,2));
    end
    Eout = zeros(size(TxData,1),size(TxData,2));
    if mod(size(TxData,1),2)
        TxData(end+1,:) = 0;
    end

    for kk=1:size(TxData,1)    
        a    = TxData(kk,:);
        pk_1 = PreBits(1,:);

        pk   = ((~a)&(~pk_1))|((a)&(pk_1));

        Eout(kk,:)   = pk;
        PreBits(1,:) = pk;
    end    
end