function [ Eout ] = DpskEncodEq( TxData,PreBits )
    if nargin<2
        PreBits = 0;
    end
    Eout = [];
    if mod(length(TxData),2)
        TxData = [TxData 0];
    end

    for kk=1:length(TxData)    
        a    = TxData(kk);
        pk_1 = PreBits(1);

        pk   = ((~a)&(~pk_1))|((a)&(pk_1));

        Eout = [Eout pk];
        PreBits(1) = pk;
    end
%     TxData(1:2) = [];
%     
%     if ~isempty(TxData)
%         [EoutIaux,EoutQaux]=DqpskEncodEq( TxData,PreBits );
%         EoutI = [EoutI EoutIaux];
%         EoutQ = [EoutQ EoutQaux];
%     end
%     Resized = (Res1.*2^1 + Res2.*2^0) + 1;
    
end