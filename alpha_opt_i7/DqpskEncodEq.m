function [ EoutI,EoutQ ] = DqpskEncodEq( TxData,PreBits )
%%                   DQPSK Enconder Through Equations
%c function [ EoutI,EoutQ ] = DqpskEncodEq( TxData,PreBits )
%c
%c This function will actualy create the I and Q components for the DQPSK 
%c modulation. This equation can be found at the book of Optical Fiber 
%c Telecommunications V B, which one of the authors is Ivan P. Kaminow at 
%c the page 144.
%c
%c
%c                                           Created by P.Marciano LG
%c                                           02/11/2017
%c                                           Last Update
%c                                           07/01/2017
%c                                           pablorafael.mcx@gmail.com
%c
%c Refences:
%c @book{kaminow2010optical,
%c   title={Optical fiber telecommunications VB: systems and networks},
%c   author={Kaminow, Ivan and Li, Tingye and Willner, Alan E},
%c   year={2010},
%c   publisher={Elsevier}
%c }
%c
%c
%c
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                    G L O B A L  V A R I A B L E S                      %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                   L I S T  O F  V A R I A B L E S                      %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c
%c    DqpskEncodEq
%c
%c  Input:
%c  TxData  : The actual data to be converted into I and Q componentes  [b]
%c  PreBits : The start point of the converssion                        [b]
%c  
%c  Output:
%c  EoutI   : The I componente created
%c  EoutQ   : The Q componente created
%c  
%c
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                             S T R U C T S                              %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c
%%                      Initializing
%This first step is to check if all variables were correctly passed as
%inputs and if they are accordingly with the expected.
    if nargin<2                                                            %If the input given was just the Txdata set the start point
        PreBits = [0 0];
    end
    EoutI = [];                                                            %Initializing the I output.
    EoutQ = [];                                                            %Initializing the Q output.
    if mod(length(TxData),2)                                               %Checking if the input data is even
        TxData = [TxData 0];                                               %Adjusting the incoming data accordingly
    end
    %The following process takes on account pair of bits, therefore the 
    %income data was split for convenience sake. 
    ResLen = linspace(1,length(TxData),length(TxData));                    %Taking the possition references
    ResB = TxData(~mod(ResLen,2));                                         %Taking the data at even positions
    ResA = TxData(mod(ResLen,2)==1);                                       %Taking the data at odd positions

    for kk=1:length(ResA)                                                  %Actualy performing the equation for each pair of bits
        %Acquiring the Current bits
        a = ResA(kk);
        b = ResB(kk);
        %Acquiring the Previous bits
        pk_1 = PreBits(1);
        qk_1 = PreBits(2);
        %Calculating the current I bit
        pk = ((~a)&(~b)&(~pk_1))|((~a)&(b)&(qk_1))|((a)&(~b)&(~qk_1))|((...
                                                            a)&(b)&(pk_1));
        %Calculating the current Q bit
        qk = ((~a)&(~b)&(~qk_1))|((~a)&(b)&(~pk_1))|((a)&(~b)&(pk_1))|((...
                                                            a)&(b)&(qk_1));
        %Construction of I and Q Components
        EoutI = [EoutI pk];
        EoutQ = [EoutQ qk];
        %Updating the Previous bits
        PreBits(1) = pk;
        PreBits(2) = qk;
    end
end