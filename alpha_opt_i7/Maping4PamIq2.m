function [Phi1,Phi2]=Maping4PamIq2(TxData,Vmin,Vmax)
%%                Creating the 4PAM electrical signals.  
%c function [PamSig]=Maping4Pam(TxData,NPPB,Polarized,MaxAmp)
%c
%c
%c This function is resposible for generating the signal wich will
%c modulate, through the MZM-I, our carrier. Basicaly it will split the 
%c ouput power level of the MZM in four parts. Remembering that each level
%c will encode two data bits. The main idea here is not to create the four
%c levels at the eletrical domain but rather at the optical domain. It can
%c be more clear if we analyze the characterist curve of an MZM
%c PowerVsBias. If we were to chose to create the 4PAM signal at the
%c eletrical domain with an even distribuited optical output we should use
%c the linear part of the curve which is very restricted. But if we were
%c willing to make use of the whole available band we needed to linearize
%c the curve. It can be achieved by dividint the optical power output in
%c four parts and finding the correspondent intervals for the polarization
%c voltage. Finaly we can create two signals that will be the input of an
%c DD-MZM that will output four levels equaly spaced.
%c
%c
%c                                           Created by P.Marciano LG
%c                                           05/01/2018
%c                                           07/01/2018
%c                                           pablorafael.mcx@gmail.com
%c
%c Refences:
%c@article{xu2017optical,
%c   title={Optical interferometric synthesis of PAM4 signals based on dual-drive Mach--Zehnder modulation},
%c   author={Xu, Jianfeng and Du, Jiangbing and Ren, Rongrong and Ruan, Zhengshang and He, Zuyuan},
%c   journal={Optics Communications},
%c   volume={402},
%c   pages={73--79},
%c   year={2017},
%c   publisher={Elsevier}
%c }
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
%c    Maping4PamMzm
%c
%c  Input:
%c  TxData : The complete tada to be transmited                         [b]
%c  Vref1  : The two levels that the first electrical signal can assume [V]
%c  Vref2  : The two levels that the second electrical signal can have  [V]
%c  
%c  Output:
%c  Phi1   : Electrical signal at one driver of the DD-MZM              [V]
%c  Phi2   : Electrical signal at other driver of the DD-MZM            [V]
%c  
%c
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                             S T R U C T S                              %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c
%%                       Initializing
%This first step is to check if all variables were correctly passed as
%inputs and if they are accordingly with the expected.
    if nargin==1
        Vmin = 0;
        Vmax = 1;
    end
%     Vref1 = [Vmin 0.4 Vmax];
%     Vref2 = [Vmin 1.2 Vmax];
    %Bom para 100 km ou mais
%     Vref1 = [Vmax/2 -0.4*Vmax -1*Vmax/2];%
%     Vref2 = [Vmax/2 1.1*Vmax -1*Vmax/2];%
    %Bom para outros KMs
    Vref1 = [Vmin 0.2105*Vmax Vmax];%
    Vref2 = [Vmin 0.6316*Vmax Vmax];%
%     Vref1 = [0 0.8];
%     Vref2 = [1.3 1.9];
%     Vref1 = [0 1.75 1.9];
%     Vref2 = [0 1.12 1.9];
    Phi1 = [];                                                             %Initializing the Signal 1 output.
    Phi2 = [];                                                             %Initializing the Signal 2 output.
    if mod(length(TxData),2)                                               %Checking if the input data is even
        TxData = [TxData 0];                                               %Adjusting the incoming data accordingly
    end
    %The following process takes on account pair of bits, therefore the 
    %income data was split for convenience sake. 
    ResLen = linspace(1,length(TxData),length(TxData));                    %Taking the possition references
    ResB = TxData(~mod(ResLen,2));                                         %Taking the data at even positions
    ResA = TxData(mod(ResLen,2)==1);                                       %Taking the data at odd positions

    for kk=1:length(ResA)                                                  %Actualy performing the signal creation for each pair of bits
        %The Most external verification looks at the data at odd positions 
        %whereas the inner verification of the data at even positions.
        if ResA(kk)%if 1
            if ResB(kk)%if 1 - Level 3
                Phi1 = [Phi1 Vref1(1)];
                Phi2 = [Phi2 Vref2(2)];
            else%if 0 - Level 4
                Phi1 = [Phi1 Vref1(1)];
                Phi2 = [Phi2 Vref2(1)];
            end
        else%if 0
            if ResB(kk)%if 1 - Level 1
                Phi1 = [Phi1 Vref1(2)];
                Phi2 = [Phi2 Vref2(1)];
            else%if 0 - Level 0
                Phi1 = [Phi1 Vref1(3)];
                Phi2 = [Phi2 Vref2(3)];
            end
        end
    end
end