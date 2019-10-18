function [Phi1,Phi2]=Maping4PamIq(TxData,Vmin,Vmax,ModSchem,MediumLength,SetCpSampZer)
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
%c  MediumLength : The length of the fiber channel                     [km]
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
    
    if ModSchem
%         Vref1 = [0.85*Vmin 1.46*Vmin 1.9*Vmin 1.71*Vmin];%[0.65*Vmin 0.885*Vmin 1.12*Vmin 1.355*Vmin];%
%         Vref2 = [0.85*Vmax 1.46*Vmax 1.9*Vmax 1.71*Vmax];%[0.65*Vmax 0.885*Vmax 1.12*Vmax 1.355*Vmax];%
        Vref1 = [0.85*Vmin 1.77*Vmin 2.85*Vmin 2.202*Vmin];%[0.65*Vmin 0.885*Vmin 1.12*Vmin 1.355*Vmin];%
        Vref2 = [0.85*Vmax 1.77*Vmax 2.85*Vmax 2.202*Vmax];%[0.65*Vmax 0.885*Vmax 1.12*Vmax 1.355*Vmax];%
%         Vref1 = [0.0*Vmin 0.82*Vmin 1.9*Vmin 1.24*Vmin];%[Vmin 0.2105*Vmax Vmax];%[Vmin 0.9*Vmax Vmax];%
%         Vref2 = [0.0*Vmin 0.82*Vmax 1.9*Vmax 1.24*Vmax];%[Vmin 0.6316*Vmax Vmax];%[Vmin 0.55*Vmax Vmax];%
    else
%         if MediumLength <= 10
%             Vref2 = [0.43*Vmax Vmin Vmin Vmax/2];
%             Vref1 = [0.78*Vmax Vmax Vmax 0.8*Vmax];
%         elseif MediumLength <= 25
%             Vref2 = [0.40*Vmax Vmin Vmin Vmax/2];
%             Vref1 = [0.77*Vmax Vmax Vmax 0.8*Vmax];
%         else
        if  MediumLength >50
%             Vref1 = [0.925 1.27 1.27];
%             Vref2 = [0.5 0 0];
            if SetCpSampZer
                Vref1 = [0.76*Vmax 0.975*Vmax 0.975*Vmax 0.76*Vmax];
                Vref2 = [0.4*Vmax 0.4*Vmax Vmin Vmin];
            else
                Vref1 = [0.74*Vmax 1*Vmax 1*Vmax 0.74*Vmax];
                Vref2 = [0.4*Vmax 0.4*Vmax Vmin Vmin];
            end
%             Vref1 = [0.74*Vmax Vmax Vmax 0.8*Vmax];
%             Vref2 = [0.4*Vmax Vmin Vmin Vmax/2];
        elseif  MediumLength >30
%             Vref1 = [0.925 1.27 1.27];
%             Vref2 = [0.5 0 0];
            if SetCpSampZer
                Vref1 = [0.765*Vmax 0.975*Vmax 0.975*Vmax 0.765*Vmax];
                Vref2 = [0.4*Vmax 0.4*Vmax Vmin Vmin];
            else
                Vref1 = [0.76*Vmax 1.01*Vmax 1.01*Vmax 0.76*Vmax];
                Vref2 = [0.4*Vmax 0.4*Vmax Vmin Vmin];
            end
        elseif  MediumLength >=10
%             Vref1 = [0.925 1.27 1.27];
%             Vref2 = [0.5 0 0];
            if SetCpSampZer
                Vref1 = [0.76*Vmax 0.975*Vmax 0.975*Vmax 0.76*Vmax];
                Vref2 = [0.4*Vmax 0.4*Vmax Vmin Vmin];
            else
                Vref1 = [0.79*Vmax 1.01*Vmax 1.01*Vmax 0.79*Vmax];
                Vref2 = [0.4*Vmax 0.4*Vmax Vmin Vmin];
            end
        else
            Vref1 = [-1.0 -2.45 -2.45 -1.0];
            Vref2 = [1.0 1.0 3.25 3.25];
%             Vref1 = [-1.35 -2.45 -2.45 -1.35];
%             Vref2 = [1.35 1.35 3.25 3.25];
%             Vref1 = [-1.3 -2.45 -2.45 -1.3];
%             Vref2 = [1.3 1.3 3.25 3.25];
        end
    end
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
                if ModSchem
                    Phi1 = [Phi1 Vref1(2)];
                    Phi2 = [Phi2 Vref2(2)];
                else
                    Phi1 = [Phi1 Vref1(2)];
                    Phi2 = [Phi2 Vref2(2)];
%                     Phi1 = [Phi1 Vref1(2)];
%                     Phi2 = [Phi2 Vref2(1)];
                end
            else%if 0 - Level 4
                Phi1 = [Phi1 Vref1(1)];
                Phi2 = [Phi2 Vref2(1)];
            end
        else%if 0
            if ResB(kk)%if 1 - Level 1
                if ModSchem
                    Phi1 = [Phi1 Vref1(4)];
                    Phi2 = [Phi2 Vref2(4)];
                else
                    Phi1 = [Phi1 Vref1(4)];
                    Phi2 = [Phi2 Vref2(4)];
%                     Phi1 = [Phi1 Vref1(1)];
%                     Phi2 = [Phi2 Vref2(2)];
                end
            else%if 0 - Level 0
                Phi1 = [Phi1 Vref1(3)];
                Phi2 = [Phi2 Vref2(3)];
            end
        end
    end
end