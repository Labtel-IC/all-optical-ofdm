function [Eout,Count,MapVet]=OpticalFFTNT(f,SymPer,StopPoint,E1,AngVec,...
                                                   ActStag,Count,E0,MapVet)
%%              All Optical Fast Furier Transform for N carriers
%c function [Eout,Count,MapVet]=OpticalFFTN(t,SymPer,StopPoint,E1,...
%c                                          AngVec,ActStag,Count,E0,MapVet)
%c  This function is responsible for implementing an all-optical FFT. It 
%c is possíble to create an FFT of any size. Albeit, one must understand 
%c it’s limitations such as the granularity problem.  This process takes on 
%c account angel from 0 to 180 degrees and this interval is divided by 2 to 
%c the power of N where N is the order of the FFT. For instance, if the 
%c user need to filter 32 carriers it will have an FFT process of order 5 
%c hence the interval will be 5.6250, therefore in practice to implement 
%c the MZI it needs to accurately create phase delays spaced of 5.6250. 
%c Another problem is, the carriers are shuffled by this FFT, and for a 
%c perfect FFT implementation, the exact phase delay must be added to the 
%c right signal at the right time. This function also returns a vector map 
%c that indicates where each carrier is, thus the user can use this vector 
%c to correctly extract the right carrier. The following illustration 
%c better demonstrates the output mapping vector for an FFT of order 5.
%c
%c               Ein : Is the input fild
%c               MZI : Is the Mach-Zehnder Interferometre
%c               Nº  : The number are each carrier from 1 to 32 within Ein
%c
%c                                                  /1
%c                                             (MZI)
%c                                            / 1   \17
%c                                       (MZI)      
%c                                      /     \ 9   /9
%c                                     / 1     (MZI)
%c                                (MZI)             \25
%c                                /    \ 5          /5
%c                               /      \      (MZI)
%c                              /        \    / 5   \21
%c                             /         (MZI)
%c                            /               \ 13  /13
%c                           / 1               (MZI)
%c                       (MZI)                      \29
%c                       /   \ 3                    /3
%c                      /     \                (MZI)
%c                     /       \              / 3   \19
%c                    /         \        (MZI)      
%c                   /           \      /     \ 11  /11
%c                  /             \    / 3     (MZI)
%c                 /              (MZI)             \27
%c                /                    \ 7          /7
%c               /                      \      (MZI)
%c              /                        \    / 7   \23
%c             /                          (MZI)
%c            /                               \ 15  /15
%c           /                                 (MZI)
%c          / 1                                     \31
%cEin---(MZI)
%c          \ 2                                     /2
%c           \                                  (MZI)
%c            \                               / 2   \18
%c             \                         (MZI)      
%c              \                       /     \ 10  /10
%c               \                     / 2     (MZI)
%c                \               (MZI)             \26
%c                 \              /    \            /6
%c                  \            /      \ 6    (MZI)
%c                   \          /        \    / 6   \22
%c                    \        /         (MZI)
%c                     \      /               \ 14  /14
%c                      \    / 2               (MZI)
%c                      (MZI)                       \30
%c                           \ 4                    /4
%c                            \                (MZI)
%c                             \              / 4   \20
%c                              \        (MZI)      
%c                               \      /     \ 12  /12
%c                                \    / 4     (MZI)
%c                                (MZI)             \28
%c                                     \ 8          /8
%c                                      \      (MZI)
%c                                       \   / 8    \24
%c                                       (MZI)
%c                                           \ 16   /16
%c                                             (MZI)
%c                                                  \32
%c Each MZI vertically aligned correspond to one Stage of the OFFT.
%c
%c
%c
%c
%c
%c                                           Created by P.Marciano LG
%c                                           11/01/2018
%c                                           Last UpDate
%c                                           12/01/2018
%c                                           pablorafael.mcx@gmail.com
%c
%c Refences:
%c@article{hillerkuss2010simple,
%c  title={Simple all-optical FFT scheme enabling Tbit/s real-time signal processing},
%c  author={Hillerkuss, D and Winter, M and Teschke, M and Marculescu, A and Li, J and Sigurdsson, G and Worms, K and Ezra, S Ben and Narkiss, N and Freude, W and others},
%c  journal={Optics express},
%c  volume={18},
%c  number={9},
%c  pages={9324--9340},
%c  year={2010},
%c  publisher={Optical Society of America}
%c}
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
%cEout,Count,MapVet]=OpticalFFTN(t,SymPer,StopPoint,E1,AngVec,...
                                                  % ActStag,Count,E0,MapVet)
%c    OpticalFFT
%c
%c  Input:
%c  t         : Time vector of the whole simulation                     [s]
%c  SymPer    : Symble Period                                           [s]
%c  StopPoint : Maximum number of stages need                           [u]
%c  E1        : Actual input field for the FFT                          [-]
%c  AngVec    : Brings the phase shifft for the current stage         [rad]
%c  ActStag   : Current interaction stage                               [u]
%c  Count     : variable to control the expansion of the phase delay    [u]
%c  E0        : Signal is equal to zero if not otherwise given          [-]
%c  MapVet    : Map vector to locate the actual carriers at the output  [u]
%c  
%c  Output:
%c  Eout      : Final and partial output signal                         [-]
%c  Count     : variable to control the expansion of the phase delay    [u]
%c  MapVet    : Map vector to locate the actual carriers at the output  [u]
%c  
%c
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                             S T R U C T S                              %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c
%%
% At the very first part of this iscript, it is necessary to be sure that
% all variable that will be used are needed to be correctly initialized.
% In addition, this function is called recursivaly. Thus, it is important
% to think in all possible steps .

    if nargin<4 || nargin==5 || nargin==6         
        %If the input argument is underspecified an error will be displayed
        error(['Input Arguments are not enougth. Please check this '
                                                        'function Call!']);

    elseif nargin == 4%The minimum input arguments required.
        AngVec  = 0;                                                       %The first angle is zero
        ActStag = 1;                                                       %For the FFT at it's first order
        Count   = ones(1,StopPoint);                                       %Control variable initialized
        E0      = zeros(size(E1,1),size(E1,2),size(E1,3));                                     %Variable for the second input of the MZI
        MapVet  = 1;                                                       %At the stage 1 the first carrier is 1
    elseif nargin<8
        E0     = zeros(size(E1,1),size(E1,2),size(E1,3));                                      %Variable for the second input of the MZI
        MapVet = 1;                                                        %At the stage 1 the first carrier is 1
    elseif nargin<9
        MapVet = 1;                                                        %At the stage 1 the first carrier is 1
    end

    %%          Delay Interferometer
    % There is where the magic of this sistme happens. Basically saying, 
    %by given the time signal a time delay and a phase delay. Those two 
    %parameters will not affect the data. But when the input signal is 
    %divided half of it passes through a light path and it is not changed. 
    %The second half will pass through another light path and through the 
    %application of a voltage it is possible to delay the signal by 
    %changing material characteristics. As a result, the second signal will 
    %have time and phase differences from the other signal. When those two 
    %components combine again at the end to DI there will be constructive 
    %and destructive results in both signals. At the end, the DI outputs 
    %two signals where one carrier the even components of the sub-carries 
    %and the other carriers the odd components of the input signal. Thus, 
    %it was implemented the first step of the Optical FFT.

    TimDel        = SymPer/(2^ActStag);                                    %Create the time delay
    PhaDel        = AngVec(Count(ActStag));                                %Acquiring the time delay         
    [Eout1,Eout2] = DelayInterfExp(f,TimDel,PhaDel,E1,E0);                 %Sends the input signal to the DI

    Count(ActStag) = Count(ActStag) + 1;                                   %Takin in account the current interaction
    ActStag        = ActStag + 1;                                          %Marking the occurence of the current stage

    %%   Checking for recursivity
 % At this point it is needed to evaluate in which interaction the program 
 %is and if it gets to an end. As stop criteria, the variable ActStag, 
 %indicates in wich stage the script is hence when the value of Actstag is 
 %equal to the StopPoint, which means that the number of stages for the 
 %program It has got to an end.
    if ActStag>StopPoint
        MapVet = [MapVet MapVet+2^(ActStag-2)];
        Eout   = cat(3,Eout1,Eout2);
    else%Otherwise move to the next stage
        n      = ActStag-1;                                                %Auxiliar variable to calculate the step
        MyStep = pi/(2^n);                                                 %Step for the phase delay
        AngVec = [AngVec AngVec+MyStep];                                   %Calculating the phase delay for this current interaction
        Vetaux = MapVet+2^(ActStag-2);                                     %Calculating the carrier of thie current interaction
        
        [EoutAux1,Count,MyVetAux1] = OpticalFFTNT(f,SymPer,StopPoint,...
                                     Eout2,AngVec,ActStag,Count,E0,Vetaux);%Calculating the next interaction of the FFT
        [EoutAux2,Count,MyVetAux2] = OpticalFFTNT(f,SymPer,StopPoint,...
                                     Eout1,AngVec,ActStag,Count,E0,MapVet);%Calculating the next interaction of the FFT
        
        MapVet = [MyVetAux2 MyVetAux1];                                    %Storing the received mapping vector.
        Eout   = cat(3,EoutAux2,EoutAux1);                                      %Storing the received Output from the FFT.
    end
end

% while Stage<=StopPoint
% 
% % else
%     n = Stage-1;
%     m = length(AngVec);
%     MyStep = pi/(2^n);
%     MyVec = zeros(1,2^n);
%     MyVec(1:m) = AngVec;
%     for kk=1:m
%         MyVec(m+kk) = AngVec(kk) + MyStep;
%         MyVec.*(180/pi)
%         a=2;
%     end
%     AngVec = MyVec;
%     Stage = Stage + 1;
% %     Eout = expfunc(Stage,AngVec,StopPoint);
% end