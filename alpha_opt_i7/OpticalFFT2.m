function [ Eout,Count ] = OpticalFFT(t,SymPer,MaxNumStag,Count,E1,E0,...
                                                            ActStag,PhaDel)
%%                        All Optical Fast Furier Transform
%c function [ Eout,Count ] = OpticalFFT(t,SymPer,MaxNumStag,Count,E1,E0,...
%c                                                          ActStag,PhaDel)
%c This function is responsible for implementing an all optical FFT. Till
%c this present moment. It is possíble to create an FFT 16 of size. Wich
%c means, it can filter just 16 subcarriers. For furter implementation it is
%c need to correctly indentify the expansion function for the phase delay
%c parameters. The proble is, the carriers are shufled by this FFt, and for
%c a perfect FFT implementation the exacly phase delay must be added to the
%c right signal at the right time. The patern of randomization of the
%c carriers was not fully understood yet.
%c
%c
%c
%c
%c
%c                                           Created by P.Marciano LG
%c                                           02/11/2017
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
%c
%c    OpticalFFT
%c
%c  Input:
%c  t          : Time vector of the whole simulation [s]
%c  SymPer     : Symble Period [s]
%c  MaxNumStag : Maximum number of stages need
%c  Count      : variable to control the expansion of the phase delay
%c  E1         : Actual input field for the FFT [-]
%c  E0         : Signal is equal to zero if not otherwise given
%c  ActStag    : Variable that marks where the FFT currently is
%c  PhaDel     : Phase Delay that will be used on the current interaction
%c  
%c  Output:
%c  Eout       : Final and partial output signal
%c  Count      : variable to control the expansion of the phase delay
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


%%   Initialization
% Firstly it will check if the Phase Delay was given as parameter. Because
% it is the mark for the first interaction of the FFT. As the example given
% outside this script. It marks the first stage of the interferometer.
%
    if nargin<5                          %If input criterias were tome 
                                         %fullfiled the code will display 
                                         %an erros
        error(['Input Arguments are not enougth. Please check this '...
                                                        'function Call!']);
    elseif nargin<6                      %Checking if the second field 
                                         %was given as an parameter
        ActStag = 1;                     
        PhaDel = pi - pi/(2^(ActStag-1));
        E0 = 0;                          
    elseif nargin<7                      %Also checking for the start point
        ActStag = 1;
        PhaDel = pi - pi/(2^(ActStag-1));
    elseif nargin<8                      %Checking if it is at the first 
        PhaDel = pi - pi/(2^(ActStag-1));%interaction by finding if the 
                                         %Delay fase was given as parameter
    end
    
%%          Delay Interferometer
% There is were the magic of this sistme hapens. Basicaly saying, by given
% the tiem signal a time delay and a phase delay. Those two parameters will
% not afect the data. But when the input signal is devided half of it
% passes throug an light path and it is not changed. The second half will
% pass through anoter light path and applying an voltage it is possible to
% delay the signal by changing material characteristics. As as result the
% second signal will have time and phase differences from the other signal.
% When those to componentes combine again at the end to DI there will be
% construtive and destrutive results in both signals. At the end, the DI
% outputs two signals where one carriers the even componentes of the
% subcarries and the other carriers the odd components of the input signal.
% Thus, it was implemented the frist step of the Optical FFT

    TimDel = SymPer/(2^ActStag);                                           %Create the time delay
    [Eout1,Eout2] = DelayInterf(t,TimDel,PhaDel,E1,E0);                    %Sends the input signal to the DI

    %%   Checking for recursivity
 % At this point it is needed to evaluate in which interaction the program
 % is and if it get to an end.
 % As an stop criteria, the variable ActStag, indicates in wich stage the
 % script is. When the value of Actstag is equal to the maximum number of
 % stages for the program It has get to an end.
    if ActStag < MaxNumStag
 % If nothing has changed, the next phasa delay is calculated
 (180/pi)*PhaDel
 Count(ActStag)
        if ~mod(Count(1),2) %Creat an phase delay for the even carriers
            PhaDel = pi - 2*Count(ActStag)*(pi/(2^ActStag));
        else%Creat an phase delay for the odd carries
            PhaDel = pi - (2*Count(ActStag)-1)*(pi/(2^ActStag));
        end
 (180/pi)*PhaDel
%%   Correction of the new interaction
% if there is more stages to be inplemented the ActStage is up dated and
% this function will be called again.
        ActStag = ActStag + 1;
        [EoutAux1,Count]= OpticalFFT(t,SymPer,MaxNumStag,Count,Eout1,0,ActStag,PhaDel);
% This will search for the last combination needed for the FFT and will
% resolve it frist. And then return for the function that called it.
% Thus it is need to calculate a new phase delay. Remember that the DI
% returned to us two signal thus it is neded to call this funtion two time
% in a row, one for each signal, with the sligtly difference for the phase
% relay.

%%      Adjust for the new phase delay
% Acording with the interaction the signal under analizes will be even or
% odd. This following lines verify which case the actual signal belongs to.
        if (ActStag - 1) == 1
            Count = ones(MaxNumStag,1);
            Count(1) = 2;
        else
            Count(ActStag) = Count(ActStag)+1;
        end
% Therefore, the new phase delay can be finded and this function can be
% called again.
 (180/pi)*PhaDel
 Count(ActStag)
        PhaDel = PhaDel - pi/2;
 (180/pi)*PhaDel
        [EoutAux2,Count] = OpticalFFT(t,SymPer,MaxNumStag,Count,Eout2,0,ActStag,PhaDel);
%%     Controling interaction and finalization of the algorithm
% As this is an recursive function the Count variable is used to control
% how many time one stage was called. It is inportant to know which will be
% the phase delay of the next stage. The problem here with the expanssion
% of this algorithm is that the singnals are shufled each time and if a new
% implementation need to be done, it is important to evaluate manualy how
% the new phase delays will be placed and their values. Then, the logic of
% this script needs to be change to mach new possible results.
        Count(ActStag) = Count(ActStag)+1;
% When all interaction were finished this function stores the results
% fields and return all of them recursively.
        Eout = [EoutAux1;EoutAux2];
    else
        Eout = [Eout1;Eout2];
    end
end

