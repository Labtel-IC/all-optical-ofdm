%c
%c                                                       ..'Ž`'..'Ž`..'Ž`..                                                   
%c       File: FindingPeaks (Detection of only the interest point for 
%c evaluation)
%c
%c       On this script is evaluated the generated optical comb as result
%c of the miriady of possible combinations of it variables and driven
%c signals. This classification is based on two parameter mainly, first is
%c planicity of the OFC, which means, the diffence in the amplitude between
%c carries. As lower is the planicity higher is the flatness of the OFC.
%c The second parameter is the number of carriers that composed the OFC.
%c The higher numbers of subcarries in the OFC the better is the result.
%c 
%c
%c
%c
%c
%c                                           by P.Marciano LG
%c                                           13/10/2017
%c                                           pablorafael.mcx@gmail.com
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
%c
%c   
%c                     
%c   peaks             :This is one variable loaded from the earlier 
%c                     simulation. In it is sotored the number of peaks 
%c                     founded that matched the minimum requirem of the 
%c                     previous simulation.                            [db]
%c
%c   locs              :This variable stores the positions in the frequency
%c                      vector for which corespond to the peaks founded.
%c
%c   num_peak          :The peaks vector have zeros within, thus this 
%c                     varible is used to count the number of peaks for
%c                     each realization on the peaks vector.
%c
%c   planicidade       :For each realization on the peaks vector this
%c                      variable stores the diference between the maximum
%c                      and the minimum that indicates how planar is the
%c                      optical COMB.
%c
%c   mais10pic         :This variable stores just the peaks with 10 or more
%c                      components 
%c                      
%c   locmais10pic      :This variable stores the position of the
%c                      planicidade for which correspond to the subset of
%c                      peaks founded.
%c
%c   plani10picmais    :This variable stored the planicity that is above a 
%c                     threshold creating a new subset from the grup of the
%c                     10 or more peaks.
%c   locpani10picmais  :This variable stores the position that corespont to
%c                     each grup of peaks with the hiest planicity from the
%c                     grup that have 10 or more peaks.
%c
%c   melhorespontos    :Here it is sotere the position of the locs vector
%c                      that is referent for the best set of the COMB
%c                      founded with hiest number of carries and the
%c                      smallest ripple.

%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                             S T R U C T S                              %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c   Here it will be added the functions that this code will call
%c   
%c
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Start of the Program                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Measure the number of peaks and the ripple
clear,close all,clc;
InputDataFindingPeaks;                                                     %Loads the inicial parameters for this algorithim
% For qualitative evaluation of the subcarriers initical classification
if Ploting 
    figure;hold all;
end
%Variable to store the main parameters for further evaluation
NumbGroup = zeros(1,size(To_Eval,2));                                      %Number of Groups composed by subcarriers
Planicity = zeros(1,size(To_Eval,2));                                      %Planicity of the respective Group
%% Main loop for each experiment realization
%   At this loop the scrip will be responsible for a first classification
% and calculate the correct parameters of planicity and number of groups
% for each combination of variables previouly stated for the generation of
% an OFC. Therefore, the script will be able to show which results brought
% to us the best OFC related to its flatness and number of subcarriers.
for kk=1:size(To_Eval,2)
    % For qualitative evaluation of the subcarriers initical classification
    if Ploting
        plot(f(To_Eval{kk}.Loc),ML*log10(To_Eval{kk}.Pea),'-d');           %It will plot the current OFC under evaluation
    end
    % Here the current OFC will be evaluated the number of groups of
    % subcarriers and the avarege planicity among them.
    for jj=1:length(To_Eval{kk}.Pea)-SubCarrPerGroup                       %select the group for evaluation
        DifAmp = abs(ML*log10(max(To_Eval{kk}.Pea(jj:jj+SubCarrPerGroup)...%Meassuring the Amplitude [dB] difference withing the group
               )) - ML*log10(min(To_Eval{kk}.Pea(jj:jj+SubCarrPerGroup))));
        if DifAmp<=Dif_Amp                                                 %If the Amplitude of the select group is relevant for the 
                                                                           %current evaluation this group will be taken in acount for 
                                                                           %further evaluation
            NumbGroup(kk)=NumbGroup(kk)+1;                                 %Taking the acount of how many groups there are withing this OFC
            Planicity(kk) = Planicity(kk) + DifAmp;                        %Measuring the avarege planicity of those selected groups
        end
    end
    Planicity(kk) = Planicity(kk)/NumbGroup(kk);
    clc;100*kk/size(To_Eval,2)                                                 %Counter to display the evolution of this script
end
pontos = 1:length(Planicity);                                              %Variable needed to recover back the index of the best 
                                                                           %selected OFC from the previusly code.
%% Ploting a graph to observe the patern acquired
% For qualitative evaluation of the subcarriers for further classification
if Ploting
    figure;
    plot(NumbGroup);%ploting the numper of peaks found for each realization
    hold all;
    plot(Planicity);%ploting the ripple for each realization selected
end
%% Finding the best Configuration
BestGroupNumb = NumbGroup(NumbGroup>=MinObsGroup);                         %Finding the sub group with the maximum number of carriers
LocBGN = NumbGroup>=MinObsGroup;                                           %Finding where those sub groups are

if ~(isempty(BestGroupNumb))                                               %Check is there are carriers to be evaluate, if not the evaluation is done
 % For qualitative evaluation of the subcarriers for further classification
    if Ploting
        figure;
        hold all;
        loc_aux = pontos(LocBGN);
        for ll=1:length(BestGroupNumb)
            plot(f(To_Eval{loc_aux(ll)}.Loc),ML*log10(To_Eval{loc_aux(ll...% Plot the selected OFC with the highest number of subcarriers.
                                                                  )}.Pea));
            a=1;
        end
    end
%% Further evaluation and classification of the selected OFCs
% Two main points must be taken in acount here. Frist, it is important to
% recover the index number of the OFC under evaluation. Which means, for
% every single configuration (combination of input variables for the
% generation of the OFC) there are one index number from which all
% characterists used can be recovered without the need of re-runn the whole
% simulation again. Therefore the auxiliar variable need to have the same
% length as the combination of variables used to generate those OFCs. The
% second parameter is the planicity. In that case the auxiliar variable
% needed to have an significant constrast with the selected parameters. As
% we are looking for the OFC with the smallest planicity hence the auxiliar
% variable needs to be as greater as possible from the actual planicity.
    aux_plan = linspace(max(Planicity),max(Planicity),size(Planicity,2));  %Creating an auxiliar variable with high constrast from the 
                                                                           %paramenter under evaluation and with the proper size to
                                                                           %keep track on the index number
    aux_plan(LocBGN) = Planicity(LocBGN);                                  %load in the auxiliar variable the planicity from the 
                                                                           %previouly select OFCs with the best number of subcarriers
    LocBestPlanicity = aux_plan<Best_Dif_Amp;                              %Finding the location of the best canditades with the 
                                                                           %lowest planicity   %(max(inve_plan)-Best_Dif_Amp);
    BestPlanicity    = Planicity(LocBestPlanicity);                        %From the grup of higher number of peaks with the smallest 
                                                                           %ripple this variable stores the ripple of the best OFCs 
                                                                           %selected      %(max(inve_plan)-Best_Dif_Amp));
                                                                          

    melhorespontos = pontos(LocBestPlanicity);                             %Return back the reference. This line is to indicate 
                                                                           %where the best canditades were located on the 
                                                                           %original family selection.
                                                                           
    %% Ploting the result founded for qualitative evaluation
    % For qualitative evaluation of the subcarriers afther classification
    % was finished.
    if Ploting
        figure;
        hold all;
        for jj=1:length(melhorespontos)
            plot(f(To_Eval{melhorespontos(jj)}.Loc),ML*log10(To_Eval{...
                                            melhorespontos(jj)}.Pea),'-o');%Observing whether the peaks actualy present to be planar
        end
    end
    save(LocalSaved,'melhorespontos','BestPlanicity');                     %Saving the results
end 
a=1;