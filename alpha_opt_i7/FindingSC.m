%c
%c                                                       ..'Ž`'..'Ž`..'Ž`..                                                   
%c       File: finding_peaks (Detection of only the interest point for 
%c evaluation)
%c
%c       This code is used to create the different types of possible
%c combinations with the variables of interest and apply then to a MZM with
%c one single arm and evaluate the number of peaks created and their
%c flatness.
%c 
%c
%c
%c
%c
%c                                           by P.Marciano LG
%c                                           08/08/2017
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
InputDataFindingPeaksRectFilter;                                           %Loads the inicial parameters for this algorithim
if Ploting
    figure;hold all;
end
for kk=1:size(To_Eval,2)                                                   %Main loop for each experiment realization
    locsaux = To_Eval{kk}.Loc;                                             %selection only the position that are diferent from zero
    peaksaux = To_Eval{kk}.Pea;                                            %selecting only the peaks from the vector
    if Ploting
        plot(ML*log10(peaksaux));
        hold all;
    end
    peaksaux = peaksaux.*Filtro_Sqr(locsaux);                              %Limit the length of the optical COMB to a maximum of carriers  chosen. if the 
                                                                           %variable Obs_Port = 21 that means the algorithim will look just the 20 
                                                                           %subcarriers around the center carier
    if Ploting
        plot(ML*log10(peaksaux));
        hold off;
    end
    peaksaux = peaksaux(peaksaux ~= 0);                                    %Eliminate the zeros to avoid error.
    sidediff(kk) = abs(ML*log10(peaksaux(1))-ML*log10(peaksaux(end)));
    centdiff(kk) = abs(ML*log10(min(peaksaux(boolean([1 0 1]))))-ML*...
                                                       log10(peaksaux(2)));
%     num_peak(kk) = length(To_Eval{kk}.Pea);                                       %making the acounting of the number of peaks
%     if ~isempty(peaksaux)
%         planicidade(kk) = max(ML*log10(peaksaux))- min(ML*log10(peaksaux));%measuring the ripple (db)
%     else
%         planicidade(kk) = inf;
%     end
    100*kk/size(To_Eval,2)
end
pontos = 1:length(centdiff);
% 
% 
%% Ploting a graph to observe the patern acquired
if Ploting
    figure;
    plot(sidediff);%ploting the numper of peaks found for each realization
    hold all;
    plot(centdiff);%ploting the ripple for each realization selected
end
%% Finding the best Configuration
BestSideDiff = sidediff(sidediff<=MaxSideDiff);
LocBSD = sidediff<=MaxSideDiff;
% mais10pic = num_peak(num_peak>=Num_Subcari);                                %Finding the sub group with the maxiimum number of carriers
% locmais10pic = num_peak>=Num_Subcari;                                       %Finding where those sub groups are

                                                                           
if ~(isempty(BestSideDiff))                                                   %Check is there are carriers to be evaluate, if not the evaluation is done
    if Ploting
        figure;
        hold all;
        loc_aux = pontos(LocBSD);
        for ll=1:length(BestSideDiff)
            plot(f(To_Eval{ll}.Loc),ML*log10(To_Eval{ll}.Pea));
            a=1;
        end
    end
    min_CenDif         = min(centdiff);
    aux_CenDif         = linspace(min_CenDif,min_CenDif,size(centdiff,2));
    aux_CenDif(LocBSD) = centdiff(LocBSD);
    LocBCD             = aux_CenDif>=MinCentDiff;
    BestCentDiff       = centdiff(LocBCD);
    BestSideDiff       = sidediff(LocBCD);
    melhorespontos     = pontos(LocBCD);
%     maxi_plan = max(planicidade);                                          %Finding the worst option for the OFCS
%     inve_plan = linspace(-maxi_plan,-maxi_plan,length(planicidade));       %create a inverse vector with the planicity vector
%     inve_plan(locmais10pic) = -1.*(planicidade(locmais10pic)-Best_Dif_Amp);%Given that this algorithm is looking for Carriers 
%                                                                            %which the Power Peak have a difference of Best_Dif_Amp, 
%                                                                            %which means the planicity of the OFC in evaluation, 
%                                                                            %the chosen planicity will be placed above 0, and for a 
%                                                                            %margin of erro it will pick the valus abouve -1 as best 
%                                                                            %canditades for the UOFCS. Those canditades will be 
%                                                                            %selected from the sub group previouly founded
%     locpani10picmais = inve_plan>-1;                                       %Finding the location of the best canditades   %(max(inve_plan)-Best_Dif_Amp);
%     plani10picmais(locpani10picmais) = inve_plan(inve_plan>-1);            %From the grup of 10 or more peaks the smallest ripple      %(max(inve_plan)-Best_Dif_Amp));
%                                                                           
%     BestPlanicity = planicidade(locpani10picmais);
%     melhorespontos = pontos(locpani10picmais);                             %Return back the reference. This line is to indicate 
                                                                           %where the best canditades were located on the 
                                                                           %original family selection.
                                                                           
    %% Ploting the result founded for qualitative evaluation
    representa = [];
    if Ploting
        figure;
        hold all;
        for jj=1:length(melhorespontos)
            plot(f(To_Eval{melhorespontos(jj)}.Loc),ML*log10(To_Eval{...
            melhorespontos(jj)}.Pea.*Filtro_Sqr(To_Eval{melhorespontos(...
                                                          jj)}.Loc)),'-o');%Observing whether the peaks actualy present to be planar
            aux_repre=ML*log10(To_Eval{melhorespontos(jj)}.Pea.*...
                              Filtro_Sqr(To_Eval{melhorespontos(jj)}.Loc));
            aux_repre = aux_repre(aux_repre~=-inf);
            representa(jj,1:length(aux_repre)) = aux_repre;
        end
    end
    save(LocalSaved,'melhorespontos','representa','BestCentDiff',...
                                                           'BestSideDiff');        %Saving the results
end
a=1;