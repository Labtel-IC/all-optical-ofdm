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
load('TwoArmLoopTestD2609T1011.mat');%Load the result of the past simulation for evaluation                                                         
for kk=1:size(peaks,1)%Main loop for each experiment realization
    locsaux = locs(kk,locs(kk,:)~=0);%selection only the position that are 
                                     %diferent from zero
    peaksaux = peaks(kk,peaks(kk,:)~=0);%selecting only the peaks from the 
                                        %vector
    num_peak(kk) = length(peaksaux);%making the acounting of the number of 
                                    %peaks
    planicidade(kk) = max(peaksaux)- min(peaksaux);%measuring the ripple
end
%% Ploting a graph to observe the patern acquired
plot(num_peak);%ploting the numper of peaks found for each realization
hold all;
plot(planicidade);%ploting the ripple for each realization selected
%% Finding the best Configuration
[mais10pic,locmais10pic] = findpeaks(num_peak,'MinPeakHeight',9);          %Finding the location of the grup
                                                                           %of the COMB that has 10 or more
                                                                           %subcaries
[plani10picmais,locpani10picmais] = findpeaks((-1.*planicidade(...         %From the grup of 10 or more peaks
                                     locmais10pic)+5),'MinPeakHeight',3);  %it is finding those that have 
                                                                           %the smallest ripple

melhorespontos = find(ismember(planicidade,planicidade(...                 %Locating in the principla vector
                                         locmais10pic(locpani10picmais))));%the position of peaks that 
                                                                           %presented the selected characteristics
%% Ploting the result founded for qualitative evaluation
figure
for jj=1:length(melhorespontos)
    plot(f(locs(melhorespontos(jj),locs(melhorespontos(jj),:)~=0)),...
           peaks(melhorespontos(jj),peaks(melhorespontos(jj),:)~=0),'-o');  %Observing whether the peaks actualy 
                                                                           %present to be planar
    hold all;
end
save('achomelhor');%Saving the results