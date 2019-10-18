function [olho_mat,tn,Olho_Aber,Aber_Lhi,Aber_Llow] = Olho(Dados,t_simb,...
                                                      n_amos_p_sim,ploting,tmax)
%%           Creating the Eye Diagram of one electrical signals.  
%c function [olho_mat,tn,Olho_Aber,Aber_Lhi,Aber_Llow] = Olho(Dados,...
%c                                             t_simb,n_amos_p_sim,ploting)
%c
%c
%c This function is resposible for generating the eye diagram of an given
%c signal.
%c
%c
%c                                           Created by P.Marciano LG
%c                                           22/11/2016
%c                                           07/01/2018
%c                                           pablorafael.mcx@gmail.com
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
%c    Olho
%c
%c  Input:
%c  Dados        : The complete tada to be evaluated                    [b]
%c  t_simb       : The time period of one symbol                        [s]
%c  n_amos_p_sim : The number of samples per symbol                  [samp]
%c  ploting      : Control variable to plot or not the eye diagram withing
%c                 this fucntion                                   [bolean]
%c  
%c  Output:
%c  olho_mat     : Matriz with the eye diagram stored                   [u]
%c  tn           : Time vector for ploting the Eye Diagram outside this
%c                 function                                             [s] 
%c  Olho_Aber    : Statical information about the Eye Opening           [u]
%c  Aber_Lhi     : Statical information about the Eye highest point     [u]
%c  Aber_Llow    : Statical information about the Eye lowest point      [u]
%c  
%c
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                             S T R U C T S                              %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c
%%                       Initializing
%% Recovering the main iformation about the income signal
if nargin<4
    error('Not enough input arguments.');
elseif nargin<5
    tmax = 1;
end
    
Dmax  = max(Dados);                                                        %Variable used to limit the ploting at the maximum point
Dmin  = min(Dados);                                                        %Variable used to limit the ploting at the minimum point
tam   = length(Dados);                                                     %Recovering the size of the income information
N_sim = tam/n_amos_p_sim;                                                  %Recovering the number of Symbols withing the information
tf    = N_sim*t_simb;                                                      %Estimating the final time of this signal
t     = linspace(0,tf,tam);                                                %Recovering the time vector of the income signal
ta    = t(2)-t(1);                                                         %Recovering the sampling ratio of this signal
%% Seting variables for controling the spliting to the income signal
% tmax  = 3;                                                                 %Seting the number eyes to be plot
if N_sim <= tmax                                                           %If the income signal just have one symbel just one eye will be plot
    n_plot = 1;
else
    %Otherwise, we need to calculate the number of samples needed to be
    %taken to create the eye diagram accordingly with the previous
    %especified.
    n_plot = tam/(tmax*n_amos_p_sim);                                      %taking the number of samples needed
    %Because of reasons, the number of samples may not be multiple of the
    %number of eyes requested for ploting hence, it must be downsizes to
    %the minimal of one eye to be ploted.
    while mod(tam,(tmax*n_amos_p_sim))                                     %Checking if the number of samples is divisible by the number of eyes to be plot         
        tmax = tmax - 1;                                                   %If not, reduz the number of eyes by one unity
        n_plot = tam/(tmax*n_amos_p_sim);                                  %Recalculate the number of samples needed to be taken
    end
end
%% Resizing the income signal for ploting
nlin     = tmax*n_amos_p_sim;                                              %Setting the number of lines of the new matriz
ncol     = n_plot;                                                         %Setting the number of columns
tn       = reshape(t,round(nlin),round(ncol));                             %Reshaping the time vector to match the signal matriz
olho_mat = reshape(Dados,round(nlin),round(ncol));                         %Reshaping the signal to an matriz

%% Resizing the income signal for measurement of the eye opening
N1     = n_amos_p_sim;                                                     %Setting the number of lines of the new matriz
N2     = N_sim;                                                            %Setting the number of columns
olho   = reshape(Dados,N1,N2);                                             %Reshaping the signal to an matriz
media  = zeros(1,N1);
hi     = zeros(1,N1);
low    = zeros(1,N1);
Lhi    = zeros(1,N1);
Llow   = zeros(1,N1);
Ab     = zeros(1,N1);
% Measuring the eye opening.
for ii = 1:N1
    %Taking the mean value for the current sample for all symbols
    media(ii) = mean(olho(ii,:));
    %Taking the points abouve the current mean value
    nhi       = find(olho(ii,:) > media(ii));
    %Taking the points under the current mean value
    nlow      = find(olho(ii,:) < media(ii));
    %Taking the mean value of the point abouve the the center of the simbol
    hi(ii)    =  mean(olho(ii,nhi));
    %Taking the mean value of the point under the the center of the simbol
    low(ii)   =  mean(olho(ii,nlow));
    %
    %Taking highest point of the eye opening by subtracting the standart
    %deviation from the mean value of the points abouve the center of the
    %symble for the current sample.
    Lhi(ii)   = mean(hi)  - std(olho(ii,nhi));
    %Taking lowest point of the eye opening by adding the standart
    %deviation from the mean value of the points under the center of the
    %symbol for the current sample.
    Llow(ii)  = mean(low) + std(olho(ii,nlow));
    %Taking the eye opening by the diffeces between the uper level and the
    %lower level.
    Ab(ii)    = Lhi(ii) - Llow(ii);
end
%The final values to be exported outside this function are based on the
%mean value calculated above.
Olho_Aber   = mean(Ab);
Aber_Lhi    = mean(Lhi);
Aber_Llow   = mean(Llow);
% if (nargout == 0),
%   maxolho = max(max(pulse));
%   minolho = min(min(pulse));
%   x       = linspace(minolho,maxolho,50);
%   y       = hist(pulse,x);
%   subplot('position',[0.10 0.09 0.15 0.85]),barh(x,y)
%   grid on,ylabel('Amplitude [u.a]')
%   subplot('position',[0.33 0.09 0.60 0.85])
%   plot(Tolho/1E-12,olho,'r','LineWidth',1)
%   grid on,xlabel('Tempo [ps]'),...
%   ylabel(sprintf('Abertura do diagrama de olho:  %2.5f [u.a]',Abertura))
% end
%% Ploting the Eye Diagram 
if ploting
    figure;
    hold all;
    maxolho = max(max(Dados));
    minolho = min(min(Dados));
    x       = linspace(minolho,maxolho,50);
    y       = hist(Dados,x);
    subplot('position',[0.10 0.09 0.15 0.85]),barh(x,y)
    grid on,ylabel('Amplitude [u.a]')
    subplot('position',[0.33 0.09 0.60 0.85])
%     for kk=1:length(olho_mat(1,:))
%        plot(tn(:,1),(olho_mat(:,kk)),'r');
       plot(tn(:,1),(olho_mat),'r');
    %    plot(tn(round(n_amos_p_sim/2),1),olho_mat(round(n_amos_p_sim/2),kk),'ob')
%     end
    grid on;
    xlabel('Tempo [s]'),ylabel(sprintf('Abertura do diagrama de olho:  %2.5f [u.a]',Olho_Aber))
    set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
end
% ylim([1.1*Dmin 1.1*Dmax]);
% warning on;
end

