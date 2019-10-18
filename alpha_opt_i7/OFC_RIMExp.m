%c
%c                                                       ..'Ž`'..'Ž`..'Ž`..                                                   
%c       File: OFC_RIM (Optical Frequency Comb Rim)
%c
%c     This frist code aim to implement the loop (or single interaction) to
%c create the multiple tones present on the Optical Comb.
%c
%c
%c
%c
%c                                           by P.Marciano LG
%c                                           18/08/2018
%c                                           pablorafael.mcx@gmail.com
%c 
%c     References:
%c       [1] Lei, Cheng, et al. "Recirculating Frequency Shifting Based 
%c           Wideband Optical Frequency Comb Generation by Phase Coherence 
%c           Control." IEEE Photonics Journal 7.1 (2015): 1-7.
%c       [2] Zhang, Junwen, et al. "Stable optical frequency-locked 
%c           multicarriers generation by double recirculating frequency
%c           shifter loops for Tb/s communication." Journal of Lightwave 
%c           Technology 30.24 (2012): 3938-3945.
%c
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                    G L O B A L  V A R I A B L E S                      %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c   
%c
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                   L I S T  O F  V A R I A B L E S                      %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c   CW              : Continuos Wave, it is the signal that represent the
%c                     lazer.
%c   Filter          : It is the filter that will be insert in the loop
%c   Eout            : It is the final signal that will get our of the loop
%c   
%c
%c
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                             S T R U C T S                              %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c   Here it will be added the functions that this code will call
%c   
%c   RFS_DSB.
%c                    .Does not have input parameters
%c   FiltroGaussiano.
%c                   Input:
%c                         .f   : Frequency
%c                         .f0  : Filter Band Width
%c                         .fc  : Central frequency
%c                         .m   : Order of super-gaussiana
%C                   Output:
%c                         .H   : Filter Transferfunctiondo filtro
%c                         .HdB : Filter Transferfunctiondo filtro [db]
%c
%c   Mach_Zehnder_Modulator_simplificado. 
%c                   Input:
%c 	                       t   : Time Vector
%c                         Ein : Electrical Field  
%c                         U1t : Eletrical input for the 1th arm of MZM  
%c                         U2t : Eletrical input for the 2nd arm of MZM
%c                         MZ_Input_File : File for configuration
%c                   Output:
%c                         Eout: Modulated Eletrical field output
%c                         H   : Transfer Function of the modulator
%c
%c
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Start of the Program                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                     Initial Signal



E1=zeros(1,length(CW));
Eout=zeros(1,length(E1));
%%                    Main Loop
for N=1:N_volta
   
    
   %Choose to use a infinit Loop or to have a output
   if InfinitLoop
       E1(S,:) = Eout(S,:) + CW;
   else      
       E1(S,:) = CuplerGain.*(Eout(S,:) + CW);
   end
   
    % Mach Zhender Modulator
  [Eout(S,:),~]=Mach_Zehnder_Modulator_simplificado(t,E1(S,:),EleSig,MzmInputInd);

   %Use the amplifyer to give a bust on the signal
   if Amplify
%        Eout(S,:) = Gain.*Eout(S,:);
%                   EDFA_SP(Eout3,pe.tb,6.0,po.NF1,po.PWsat1,pe.fsmp_b./2);
        Eout(S,:) = EDFA_SP(Eout(S,:),t,Gain,NF,PWsat,fsample/2,c,lambdac);
   end
   
   %Use the filter on the Loop signal to limit the Nmb of Carriers
   if Filtering
       Eout(S,:) = ifft(fft(Eout(S,:)).*Filtro);    
   end
   
   % add the noise to the Loop signal before re-enter into the MZM
%    if AddNoise
%       EsN0 = EbN0 + 3 + 10*log10(k) - 10*log10(nsamp);%Calculate the 
%                                                      %signal to noise ratio
%       Eout(S,:) = Eout(S,:) + wgn(1,length(t),-1*EsN0);
%    end
%    if N>80
   if Ploting
     
%        figure(1);
%        subplot(3,1,1),plot(f,db(abs((fftshift(fft(E1))))/length(E1)));
%        title('Input Signal, Ein','FontSize',14);
%        ylabel('Magnitude [db]','FontSize',12);
%        xlabel('Frequency [Hz]','FontSize',12);
%        grid on;
%        if AxisOn
%            axis(AxisFig1);
%        end
       figure(1);
%        subplot(3,1,2)
       plot(f,db(abs((fftshift(fft(Eout))))/length(Eout)));
       title('Output Signal, Eout','FontSize',14);
       ylabel('Magnitude [db]','FontSize',12);
       xlabel('Frequency [Hz]','FontSize',12);
       grid on;
       axis([0.9e10 4.55e11 -26 -17])
%        if AxisOn
%            axis(AxisFig2);
%        end
%        figure(1);
%        subplot(3,1,3),plot(f,fftshift(Filtro));
%        title('Optical Ring Filter','FontSize',14);
%        ylabel('Magnitude [db]','FontSize',12);
%        xlabel('Frequency [Hz]','FontSize',12);
%        grid on; 
%        if AxisOn
%            axis(AxisFig3);
%        end
%        N
       
       drawnow;
%        end
   end
end
