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

   %Use the filter on the Loop signal to limit the Nmb of Carriers
   if Filtering
       Eout(S,:) = ifft(fft(Eout(S,:)).*Filtro);    
   end
   
   %Use the amplifyer to give a bust on the signal
   if Amplify
       N
       [~,EoutPdBm] =  MeasPower(Eout(S,:))
       SatGain
       EoutPdBm+Gain
       if EoutPdBm<SatGain
           if (EoutPdBm+Gain)>SatGain
               AuxGain = SatGain-EoutPdBm;
               Eout(S,:) = (10^((AuxGain-30)/10)).*Eout(S,:);
           else
               Eout(S,:) = (10^((Gain-30)/10)).*Eout(S,:);
           end
       end
   end
   
   % add the noise to the Loop signal before re-enter into the MZM
   if AddNoise
      EsN0 = 10*log10(AsePower);% + EbN0 + 3 + 10*log10(k) - 10*log10(nsamp);%Calculate the %signal to noise ratio
      Eout(S,:) = Eout(S,:) + wgn(1,length(t),EsN0);
   end
   
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
%        subplot(3,1,2),
       plot(f,20*log10(abs((fftshift(fft(Eout))))/length(Eout)),f,20*log10(fftshift(Filtro)));
%        axis([0 4.7e11 -5 -3]);
       xmax = max(20*log10(abs(fftshift(fft(Eout)./length(Eout)))));
       set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%        axis([-3e12 4e12 -300 100]);
        if AddNoise
%             axis([0 200*12.5e9 -100 0]);
            axis([-1e11 1.8e12 xmax-100 xmax+0.1]);
        else
            if fc > 12.5e9
                axis([-1e11 0.5e13 xmax-1 xmax+0.1]);
            else
                axis([-1e11 1.8e12 xmax-1 xmax+0.1]);
            end
        end
%         axis([-12.5*5e9 max(f) -40 0])
%        title('Output Signal, Eout','FontSize',14);
%        ylabel('Magnitude [db]','FontSize',12);
%        xlabel('Frequency [Hz]','FontSize',12);
%        grid on;
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
       drawnow;
   end
end
close all;