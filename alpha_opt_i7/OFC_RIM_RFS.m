%c
%c                                                       ..'Ž`'..'Ž`..'Ž`..                                                   
%c       File: RFS_DSB (Recirculating Frequency Shifting Double Sided Band)
%c
%c     This frist code aim to explore how it is possible to generate a comb
%c of multy carriers. Those carriers are criated from one single Lazerbeam.
%c This beam passes through a loop in which it will be modulated in a  such 
%c that the new frequency will be shifted by fc. Then this shifted  carrier
%c will be added to the main signal again. Then, once more this new  signal
%c pass again through the loop. The number of times that the signal passes
%c through the loop will determinate the number of multiple carriersshifted
%c by fr. In the loop a filter will provide a cut off limit to avoid   that
%c infinit carriers be generated. The modulation is done by a  Mach-Zeander
%c optical device. In the loop there are also a optical amplifier        to
%c compensate loses from the loop. The folowing equation describes who does
%c this process works:
%c
%c
%c
%c      Eout = Ein*cos( (phi_1-phi_2)/2 )*e^[j*( (phi_1+phi_2)/2 )]
%c                              
%c
%c          phi_1 = (pi/V_pi)*V*sin(2*pi*f*t)
%c
%c          phi_1 = (pi/V_pi)*V*sin(2*pi*f*t)
%c
%c
%c
%c                                           by P.Marciano LG
%c                                           02/08/2016
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

%RFS_DSB_INPUT_DATA;%            
% 
% EleSig1           = 0.707*sin(2*pi*fc*t);
% EleSig2           = 0.707*sin(2*pi*fc*t + phase_is);
% if length(DoubleArm) == 2
%     EleSig.U1t=EleSig1;
%     EleSig.U2t=EleSig2;
% else
%     EleSig = EleSig1;
% end
%     
% [Filtro,FiltroDB] = FiltroGaussiano(f,FilterBand,Filter_Center_Freq,m);
% Filtro = fftshift(Filtro);
% E1(S,:)=CW; 
% mCW = max(abs(CW));
E1(S,:)=zeros(1,length(CW));
Eout(S,:)=zeros(1,length(E1(S,:)));

Nvoltas = 1;
%%                    Main Loop
for N=1:N_volta
   
    
   %Choose to use a infinit Loop or to have a output
   if InfinitLoop
       E1(S,:) = Eout(S,:) + CW;
   else      
       E1(S,:) = CuplerGain*(Eout(S,:) + CW);
   end
   
    % Mach Zhender Modulator
   [Eout(S,:),H(S,:)]=Mach_Zehnder_Modulator_simplificado(t,E1(S,:),EleSig,1);
   
   
   
   %Use the amplifyer to give a bust on the signal
   if Amplify
       Eout(S,:) = Gain*Eout(S,:);
   end
   
   %Use the filter on the Loop signal to limit the Nmb of Carriers
   if Filtering
       Eout(S,:) = ifft(fft(Eout(S,:)).*Filtro);    
   end
   
   %
%    if length(DoAmpliOne)== ON
% %        Eout(S,:) = (18.6/max(abs(fft(Eout(S,:))/length(Eout(S,:)))))*Eout(S,:);
%    end
   
   % add the CW to the Loop signal before re-enter into the MZM
%    if length(AddCWToLoop) == ON
% %       Eout(S,:) = Eout(S,:) + CW;
%    end
   
   % add the noise to the Loop signal before re-enter into the MZM
   if AddNoise
      Eout(S,:) = Eout(S,:) + wgn(1,length(t),-40);
   end
   
  %To plot the grphic relevant for this simualtion
   if Ploting
       figure(1);subplot(3,1,1),plot(f,db(abs((fftshift(fft(E1(S,:)))))/length(E1(S,:))));
       title('Input Signal, Ein','FontSize',14);
       ylabel('Magnitude [db]','FontSize',12);
       xlabel('Frequency [Hz]','FontSize',12);
       grid on;%(GridIs);
       if AxisOn
           axis(AxisFig1);
       end
       figure(1);subplot(3,1,2),plot(f,db(abs((fftshift(fft(Eout(S,:)))))/length(Eout(S,:))));
       title('Output Signal, Eout','FontSize',14);
       ylabel('Magnitude [db]','FontSize',12);
       xlabel('Frequency [Hz]','FontSize',12);
       grid on;%(GridIs);
       if AxisOn
           axis(AxisFig2);
       end
       figure(1);subplot(3,1,3),plot(f,fftshift(Filtro));
       title('Optical Ring Filter','FontSize',14);
       ylabel('Magnitude [db]','FontSize',12);
       xlabel('Frequency [Hz]','FontSize',12);
       grid on;
       if AxisOn
           axis(AxisFig3);
       end
%        set (gcf, 'Units', 'normalized', 'Position', [0,0,1,1]);
       drawnow;
   end
end
%E2(S,:) = ifft(fft(Eout(S,:)).*(fft(SingData)*10^(-4.2)));
