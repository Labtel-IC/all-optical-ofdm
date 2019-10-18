
        if AddingNoise
            EsN0 = 10*log10(AsePower);
            Eout = Eout + wgn(1,length(t),EsN0);
            AddingNoise = 0;
        end
        
       
                %% Finding the OSNR per carrier
        EoutF   = 20*log10(abs(fftshift(fft(Eout)./length(...
                                                      Eout))));%Measurements will be done at frequency domain
        [peaks_aux,locs_aux] = findpeaks(EoutF,'MinPeakHeight',max(EoutF)-5);%Finding the carriers;
        PortCon = 1;
        for kk=1:length(locs_aux)                              %Loop for calculating the OSNR for each carrier
            FreqLoc = locs_aux(kk);                            %Selecting the actual carrier
            CarFreq = f(FreqLoc);
            BaseLoc = [find((f>=(CarFreq-fc/2))&(f<=(CarFreq...
            -fc/4))) find((f<=(CarFreq+fc/2))&(f>=(CarFreq+...
                                                      fc/4)))];%Selecting the location of the noise
            [ToCalc,~]=findpeaks(EoutF(BaseLoc));              %Finding the top flor
            BasePow = mean(ToCalc);                            %Calculating the mean value in dB
            StorPoss = round(f(locs_aux(kk))/fc);          %Mapping the actual carrier for one possition on the OSNR vector
            if StorPoss>0                                      %Ignoring the carrier at f=0
                OSNRPC2(StorPoss) = abs(abs(peaks_aux(kk))-...
                                                 abs(BasePow));%Measuring and storing the OSNR
                Port2(PortCon)=StorPoss;
                PeakPower(StorPoss) = peaks_aux(kk);
                PeakPowerL(StorPoss) = locs_aux(kk);
                PortCon = PortCon +1;
                
%                 NoSpDe = 
                
%                 NoFi = ((2*abs(BasePow))/(2.0445*planck*(LightSpeed/lambdac)*(LightSpeed/(lambdac^2))*Bref)) + 1/2.0445;
%                 10*log10(NoFi)
%                 OSNRPC1(StorPoss) = abs(peaks_aux(kk))/(NoFi*planck*(LightSpeed/lambdac)*(LightSpeed/(lambdac^2))*Bref);
%                 OSNRPC1(StorPoss) = abs(peaks_aux(kk)*(lambdac^2))/(abs(BasePow)*LightSpeed*Bref);
%                 10*log10(OSNRPC1(StorPoss))
                
            end
        end
        hold on;
        plot(f(PeakPowerL(1:128)),PeakPower(1:128),'x');
        Planicity = abs(abs(max(PeakPower(1:128))) - abs(min(PeakPower(1:128))))
        fcurv2 = fit(Port2(1:128)',OSNRPC2(1:128)','poly1');
        figure;plot(fcurv2,Port2(1:128),OSNRPC2(1:128),'b*');
        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
        a=9;