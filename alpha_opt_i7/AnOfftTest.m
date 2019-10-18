close all;clear;clc;
for mm=1:1
    IfftOrSum = mm;
    UsingAllCarriers = 1;
    SavingAT         = 'Test_0';
    SavingAT2         = ['Test_OOK_car_32_offt_' num2str(mm)];
    NumCarr          = 32;                                                      %Total number of carriers (DownStream and UpStream) - This variable must be even
    NumCarrUsed      = 32;
    lambdac          = 1550e-9;                                                %Central wave length of the signal throughout the fiber
    FiberLength      = 10;
    RefCarr          = 60;
    
    if ~mod(RefCarr,NumCarr)
        ObsCarrPos = 1:NumCarr;
    else
        ObsCarrPos = [mod(RefCarr,NumCarr):NumCarr 1:mod(RefCarr,NumCarr)-1];
    end
    
    CarrPosInt       = 1:NumCarrUsed;
    CarrOddPos       = 1:2:NumCarrUsed;
    CarrEvenPos      = 2:2:NumCarrUsed;
for zz=1:1
    %%
    SendingCarr = 0;
    if SendingCarr
        CarrToSend       = zeros(1,NumCarrUsed);
        for kk=0:length(CarrOddPos)-1
            CarrToSend(kk+1,CarrOddPos(1:1+kk)) = 1;
        end
        [~,ObsCarrLocal] = max(sum(CarrToSend));
        ObsCarrUsed      = zeros(1,NumCarrUsed);
        ObsCarr          = ObsCarrPos(ObsCarrLocal);
        
        for ll=1:size(CarrToSend,1)
            CarrUsed = CarrToSend(ll,:);
            ObsCarrUsed(logical(CarrUsed)) = ObsCarrPos(logical(CarrUsed));
            MAIN;
            if CurrentModula == 1
                VetOpen(ll)   = AberLev(logical(CarrToSend(1,:)));
                VetNumCar(ll) = sum(CarrUsed);
                figure(661);
                hold all;
                grid on;
                plot(VetNumCar(ll),VetOpen(ll),'x','LineWidth',2);
                set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                drawnow;
            elseif CurrentModula == 2
            end
        end
        plot(VetNumCar,VetOpen,'b:','LineWidth',2);
    end
    %%
    SendingCarr = 0;
    if SendingCarr
        CarrToSendInt       = zeros(1,NumCarrUsed);
        for kk=0:NumCarrUsed-1
            CarrToSendInt(kk+1,CarrPosInt(1:1+kk)) = 1;
        end
        [~,ObsCarrLocal] = max(sum(CarrToSendInt));
        ObsCarr          = ObsCarrPos(ObsCarrLocal);
        
        for ll=1:size(CarrToSendInt,1)
            CarrUsed = CarrToSendInt(ll,:);
            ObsCarrUsed(logical(CarrUsed)) = ObsCarrPos(logical(CarrUsed));
            MAIN;
            if CurrentModula == 1
                VetOpen(ll)   = AberLev(logical(CarrToSendInt(1,:)));
                VetNumCar(ll) = sum(CarrUsed);
                figure(661);
                hold all;
                grid on;
                plot(VetNumCar(ll),VetOpen(ll),'x','LineWidth',2);
                set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                drawnow;
            elseif CurrentModula == 2
            end
        end
        plot(VetNumCar,VetOpen,'r:','LineWidth',2);
    end
    %%
    SendingCarr = 0;
    if SendingCarr
        CarrToSend2       = zeros(1,NumCarrUsed);
        for kk=0:length(CarrOddPos)-1
            CarrToSend2(kk+1,CarrOddPos(end-kk:end)) = 1;
        end
        [~,ObsCarrLocal] = max(sum(CarrToSend2));
        ObsCarr          = ObsCarrPos(ObsCarrLocal);
        
        for ll=1:size(CarrToSend2,1)
            CarrUsed = CarrToSend2(ll,:);
            ObsCarrUsed(logical(CarrUsed)) = ObsCarrPos(logical(CarrUsed));
            MAIN;
            if CurrentModula == 1
                VetOpen(ll)   = AberLev(logical(CarrToSend2(1,:)));
                VetNumCar(ll) = sum(CarrUsed);
                figure(661);
                hold all;
                grid on;
                plot(VetNumCar(ll),VetOpen(ll),'x','LineWidth',2);
                set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                drawnow;
            elseif CurrentModula == 2
            end
        end
        plot(VetNumCar,VetOpen,'g:','LineWidth',2);
    end
    %%
    SendingCarr = 0;
    if SendingCarr
        CarrToSendInt2       = zeros(1,NumCarrUsed);
        for kk=0:NumCarrUsed-1
            CarrToSendInt2(kk+1,CarrPosInt(end-kk:end)) = 1;
        end
        [~,ObsCarrLocal] = max(sum(CarrToSendInt2));
        ObsCarr          = ObsCarrPos(ObsCarrLocal);
        
        for ll=1:size(CarrToSendInt2,1)
            CarrUsed = CarrToSendInt2(ll,:);
            ObsCarrUsed(logical(CarrUsed)) = ObsCarrPos(logical(CarrUsed));
            MAIN;
            if CurrentModula == 1
                VetOpen(ll)   = AberLev(logical(CarrToSendInt2(1,:)));
                VetNumCar(ll) = sum(CarrUsed);
                figure(661);
                hold all;
                grid on;
                plot(VetNumCar(ll),VetOpen(ll),'x','LineWidth',2);
                set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                drawnow;
            elseif CurrentModula == 2
            end
        end
        plot(VetNumCar,VetOpen,'y:','LineWidth',2);
    end
    %%
    SendingCarr = 0;
    if SendingCarr
        CarrToSend3       = zeros(1,NumCarrUsed);
        for kk=0:length(CarrEvenPos)-1
            CarrToSend3(kk+1,CarrEvenPos(1:1+kk)) = 1;
        end
        [~,ObsCarrLocal] = max(sum(CarrToSend3));
        ObsCarr          = ObsCarrPos(ObsCarrLocal);
        
        for ll=1:size(CarrToSend3,1)
            CarrUsed = CarrToSend3(ll,:);
            ObsCarrUsed(logical(CarrUsed)) = ObsCarrPos(logical(CarrUsed));
            MAIN;
            if CurrentModula == 1
                VetOpen(ll)   = AberLev(logical(CarrToSend3(1,:)));
                VetNumCar(ll) = sum(CarrUsed);
                figure(661);
                hold all;
                grid on;
                plot(VetNumCar(ll),VetOpen(ll),'x','LineWidth',2);
                set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                drawnow;
            elseif CurrentModula == 2
            end
        end
        plot(VetNumCar,VetOpen,'c:','LineWidth',2);
    end
    %%
    SendingCarr = 0;
    if SendingCarr
        CarrToSend4       = zeros(1,NumCarrUsed);
        for kk=0:length(CarrEvenPos)-1
            CarrToSend4(kk+1,CarrEvenPos(end-kk:end)) = 1;
        end
        [~,ObsCarrLocal] = max(sum(CarrToSend4));
        ObsCarr          = ObsCarrPos(ObsCarrLocal);
        
        for ll=1:size(CarrToSend4,1)
            CarrUsed = CarrToSend4(ll,:);
            ObsCarrUsed(logical(CarrUsed)) = ObsCarrPos(logical(CarrUsed));
            MAIN;
            if CurrentModula == 1
                VetOpen(ll)   = AberLev(logical(CarrToSend4(1,:)));
                VetNumCar(ll) = sum(CarrUsed);
                figure(661);
                hold all;
                grid on;
                plot(VetNumCar(ll),VetOpen(ll),'x','LineWidth',2);
                set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                drawnow;
            elseif CurrentModula == 2
            end
        end
        plot(VetNumCar,VetOpen,'k:','LineWidth',2);
    end
    
    %%
    SendingCarr = 0;
    if SendingCarr
        CarrToSend6       = zeros(1,NumCarrUsed);
        for kk=0:(length(CarrEvenPos))-17
            CarrToSend6(kk+1,CarrEvenPos((end/2)-kk:(end/2)+kk)) = 1;
        end
        [~,ObsCarrLocal] = max(sum(CarrToSend6));
        ObsCarr          = ObsCarrPos(ObsCarrLocal);
        
        for ll=1:size(CarrToSend6,1)
            CarrUsed = CarrToSend6(ll,:);
            ObsCarrUsed(logical(CarrUsed)) = ObsCarrPos(logical(CarrUsed));
            MAIN;
            if CurrentModula == 1
                VetOpen(ll)   = AberLev(logical(CarrToSend6(1,:)));
                VetNumCar(ll) = sum(CarrUsed);
                figure(661);
                hold all;
                grid on;
                plot(VetNumCar(ll),VetOpen(ll),'x','LineWidth',2);
                set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                drawnow;
            elseif CurrentModula == 2
            end
        end
        plot(VetNumCar,VetOpen,'b:','LineWidth',2);
        drawnow;
    end
end
    %%
    SendingCarr = 0;
    if SendingCarr
        CarrToSend5       = zeros(1,NumCarrUsed);
        for kk=0:(length(CarrOddPos))-9%for 16 set this to 5; for 32 set this to 9; for 64 set this to 17
            CarrToSend5(kk+1,CarrOddPos((end/2)-kk:(end/2)+kk)) = 1;
        end
        [~,ObsCarrLocal] = max(sum(CarrToSend5));
        ObsCarr          = ObsCarrPos(ObsCarrLocal);
        for nn=1:10
            nn
            for ll=1:size(CarrToSend5,1)
                ll
                CarrUsed = CarrToSend5(ll,:);
                ObsCarrUsed(logical(CarrUsed)) = ObsCarrPos(logical(CarrUsed));
                MAIN;
                if CurrentModula == 1
                    BERtOOK(nn,ll)    = BerOOK(logical(CarrToSend5(1,:)));
                    VetOpen2(nn,ll)   = AberLev(logical(CarrToSend5(1,:)));
                    VetNumCar2(nn,ll) = sum(CarrUsed);
                    if IfftOrSum
                        figure(661);
                        hold all;
                        grid on;
                        plot(VetNumCar2(nn,ll),VetOpen2(nn,ll)/VetOpen2(nn,1),'gd','LineWidth',2);
                        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                        drawnow;
                    else
                        figure(662);
                        hold all;
                        grid on;
                        plot(VetNumCar2(nn,ll),VetOpen2(nn,ll)/VetOpen2(nn,1),'bx','LineWidth',2);
                        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                        drawnow;
                    end
                elseif CurrentModula == 2
                    BERt4PAM(nn,ll)    = Ber4PAM(logical(CarrToSend5(1,:)));
                    VetOpenO1(nn,ll)  = AberLev(1,logical(CarrToSend5(1,:)),CurrentTest);
                    VetOpenO2(nn,ll)  = AberLev(2,logical(CarrToSend5(1,:)),CurrentTest);
                    VetOpenO3(nn,ll)  = AberLev(3,logical(CarrToSend5(1,:)),CurrentTest);
                    VetNumCar2(nn,ll) = sum(CarrUsed);
                    if IfftOrSum
                        figure(663);
                        hold all;
                        grid on;
                        plot(VetNumCar2(nn,ll),VetOpenO1(nn,ll)/VetOpenO1(nn,1),'gd','LineWidth',2);
                        plot(VetNumCar2(nn,ll),VetOpenO2(nn,ll)/VetOpenO2(nn,1),'yd','LineWidth',2);
                        plot(VetNumCar2(nn,ll),VetOpenO3(nn,ll)/VetOpenO3(nn,1),'md','LineWidth',2);
                        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                        drawnow;
                    else
                        figure(664);
                        hold all;
                        grid on;
                        plot(VetNumCar2(nn,ll),VetOpenO1(nn,ll)/VetOpenO1(nn,1),'bd','LineWidth',2);
                        plot(VetNumCar2(nn,ll),VetOpenO2(nn,ll)/VetOpenO2(nn,1),'rd','LineWidth',2);
                        plot(VetNumCar2(nn,ll),VetOpenO3(nn,ll)/VetOpenO3(nn,1),'kd','LineWidth',2);
                        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                        drawnow;
                    end
                end
            end
            if CurrentModula == 1
                if IfftOrSum
                    plot(VetNumCar2(nn,:),VetOpen2(nn,:)./VetOpen2(nn,1),'g:','LineWidth',2);
                    drawnow;
                else
                    plot(VetNumCar2(nn,:),VetOpen2(nn,:)./VetOpen2(nn,1),'b:','LineWidth',2);
                    drawnow;
                end
                save(SavingAT2,'VetNumCar2','VetOpen2','BERtOOK');
            elseif CurrentModula == 2
                if IfftOrSum
                    plot(VetNumCar2(nn,:),VetOpenO1(nn,:)/VetOpenO1(nn,1),'g:','LineWidth',2);
                    plot(VetNumCar2(nn,:),VetOpenO2(nn,:)/VetOpenO2(nn,1),'y:','LineWidth',2);
                    plot(VetNumCar2(nn,:),VetOpenO3(nn,:)/VetOpenO3(nn,1),'m:','LineWidth',2);
                    drawnow;
                else
                    plot(VetNumCar2(nn,:),VetOpenO1(nn,:)/VetOpenO1(nn,1),'b:','LineWidth',2);
                    plot(VetNumCar2(nn,:),VetOpenO2(nn,:)/VetOpenO2(nn,1),'r:','LineWidth',2);
                    plot(VetNumCar2(nn,:),VetOpenO3(nn,:)/VetOpenO3(nn,1),'k:','LineWidth',2);
                    drawnow;
                end
                save(SavingAT2,'VetNumCar2','VetOpenO1','VetOpenO2','VetOpenO3','BERt4PAM');
            end
        end
        
    end
    %%
    SendingCarr = 0;
    if SendingCarr
        CarrToSendInt3       = zeros(1,NumCarrUsed);
        for kk=0:(NumCarrUsed)-17%for 16 set this to 9; for 32 set this to 17; for 64 set this to 33
            CarrToSendInt3(kk+1,CarrPosInt((end/2)-kk:(end/2)+kk)) = 1;
        end
        [~,ObsCarrLocal] = max(sum(CarrToSendInt3));
        ObsCarr          = ObsCarrPos(ObsCarrLocal);
        for nn=1:10
            nn
            for ll=1:size(CarrToSendInt3,1)
                CarrUsed = CarrToSendInt3(ll,:);
                ObsCarrUsed(logical(CarrUsed)) = ObsCarrPos(logical(CarrUsed));
                MAIN;
                if CurrentModula == 1
                    BERtOOK(nn,ll)    = BerOOK(logical(CarrToSend5(1,:)));
                    VetOpen(nn,ll)   = AberLev(logical(CarrToSendInt3(1,:)));
                    VetNumCar(nn,ll) = sum(CarrUsed);
                    if IfftOrSum
                        figure(661);
                        hold all;
                        grid on;
                        plot(VetNumCar(nn,ll),VetOpen(nn,ll)/VetOpen(nn,1),'kd','LineWidth',2);
                        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                        drawnow;
                    else
                        figure(662);
                        hold all;
                        grid on;
                        plot(VetNumCar(nn,ll),VetOpen(nn,ll)/VetOpen(nn,1),'rx','LineWidth',2);
                        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                        drawnow;
                    end
                elseif CurrentModula == 2
                    BERt4PAM(nn,ll)    = Ber4PAM(logical(CarrToSend5(1,:)));
                    VetOpenA1(nn,ll)  = AberLev(1,logical(CarrToSend5(1,:)),CurrentTest);
                    VetOpenA2(nn,ll)  = AberLev(2,logical(CarrToSend5(1,:)),CurrentTest);
                    VetOpenA3(nn,ll)  = AberLev(3,logical(CarrToSend5(1,:)),CurrentTest);
                    VetNumCar(nn,ll) = sum(CarrUsed);
                    if IfftOrSum
                        figure(663);
                        hold all;
                        grid on;
                        plot(VetNumCar(nn,ll),VetOpenA1(nn,ll)/VetOpenA1(nn,1),'gd','LineWidth',2);
                        plot(VetNumCar(nn,ll),VetOpenA2(nn,ll)/VetOpenA2(nn,1),'yd','LineWidth',2);
                        plot(VetNumCar(nn,ll),VetOpenA3(nn,ll)/VetOpenA3(nn,1),'md','LineWidth',2);
                        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                        drawnow;
                    else
                        figure(664);
                        hold all;
                        grid on;
                        plot(VetNumCar(nn,ll),VetOpenA1(nn,ll)/VetOpenA1(nn,1),'bd','LineWidth',2);
                        plot(VetNumCar(nn,ll),VetOpenA2(nn,ll)/VetOpenA2(nn,1),'rd','LineWidth',2);
                        plot(VetNumCar(nn,ll),VetOpenA3(nn,ll)/VetOpenA3(nn,1),'kd','LineWidth',2);
                        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
                        drawnow;
                    end
                end
            end
            if CurrentModula == 1
                if IfftOrSum
                    plot(VetNumCar(nn,:),VetOpen2(nn,:)./VetOpen2(nn,1),'g:','LineWidth',2);
                    drawnow;
                else
                    plot(VetNumCar(nn,:),VetOpen2(nn,:)./VetOpen2(nn,1),'b:','LineWidth',2);
                    drawnow;
                end
                save(SavingAT2,'VetNumCar2','VetOpen2','BERtOOK');
            elseif CurrentModula == 2
                if IfftOrSum
                    plot(VetNumCar(nn,:),VetOpenA1(nn,:)/VetOpenA1(nn,1),'g:','LineWidth',2);
                    plot(VetNumCar(nn,:),VetOpenA2(nn,:)/VetOpenA2(nn,1),'y:','LineWidth',2);
                    plot(VetNumCar(nn,:),VetOpenA3(nn,:)/VetOpenA3(nn,1),'m:','LineWidth',2);
                    drawnow;
                else
                    plot(VetNumCar(nn,:),VetOpenA1(nn,:)/VetOpenA1(nn,1),'b:','LineWidth',2);
                    plot(VetNumCar(nn,:),VetOpenA2(nn,:)/VetOpenA2(nn,1),'r:','LineWidth',2);
                    plot(VetNumCar(nn,:),VetOpenA3(nn,:)/VetOpenA3(nn,1),'k:','LineWidth',2);
                    drawnow;
                end
                save(SavingAT2,'VetNumCar2','VetOpenO1','VetOpenO2','VetOpenO3','BERt4PAM');
            end
        end
    end
%     if IfftOrSum
%         save(SavingAT2,'VetNumCar2','VetOpenO1','VetOpenO2','VetOpenO3');
%     else
%         save(SavingAT2,'VetNumCar2','VetOpenA1','VetOpenA2','VetOpenA3');
%     end
    
end
a=1;
% for kk=0:length(CarrOddPos)
%     CarrToSend2(kk+1,CarrOddPos(end-kk:end)) = 1;
% end
% CarrToSend2

%%
% MAIN;