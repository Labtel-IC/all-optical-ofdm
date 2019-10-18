function [ Eout ] = OpticalFFT(t,SymPer,MaxNumStag,E1,E0,ActStag...
                                                                   ,PhaDel)
    if nargin<7
        PhaDel = pi - pi/(2^(ActStag-1));
    elseif nargin<6
        ActStag = 1;
        PhaDel = pi - pi/(2^(ActStag-1));
    elseif nargin<5
        ActStag = 1;
        PhaDel = pi - pi/(2^(ActStag-1));
        E0 = 0;
    elseif nargin<4
        error(['Input Arguments are not enougth. Please check this function'...
                                                                    ' Call!']);
    end
    TimDel = SymPer/(2^ActStag);
    [Eout1,Eout2] = DelayInterf(t,TimDel,PhaDel,E1,E0);
    
%     if Ploting
%         figure;
%         hold all;
%         plot(time2freq(t),db(abs((fftshift(fft(Eout1))))/length(Eout1)));
%         plot(time2freq(t),db(abs((fftshift(fft(Eout1))))/length(Eout1)));
%         axis([1e11 3e11 -370 10]);
%     end
    if ActStag < MaxNumStag
        if (PhaDel < pi/2)&&(ActStag>1)
            ModPhaDel = pi/2;
        else
            ModPhaDel = pi/(2^ActStag);
        end
        PhaDel = PhaDel + ModPhaDel;
        ActStag = ActStag + 1;
        [EoutAux1] = OpticalFFT(t,SymPer,MaxNumStag,Eout1,0,ActStag,PhaDel);
        PhaDel = PhaDel - pi/2;
        [EoutAux2] = OpticalFFT(t,SymPer,MaxNumStag,Eout2,0,ActStag,PhaDel);
        Eout = [EoutAux1;EoutAux2];
    else
        Eout = [Eout1;Eout2];
    end
end

