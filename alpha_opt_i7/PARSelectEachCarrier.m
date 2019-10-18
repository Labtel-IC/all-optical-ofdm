function [Eout]=PARSelectEachCarrier(Ein,NumCarr,f,fin,FBWD,Order)
% figure;
% hold all;
% grid on;
% plot(f,20.*log10(abs(fftshift(fft(Ein./length(Ein))))));
parfor kk=1:NumCarr
    [FiltEin,~]=FiltroGaussiano(f,FBWD,kk*fin,Order);
    FiltEin = fftshift(FiltEin);
    Eout(kk,:) = ifft(fft(Ein).*FiltEin);
%     plot(f,20.*log10(fftshift(FiltEin)),f,20.*log10(abs(fftshift(fft(Eout(kk,:)./length(Eout(kk,:)))))));
end
% ThisFigure = 1;
% CriandoFiguras;
