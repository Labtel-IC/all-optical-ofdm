function [Eout]=SelectEachCarrier(Ein,NumCarr,f,fin,FBWD,Order,fc)
% figure;
% hold all;
% grid on;
% plot(f,20.*log10(abs(fftshift(fft(Ein./length(Ein))))));
for kk=0:NumCarr-1
    [FiltEin,~]=FiltroGaussiano(f,FBWD,kk*fc+fin,Order);
    FiltEin = fftshift(FiltEin);
    Eout(1+kk,:) = ifft(fft(Ein).*FiltEin);
%     plot(f,20.*log10(fftshift(FiltEin)),f,20.*log10(abs(fftshift(fft(Eout(kk,:)./length(Eout(kk,:)))))));
end
% ThisFigure = 1;
% CriandoFiguras;
