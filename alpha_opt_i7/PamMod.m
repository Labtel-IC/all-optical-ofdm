function [Eout] = PamMod (Ein,Sig1,Sig2,U_pi2,Vbias,t)

%This test will implemente the IQ modulator as theoretically as possible.
%The CW signal will pass through an optical coupler that will divided the
%signal in two. Each part will then pass through a MZM. At the end of each
%MZM the signal will be combined with an optical couples as well. One of
%thos signals will be out phased by 90 degrees. Finaly anoter optical
%coupler will combine those two signal to form the  output.

%%
U0    = Vbias;
U1t   = Sig1;
U2t   = Sig2;
U_pi  = U_pi2;
Vphi  = (U1t - U2t - U0);
Eout = Ein.*(1/2).*cos(Vphi.*(pi/(2*U_pi))).*exp(-1j*(pi/2).*((U1t + U2t)/U_pi));

%%
if nargin==6
    f     = time2freq(t);
    figure;
    plot(f,20*log10(abs(fftshift(fft(Eout)./length(Eout))))); 

    title('Spectrum of the Optical Field from the IQ-MZM','FontSize',16,...
                                                      'FontWeight','bold');
    xlabel('Frequency [Hz]','FontSize',14);%,'FontWeight','bold');
    ylabel('Amplitude [dB]','FontSize',14);%,'FontWeight','bold');
    legend({'Eout'},'FontSize',12,'Location','best','FontWeight','bold',...
                                                              'Box','off');
    grid on;
    figure;
    plot(t,abs(Eout),t,real(Eout),t,imag(Eout));
%     axis([0 16*1/(12.5e9) 1.1*min(real(Eout)) 1.1*max(real(Eout))]);

    title('Comparison of Optical Field Components IQ-MZM','FontSize',16,...
                                                      'FontWeight','bold');
    xlabel('Time [s]','FontSize',14);%,'FontWeight','bold');
    ylabel('Amplitude [u]','FontSize',14);%,'FontWeight','bold');
    legend({'Eout','Eout-Re','Eout-Im'},'FontSize',12,'Location','best',...
                                          'FontWeight','bold','Box','off');
    grid on;
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    a=1;
end
