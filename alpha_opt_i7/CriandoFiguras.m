switch ThisFigure
    case 1
        figure(ThisFigure);
        hold all;
        UpLink = 0;
        DownLink = 0;
        for kk=1:NumCarr
            if ~mod(kk,2)
                UpLink = UpLink + Eout(kk,:);
            else
                DownLink = DownLink + Eout(kk,:);
            end
        end
        plot(f,20*log10(abs(fftshift(fft(UpLink./length(UpLink))))),'r');
        plot(f,20*log10(abs(fftshift(fft(DownLink./length(DownLink))))),'b');
        title('DownLink and UpLink Carriers','FontSize',16,'FontWeight','bold');
        xlabel('Frequency [Hz]','FontSize',14);
        ylabel('Amplitude [u]','FontSize',14);
        legend({'UpLink','DownLink'},'FontSize',12,'Location','best',...
                                  'FontWeight','bold','Box','off');
        grid on;
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
        a=1;
        close all;
    case 2
        figure(ThisFigure);
        plot(f,db(abs((fftshift(fft(Eout))))/length(Eout)));
        title('Output Signal, Eout','FontSize',14);
        ylabel('Magnitude [db]','FontSize',12);
        xlabel('Frequency [Hz]','FontSize',12);
        grid on;
        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
        
        a=2;
        close all;
    case 3
        figure(ThisFigure);
        plot(f,db(abs((fftshift(fft(CW))))/length(CW)),'x');
        title('Output Signal, Eout','FontSize',14);
        ylabel('Magnitude [db]','FontSize',12);
        xlabel('Frequency [Hz]','FontSize',12);
        grid on;
        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
        
        a=3;
        close all;
    case 4
        figure(ThisFigure);
        hold all;
        for jj=1:NPPB:length(EoutModAux)
            plot(real(EoutModAux(jj+NPPB/2))./max(real(EoutModAux)),imag(EoutModAux(jj+NPPB/2))./max(imag(EoutModAux)),'xr','LineWidth',2);
        end
        title('Scatter Plor DPSK','FontSize',14);
        ylabel('Imaginary [u]','FontSize',12);
        xlabel('Real [u]','FontSize',12);
        grid on;
        axis([-1 1 -1 1]);
        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
        a=4;
        close all;
    case 5
        figure(ThisFigure);
        hold all;
        for jj=1:NPPB:length(EoutModAux)
            plot(real(EoutModAux(jj+NPPB/2))./max(real(EoutModAux)),imag(EoutModAux(jj+NPPB/2))./max(imag(EoutModAux)),'xr','LineWidth',2);
%             rex = (abs(EoutModAux(jj+NPPB/2)).^2)/max(abs(EoutModAux).^2);
%             rex = rex-min(rex);
%                 plot(rex,rex,'xr','LineWidth',2);
        end
        title('Scatter Plor OOK','FontSize',14);
        ylabel('Imaginary [u]','FontSize',12);
        xlabel('Real [u]','FontSize',12);
        grid on;
%         axis([-1 1 -1 1]);
        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
        a=5;
        close all;
    case 6
        figure(ThisFigure);
        hold all;
        for jj=1:NPPB:length(EoutModAux)
%             plot(real(EoutModAux(jj+NPPB/2))./max(real(EoutModAux)),imag(EoutModAux(jj+NPPB/2))./max(imag(EoutModAux)),'xr','LineWidth',2);
%             plot(abs(EoutModAux(jj+NPPB/2))./max(abs(EoutModAux)),abs(EoutModAux(jj+NPPB/2))./max(abs(EoutModAux)),'xr','LineWidth',2);
            rex = (abs(EoutModAux(jj+NPPB/2)).^2)/max(abs(EoutModAux).^2);
%             rex = rex-min(rex);
                plot(rex,rex,'xr','LineWidth',2);
        end
        title('Scatter Plor 4PAM','FontSize',14);
        ylabel('Imaginary [u]','FontSize',12);
        xlabel('Real [u]','FontSize',12);
        grid on;
%         axis([-1 1 -1 1]);
        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
        a=5;
        close all;
    case 7
        figure(ThisFigure);
        plot(f,db(abs((fftshift(fft(Eout))))/length(Eout)));
       title('Output Signal, Eout','FontSize',14);
       ylabel('Magnitude [db]','FontSize',12);
       xlabel('Frequency [Hz]','FontSize',12);
       grid on;
%         axis([-1 1 -1 1]);
        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
        a=5;
        close all;
    case 8
        figure(ThisFigure);
        hold all;
        for jj=1:NPPB:length(EoutModAux)
            plot(real(EoutModAux(jj+NPPB/2))./max(real(EoutModAux)),imag(EoutModAux(jj+NPPB/2))./max(imag(EoutModAux)),'xr','LineWidth',2);
        end
        title('Scatter Plor DQPSK','FontSize',14);
        ylabel('Imaginary [u]','FontSize',12);
        xlabel('Real [u]','FontSize',12);
        grid on;
        axis([-1 1 -1 1]);
        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
        a=5;
        close all;
    case 9
        Olho(abs(EoutModAux),T,NPPB,1);
        title('Eye Diagram of DQPSK','FontSize',18);
        ylabel('Amplitude [u]','FontSize',16);
        xlabel('Time [s]','FontSize',16);
        grid on;
%         axis([-1 1 -1 1]);
        a=5;
        close all;
    case 10
        Olho(abs(EoutModAux),T,NPPB,1);
        title('Eye Diagram of DPSK','FontSize',18);
        ylabel('Amplitude [u]','FontSize',16);
        xlabel('Time [s]','FontSize',16);
        grid on;
%         axis([-1 1 -1 1]);
        a=5;
        close all;
    case 11
        Olho(abs(EoutModAux).^2,T,NPPB,1);
        title('Eye Diagram of 4PAM','FontSize',18);
        ylabel('Amplitude [u]','FontSize',16);
        xlabel('Time [s]','FontSize',16);
        grid on;
%         axis([-1 1 -1 1]);
        a=5;
        close all;
    case 12
        Olho(abs(EoutModAux),T,NPPB,1);
        title('Eye Diagram of OOK','FontSize',18);
        ylabel('Amplitude [u]','FontSize',16);
        xlabel('Time [s]','FontSize',16);
        grid on;
%         axis([-1 1 -1 1]);
        a=5;
        close all;
    case 13
        Olho(Ix,T,NPPB,1);
        title('Eye Diagram of 4PAM','FontSize',18);
        ylabel('Amplitude [u]','FontSize',16);
        xlabel('Time [s]','FontSize',16);
        grid on;
%         axis([-1 1 -1 1]);
        a=5;
        close all;
    otherwise
end