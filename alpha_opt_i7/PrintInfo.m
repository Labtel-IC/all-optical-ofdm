function []=PrintInfo(ThisCase,A,B,C,D,E,F,G,H,I,J,K)
%%
if nargin <1
    error('Not enough input arguments. Please check this function call');
end

switch ThisCase
    case 0
        
    case 1%The main input signal
        figure;
        plot(A,db(abs((fftshift(fft(B))))/length(B)));
        title('Output Signal, Eout');
        ylabel('Magnitude [db]');
        xlabel('Frequency [Hz]');
        ThisFig = gca;
        ThisFig.FontSize = 40;
        ThisFig.FontName = 'Times New Roman';
        ThisFig.Box = 'on';
%         ThisFig.LineWidth = 2;
        grid on;
        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
    %     axis([-4e11 4e11 -440 0]);
    case 2%Eye of the DPSK input signal
        Olho(A,B,C,1);

    case 3%Eye of the DQPSK input signal
        Olho(A,B,C,1);
        Olho(D,B,C,1);
    case 4%Eye of the 4PAM input signal
        Olho(A,B,C,1);
        Olho(D,B,C,1);
    case 5%Comparison between data and OOK input signal
        figure;
        hold all;
        plot(A,B);
        plot(A,C);
        title('Checking Modulating Signal','FontSize',14);
        ylabel('Amplitude [A]','FontSize',12);
        xlabel('Time [s]','FontSize',12);
%             axis([-0.5*T 0.5e2*(T) -1.2 1.2]);
        grid on;
        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
        a=1;
    case 6%Eye diagram of DPSK modulated signal
        Olho(A,B,C,1);
    case 7%Eye diagram of DQPSK modulated signal
        Olho(A,B,C,1);
    case 8%Eye diagram of 4PAM modulated signal after photo diodo
        Olho(A,B,C,1);
    case 9%Eye diagram of OOK modulated signal
        Olho(A,B,C,1);
    case 10
        figure;
        hold all
        plot(A,B);
        plot(A,C);
        plot(A,D);
        plot(A,E);
        a=1;
        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
%     axis([1e11 2e11 -120 20]);
    case 11%Eye diagram of OOK modulated signal
        Olho(A,B,C,1);
    case 12
        figure;
        hold all;
        plot(A,B);
        for kk=1:size(C,1)
            plot(A,db(abs((fftshift(fft(C(kk,:)))))/length(C(kk,:))));
        end
        axis([0.9*D D+(size(C,1))*1.1*12.5e9 -250 10]);
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
    case 13
        figure;
        plot(A,B);
        title('Received Optical signal at time domain',...
                                'FontSize',16,'FontWeight','bold');
        xlabel('Time [s]','FontSize',14);
        ylabel('Amplitude [u]','FontSize',14);
        legend({'Erec'},'FontSize',12,'Location','best',...
                                  'FontWeight','bold','Box','off');
        grid on;
        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
    case 14
        figure;
        plot(A,B);
        title('Spectrum of the Received Optical signal',...
                                'FontSize',16,'FontWeight','bold');
        xlabel('Frequency [Hz]','FontSize',14);
        ylabel('Amplitude [dB]','FontSize',14);
        legend({'Erec'},'FontSize',12,'Location','best',...
                                  'FontWeight','bold','Box','off');
        grid on;
        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
    case 15
        figure;
        plot(A,B,A,C,A,D,A,E);
        title('Comparing received symbol and synchronism symbol'...
                               ,'FontSize',16,'FontWeight','bold');
        xlabel('Time [s]','FontSize',14);
        ylabel('Amplitude [u]','FontSize',14);
        legend({'Rec-Sym','Sycn','Rec-Sym1','Rec-Sym2'},...
        'FontSize',12,'Location','best','FontWeight','bold',...
                                                      'Box','off');
        grid on;
        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
    case 16
        figure;plot(A,B,A,C,A,D,A,E);
        title('Comparing received symbol and synchronism symbol',...
                                        'FontSize',16,'FontWeight','bold');
        xlabel('Time [s]','FontSize',14);
        ylabel('Amplitude [u]','FontSize',14);
        legend({'Rec-Sym','Sycn','Rec-Sym1','Rec-Sym2'},'FontSize',12,...
                        'Location','best','FontWeight','bold','Box','off');
        grid on;
        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
    case 17
        Olho(A,B,C,1);
        title('Eye diagram of the received signal');
        xlabel('Time [s]');
        ylabel('Amplitude [u]');
        subplot('position',[0.10 0.09 0.15 0.85]);
        ThisFig = gca;
        ThisFig.FontSize = 20;
        ThisFig.FontName = 'Times New Roman';
        ThisFig.FontWeight = 'bold';
        ThisFig.Box = 'on';
        subplot('position',[0.33 0.09 0.60 0.85]);
        ThisFig = gca;
        ThisFig.FontSize = 20;
        ThisFig.FontName = 'Times New Roman';
        ThisFig.FontWeight = 'bold';
        ThisFig.Box = 'on';
%         ThisFig.LineWidth = 2;
        legend({'Erec'},'FontSize',12,'Location','best',...
                                  'FontWeight','bold','Box','off');
        grid on;
    case 18
        figure;
        plot(A,B,A,C);
        axis([0 12.8e-10 -1.1 1.1]);
        title('Comparing Transmited and Received data',...
                                'FontSize',16,'FontWeight','bold');
        xlabel('Time [s]','FontSize',14);
        ylabel('Amplitude [u]','FontSize',14);
        legend({'Tx','Rec-Odd'},'FontSize',12,'Location','best',...
                                  'FontWeight','bold','Box','off');
        grid on;
        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
    case 19
        figure;hold all;
        plot(A);
        plot(B);
        title('Final Comparison of Tx and RX data','FontSize',16,...
                                                      'FontWeight','bold');
        xlabel('Time [s]','FontSize',14);
        ylabel('Amplitude [u]','FontSize',14);
        legend({'Tx','Rx'},'FontSize',12,'Location',...
                           'best','FontWeight','bold','Box','off');
        grid on;
        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
    case 20
        figure;
        plot(A,B);
        title('Received Optical signal at time domain',...
                                'FontSize',16,'FontWeight','bold');
        xlabel('Time [s]','FontSize',14);
        ylabel('Amplitude [u]','FontSize',14);
        legend({'Erec'},'FontSize',12,'Location','best',...
                                  'FontWeight','bold','Box','off');
        grid on;
        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
    case 21
        figure;
        plot(A,B);
        title('Spectrum of the Received Optical signal','FontSize',16,...
                                                      'FontWeight','bold');
        xlabel('Frequency [Hz]','FontSize',14);
        ylabel('Amplitude [dB]','FontSize',14);
        legend({'Erec'},'FontSize',12,'Location','best','FontWeight',...
                                                       'bold','Box','off');
        grid on;
        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
    case 22
        figure;
        plot(A,B,A,C,A,D,A,E);
        title('Comparing received symbol and synchronism symbol'...
                               ,'FontSize',16,'FontWeight','bold');
        xlabel('Time [s]','FontSize',14);
        ylabel('Amplitude [u]','FontSize',14);
        legend({'Rec-Sym','Sycn','Rec-Sym1','Rec-Sym2'},...
        'FontSize',12,'Location','best','FontWeight','bold',...
                                                      'Box','off');
        grid on;
        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
    case 23
        figure;
        plot(A,B,A,C,A,D,A,E);
        title('Comparing received symbol and synchronism symbol'...
                               ,'FontSize',16,'FontWeight','bold');
        xlabel('Time [s]','FontSize',14);
        ylabel('Amplitude [u]','FontSize',14);
        legend({'Rec-Sym','Sycn','Rec-Sym1','Rec-Sym2'},...
        'FontSize',12,'Location','best','FontWeight','bold',...
                                                      'Box','off');
        grid on;
        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
    case 24
        Olho(A,B,C,1);
        title('Eye diagram of the received signal - DataOdd',...
                                'FontSize',16,'FontWeight','bold');
        xlabel('Time [s]','FontSize',14);
        ylabel('Amplitude [u]','FontSize',14);
        legend({'Icoded'},'FontSize',12,'Location','best',...
                                  'FontWeight','bold','Box','off');
        grid on;
    case 25
        Olho(A,B,C,1);
        title('Eye diagram of the Theoretical signal - DataOdd',...
                                'FontSize',16,'FontWeight','bold');
        xlabel('Time [s]','FontSize',14);
        ylabel('Amplitude [u]','FontSize',14);
        legend({'Icoded'},'FontSize',12,'Location','best',...
                                  'FontWeight','bold','Box','off');
        grid on;
    case 26
        Olho(A,B,C,1);
        title('Eye diagram of the received signal - DataEven',...
                                'FontSize',16,'FontWeight','bold');
        xlabel('Time [s]','FontSize',14);
        ylabel('Amplitude [u]','FontSize',14);
        legend({'Qcoded'},'FontSize',12,'Location','best',...
                                  'FontWeight','bold','Box','off');
        grid on;
    case 27
        Olho(A,B,C,1);
        title('Eye diagram of the Theoretical signal - DataEven'...
                               ,'FontSize',16,'FontWeight','bold');
        xlabel('Time [s]','FontSize',14);
        ylabel('Amplitude [u]','FontSize',14);
        legend({'Qcoded'},'FontSize',12,'Location','best',...
                                  'FontWeight','bold','Box','off');
        grid on;
    case 28
        figure;
        plot(A,B,A,C,A,D,A,E,A,F,A,G,A,H,A,I,A,J,A,H,A,K);
        axis([0 6.72e-9 -1.1 1.1]);
        title(['Comparing Transmited, Received and Theoretical '...
                        'data'],'FontSize',16,'FontWeight','bold');
        xlabel('Time [s]','FontSize',14);
        ylabel('Amplitude [u]','FontSize',14);
        legend({'Tx-Odd','Tx-Even','Rec-Odd','Rec-OddCon',...
        'Rec-OddDes','Rec-Even','Rec-EvenCon','Rec-EvenDes',...
        'Theo-Odd','Theo-Even'},'FontSize',12,'Location','best',...
                                  'FontWeight','bold','Box','off');
        grid on;
        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
    case 29
        figure;hold all;
        plot(A);
        plot(B);
        title('Final Comparison of Tx and RX data - Odd',...
                                'FontSize',16,'FontWeight','bold');
        xlabel('Time [s]','FontSize',14);
        ylabel('Amplitude [u]','FontSize',14);
        legend({'Tx-Odd','Rx-Odd'},'FontSize',12,'Location',...
                           'best','FontWeight','bold','Box','off');
        grid on;
        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
    case 30
        figure;hold all;
        plot(A);
        plot(B);
        title('Final Comparison of Tx and RX data - Even',...
                                'FontSize',16,'FontWeight','bold');
        xlabel('Time [s]','FontSize',14);
        ylabel('Amplitude [u]','FontSize',14);
        legend({'Tx-Even','Rx-Even'},'FontSize',12,'Location',...
                           'best','FontWeight','bold','Box','off');
        grid on;
        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
    case 31
        figure;hold all;
        plot(A,B);
        title('Spectrum of the Received eletrical signal',...
                                'FontSize',16,'FontWeight','bold');
        xlabel('Frequency [Hz]','FontSize',14);
        ylabel('Amplitude [dB]','FontSize',14);
        legend({'Ix'},'FontSize',12,'Location','best',...
                                  'FontWeight','bold','Box','off');
        grid on;
        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
    case 32
        figure;hold all;grid on;
        plot(A,B);
        title(['Comparison between transmited and received '...
                      'signal'],'FontSize',16,'FontWeight','bold');
        xlabel('Time [s]','FontSize',14);
        ylabel('Amplitude [u]','FontSize',14);
        legend({'Rx'},'FontSize',12,'Location','best',...
                                  'FontWeight','bold','Box','off');
        grid on;
        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
    case 33
        Olho(A,B,C,1);
        title('Eye diagram of the received signal','FontSize',16,...
                                                      'FontWeight','bold');
        xlabel('Time [s]','FontSize',14);
        ylabel('Amplitude [u]','FontSize',14);
        legend({'Ix'},'FontSize',12,'Location','best','FontWeight',...
                                                       'bold','Box','off');
        grid on;
    case 34
        figure;
        plot(A,B,A,C);
        title('Comparing received symbol and synchronism symbol'...
                               ,'FontSize',16,'FontWeight','bold');
        xlabel('Time [s]','FontSize',14);
        ylabel('Amplitude [u]','FontSize',14);
        legend({'Rec-Sym','Sycn'},'FontSize',12,'Location',...
                           'best','FontWeight','bold','Box','off');
        grid on;
        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
    case 35
        figure;plot(A,B,A,C);
        title(['Comparing received symbol and synchronism ' ...
                  'symbol'],'FontSize',16,'FontWeight','bold');
        xlabel('Time [s]','FontSize',14);
        ylabel('Amplitude [u]','FontSize',14);
        legend({'Rec-Sym','Sycn'},'FontSize',12,'Location',...
                       'best','FontWeight','bold','Box','off');
        grid on;
        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
    case 36
        Olho(A,B,C,1);
        title('Eye diagram of the received signal','FontSize',...
                                           16,'FontWeight','bold');
        xlabel('Time [s]','FontSize',14);
        ylabel('Amplitude [u]','FontSize',14);
        legend({'Ix'},'FontSize',12,'Location','best',...
                                  'FontWeight','bold','Box','off');
        grid on;
    case 37
        figure;
        bar(A,B);
        title('Histogram of the income signal - Decission Level'...
                               ,'FontSize',16,'FontWeight','bold');
        xlabel('Countble Interval [u]','FontSize',14);
        ylabel('Number of Occurrence [u]','FontSize',14);
        grid on;
        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
    case 38
        figure;hold all;
        plot(C,A,'LineWidth',2);
        plot(C(B),A(B),'xr','LineWidth',2);
        title('Profile of the Eye Diagram','FontSize',16,...
                                              'FontWeight','bold');
        xlabel('Countble Interval [u]','FontSize',14);
        ylabel('Number of Occurrence [u]','FontSize',14);
        legend({'Hist','PeakFound'},'FontSize',12,'Location',...
                           'best','FontWeight','bold','Box','off');
        grid on;
        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
    case 39
        PossX = linspace(1,A,A);                                 %Auxiliar variable for ploting
        figure;
        bar(PossX,B);
        title('Histogram of the income signal - Symmetry',...
                                'FontSize',16,'FontWeight','bold');
        xlabel('Countble Interval [u]','FontSize',14);
        ylabel('Number of Occurrence [u]','FontSize',14);
        grid on;
        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
    case 40
        figure;hold all;
        plot(A,'LineWidth',2);
        plot(B,A(B),'xr','LineWidth',2);
        title('Profile of the Eye Diagram','FontSize',16,'FontWeight',...
                                                                   'bold');
        xlabel('Countble Interval [u]','FontSize',14);
        ylabel('Number of Occurrence [u]','FontSize',14);
        legend({'Hist','PeakFound'},'FontSize',12,'Location','best',...
                                          'FontWeight','bold','Box','off');
        grid on;
        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
    case 41
        figure;
        plot(linspace(0,A,2*B),C,linspace(0,A,2*B),D);
        title('Comparison between transmited and received Data',...
                                'FontSize',16,'FontWeight','bold');
        xlabel('Time [s]','FontSize',14);
        ylabel('Amplitude [u]','FontSize',14);
        legend({'TxData','RxData'},'FontSize',12,'Location',...
                           'best','FontWeight','bold','Box','off');
        grid on;
        set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
    case 42
        figure;hold all;
        plot(A,B);
        title('Spectrum of the Received eletrical signal',...
                                'FontSize',16,'FontWeight','bold');
        xlabel('Frequency [Hz]','FontSize',14);
        ylabel('Amplitude [dB]','FontSize',14);
        legend({'Ix'},'FontSize',12,'Location','best',...
                                  'FontWeight','bold','Box','off');
        grid on;
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
    case 43
        figure;hold all;
        plot(A,B,A,C);
        title(['Comparison between transmited and received ' ...
                      'signal'],'FontSize',16,'FontWeight','bold');
        xlabel('Time [s]','FontSize',14);
        ylabel('Amplitude [u]','FontSize',14);
        legend({'Tx','Rx'},'FontSize',12,'Location','best',...
                                  'FontWeight','bold','Box','off');
        grid on;
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
    case 44
        Olho(A,B,C,1);
        title('Eye diagram of the received signal','FontSize',...
                                           16,'FontWeight','bold');
        xlabel('Time [s]','FontSize',14);
        ylabel('Amplitude [u]','FontSize',14);
        legend({'Ix'},'FontSize',12,'Location','best',...   
                                  'FontWeight','bold','Box','off');
        grid on;
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
    case 45
        figure;
        plot(A,B,A,C);
        title('Comparing received symbol and synchronism symbol'...
                               ,'FontSize',16,'FontWeight','bold');
        xlabel('Time [s]','FontSize',14);
        ylabel('Amplitude [u]','FontSize',14);
        legend({'Rec-Sym','Sycn'},'FontSize',12,'Location',...
                           'best','FontWeight','bold','Box','off');
        grid on;
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
    case 46
        figure;hold all;
        plot(A,B);
        title(['Received Signal signal'],'FontSize',16,'FontWeight',...
                                                                   'bold');
        xlabel('Time [s]','FontSize',14);
        ylabel('Amplitude [u]','FontSize',14);
        legend({'Rx'},'FontSize',12,'Location','best','FontWeight',...
                                                       'bold','Box','off');
        grid on;
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
    case 47
        Olho(A,B,C,1);
        title('Eye diagram of the synchronized signal',...
                                'FontSize',16,'FontWeight','bold');
        xlabel('Time [s]','FontSize',14);
        ylabel('Amplitude [u]','FontSize',14);
        legend({'Ix'},'FontSize',12,'Location','best',...
                                  'FontWeight','bold','Box','off');
        grid on;
    case 48
        PossX = linspace(1,A,A);                                 %Auxiliar variable for ploting
        figure;grid on;
        bar(PossX,B);
        title('Histogram of the income signal - Symmetry',...
                                'FontSize',16,'FontWeight','bold');
        xlabel('Countble Interval [u]','FontSize',14);
        ylabel('Number of Occurrence [u]','FontSize',14);
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
    case 49
        figure;hold all;
        plot(A,'LineWidth',2);
        plot(B,A(B),'xr','LineWidth',2);
        title('Profile of the Eye Diagram','FontSize',16,...
                                              'FontWeight','bold');
        xlabel('Countble Interval [u]','FontSize',14);
        ylabel('Number of Occurrence [u]','FontSize',14);
        legend({'Hist','SymPoi'},'FontSize',12,'Location',...
                           'best','FontWeight','bold','Box','off');
        grid on;
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
    case 50
        figure;
        plot(A,B,A,C);
        title('Comparison between transmited and received Data',...
                                'FontSize',16,'FontWeight','bold');
        xlabel('Time [s]','FontSize',14);
        ylabel('Amplitude [u]','FontSize',14);
        legend({'TxData','RxData'},'FontSize',12,'Location',...
                           'best','FontWeight','bold','Box','off');
        grid on;
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
    case 51
        figure;
        hold all;
        for kk=1:size(A,1)
            plot(B,20*log10(abs(fftshift(fft(A(kk,:)./length(A(kk,:)))))));
        end
        [~,locmin] = findpeaks(20*log10(abs(fftshift(fft(A(1,:)./length(A(1,:)))))),'SortStr','descend','NPeaks',1);
        [~,locmax] = findpeaks(20*log10(abs(fftshift(fft(A(end,:)./length(A(end,:)))))),'SortStr','descend','NPeaks',1);
%         axis([0.9*B(locmin) 1.1*B(locmin) -200 0]);
        title('Splited COMB','FontSize',16,'FontWeight','bold');
        xlabel('Frequency [Hz]','FontSize',14);
        ylabel('Magnitude [dB]','FontSize',14);
        grid on;
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
    otherwise
end