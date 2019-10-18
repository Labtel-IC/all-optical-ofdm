%%
NoiFig = figure(666);
hold all;

if ~exist('CompEbn0','var')
    CompEbn0 = VetSnrIni:VetSnrPass:VerSnrEnd;
    if UseOfdmElect
        VerThisMod = OfdMod;
        ThisM      = max(DmtMve); 
    else
        VerThisMod = Modulation;
    end
    switch VerThisMod
        case 'qam'
            BerTheo    = berawgn(CompEbn0, 'qam',  ThisM);
            LegendThis = [{[num2str(ThisM) 'QAM-Theo']} {[num2str(ThisM) 'QAM-Simu']}];
        case 'dpsk'
            BerTheo    = berawgn(CompEbn0, 'dpsk', ThisM);
            LegendThis = [{[num2str(ThisM) 'DPSK-Theo']} {[num2str(ThisM) 'DPSK-Simu']}];
        case '4PAM'
            ThisM      = 4;
            BerTheo    = berawgn(CompEbn0, 'pam',  ThisM);
            LegendThis = [{[num2str(ThisM) 'PAM-Theo']} {[num2str(ThisM) 'PAM-Simu']}];
        case 'DQPSK'
            ThisM      = 4;
            BerTheo    = berawgn(CompEbn0, 'dpsk', ThisM);
            LegendThis = [{[num2str(ThisM) 'DQPSK-Theo']} {[num2str(ThisM) 'DQPSK-Simu']}];
        case 'DPSK'
            ThisM      = 2;
            BerTheo    = berawgn(CompEbn0, 'dpsk', ThisM);
            LegendThis = [{'DPSK-Theo'} {'DPSK-Simu'}];
        otherwise
            ThisM      = 2;
            BerTheo    = berawgn(CompEbn0, 'pam',  ThisM);
            LegendThis = [{'OOK-Theo'} {'OOK-Simu'}];
    end
end

b = b(berpos1);
ThisSimuBer(ThisPlotCont)= mean(b(:,ceil(size(b,2)/2)));
TheoColor      = [0 0 0];
PlotColor      = [1 0.6 0];
MarkerColor    = PlotColor;
ThisLineWidth  = 2;
ThisMarkerSize = 6;
PlotMarker     = 'o';

plot(CompEbn0,BerTheo,'-','color',TheoColor,'LineWidth',ThisLineWidth);
plot(CompEbn0(1:length(ThisSimuBer)),ThisSimuBer,PlotMarker,'color',...
PlotColor,'LineWidth',ThisLineWidth,'MarkerFaceColor',MarkerColor,...
'MarkerSize',ThisMarkerSize);

xlabel('EbN0','FontSize',20);
ylabel('BER','FontSize',20);
ThisFig = gca;
ThisFig.FontSize = 20;
ThisFig.FontName = 'CMR10';
ThisFig.Box = 'on';
ThisFig.LineWidth = 2;
ThisFig.YScale = 'log';
ThisFig.XGrid = 'on';
ThisFig.YGrid = 'on';
ThisFig.XMinorGrid = 'off';
ThisFig.YMinorGrid = 'off';
axis([0 length(CompEbn0)-1 1e-7 1e0]);
legend(LegendThis,'FontName','CMR10','FontSize',20,'Box','off','Location','Best');

set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
drawnow;
ThisPlotCont = ThisPlotCont + 1;
a=1;