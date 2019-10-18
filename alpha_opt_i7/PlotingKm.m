%%
NoiFig = figure(666);
hold all;

if ~exist('FiberLengthKm','var')
    FiberLengthKm = VetSnrIni:VetSnrPass:VerSnrEnd;
    if UseOfdmElect
        VerThisMod = OfdMod;
        ThisM      = max(DmtMve); 
    else
        VerThisMod = Modulation;
    end
    switch VerThisMod
        case 'qam'
            LegendThis = [{'FEC-LIMIT: 3.8e-3'} {'OFDM-QAM'}];
        case 'dpsk'
            LegendThis = [{'FEC-LIMIT: 3.8e-3'} {'OFDM-DPSK'}];
        case '4PAM'
            LegendThis = [{'FEC-LIMIT: 3.8e-3'} {'PAM-4'}];
        case 'DQPSK'
            LegendThis = [{'FEC-LIMIT: 3.8e-3'} {'DQPSK'}];
        case 'DPSK'
            LegendThis = [{'FEC-LIMIT: 3.8e-3'} {'DPSK'}];
        otherwise
            LegendThis = [{'FEC-LIMIT: 3.8e-3'} {'OOK'}];
    end
end
    
FecLimit = linspace(3.8e-3,3.8e-3,length(FiberLengthKm));
ThisSimuBer(ThisPlotCont)= mean(b(:,ceil(size(b,2)/2)));
TheoColor      = [0 0 0];
PlotColor      = [1 0.6 0];
MarkerColor    = PlotColor;
ThisLineWidth  = 2;
ThisMarkerSize = 6;
PlotMarker     = 'd';

plot(FiberLengthKm,FecLimit,'-','color',TheoColor,'LineWidth',ThisLineWidth);
plot(FiberLengthKm(1:length(ThisSimuBer)),ThisSimuBer,PlotMarker,'color',...
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
axis([VetSnrIni VerSnrEnd 1e-7 1e0]);
legend(LegendThis,'FontName','CMR10','FontSize',20,'Box','off','Location','Best');

set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
drawnow;
ThisPlotCont = ThisPlotCont + 1;
a=1;