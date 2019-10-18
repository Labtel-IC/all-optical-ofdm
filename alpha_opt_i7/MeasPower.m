function [ThisPower,ThisPowerdBm] = MeasPower(Ein,t,InfLim,UpeLim)
%%
    if nargin<2
        t = linspace(0,length(Ein),length(Ein));
        InfLim = t(1);
        UpeLim = t(end);
    elseif nargin<3
        InfLim = t(1);
        UpeLim = t(end);
    end

    ThisPower = sum(abs(Ein(find(t==InfLim):find(t==UpeLim))).^2)/(...
                             length(Ein(find(t==InfLim):find(t==UpeLim))));
    if nargout == 2
        ThisPowerdBm = 30 + 10*log10(ThisPower);
    end

end