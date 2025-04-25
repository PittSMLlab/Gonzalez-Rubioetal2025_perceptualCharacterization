function [fh]=rtPlots(trialData,goodOnly, flagWF)
if nargin<2 || isempty(goodOnly)
    goodOnly=0;
end

if nargin<3 || isempty(flagWF)
    flagWF=0;  %default
end

%% Figure Features
[cmap,unsignedMap]=probeColorMap(23);
sSize=25;

%% reaction times vs. pert size
rt=trialData.reactionTime;
fun=@nanmean;

B2=findgroups(trialData.pertSize); %pertSize>0 means vR>vL
pp2=unique(trialData.pertSize);

subplot(2,3,[3,6])
RT=splitapply(fun,rt,B2); 
eRT=splitapply(@(x) 1.96*(nanstd(x)/sqrt(sum(~isnan(x)))),rt,B2); 
ss=scatter(pp2,RT,sSize,pp2,'filled','MarkerEdgeColor','w');

set(gca,'Colormap',unsignedMap)
hold on
grid on
ylabel('mean RT (s)') 
if flagWF==1
    xlabel('|\Delta V / V|')
else
    xlabel('|\Delta V| (mm/s)')
end

%Adding vR>vL and vL<vR overlapping
if flagWF==1
    set(gca, 'XLim', [-0.31 0.31]) 
else
    set(gca, 'XLim', [-310 310])
end

set(gca, 'YLim', [1.5 5])
e=errorbar(pp2,RT,eRT,'k');
e.LineStyle='none';
uistack(ss,'top')
end

