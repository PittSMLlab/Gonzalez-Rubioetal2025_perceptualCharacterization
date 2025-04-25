function [fh,f2]=accrtPlot(trialData, flagWF)

if nargin<2 || isempty(flagWF)
    flagWF=0;  %default
end

if flagWF == 1
    %Perturbation size as a proportion of mean speed between the legs
    %(Weber Fraction)
    trialData.pertSize=trialData.pertSize./1050;
    trialData.prevSize=trialData.prevSize./1050;
end

%% Creating a (modified) subject ID field:

aux=trialData.subID;
trialData.ID=categorical(aux);

%% Remove null and no response trials (for accurate counting of DOF in stats)

trialData=trialData(~trialData.noResponse,:);

%% Get probe sizes

B=findgroups(trialData.pertSize); %pertSize>0 means vR>vL
pp=unique(trialData.pertSize);

%% First figure: Logistic fits, accuracy and reaction times

fh=figure('Units','pixels','InnerPosition',[100 100 3*300 1*300]);
sSize=25;
[cmap,unsignedMap]=probeColorMap(23);

%% Backward selection: test the effect of laterality, task learning, and habituation using Mixed effects model

X=trialData;
X.pertSign=sign(X.pertSize);
X.ID = categorical(X.ID);
X.blockNo = categorical(X.blockNo);

% Model
frml='leftResponse~1+pertSize+prevSize+pertSize:blockNo+pertSize:pertSign+(pertSize+prevSize+pertSize:blockNo+pertSize:pertSign|ID)';
mm0=fitglme(X,frml,'Distribution','binomial','Link','logit','FitMethod','Laplace')

% % Drop the biggest non significant term, blockNo
% frml='leftResponse~1+pertSize+prevSize+pertSize:pertSign+(pertSize+prevSize+pertSize:pertSign|ID)';
% mm0=fitglme(X,frml,'Distribution','binomial','Link','logit','FitMethod','Laplace');

% % Drop the biggest non significant term, prevSize
% frml='leftResponse~1+pertSize+pertSize:pertSign+(pertSize+pertSize:pertSign|ID)';
% mm0=fitglme(X,frml,'Distribution','binomial','Link','logit','FitMethod','Laplace');

% Drop the biggest non significant term, pertSign
frml='leftResponse~1+pertSize+(pertSize|ID)';
mm0=fitglme(X,frml,'Distribution','binomial','Link','logit','FitMethod','Laplace')


%% First sub-plot Figure 1: proportion of left choices as function of probe size/weber fraction

subplot(2,3,[1,4])
hold on
set(gca,'Colormap',cmap);
S=splitapply(@(x) sum(x==-1)/sum(~isnan(x)),trialData.initialResponse,B); %Not counting NR responses
E=splitapply(@(x) 1.96*(nanstd(x==-1)/sqrt(sum(~isnan(x)))),trialData.initialResponse,B); %95% CI
% E=splitapply(@(x) nanstd(x==-1)/sqrt(sum(~isnan(x))),trialData.initialResponse,B); %Not counting NR responses
ss=scatter(pp,S,sSize,pp,'filled','MarkerEdgeColor','w');
grid on;
ylabel('proportion of left choices')
if flagWF==1
    axis([-0.310 0.310 0 1.01])
else
    axis([-310 310 0 1.01])
end

%Add fits:
hold on
errorbar(pp,S,E,'k','LineStyle','none')
set(gca,'Colormap',unsignedMap);

Nsubs=unique(trialData.ID);
if mm0.Formula.FELinearFormula.HasIntercept  %mm0.Formula.HasIntercept
    hi='';
else
    hi='-1';
end

if flagWF==1
    xx=[-0.310:0.001:0.310];
else
    xx=[-310:310];
end

%Add group fit:
betas=[]; 
y_all=[];
for i=1:length(Nsubs)
    mm{i}=fitglm(X(X.ID==Nsubs(i),:),[char(mm0.Formula.FELinearFormula) hi],'Distribution','binomial');
    mm{i}.plotPartialDependence('pertSize')
    betas=[betas; mm{i}.Coefficients.Estimate(1) mm{i}.Coefficients.Estimate(2)]; 
    y_all=[y_all; 1./(1+exp(-1*(betas(i,1)+betas(i,2)*xx)))]; % Individual participants' curves
end

plot(xx,nanmean(y_all), 'k', 'LineWidth',2); % Average individual curves

ll=findobj(gca,'Type','Line');
set(ll(2:end),'Color',.7*ones(1,3));
uistack(ll(1:end-1),'bottom')
uistack(ss,'top')
ylabel('proportion of "left" choices')
xlabel('R slower       same        L slower')

if flagWF==1
    set(gca,'XLim',[-0.310 0.310])
else
    set(gca,'XLim',[-310 310])
end


%% Second sub-plot Figure 1: Accuracy as a function of probe size/weber fraction

subplot(2,3,[2,5])
set(gca,'Colormap',unsignedMap);
hold on;
grid on;
trialData.correctResponses=double(trialData.correctResponses);
trialData.correctResponses(isnan(trialData.initialResponse))=nan;
B2=findgroups(abs(trialData.pertSize));
S2=splitapply(@(x) nansum(x)/sum(~isnan(x)),trialData.correctResponses,B2); %Not counting NR responses
S2(S2==0)=NaN;
ap=sort(unique(abs(pp)));
ss=scatter(ap,S2,sSize,ap,'filled','MarkerEdgeColor','w');
E2=splitapply(@(x) 1.96*(nanstd(x==1)/sqrt(sum(~isnan(x)))),trialData.correctResponses,B2); %95% CI
grid on
errorbar(ap,S2,E2,'k','LineStyle','none')
ylabel('accuracy')
if flagWF==1
    axis([0 0.31 .5 1]) % For dataset 2 axis([0 0.36 .5 1])
else
    axis([0 310 .5 1]) % For dataset 2 axis([0 0.36 .5 1])
end

%Run binomial tests and BH on accuracy results:
clear pval h
for i=2:length(S2)
    Ntrials(i)=sum(~isnan(trialData.initialResponse(abs(trialData.pertSize)==ap(i))));
    correctTrials=Ntrials(i)*S2(i);
    pval(i-1)=binocdf(correctTrials-1,Ntrials(i),.5,'upper');
end
disp('Significance testing on accuracy:')
for i=1:length(pval)
    disp(['\Delta v=' num2str(ap(i+1),3) ', p=' num2str(pval(i),4) ', hits=' num2str(Ntrials(i+1)*S2(i+1),2) '/' num2str(Ntrials(i+1),2)])
end
h=BenjaminiHochberg(pval, .05); %Two-stage BKY procedure
if all(h)
    [mp,mpi]=max(pval);
    disp(['All test were significant with BKY. Largest p=' num2str(mp) ' for probe size=' num2str(ap(mpi+1)) 'mm/s'])
else
    disp('Some non-sig tests!')
    disp(num2str(ap(logical([0 h]))))
end

if flagWF==1
    xx=[0:0.001:0.310];
else
    xx=[0:310];
end

%Add group fit:
b=mm0.Coefficients.Estimate;
if ~mm0.Formula.FELinearFormula.HasIntercept
    y=1-1./(1+exp(b(1)*xx));
else
    y=.5*(1-1./(1+exp(b(1)+b(2)*xx))+1./(1+exp(b(1)+b(2)*-xx)));
end
thm=find(y>.75,1,'first');

%Add individual fits:
for i=1:length(Nsubs)
    XX=X(X.ID==Nsubs(i),:);
    b=mm{i}.Coefficients.Estimate;
    if ~mm{i}.Formula.HasIntercept
        y=1-1./(1+exp(b(1)*xx));
    else
        y=.5*(1-1./(1+exp(b(1)+b(2)*xx))+1./(1+exp(b(1)+b(2)*-xx)));
    end

    if isempty(find(y>.75,1,'first'));
        th(i) = nan;
    else
        th(i)=find(y>.75,1,'first');
    end
end

ll=findobj(gca,'Type','Line','LineWidth',1');
uistack(ll,'bottom')
uistack(ss,'top')

yline(0.75,'--')

%% Third sub-plot Figure 1: information on reaction times vs probe size/weber fraction

rtPlots(trialData,[],flagWF);

%% Figure 3

sSize=40;

k=1.0986;
X=trialData; %Excludes no response trials already
X.acc=X.correctResponses+.5*trialData.nullTrials;
b0=cellfun(@(x) x.Coefficients.Estimate(1),mm);
b1=cellfun(@(x) x.Coefficients.Estimate(2),mm);
[~,idx]=sort(b1); % Sorting the bar plots depending on the individual's JND
CI=cell2mat(cellfun(@(x) x.coefCI,mm,'UniformOutput',false));
nSubs=length(unique(X.ID));
b0CI=reshape(CI(1,:),2,nSubs);
b1CI=reshape(CI(2,:),2,nSubs);

% Start plotting
f2=figure('Units','pixels','InnerPosition',[100 100 3*300 1*300]);
subplot(1,3,1) %Bias
hold on

biases=-b0(idx)./b1(idx); %PSE
biasCI=-b0CI(:,idx)./b1(:,idx);
gb=mean(biases);
gSE=1.96*std(biases)./sqrt(length(biases)); % SE of the biases
gb_mixed=-mm0.Coefficients.Estimate(1)/mm0.Coefficients.Estimate(2);
gCI_mixed=-mm0.coefCI./mm0.Coefficients.Estimate(2);
if flagWF==1
    ylabel('PSE [-\beta_0/\beta_1]')
else
    ylabel('PSE [-\beta_0/\beta_1] (mm/s)')
end

%plot:
bb=bar(1:nSubs,biases,'FaceColor','flat','EdgeColor','none');
bar(nSubs+2,gb,'FaceColor',.2*ones(1,3),'EdgeColor','none'); 
errorbar(nSubs+2,gb,gSE,'k','LineStyle','none','LineWidth',1);
yline(gb,'--k');
ee=errorbar(1:nSubs,biases,biases-biasCI(2,:),biasCI(1,:)-biases,'k','LineStyle','none','LineWidth',1);

title('bias')
xlabel('subject')

ticks=cellstr(string([1:nSubs]));
ticks{end+1}='Group';
set(gca,'XTick',[1:nSubs,nSubs+2],'XTickLabel',ticks)

disp('-------------Stats PSE:------------')
disp(['Group=' num2str(gb) ', mean=' num2str(mean(biases))  ', std=' num2str(std(biases)) ', range=[' num2str(min(biases)) ',' num2str(max(biases)) ']']);


subplot(1,3,2)
hold on
%Alt: plot 1.1/beta_1
slopes=k./b1(idx);
slopeCI=k./b1CI(:,idx);
bs=bar(slopes,'FaceColor','flat','EdgeColor','none');
ee=errorbar(1:nSubs,slopes,slopes-slopeCI(2,:),slopeCI(1,:)-slopes,'k','LineStyle','none','LineWidth',1); % 1:nSubs
gb=mean(slopes);
bar(nSubs+2,gb,'FaceColor',.2*ones(1,3),'EdgeColor','none');
yline(gb,'--k');
gSE=1.96*std(slopes)./sqrt(length(slopes));
errorbar(nSubs+2,gb,gSE,'k','LineStyle','none','LineWidth',1);

title('probe size effect')
if flagWF==1
    ylabel('JND [1.1\beta_1^{-1}]')
else
    ylabel('JND [1.1\beta_1^{-1}] (mm/s)')
end
set(gca,'XTick',[1:nSubs,nSubs+2],'XTickLabel',ticks)
xlabel('subject')

% Color map

ex2=[0.8500 0.3250 0.0980]; 
ex1=[0 0.4470 0.7410];
mid1=0.5*ex1; 
mid2=0.5*ex2;
gamma=0.5;
N=round(length(Nsubs)/2);
mapRB=[flipud((mid1 + (ex1-mid1).*([1:N]'/(N)).^gamma)); (mid2+ (ex2-mid2).*([1:N]'/(N)).^gamma)];

bs.CData=mapRB([1:length(Nsubs)],:);
bb.CData=mapRB([1:length(Nsubs)],:);

colormap(mapRB);
cb = colorbar('Location','north', 'TickLabels', {num2str(round(min(slopes))), num2str((round(max(slopes(idx(13)))))), num2str(round(max(slopes)))});

subplot(1,3,3) %Avg. accuracy vs. avg. RT
hold on
G=findgroups(X.ID);
acc=splitapply(@nanmean,X.acc,G); acc=acc(idx);
eacc=splitapply(@(x) 1.96*(nanstd(x)/sqrt(numel(x))),X.acc,G); eacc=eacc(idx);
RT=splitapply(@nanmean,X.reactionTime,G); RT=RT(idx);
eRT=splitapply(@(x) 1.96*(nanstd(x)/sqrt(numel(x))),X.reactionTime,G); eRT=eRT(idx);
errorbar(acc,RT,eRT,'k','LineStyle','none')
ee=errorbar(acc,RT,eacc,'k','Horizontal','LineStyle','none','DisplayName','ste');
ss=scatter(acc,RT,sSize,mapRB(1:length(Nsubs),:),'filled','MarkerEdgeColor','w','DisplayName','Individual subject data');
rta=nanmean(X.reactionTime);
acca=nanmean(X.acc);
scatter(acca, rta,sSize,.2*ones(1,3),'filled','MarkerEdgeColor','w');
eRTa=nanstd(X.reactionTime)/sqrt(sum(~isnan(X.reactionTime)));
eacca=nanstd(X.acc)/sqrt(sum(~isnan(X.reactionTime)));
errorbar(acca,rta,eacca,'k','Horizontal','LineStyle','none','DisplayName','ste');
errorbar(acca,rta,eRTa,'k','LineStyle','none');
xlabel('accuracy');
ylabel('mean RT (s)');
title('indiv. accuracy vs. RT');

% Correlation Analysis
[rho,pval]=corr(acc,RT); disp(['Accuracy vs. Mean RT: rho=' num2str(rho) ', p=' num2str(pval)]);
[rho,pval]=corr(slopes',acc); disp(['JND vs. Mean Accuracy: rho=' num2str(rho) ', p=' num2str(pval)]);
[rho,pval]=corr(slopes',RT); disp(['JND vs. Mean RT: rho=' num2str(rho) ', p=' num2str(pval)]);

% extendedPanelWidth(f2,.1);

hold off;

disp('-------------Stats JND:------------')
disp(['Group=' num2str(gb) ', mean=' num2str(mean(slopes))  ', std=' num2str(std(slopes)) ', range=[' num2str(min(slopes)) ',' num2str(max(slopes)) ']']);

end



