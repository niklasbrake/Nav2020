function plotfigure4

warning('off','curvefit:prepareFittingData:removingNaNAndInf')
basePath = fileparts(fileparts(mfilename('fullpath')));
dataPath = fullfile(basePath,'figures','dependencies','data');
addpath(fullfile(basePath,'data_processing'));
addpath(fullfile(basePath,'figures','dependencies','functions'));

load(fullfile(dataPath,'VCF','VCF_SummaryData.mat'));
F(:,3) = F(:,3)/mean(F(vF>-50,3));
D(1) = load(fullfile(dataPath,'VCF','D1_20180510c2.mat'));
baseline(1) = 1.8849;
D(2) = load(fullfile(dataPath,'VCF','D2_20180512c1.mat'));
baseline(2) = 1.8307;
D(3) = load(fullfile(dataPath,'VCF','D3_20180518c4.mat'));
baseline(3) = 1.8307;
D(4) = load(fullfile(dataPath,'VCF','D4_20180514c5.mat'));
baseline(4) = 2.0698;

WT = load(fullfile(dataPath,'ephys_SummaryData.mat'));

fig = figure('Units','Centimeters','Color','w');
% fig.Position(3) = 8.5;
fig.Position(3) = 8.9;
fig.Position(4) = 15.25;
fig.Position(2) = fig.Position(2)-4;
figh = fig.Position(4);


ax(1) = axes('units','Centimeters','fontsize',7,'LineWidth',1,'box','off','tickdir','out');
ax(1).Position = [4.7,figh-2.25,3,2];
ax(2) = axes('units','Centimeters','fontsize',7,'LineWidth',1,'box','off','tickdir','out');
ax(2).Position = [4.7,figh-6+0.75,3,2];
ax(3) = axes('units','Centimeters','fontsize',7,'LineWidth',1,'box','off','tickdir','out');
ax(3).Position = [4.7,figh-9+0.75,3,2];
ax(4) = axes('units','Centimeters','fontsize',7,'LineWidth',1,'box','off','tickdir','out');
ax(4).Position = [4.7,figh-12+0.75,3,2];

for i = 1:4
	axes(ax(i)); hold on;
	FB = FitBoltzmanCurve2(WT.Nav15e.activation(:,1),WT.Nav15e.activation(:,2),-50,-5);
	EB(2) = errorbar(WT.Nav15e.activation(:,1),WT.Nav15e.activation(:,2),WT.Nav15e.activation(:,3),'LineWidth',0.25,'LineStyle','none','Marker','v','MarkerSize',3,'color',ones(1,3)*0.4+0.6*[1,0,0]); hold on;
	EB(2).MarkerFaceColor = EB(2).Color;
	plot(vF(1):vF(end),FB(vF(1):vF(end)),'Color',EB(2).Color);
	FB = FitBoltzmanCurve2(vF,F(:,i),-50,-5);
	EB(1) = errorbar(vF,F(:,i),sF(:,i),'LineWidth',0.25,'LineStyle','none','Marker','^','MarkerSize',3,'color',[0,1,0]); hold on;
	EB(1).MarkerFaceColor = EB(1).Color;
	plot(vF(1):vF(end),FB(vF(1):vF(end)),'Color',EB(1).Color);
	FB = FitBoltzmanCurve2(WT.Nav15e.inactivation(3:end,1),WT.Nav15e.inactivation(3:end,2),-100,5);
	EB(3) = errorbar(WT.Nav15e.inactivation(3:end,1),WT.Nav15e.inactivation(3:end,2),WT.Nav15e.inactivation(3:end,3),'LineWidth',0.25,'LineStyle','none','Marker','o','MarkerSize',3,'color',ones(1,3)*0.4+0.6*[0,0.3,1]); hold on;
	EB(3).MarkerFaceColor = EB(3).Color;
	plot(vF(1):vF(end),FB(vF(1):vF(end)),'Color',EB(3).Color);
	xlim([-200 50]);
	ylim([0,1.2]);
	if(i==4)
		xlabel('Voltage (mV)')
	end
	if(i==3)
		FB = FitBoltzmanCurve2(vF,F(:,4),-50,-5);
		FB2 = @(V) FB.v50/FB.k.*exp((FB.v50-V)/FB.k)./(1+exp((FB.v50-V)/FB.k)).^2;
		FB = FitBoltzmanCurve2(vF,F(:,3),-50,-5);
		bmp = FB2(vF(1):vF(end));
		bmp = bmp(:)/max(bmp(:))*(max(F(:,3)-1));
		plot(vF(1):vF(end),FB(vF(1):vF(end))+bmp,'--k');
		% FB = FitBoltzmanCurve2(WT.minusD4.inactivation(3:end,1),WT.minusD4.inactivation(3:end,2),-100,5);
		% EB(4) = errorbar(WT.minusD4.inactivation(3:end,1),WT.minusD4.inactivation(3:end,2),WT.minusD4.inactivation(3:end,3),'LineWidth',0.25,'LineStyle','none','Marker','^','MarkerSize',3); hold on;
		% EB(4).Color = EB(3).Color*0.6+(ones(1,3)-0.6);
		% EB(4).MarkerFaceColor = EB(4).Color;
		% plot(vF(1):vF(end),FB(vF(1):vF(end)),'Color',EB(4).Color);
	end

end
L = legend(EB(1:3),{'FV','GV (WT)','SSI (WT)'});
L.Position = [0.8308 0.2682 0.1760 0.0649];
L.ItemTokenSize = [10,5];
L.Box = 'off';

C1 = EB(1).Color;
CM = [linspace(0,C1(1),5)' linspace(0,C1(2),5)' linspace(0,C1(3),5)'];

axs(1) = axes('units','Centimeters');
axs(1).Position = [1.15,figh-2.5,2.75,2.5]; axs(1).YLim = [-0.025 0.01];
axs(2) = axes('units','Centimeters');
axs(2).Position = [1.15,figh-5.5,2.75,2.5]; axs(2).YLim = [-0.03 0.01];
axs(3) = axes('units','Centimeters');
axs(3).Position = [1.15,figh-8.5,2.75,2.5]; axs(3).YLim = [-0.015 0.02];
axs(4) = axes('units','Centimeters');
axs(4).Position = [1.15,figh-11.5,2.75,2.5]; axs(4).YLim = [-0.035 0.01];

T = 1.7e4-1.2e4:1.7e4+1.2e4;
tLine = 0.8e4;
for i = 1:4
	axes(axs(i)); hold on;
	yl = axs(i).YLim;
	for j = 1:5
		plot(T,D(i).VCF(T,j),'LineWidth',1,'Color',CM(j,:));
	end
	line([T(1) T(end)],[0,0],'Color','k','LineWidth',1,'LineStyle','--');
	xlim([T(1) T(end)]);
	line([tLine tLine+0.5e4],[yl(1),yl(1)],'Color','k','LineWidth',1);
			text(tLine+0.25e4,yl(1), '50 ms', 'FontSize',6,'HorizontalAlignment','center',...
				'VerticalAlignment','top');
	if(sum(i==[1,3]))
		line([tLine tLine],[yl(1),yl(1)+0.003*baseline(i)],'Color','k','LineWidth',1);
		text(tLine-5e2,yl(1)+(0.003*baseline(i))/2, '0.3%', 'FontSize',6,'HorizontalAlignment','right');
	else
		line([tLine tLine],[yl(1),yl(1)+0.004*baseline(i)],'Color','k','LineWidth',1);
		text(tLine-5e2,yl(1)+(0.004*baseline(i))/2, '0.4%', 'FontSize',6,'HorizontalAlignment','right');
		if(i==4)
			
		end
	end
			
	set(get(gca,'xaxis'),'visible','off');
	set(get(gca,'yaxis'),'visible','off');
end


A=labelpanel(0.015,0.93,'a');
A=labelpanel(0.015,0.735,'b');
A=labelpanel(0.015,0.54,'c');
A=labelpanel(0.015,0.35,'d');


A=labelpanel(0.06,0.93,'DI*','Cell'); A.FontWeight = 'normal';
A=labelpanel(0.06,0.735,'DII*','Cell'); A.FontWeight = 'normal';
A=labelpanel(0.06,0.54,'DIII*','Cell'); A.FontWeight = 'normal';
A=labelpanel(0.06,0.35,'DIV*','Cell'); A.FontWeight = 'normal';


clrs = lines(6);
clrs=clrs([2,3,4,6],:);
mrks = {'^','v','d','o'};

iD4 = WT.minusD4.inactivation;
[V,iV1,iV2] = intersect(vF,iD4(:,1));
FF = F(iV1,3);
FF(FF>1.05) = nan;
SSI = iD4(iV2,2);

ax(5) = axes('units','Centimeters','fontsize',7,'LineWidth',1,'box','off','tickdir','out'); hold on;
ax(5).Position = [4.7,figh-15+0.7,2,2];
	h(1) = errorbar(SSI,F(iV1,1),sF(iV1,1),sF(iV1,1),iD4(iV2,3),iD4(iV2,3),mrks{1},'Color',clrs(1,:),'LineWidth',0.25,'MarkerSize',3); h(1).MarkerFaceColor = h(1).Color;
	hold on;
	h(2) = errorbar(SSI,F(iV1,2),sF(iV1,2),sF(iV1,2),iD4(iV2,3),iD4(iV2,3),mrks{2},'Color',clrs(2,:),'LineWidth',0.25,'MarkerSize',3); h(2).MarkerFaceColor = h(2).Color;
	h(3) = errorbar(SSI,F(iV1,3),sF(iV1,3),sF(iV1,3),iD4(iV2,3),iD4(iV2,3),mrks{3},'Color',clrs(3,:),'LineWidth',0.25,'MarkerSize',3); h(3).MarkerFaceColor = h(3).Color;
	FT = fitlm(SSI,FF); hold on;
	plot(SSI,FT.predict(SSI),'--k','LineWidth',1);
	xlim([-0.05,1])
	ylim([-0.1,1.2])
	box off
	set(gca,'TickDir','out')
	xlabel('SSI (DIV-CN)')
	axis square;


ax(6) = axes('units','Centimeters','fontsize',7,'LineWidth',1,'box','off','tickdir','out'); hold on;
ax(6).Position = [1.15,figh-15+0.7,2.75,2];
	FB = FitBoltzmanCurve2(WT.minusD4.inactivation(:,1),WT.minusD4.inactivation(:,2),-100,5);
	EB(4) = errorbar(WT.minusD4.inactivation(:,1),WT.minusD4.inactivation(:,2),WT.minusD4.inactivation(:,3),'LineWidth',0.25,'LineStyle','none','Marker','o','MarkerSize',3); hold on;
	EB(4).Color = EB(3).Color*0.6+(ones(1,3)-0.6);
	EB(4).MarkerFaceColor = EB(4).Color;
	plot(vF(1):vF(end),FB(vF(1):vF(end)),'Color',EB(4).Color);
	for i = 1:3
		FB = FitBoltzmanCurve2(vF,F(:,i),-50,-5);
		EB(i) = errorbar(vF,F(:,i),sF(:,i),'LineWidth',0.25,'LineStyle','none','Marker',mrks{i},'MarkerSize',3,'color',clrs(i,:)); hold on;
		EB(i).MarkerFaceColor = EB(i).Color;
		plot(vF(1):vF(end),FB(vF(1):vF(end)),'Color',EB(i).Color);
		end
	xlim([-200 20]);
	ylim([-0.1,1.2]);
	xlabel('Voltage (mV)')
	% ylabel('Normalized FV');

L = legend(EB(1:4),{'FV (DI*)','FV (DII*)','FV (DIII*)','SSI (DIV-CN)'});
L.ItemTokenSize = [10,5];
L.Box = 'off';
L.Position = [0.8 0.1 0.1760 0.0649];


labelpanel(0.06,0.16,'e');