function plotfigure6

basePath = fileparts(fileparts(mfilename('fullpath')));
dataPath = fullfile(basePath,'figures','dependencies','data');
addpath(fullfile(basePath,'data_processing'));
addpath(fullfile(basePath,'figures','dependencies','functions'));


fig = figure('color','w','units','centimeters');
fig.Position(3) = 8.5; fig.Position(4) = 12.75;
fig.Position(2) = 3;

plotdata(dataPath);
plotWTCSI(dataPath);
plotcsidata(dataPath);

labelpanel(0.0,0.95,'a');
labelpanel(0.51,0.95,'b');
labelpanel(0.0,0.525,'c');
labelpanel(0.49,0.525,'d');
labelpanel(0.0,0.25,'e');
labelpanel(0.49,0.25,'f');

function plotdata(dataPath)
	figPos = get(gcf,'position');
	fy =figPos(4);

	Nav15 = load(fullfile(dataPath,'Nav1.5e','GV_Curves.mat'));
	Nav14 = load(fullfile(dataPath,'Nav1.4','GV_Curves.mat'));
	Nav16 = load(fullfile(dataPath,'Nav1.6','GV_Curves.mat'));


	clrs(1,:) = [0,0,0];
	clrs(2,:) = [1,0,0];
	clrs(3,:) = [0,0.3,1];

	ax(1) = axes('units','centimeters','FontSize',7); hold on;
	ax(1).Position = [0.6,0.25+fy-4.6,3.4,2.7];
		V=mean(Nav15.A_V(:,1:30));
		G = Nav15.ActFit(:,2).*Nav15.A_I(:,1:30)./(V-Nav15.ActFit(:,1));
		outliers = find(or(Nav15.ActFit(:,1)<45,Nav15.ActFit(:,1)>75));
		outliers = [outliers;find(isoutlier(Nav15.ActFit(:,3)))];
		G(outliers,:)=[];
		G = G./max(G,[],2);
		Nav15.Av50 = Nav15.ActFit(:,4);
		Nav15.Av50(outliers) = [];
		plot(V,G','.','Color',clrs(1,:),'MarkerSize',2);
		hold on;
		FT=FitBoltzmanCurve2(V,mean(G),-50,-10);
		V = linspace(V(1),V(end),1e2);
		plot(V,FT(V),'color',clrs(1,:),'LineWidth',1);

		V=mean(Nav14.A_V(:,1:30));
		G = Nav14.ActFit(:,2).*Nav14.A_I(:,1:30)./(V-Nav14.ActFit(:,1));
		outliers = find(or(Nav14.ActFit(:,1)<45,Nav14.ActFit(:,1)>75));
		outliers = [outliers;find(isoutlier(Nav14.ActFit(:,3)))];
		G(outliers,:)=[];
		G = G./max(G,[],2);
		Nav14.Av50 = Nav14.ActFit(:,4);
		plot(V,G','.','Color',clrs(2,:),'MarkerSize',2);
		hold on;
		FT=FitBoltzmanCurve2(V,mean(G),-50,-10);
		V = linspace(V(1),V(end),1e2);
		plot(V,FT(V),'color',clrs(2,:),'LineWidth',1);

		V=mean(Nav16.A_V(:,1:30));
		G = Nav16.ActFit(:,2).*Nav16.A_I(:,1:30)./(V-Nav16.ActFit(:,1));
		outliers = find(or(Nav16.ActFit(:,1)<45,Nav16.ActFit(:,1)>75));
		outliers = [outliers;find(isoutlier(Nav16.ActFit(:,3)))];
		G(outliers,:)=[];
		G = G./max(G,[],2);
		Nav16.Av50 = Nav16.ActFit(:,4);
		Nav16.Av50(outliers) = [];
		plot(V,G','.','Color',clrs(3,:),'MarkerSize',2);
		hold on;
		FT=FitBoltzmanCurve2(V,mean(G),-50,-10);
		V = linspace(V(1),V(end),1e2);
		plot(V,FT(V),'color',clrs(3,:),'LineWidth',1);

		xlim([-60,35])
		ylim([0,1.05]);
		set(gca,'tickdir','out');
		set(gca,'LineWidth',0.75);
		set(gca,'xtick',[-120:30:30]);
		% set(get(gca,'yaxis'),'visible','off');
		XL = xlabel('Test pulse (mV)');

		h(1)=plot(nan,nan,'.','color',clrs(1,:),'MarkerSize',10);
		h(2)=plot(nan,nan,'.','color',clrs(2,:),'MarkerSize',10);
		h(3)=plot(nan,nan,'.','color',clrs(3,:),'MarkerSize',10);
		L=legend(h,{'Nav1.5e','Nav1.4','Nav1.6'},'box','off',...
					'fontsize',7,'position',[0.28,0.68,0.18,0.07]);
		L.ItemTokenSize=[10,5];

	ax(2) = axes('units','centimeters','FontSize',7); hold on;
	ax(2).Position = [4.8,0.25+fy-4.6,3.4,2.7];
		V = mean(Nav15.I_V);
		I = Nav15.I_I.*Nav15.InactFit(:,2);
		I(max(I,[],2)>1.5,:) = [];
		Nav15.Iv50 = Nav15.InactFit(:,4);
		plot(V,I','.','Color',clrs(1,:),'MarkerSize',2);
		hold on;
		FT=FitBoltzmanCurve2(V,mean(I),-50,10);
		V = linspace(V(1),V(end),1e2);
		plot(V,FT(V),'color',clrs(1,:),'LineWidth',1);

		V = mean(Nav14.I_V);
		I = Nav14.I_I.*Nav14.InactFit(:,2);
		Nav14.Iv50 = Nav14.InactFit(:,4);
		plot(V,I','.','Color',clrs(2,:),'MarkerSize',2);
		FT=FitBoltzmanCurve2(V,mean(I),-50,10);
		V = linspace(V(1),V(end),1e2);
		plot(V,FT(V),'color',clrs(2,:),'LineWidth',1);

		V = mean(Nav16.I_V);
		I = Nav16.I_I.*Nav16.InactFit(:,2);
		Nav16.Iv50 = Nav16.InactFit(:,4);
		plot(V,I','.','Color',clrs(3,:),'MarkerSize',2);
		FT=FitBoltzmanCurve2(V,mean(I),-50,10);
		V = linspace(V(1),V(end),1e2);
		plot(V,FT(V),'color',clrs(3,:),'LineWidth',1);

		xlim([-120,-20])
		ylim([0,1.05]);
		set(gca,'tickdir','out');
		set(gca,'LineWidth',0.75);
		set(gca,'xtick',[-120:30:30]);
		% set(get(gca,'yaxis'),'visible','off');
		XL = xlabel('Conditioning pulse (mV)');



	a14 = load(fullfile(dataPath,'Nav1.4','20181205c3','activation.mat'));
	ax(1) = axes('units','centimeters'); hold on;
	ax(1).Position = [0.7,0.25+fy-1.8,1.4,1.75];
		axis off;
			V = mean(a14.Voltage(:,a14.Epochs(4):a14.Epochs(4)+500)');
			a14.Current = activationleakcorrection(a14.Voltage,a14.Current,a14.Epochs);
			X = a14.Current(:,a14.Epochs(4):a14.Epochs(4)+500)';
			[M,I] = max(abs(X));
			sig = sign(X(I+[0:length(I)-1]*501));
			FT = FitBoltzman(V,M.*sig/max(M(:)),-50,-10,60,1);
			X = X(:,find(V<FT.ERev-15));
			mX = max(abs(X(:)));
		plot(X,'LineWidth',0.25,'color',clrs(2,:));
		xlim([1,500]); ylim([-1.1,0.2]*mX);
		fr = 2e3/mX;
		line([300,300],[-0.9,-0.9+fr]*mX,'LineWidth',1,'color',clrs(1,:));
		text(340,mX*(-0.9+fr/2),['2 nA'],'FontSize',6,'HorizontalAlignment','left');

	a14 = load(fullfile(dataPath,'Nav1.4','20181205c3','inactivation.mat'));
	ax(1) = axes('units','centimeters'); hold on;
	ax(1).Position = [5,0.25+fy-1.8,1.4,1.75];
		axis off;
			a14.Current = inactivationleakcorrection(a14.Voltage,a14.Current,a14.Epochs);
			X = a14.Current(:,a14.Epochs(5)+21:a14.Epochs(5)+520)';
		plot(X,'LineWidth',0.25,'color',clrs(2,:));
		mX = max(abs(X(:)));
		xlim([1,500]); ylim([-1.1,0.2]*mX);
		fr = 2e3/mX;
		line([300,300],[-0.9,-0.9+fr]*mX,'LineWidth',1,'color',clrs(1,:));
		text(340,mX*(-0.9+fr/2),['2 nA'],'FontSize',6,'HorizontalAlignment','left');
		
	a16 = load(fullfile(dataPath,'Nav1.6','20170419c4','activation.mat'));
	ax(1) = axes('units','centimeters'); hold on;
	ax(1).Position = [2.5,0.25+fy-1.8,1.4,1.75];
		axis off;
			V = mean(a16.Voltage(:,a16.Epochs(4):a16.Epochs(4)+500)');
			a16.Current = activationleakcorrection(a16.Voltage,a16.Current,a16.Epochs);
			X = a16.Current(:,a16.Epochs(4):a16.Epochs(4)+500)';
			[M,I] = max(abs(X));
			sig = sign(X(I+[0:length(I)-1]*501));
			FT = FitBoltzman(V,M.*sig/max(M(:)),-50,-10,60,1);
			X = X(:,find(V<FT.ERev-15));
			mX = max(abs(X(:)));
		plot(X,'LineWidth',0.25,'color',clrs(3,:));
		xlim([1,500]); ylim([-1.1,0.2]*mX);
		fr = 1e3/mX;
		line([300,300],[-0.9,-0.9+fr]*mX,'LineWidth',1,'color',clrs(1,:));
		text(340,mX*(-0.9+fr/2),['1 nA'],'FontSize',6,'HorizontalAlignment','left');
		plot([300,300,500],[-0.9+fr,-0.9,-0.9]*mX,'LineWidth',1,'Color',clrs(1,:));
		text(400,-0.9*mX,['2 ms'],'FontSize',6,'HorizontalAlignment','center',...
				'VerticalAlignment','top');

	a16 = load(fullfile(dataPath,'Nav1.6','20170419c4','inactivation.mat'));
	ax(1) = axes('units','centimeters'); hold on;
	ax(1).Position = [6.8,0.25+fy-1.8,1.4,1.75];
		axis off;
			a16.Current = inactivationleakcorrection(a16.Voltage,a16.Current,a16.Epochs);
			X = a16.Current(:,a16.Epochs(5)+21:a16.Epochs(5)+520)';
		plot(X,'LineWidth',0.25,'color',clrs(3,:));
		mX = max(abs(X(:)));
		xlim([1,500]); ylim([-1.1,0.2]*mX);
		fr = 1e3/mX;
		plot([300,300,500],[-0.9+fr,-0.9,-0.9]*mX,'LineWidth',1,'Color',clrs(1,:));
		text(340,mX*(-0.9+fr/2),['1 nA'],'FontSize',6,'HorizontalAlignment','left');
		text(400,-0.9*mX,['2 ms'],'FontSize',6,'HorizontalAlignment','center',...
				'VerticalAlignment','top');

function plotWTCSI(dataPath)
	Folder = fullfile(dataPath,'Nav1.5e','20170418c3');

	load(fullfile(Folder,'activation.mat')); % Load inactivation data
	V = mean(Voltage(:,Epochs(4):Epochs(4)+500)');
	Current = activationleakcorrection(Voltage,Current,Epochs);
	X = Current(:,Epochs(4)+30:Epochs(4)+530)';
	[~,idxMax] = max(abs(X));
	for j = 1:size(X,2)
		Imax(j) = X(idxMax(j),j);
	end
	[~,idx] = min(Imax);
	Vlo = V(idx)+10;
	Vhi = 65;
	idcsLinear = and(V>Vlo,V<Vhi);
	FT = fitlm(V(idcsLinear),Imax(idcsLinear));
	ERev = -FT.Coefficients.Estimate(1)/FT.Coefficients.Estimate(2);

	load(fullfile(Folder,'inactivation.mat')); % Load inactivation data
	Current = inactivationleakcorrection(Voltage,Current,Epochs); % Correct for leak current
	V = mean(Voltage(:,Epochs(4):Epochs(4)+500),2); % Get the voltage of pre pulse
	preCurrent = Current(:,Epochs(4)+42:Epochs(5)-50); % Get the current in response to the pre-pulse (42 used to get rid of capactive transient)
	preCurrent = preCurrent-median(preCurrent,2); % Responses are baselined to return to 0
	preG = preCurrent ./(V-ERev); % Correct for reversal potential to estimate open probability of channels
	% peakG = max(abs(preG(end,:))); % Get peak conductance
	% O = preCurrent./peakG; % Normalize to get a probablility
	OSI = sum(preG,2);%/max(sum(preG,2)); % OSI is proportional to the total occupancy of the open state
	FT = FitBoltzmanCurve(V,OSI,-10,-5,max(OSI));
	OSI=OSI/FT.Gmx;


	I_V = mean(Voltage(:,Epochs(4):Epochs(5)),2); % Get the voltages of the test-pulses
	I_V = round(I_V/5)*5; % Round 
	postCurrent = Current(:,Epochs(5)+42:Epochs(6)); % Get currents in response to test-pulse
	I_I = min(postCurrent,[],2); % get peak current minimum (they should all be negative)
	Inact = 1 - I_I/min(I_I); % Get the fraction of unavailable channels.

	CSI = Inact - OSI; % CSI are channels is amount of unavailble channels that can't be accounted 
												% for by the pre-pulse current


	ax = axes('units','centimeters');
	ax.Position = [0.2,4,4,3];
	CM = parula(15);
	CM = CM(1:13,:);
	for i = 2:2:24
		plot(2201:3200,(postCurrent(i+8,1:1000)')/max(abs(preCurrent(:))),'color',CM(i/2,:)); hold on;
	end

	for i = 8:2:24
		plot(cumsum(preG(8+i,1:1800))/FT.Gmx,'color',CM(i/2,:));
		plot(preCurrent(8+i,1:1800)'/max(abs(preCurrent(:))),'color',CM(i/2,:));
	end
	axis off;
	line([0,3200],[0,0],'color','k');
	text(1900,0.02,'...','fontsize',7);
	ylim([-1,1.05])
	xlim([0,3200]);
	scatter(2201,0.1,10,'vk','filled')
	text(2300,0.11,'-10 mV','HorizontalAlignment','left','VerticalAlignment','middle','FontSize',6)

	colormap(CM);
	C = colorbar('location','south');
	set(gca,'CLim',[-120,0])
	C.Units = 'centimeters';
	C.Limits = [-120,0];
	C.Position = [ax.Position(1)+0.71,ax.Position(2)+0.85,1.6,0.3];
	C.FontSize = 6;
	C.Ticks = [-120,-60,0];
	C.TickLabels = {'-120','-60','0 mV'};
	C.TickDirection = 'out';
	C.LineWidth = 0.25;
	C.Label.String = 'Conditioning pulse';
	C.Label.Position(2) = -1.2;

	drawbrace([1825,1.05],[1825,0.01],0.75,'color','k');

	str = '$$\int P$$';
	text(1975,0.5,str,'Interpreter','latex','FontSize',6)
	text(2250,0.5,'_{Open}','FontSize',6)
	str = '$$\propto$$';
	text(2600,0.5,str,'Interpreter','latex','FontSize',6)
	text(2800,0.5,'OSI','FontSize',6);

	axis off;

	text(1600,1.15,['CSI = SSI ' char(8211) ' OSI'],'FontSize',7,'HorizontalAlignment','center','VerticalAlignment','bottom')

function plotcsidata(dataPath)
	clrs(1,:) = [0,0,0];
	clrs(2,:) = [1,0,0];
	clrs(3,:) = [0,0.3,1];

	Nav14 = load(fullfile(dataPath,'Nav1.4','CSI.mat'));
	Nav15 = load(fullfile(dataPath,'Nav1.5e','CSI.mat'));
	Nav16 = load(fullfile(dataPath,'Nav1.6','CSI.mat'));

	ax = axes('units','centimeters'); hold on;
	ax.Position = [4.8,4.6,3.4,2.7];
		CM = parula(15);
		CM = CM(1:13,:);
		errorbar(Nav15.V,mean(Nav15.CSI),std(Nav15.CSI)/sqrt(size(Nav15.CSI,1)),'LineWidth',0.75,'Color','r')
		hold on;
		errorbar(Nav15.V,mean(Nav15.OSI),std(Nav15.OSI)/sqrt(size(Nav15.CSI,1)),'LineWidth',0.75,'Color',CM(end,:))
		errorbar(Nav15.V,mean(Nav15.Inact),std(Nav15.Inact)/sqrt(size(Nav15.CSI,1)),'LineWidth',0.75,'color',CM(1,:));
		box off;
		xlim([-120,0])
		xticks([-120:60:0]);
		ylim([-0.1,1.05]);
		set(gca,'TickDir','out')
		set(gca,'FontSize',7);
		set(gca,'LineWidth',0.75);
		xlabel('Conditioning Pulse (mV)');

		text(-85,0.2,'SSI','color',CM(1,:),'FontSize',7);
		text(-20,0.6,'OSI','color',CM(end,:),'FontSize',7);
		text(-60,0.75,'CSI','color','r','FontSize',7);


	ax = axes('Units','centimeters');
	ax.Position = [0.92,0.85,2.9,2.75];
		FT14 = FitBoltzmanCurve4(repmat(Nav14.V',[size(Nav14.CSI,1),1]),Nav14.CSI);
		FT15 = FitBoltzmanCurve4(repmat(Nav15.V',[size(Nav15.CSI,1),1]),Nav15.CSI);
		FT16 = FitBoltzmanCurve4(repmat(Nav16.V',[size(Nav16.CSI,1),1]),Nav16.CSI);
		V = linspace(-120,0,1e2);

		plot(Nav15.V,Nav15.CSI,'.','Color',clrs(1,:),'MarkerSize',2); hold on;
		plot(Nav14.V,Nav14.CSI,'.','Color',clrs(2,:),'MarkerSize',2);
		plot(Nav16.V,Nav16.CSI,'.','Color',clrs(3,:),'MarkerSize',2);
		plot(V,FT15(V),'Color',clrs(1,:),'LineWidth',1);
		plot(V,FT14(V),'Color',clrs(2,:),'LineWidth',1);
		plot(V,FT16(V),'Color',clrs(3,:),'LineWidth',1);
		xlim([-120,0])
		xticks([-120:60:0]);
		ylim([-0.1,1.05]);
		box off;
		set(gca,'TickDir','out')
		set(gca,'FontSize',7);
		set(gca,'LineWidth',0.75);
		xlabel('Conditioning Pulse (mV)');
		ylabel('% CSI');


	for i =1 :size(Nav14.CSI,1)
		FT = FitBoltzmanCurve4(Nav14.V',Nav14.CSI(i,:));
		total14(i) = sum(FT(-200:0.1:200)*0.1);
	end
	for i =1 :size(Nav15.CSI,1)
		FT = FitBoltzmanCurve4(Nav15.V',Nav15.CSI(i,:));
		total15(i) = sum(FT(-200:0.1:200)*0.1);
	end
	for i =1 :size(Nav16.CSI,1)
		FT = FitBoltzmanCurve4(Nav16.V',Nav16.CSI(i,:));
		total16(i) = sum(FT(-200:0.1:200)*0.1);
	end
	total14 = total14(:);
	total15 = total15(:);
	total16 = total16(:);

	Y = [total15;total14;total16];
	G = [ones(size(total15));2*ones(size(total14));3*ones(size(total16))];
	[p,tbl,stats] = anova1(Y,G,'off');
	k = length(unique(G))
	df = length(G)-k
	c = multcompare(stats,'display','off')
	c(2:3,end)

	ax = axes('Units','centimeters');
	ax.Position = [4.8,0.85,1.7,2.75];
		plot(1+(mod(1:length(total15),5)-1.5)/15,total15,'square','color',clrs(1,:), ...
			'MarkerSize',1.5,'MarkerFaceColor',clrs(1,:)); hold on;
		plot(2+(mod(1:length(total14),5)-1.5)/15,total14,'square','color',clrs(2,:), ...
			'MarkerSize',1.5,'MarkerFaceColor',clrs(2,:));
		plot(3+(mod(1:length(total16),5)-1.5)/15,total16,'square','color',clrs(3,:), ...
			'MarkerSize',1.5,'MarkerFaceColor',clrs(3,:));
		xticks([1,2,3]);
		xticklabels({'Nav1.5e','Nav1.4','Nav1.6'});
		ylim([25,67]);
		xlim([0.5,3.5]);
		box off;
		set(gca,'TickDir','out')
		set(gca,'FontSize',7);
		set(gca,'LineWidth',0.75);
		ylabel('Total CSI (mV)');
		ax = get(gca,'xaxis');
		
		ax.TickLabelRotation = -30;
		line([1,3],[64,64],'color','k')
		text(2,65,'***','FontSize',8,'VerticalAlignment','middle','HorizontalALignment','center')
		line([1,2],[60,60],'color','k')
		text(1.5,61,'**','FontSize',8,'VerticalAlignment','middle','HorizontalALignment','center')
		line([2,3],[53,53],'color','k')
		text(2.5,54,'***','FontSize',8,'VerticalAlignment','middle','HorizontalALignment','center')

		% xticklabels({'1.5e','1.4','1.6'});
		% xticklabels({'Nav1.5e','Nav1.4','Nav1.6'});

	return;
	ax = axes('Units','centimeters');
	ax.Position = [6.6,0.85,1.7,2.75];

		bar(1,mean(total15),'FaceColor',clrs(1,:),'LineStyle','none','BarWidth',0.3)
		hold on;
		plot([1,1],[-stderror(total15) stderror(total15)] + mean(total15),'color',clrs(1,:),'LineWidth',1); 

		bar(2,mean(total14),'FaceColor',clrs(2,:),'LineStyle','none','BarWidth',0.3)
		plot([2,2],[-stderror(total14) stderror(total14)] + mean(total14),'color',clrs(2,:),'LineWidth',1);


		bar(3,mean(total16),'FaceColor',clrs(3,:),'LineStyle','none','BarWidth',0.3)
		plot([3,3],[-stderror(total16) stderror(total16)] + mean(total16),'color',clrs(3,:),'LineWidth',1);

		xticks([1,2,3]);
		xticklabels({'1.5e','1.4','1.6'});
		ylim([25,67]);
		xlim([0.5,3.5]);
		box off;
		set(gca,'TickDir','out')
		set(gca,'FontSize',7);
		set(gca,'LineWidth',0.75);
		set(get(gca,'yaxis'),'visible','off');

		T = text(-0.9,14.3,'Nav Isoform','FontSize',7.7);

		line([1,3],[60,60],'color','k')
		text(2,61,'***','FontSize',8,'VerticalAlignment','middle','HorizontalALignment','center')
		line([1,2],[55,55],'color','k')
		text(1.5,56,'**','FontSize',8,'VerticalAlignment','middle','HorizontalALignment','center')
		line([2,3],[50,50],'color','k')
		text(2.5,51,'***','FontSize',8,'VerticalAlignment','middle','HorizontalALignment','center')
