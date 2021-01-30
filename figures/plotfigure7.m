function plotfigure7

	basePath = fileparts(fileparts(mfilename('fullpath')));
	dataPath = fullfile(basePath,'figures','dependencies','data');
	addpath(fullfile(basePath,'data_processing'));
	addpath(fullfile(basePath,'figures','dependencies','functions'));

	fig = figure('color','w','units','centimeters');
	fig.Position(3) = 8.5; fig.Position(4) = 12.75;
	fig.Position(2) = 3;

	WT = load(fullfile(dataPath,'Nav1.5e','GV_Curves.mat'));
	B1 = load(fullfile(dataPath,'Nav1.5e+B1','GV_Curves.mat'));
	B3 = load(fullfile(dataPath,'Nav1.5e+B3','GV_Curves.mat'));
	csiWT = load(fullfile(dataPath,'Nav1.5e','CSI.mat'));
	csiB1 = load(fullfile(dataPath,'Nav1.5e+B1','CSI.mat'));
	csiB3 = load(fullfile(dataPath,'Nav1.5e+B3','CSI.mat'));
	ax=plotbeta(WT,B1,B3,csiWT,csiB1,csiB3,2.45);
	axes(ax(1)); text(-55,1,'Nav1.5e','FontSize',8)
	axes(ax(3));
	ylim([37,62]);
	line([1,2],[59,59],'color','k'); text(1.5,59.5,'***','FontSize',8,'VerticalAlignment','middle','HorizontalALignment','center')
	line([1,3],[61.5,61.5],'color','k'); text(2,62,'***','FontSize',8,'VerticalAlignment','middle','HorizontalALignment','center')

	WT = load(fullfile(dataPath,'Nav1.4','GV_Curves.mat'));
	B1 = load(fullfile(dataPath,'Nav1.4+B1','GV_Curves.mat'));
	B3 = load(fullfile(dataPath,'Nav1.4+B3','GV_Curves.mat'));
	csiWT = load(fullfile(dataPath,'Nav1.4','CSI.mat'));
	csiB1 = load(fullfile(dataPath,'Nav1.4+B1','CSI.mat'));
	csiB3 = load(fullfile(dataPath,'Nav1.4+B3','CSI.mat'));	
	[ax,T]=plotbeta(WT,B1,B3,csiWT,csiB1,csiB3,5.5);
	axes(ax(1)); text(-55,1,'Nav1.4','FontSize',8)
	ax(2).XLim = [-110,-30];
	ax(2).XTick = [-100:30:-40];
	ax(2).XTickLabels = [-100:30:-40];
	T.Position(1) = -31.3;
	axes(ax(3));
	ax(3).XTickLabels{1} = [char(945) '1.4'];
	ylim([29,54]);

	line([1,2],[49.5,49.5],'color','k'); text(1.5,50,'**','FontSize',8,'VerticalAlignment','middle','HorizontalALignment','center')
	line([1,3],[52,52],'color','k'); text(2,52.5,'***','FontSize',8,'VerticalAlignment','middle','HorizontalALignment','center')


	WT = load(fullfile(dataPath,'Nav1.6','GV_Curves.mat'));
	B1 = load(fullfile(dataPath,'Nav1.6+B1','GV_Curves.mat'));
	B3 = load(fullfile(dataPath,'Nav1.6+B3','GV_Curves.mat'));
	csiWT = load(fullfile(dataPath,'Nav1.6','CSI.mat'));
	csiB1 = load(fullfile(dataPath,'Nav1.6+B1','CSI.mat'));
	csiB3 = load(fullfile(dataPath,'Nav1.6+B3','CSI.mat'));
	[ax,T] = plotbeta(WT,B1,B3,csiWT,csiB1,csiB3,8.6);
	axes(ax(1)); text(-55,1,'Nav1.6','FontSize',8)
	ax(2).XLim = [-100,-20];
	ax(2).XTick = [-90:30:-30];
	ax(2).XTickLabels = [-90:30:-30];
	T.Position(1) = -21.3;
	ax(3).XTickLabels{1} = [char(945) '1.6'];
	ax(3).YLim = [17,42];

	clrs(1,:) = [0,0,0];
	clrs(2,:) = [1,0,0];
	clrs(3,:) = [0,0.3,1];
	h(1) = plot(nan,nan,'o','color',clrs(2,:),'MarkerSize',4,'MarkerFaceColor',clrs(2,:)); hold on;
	h(2) = plot(nan,nan,'square','color',clrs(3,:),'MarkerSize',4,'MarkerFaceColor',clrs(3,:)); hold on;
	L=legend(h,{[' + ' char(946) '1'],[' + ' char(946) '3']},'box','off',...
						'fontsize',7,'position',[0.19,0.83,0.18,0.07]);
	L.ItemTokenSize=[10,5];

	plotlitdata(12.2);

	labelpanel(0.01,0.95,'a');
	labelpanel(0.72,0.95,'b');
	labelpanel(0.01,0.715,'c');
	labelpanel(0.72,0.715,'d');
	labelpanel(0.01,0.47,'e');
	labelpanel(0.72,0.47,'f');
	labelpanel(0.01,0.23,'g');

function [ax,T] = plotbeta(WT,B1,B3,csiWT,csiB1,csiB3,offset)
	fig = gcf;
	fy = fig.Position(4);

	clrs(1,:) = [0,0,0];
	clrs(2,:) = [1,0,0];
	clrs(3,:) = [0,0.3,1];

	ax(1) = axes('units','centimeters','FontSize',7); hold on;
	ax(1).Position = [0.6,0.25+fy-offset,2.3,2];
	ax(2) = axes('units','centimeters','FontSize',7); hold on;
	ax(2).Position = [3.55,0.25+fy-offset,2.3,2];
		V=mean(WT.A_V(:,1:30));
		G = WT.ActFit(:,2).*WT.A_I(:,1:30)./(V-WT.ActFit(:,1));
		outliers = find(or(WT.ActFit(:,1)<45,WT.ActFit(:,1)>75));
		outliers = [outliers;find(isoutlier(WT.ActFit(:,3)))];
		G(outliers,:)=[];
		G = G./max(G,[],2);
		WT.Av50 = WT.ActFit(:,4);
		WT.Av50(outliers) = [];
		axes(ax(1));
		h(1)=errorbar(V,mean(G),stderror(G),'LineWidth',0.5,'Marker','v','color',clrs(1,:),...
				'MarkerFaceColor',clrs(1,:),'LineStyle','none','MarkerSize',3);
		hold on;
		FT=FitBoltzmanCurve2(V,mean(G),-50,-10);
		plot(V,FT(V),'color',clrs(1,:));
		V = mean(WT.I_V);
		I = WT.I_I.*WT.InactFit(:,2);
		I(max(I,[],2)>1.5,:) = [];
		WT.Iv50 = WT.InactFit(:,4);
		axes(ax(2));
		errorbar(V,mean(I),stderror(I),'LineWidth',0.5,'Marker','v','color',clrs(1,:),...
				'MarkerFaceColor',clrs(1,:),'LineStyle','none','MarkerSize',3);
		FT=FitBoltzmanCurve2(V,mean(I),-50,10);
		plot(V,FT(V),'color',clrs(1,:));

		V=mean(B1.A_V(:,1:30));
		G = B1.A_I(:,1:30)./(V-B1.ActFit(:,1));
		G = G.*B1.ActFit(:,2);
		G = G./max(G,[],2);
		outliers = find(or(B1.ActFit(:,1)<45,B1.ActFit(:,1)>75));
		outliers = [outliers;find(isoutlier(B1.ActFit(:,3)))];
		G(outliers,:)=[];
		B1.Av50 = B1.ActFit(:,4);
		B1.Av50(outliers,:) = [];
		axes(ax(1));
		h(2)=errorbar(V,mean(G),stderror(G),'LineWidth',0.5,'Marker','o','color',clrs(2,:),...
				'MarkerFaceColor',clrs(2,:),'LineStyle','none','MarkerSize',3);
		FT=FitBoltzmanCurve2(V,mean(G),-50,-10);
		plot(V,FT(V),'color',clrs(2,:));
		V = mean(B1.I_V);
		I = B1.I_I.*B1.InactFit(:,2);
		I(max(I,[],2)>1.5,:) = [];
		B1.Iv50 = B1.InactFit(:,4);
		axes(ax(2));
		errorbar(V,mean(I),stderror(I),'LineWidth',0.5,'Marker','o','color',clrs(2,:),...
				'MarkerFaceColor',clrs(2,:),'LineStyle','none','MarkerSize',3);
		FT=FitBoltzmanCurve2(V,mean(I),-50,10);
		plot(V,FT(V),'color',clrs(2,:));

		V=mean(B3.A_V(:,1:30));
		G = B3.ActFit(:,2).*B3.A_I(:,1:30)./(V-B3.ActFit(:,1));
		G = G./max(G,[],2);
		outliers = find(or(B3.ActFit(:,1)<45,B3.ActFit(:,1)>75));
		outliers = [outliers;find(isoutlier(B3.ActFit(:,3)))];
		G(outliers,:)=[];
		B3.Av50 = B3.ActFit(:,4);
		B3.Av50(outliers,:) = [];
		axes(ax(1));
		h(3)=errorbar(V,mean(G),stderror(G),'LineWidth',0.5,'Marker','square','color',clrs(3,:),...
				'MarkerFaceColor',clrs(3,:),'LineStyle','none','MarkerSize',3);
		FT=FitBoltzmanCurve2(V,mean(G),-50,-10);
		plot(V,FT(V),'color',clrs(3,:));
		V = mean(B3.I_V);
		I = B3.I_I.*B3.InactFit(:,2);
		outliers = 1;
		B3.InactFit(outliers,:) = [];
		I(max(I,[],2)>1.5,:) = [];
		B3.Iv50 = B3.ActFit(:,4);
		axes(ax(2));
		errorbar(V,mean(I),stderror(I),'LineWidth',0.5,'Marker','square','color',clrs(3,:),...
				'MarkerFaceColor',clrs(3,:),'LineStyle','none','MarkerSize',3);
		FT=FitBoltzmanCurve2(V,mean(I),-50,10);
		plot(V,FT(V),'color',clrs(3,:));

		axes(ax(1));
		xlim([-60,35])
		ylim([0,1.05]);
		set(gca,'tickdir','out');
		set(gca,'LineWidth',0.75);
		set(gca,'xtick',[-50:25:25]);
		text(32,-0.156,'mV','FontSize',7);

		axes(ax(2));
		xlim([-120,-40])
		ylim([0,1.05]);
		set(gca,'tickdir','out');
		set(gca,'LineWidth',0.75);
		set(gca,'xtick',[-110:30:-50]);
		T=text(-41.3,-0.156,'mV','FontSize',7);

	ax(3) = axes('units','centimeters','FontSize',7); hold on;
	ax(3).Position = [6.9,0.25+fy-offset,1.4,2];

		for i =1 :size(csiWT.CSI,1)
			FT = FitBoltzmanCurve4(csiWT.V',csiWT.CSI(i,:));
			totalWT(i) = sum(FT(-200:0.1:200)*0.1);
		end
		for i =1 :size(csiB1.CSI,1)
			FT = FitBoltzmanCurve4(csiB1.V',csiB1.CSI(i,:));
			totalB1(i) = sum(FT(-200:0.1:200)*0.1);
		end
		for i =1 :size(csiB3.CSI,1)
			FT = FitBoltzmanCurve4(csiB3.V',csiB3.CSI(i,:));
			totalB3(i) = sum(FT(-200:0.1:200)*0.1);
		end
		totalWT = sort(totalWT(:));
		totalB1 = sort(totalB1(:));
		totalB3 = sort(totalB3(:));

		Y = [totalWT;totalB1;totalB3];
		G = [ones(size(totalWT));2*ones(size(totalB1));3*ones(size(totalB3))];
		[p,tbl,stats] = anova1(Y,G,'off');
		c = multcompare(stats,'display','off')
		plot(1+(mod(1:length(totalWT),5)-1.5)/15,totalWT,'.','color',clrs(1,:)); hold on;
		plot(2+(mod(1:length(totalB1),5)-1.5)/15,totalB1,'.','color',clrs(2,:));
		plot(3+(mod(1:length(totalB3),5)-1.5)/15,totalB3,'.','color',clrs(3,:));

		xlim([0.5,3.5]);
		xticks([1,2,3]);
		xticklabels({[char(945) '1.5e'],['+ ' char(946) '1'],['+ ' char(946) '3']});
		set(gca,'FontSize',7);
		box off;
		set(gca,'TickDir','out');
		ylabel('Total CSI (mV)')
		xax = get(gca,'xaxis'); xax.TickLabelRotation = -45;

function plotlitdata(offset)

	fig = gcf;
	fy = fig.Position(4);

	T = [1,-55.5,1.3,-54.2,0.8,-18.3,1.5,-19.0,1.8,1;
	1,-35.0,1.0,-43.0,1.0,-17.0,1.0,-18.0,1.0,1;
	2,-51.7,2.6,-58.1,1.3,-18.9,1.4,-20.1,3.1,1;
	2,-42.0,2.0,-52.0,1,-18.0,2,-21.0,2.0,1;
	3,-64.9,1.5,-59.9,0,-25.5,1.6,-25.5,0,0;
	4,-54.0,1.0,-60.0,1.0,-27.0,1.0,-26.0,1.0,1;
	4,-54.0,0.4,-62.0,0.1,-25.0,0.5,-33.0,0.3,1;
	4,-66.1,0.6,-66.5,0.5,-22.8,0.9,-23.9,1.2,0;
	4,-74.2,1.9,-65.3,1.6,-24.8,1.2,-24.7,2.3,0;
	4,-67.7,1.3,-66.0,2.2,-19.2,2.0,-19.1,4.7,0;
	4,-64.0,0.7,-59.9,0.5,-8.0,0.7,-6.1,0.56,0;
	5,-76.0,1.0,-78.0,1.0,-36.0,1.0,-37.0,1.0,1;
	5,-52.0,4.2,-56.1,3.5,-27.8,2.3,-28.7,1.6,1; % Qu two-microelectrode
	5,-65.3,0.9,-65.9,0.8,-18.6,4.2,-23.3,2.0,1;
	5,-84.8,2.5,-74.0,2.4,-35.8,1.4,-34.6,1.9,1;
	5,-77.1,0.5,-72.9,1.0,-25.6,1.0,-24.9,1.0,0;
	5,-82.0,0.5,-73.5,1.0,-16.8,0.5,-14.8,0.7,0;
	6,-74.3,2.3,-72.2,0.6,-36.7,1.1,-34.8,1.7,0;
	6,-51.5,0.4,-50.8,0.3,-13.4,0.8,-14.2,0.3,1;
	6,-54.2,0.7,-50.7,1.3,-11.5,0.5,-11.6,0.6,1;
	7,-70.9,0.5,-65.7,0.5,-18.6,0.4,-17.4,1.8,0;
	7,-68.2,0.4,-69.2,0.4,-22.0,2.7,-27.7,1.3,1;
	% 8,-54.8,1.7,-62.6,2.3,4.7,0.7,-3.3,1.0,1; % Outlier (see Methods)
	8,-43.2,2.0,-47.8,1.5,-12.5,1.7,-16.5,1.4,0];

	mrks = {'o','v','square','square'};
	mrks = {'o','v','o'};

	clrs = lines(8);
	clrs(8,:) = clrs(4,:);
	clrs(4,:) = [0,1,1];
	clrs(5,:) = [0,0,1];
	clrs(6,:) = [0.5,1,0];

	ax(1) = axes('units','centimeters','FontSize',7); hold on;
	ax(1).Position = [1,0.3+fy-offset,3.4,2.5];
	for i = [1,2,3,7,8,6,4,5]
		tempT = T(T(:,1)==i,:);
		dCSI = tempT(:,4)-tempT(:,2);
		sdCSI = sqrt(tempT(:,5).^2+tempT(:,3).^2);
		L = size(tempT,1);
		locs = ([1:L]/2)-mean(1:L)/2;
		if(L>1)
			locs = locs*0.2;
		end
		for k = 1:L
			x = i+locs(k);
			y = dCSI(k);
			sy = sdCSI(k);
			h = plot([x,x],[-sy sy]+y,'color',clrs(i,:)); hold on;
			h.Color = [h.Color 0.6];
			scatter(x,y,20,clrs(i,:),'filled','Marker',mrks{1+tempT(k,end)},'MarkerFaceAlpha',0.6);
		end
	end
	xlim([0.5,8.5]);
	xticks([1:8]);
	xlabel('Nav1.x');
	ylabel(['\Delta SSI (mV)']);
	box off
	set(gca,'TickDir','out')
	set(gca,'LineWidth',1)
	line([0.5,8.5],[0,0],'color','k','LineStyle','--','LineWidth',1)
	set(gca,'FontSize',7);
	ylim([-15,15])
	h(1) = plot(nan,nan,'ok','MarkerSize',4,'MarkerFaceColor','k');
	h(2) = plot(nan,nan,'vk','MarkerSize',4,'MarkerFaceColor','k');
	L = legend(h,'HEK293','X. Oocyte');
	L.ItemTokenSize = [10,10];
	L.Position = [0.1107 0.2261 0.1916 0.0541];
	L.Box = 'off'

	ax(2) = axes('units','centimeters','FontSize',7); hold on;
	ax(2).Position = [4.9,0.3+fy-offset,3.4,2.5];
	for i = [1,2,3,7,8,6,4,5]
		tempT = T(T(:,1)==i,:);
		sCSI = sqrt(tempT(:,7).^2+tempT(:,3).^2);
		CSI = tempT(:,6)-tempT(:,2);
		dCSI = tempT(:,4)-tempT(:,2);
		sdCSI = sqrt(tempT(:,5).^2+tempT(:,3).^2);

		L = size(tempT,1);
		for k = 1:L
			x = CSI(k); sx = sCSI(k); y = dCSI(k); sy = sdCSI(k);
			h = plot([-sx sx]+x,[y, y],'color',clrs(i,:)); hold on;
			h.Color = [h.Color,0.6];
			h = plot([x x],[-sy, sy]+y,'color',clrs(i,:));
			h.Color = [h.Color,0.6];

			scatter(x,y,20,clrs(i,:),'filled','Marker',mrks{1+tempT(k,end)},'MarkerFaceAlpha',0.6);
		end
	end

	xlabel(['CSI* (mV)']);
	box off
	set(gca,'TickDir','out')
	set(gca,'LineWidth',1)
	xlim([10,80])
	line(get(gca,'xlim'),[0,0],'color','k','LineStyle','--','LineWidth',1)
	ylim([-15,15])

	FT = fitlm((T(:,6)-T(:,2)),T(:,4)-T(:,2))
	plot(10:80,FT.predict((10:80)'),'color','k')
	set(gca,'FontSize',7);

	R2 = FT.Rsquared.Ordinary;
	pV = fix(log(FT.anova{1,end})/log(10));
	text(60,-8,['R^2 = ' num2str(R2,2)],'FontSize',6);
	text(60,-12,sprintf(['p < 10^{' int2str(pV) '}']),'FontSize',6);