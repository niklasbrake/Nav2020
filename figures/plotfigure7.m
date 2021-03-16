function plotfigure7

	basePath = fileparts(fileparts(mfilename('fullpath')));
	dataPath = fullfile(basePath,'figures','dependencies','data');
	addpath(fullfile(basePath,'data_processing'));
	addpath(fullfile(basePath,'figures','dependencies','functions'));

	fig = figure('color','w','units','centimeters');
	fig.Position(3) = 8.5; fig.Position(4) = 12;
	fig.Position(2) = 3;

	simulateDIVCN(dataPath);

	WT = load(fullfile(dataPath,'Nav1.5e','GV_Curves.mat'));
	B1 = load(fullfile(dataPath,'Nav1.5e+B1','GV_Curves.mat'));
	B3 = load(fullfile(dataPath,'Nav1.5e+B3','GV_Curves.mat'));
	ax=plotbeta(WT,B1,B3,5.3);
	% T.Position(1) = 15;
	axes(ax(1)); title(['Nav1.5e + ' char(946) '1'],'FontWeight','normal');
	axes(ax(2)); title(['Nav1.5e + ' char(946) '3'],'FontWeight','normal');
	axes(ax(3));
	ylim([-100,-60]);
	YL=ylabel('V_{1/2} SSI (mV)');		YL.Position(1) = -0.75;
	line([1,2],[-65,-65],'color','k'); text(1.5,-64,'***','FontSize',8,'VerticalAlignment','middle','HorizontalALignment','center')
	line([1,3],[-61,-61],'color','k'); text(2,-60,'***','FontSize',8,'VerticalAlignment','middle','HorizontalALignment','center')

	WT = load(fullfile(dataPath,'Nav1.4','GV_Curves.mat'));
	B1 = load(fullfile(dataPath,'Nav1.4+B1','GV_Curves.mat'));
	B3 = load(fullfile(dataPath,'Nav1.4+B3','GV_Curves.mat'));
	[ax,T]=plotbeta(WT,B1,B3,8.4);
	% T.Position(1) = 15;
	axes(ax(1)); title(['Nav1.4 + ' char(946) '1'],'FontWeight','normal');
	axes(ax(2)); title(['Nav1.4 + ' char(946) '3'],'FontWeight','normal');
	axes(ax(3));
	ax(3).YLim = [-80,-40];
	YL=ylabel('V_{1/2} SSI (mV)');		YL.Position(1) = -0.75;
	ax(3).XTickLabels{1} = [char(945) '1.4'];
	line([1,2],[-45,-45],'color','k'); text(1.5,-44,'***','FontSize',8,'VerticalAlignment','middle','HorizontalALignment','center')
	line([1,3],[-41,-41],'color','k'); text(2,-40,'***','FontSize',8,'VerticalAlignment','middle','HorizontalALignment','center')


	WT = load(fullfile(dataPath,'Nav1.6','GV_Curves.mat'));
	B1 = load(fullfile(dataPath,'Nav1.6+B1','GV_Curves.mat'));
	B3 = load(fullfile(dataPath,'Nav1.6+B3','GV_Curves.mat'));
	[ax,T] = plotbeta(WT,B1,B3,11.5);
	% T.Position(1) = 15;
	axes(ax(1)); title(['Nav1.6 + ' char(946) '1'],'FontWeight','normal');
	axes(ax(2)); title(['Nav1.6 + ' char(946) '3'],'FontWeight','normal');
	axes(ax(3));
	ax(3).YLim = [-70,-30];
	YL=ylabel('V_{1/2} SSI (mV)');		YL.Position(1) = -0.75;
	ax(3).XTickLabels{1} = [char(945) '1.6'];
	line([1,2],[-42,-42],'color','k'); text(1.7,-41,'p = 0.057','FontSize',6,'VerticalAlignment','bottom','HorizontalALignment','center')
	line([1,3],[-34,-34],'color','k'); text(2,-33,'p = 0.35','FontSize',6,'VerticalAlignment','bottom','HorizontalALignment','center')


	% clrs(1,:) = [0,0,0];
	% clrs(2,:) = [1,0,0];
	% clrs(3,:) = [0,0.3,1];
	% h(1) = plot(nan,nan,'o','color',clrs(2,:),'MarkerSize',4,'MarkerFaceColor',clrs(2,:)); hold on;
	% h(2) = plot(nan,nan,'square','color',clrs(3,:),'MarkerSize',4,'MarkerFaceColor',clrs(3,:)); hold on;
	% L=legend(h,{[' + ' char(946) '1'],[' + ' char(946) '3']},'box','off',...
	% 					'fontsize',7,'position',[0.19,0.83,0.18,0.07]);
	% L.ItemTokenSize=[10,5];

	labelpanel(0.01,0.95,'a');
	labelpanel(0.01,0.72,'b');
	labelpanel(0.72,0.72,'c');
	labelpanel(0.01,0.46,'d');
	labelpanel(0.72,0.46,'e');
	labelpanel(0.01,0.22,'f');
	labelpanel(0.72,0.22,'g');

function [ax,T] = plotbeta(WT,B1,B3,offset)
	fig = gcf;
	fy = fig.Position(4);

	clrs(1,:) = [0,0,0];
	clrs(2,:) = [1,0,0];
	clrs(3,:) = [0,0.3,1];

	ax(1) = axes('units','centimeters','FontSize',7); hold on;
	ax(1).Position = [0.6,0.25+fy-offset,2.3,2];
	ax(2) = axes('units','centimeters','FontSize',7); hold on;
	ax(2).Position = [3.55,0.25+fy-offset,2.3,2];
		
	axes(ax(1));
		V=mean(WT.A_V(:,1:30));
		% for i = 1:size(WT.A_V,1)
		% 	FT = FitBoltzman(WT.A_V(i,1:33),WT.A_I(i,1:33),-20,-10,60,0.1);
		% 	WT.ActFit(i,:) = coeffvalues(FT);
		% end
		G = WT.ActFit(:,2).*WT.A_I(:,1:30)./(V-WT.ActFit(:,1));
		outliers = find(or(WT.ActFit(:,1)<45,WT.ActFit(:,1)>75));
		outliers = [outliers;find(isoutlier(WT.ActFit(:,3)))];
		G(outliers,:)=[];
		G = G./max(G,[],2);
		WT.Av50 = WT.ActFit(:,4);
		WT.Av50(outliers) = [];
		h(1)=errorbar(V,mean(G),stderror(G),'LineWidth',0.5,'Marker','v','color',clrs(1,:),...
				'MarkerFaceColor',clrs(1,:),'LineStyle','none','MarkerSize',3);
		hold on;
		FT=FitBoltzmanCurve2(V,mean(G),-50,-10);
		plot(V,FT(V),'color',clrs(1,:));
	axes(ax(2));
		h(1)=errorbar(V,mean(G),stderror(G),'LineWidth',0.5,'Marker','v','color',clrs(1,:),...
				'MarkerFaceColor',clrs(1,:),'LineStyle','none','MarkerSize',3);
		hold on;
		FT=FitBoltzmanCurve2(V,mean(G),-50,-10);
		plot(V,FT(V),'color',clrs(1,:));

	axes(ax(1));
		V = mean(WT.I_V);
		for i = 1:size(WT.InactFit,1)
			FT = FitBoltzman2(V,WT.I_I(i,:),-70,10,-1);
			WT.InactFit(i,2:end) = coeffvalues(FT);
		end
		WT.Iv50 = WT.InactFit(:,4);
		outliers = find(isoutlier(WT.Iv50,'grubbs'));
		WT.InactFit(outliers,:) = [];
		WT.I_I(outliers,:) = [];
		I = WT.I_I.*WT.InactFit(:,2);
		WT.Iv50 = WT.InactFit(:,4);
		errorbar(V,mean(I),stderror(I),'LineWidth',0.5,'Marker','v','color',clrs(1,:),...
				'MarkerFaceColor',clrs(1,:),'LineStyle','none','MarkerSize',3);
		FT=FitBoltzmanCurve2(V,mean(I),-50,10);
		plot(V,FT(V),'color',clrs(1,:));
	axes(ax(2));
		errorbar(V,mean(I),stderror(I),'LineWidth',0.5,'Marker','v','color',clrs(1,:),...
				'MarkerFaceColor',clrs(1,:),'LineStyle','none','MarkerSize',3);
		FT=FitBoltzmanCurve2(V,mean(I),-50,10);
		plot(V,FT(V),'color',clrs(1,:));

	axes(ax(1));
		V=mean(B1.A_V(:,1:30));
		G = B1.A_I(:,1:30)./(V-B1.ActFit(:,1));
		G = G.*B1.ActFit(:,2);
		G = G./max(G,[],2);
		outliers = find(or(B1.ActFit(:,1)<45,B1.ActFit(:,1)>75));
		outliers = [outliers;find(isoutlier(B1.ActFit(:,3)))];
		G(outliers,:)=[];
		B1.Av50 = B1.ActFit(:,4);
		B1.Av50(outliers,:) = [];
		h(2)=errorbar(V,mean(G),stderror(G),'LineWidth',0.5,'Marker','o','color',clrs(2,:),...
				'MarkerFaceColor',clrs(2,:),'LineStyle','none','MarkerSize',3);
		FT=FitBoltzmanCurve2(V,mean(G),-50,-10);
		plot(V,FT(V),'color',clrs(2,:));

	axes(ax(1));
		V = mean(B1.I_V);
		for i = 1:size(B1.InactFit,1)
			FT = FitBoltzman2(V,B1.I_I(i,:),-70,10,-1);
			B1.InactFit(i,2:end) = coeffvalues(FT);
		end
		B1.Iv50 = B1.InactFit(:,4);
		outliers = find(isoutlier(B1.Iv50,'grubbs'));
		B1.InactFit(outliers,:) = [];
		B1.I_I(outliers,:) = [];
		I = B1.I_I.*B1.InactFit(:,2);
		B1.Iv50 = B1.InactFit(:,4);
		errorbar(V,mean(I),stderror(I),'LineWidth',0.5,'Marker','o','color',clrs(2,:),...
				'MarkerFaceColor',clrs(2,:),'LineStyle','none','MarkerSize',3);
		FT=FitBoltzmanCurve2(V,mean(I),-50,10);
		plot(V,FT(V),'color',clrs(2,:));

	axes(ax(2));
		V=mean(B3.A_V(:,1:30));
		G = B3.ActFit(:,2).*B3.A_I(:,1:30)./(V-B3.ActFit(:,1));
		G = G./max(G,[],2);
		outliers = find(or(B3.ActFit(:,1)<45,B3.ActFit(:,1)>75));
		outliers = [outliers;find(isoutlier(B3.ActFit(:,3)))];
		G(outliers,:)=[];
		B3.Av50 = B3.ActFit(:,4);
		B3.Av50(outliers,:) = [];
		h(3)=errorbar(V,mean(G),stderror(G),'LineWidth',0.5,'Marker','square','color',clrs(3,:),...
				'MarkerFaceColor',clrs(3,:),'LineStyle','none','MarkerSize',3);
		FT=FitBoltzmanCurve2(V,mean(G),-50,-10);
		plot(V,FT(V),'color',clrs(3,:));
	axes(ax(2));
		V = mean(B3.I_V);
		for i = 1:size(B3.InactFit,1)
			FT = FitBoltzman2(V,B3.I_I(i,:),-70,10,-1);
			B3.InactFit(i,2:end) = coeffvalues(FT);
		end
		B3.Iv50 = B3.InactFit(:,4);
		outliers = find(isoutlier(B3.Iv50,'grubbs'));
		B3.InactFit(outliers,:) = [];
		B3.I_I(outliers,:) = [];
		I = B3.I_I.*B3.InactFit(:,2);
		B3.Iv50 = B3.InactFit(:,4);
		errorbar(V,mean(I),stderror(I),'LineWidth',0.5,'Marker','square','color',clrs(3,:),...
				'MarkerFaceColor',clrs(3,:),'LineStyle','none','MarkerSize',3);
		FT=FitBoltzmanCurve2(V,mean(I),-50,10);
		plot(V,FT(V),'color',clrs(3,:));

	axes(ax(1));
		xlim([-120,35])
		ylim([-0.05,1.05]);
		set(gca,'tickdir','out');
		set(gca,'LineWidth',0.75);
		% set(gca,'xtick',[-50:25:25]);
		text(14,-0.2142,'mV','FontSize',7);

	axes(ax(2));
		xlim([-120,35])
		ylim([-0.05,1.05]);
		set(gca,'tickdir','out');
		set(gca,'LineWidth',0.75);
		% set(gca,'xtick',[-110:30:-50]);
		T=text(14,-0.2142,'mV','FontSize',7);

	ax(3) = axes('units','centimeters','FontSize',7); hold on;
	ax(3).Position = [7,0.25+fy-offset,1.4,2];

		I50_WT = WT.Iv50;
		I50_B1 = B1.Iv50;
		I50_B3 = B3.Iv50;

		Y = [I50_WT;I50_B1;I50_B3];
		G = [ones(size(I50_WT));2*ones(size(I50_B1));3*ones(size(I50_B3))];

		k = length(unique(G))
		df = length(G)-k
		[p,tbl,stats] = anova1(Y,G,'off');
		c = multcompare(stats,'display','off','CType','hsd','Alpha',0.01)

		plot(1+(mod(1:length(WT.Iv50),5)-1.5)/15,WT.Iv50,'.','color',clrs(1,:)); hold on;
		plot(2+(mod(1:length(B1.Iv50),5)-1.5)/15,B1.Iv50,'.','color',clrs(2,:));
		plot(3+(mod(1:length(B3.Iv50),5)-1.5)/15,B3.Iv50,'.','color',clrs(3,:));

		xlim([0.5,3.5]);
		xticks([1,2,3]);
		xticklabels({[char(945) '1.5e'],['+ ' char(946) '1'],['+ ' char(946) '3']});
		set(gca,'FontSize',7);
		box off;
		set(gca,'TickDir','out');
		xax = get(gca,'xaxis'); xax.TickLabelRotation = -45;


function simulateDIVCN(dataPath)
	fig = gcf;
	figh = fig.Position(4);

	clrs = lines(6);
	clrs(1,:) = [0,0,0];
	clrs(2,:) = clrs(6,:);

	% Nav1.5 model
	load(fullfile(dataPath,'Nav15ParsNB'));
	parsNav15 = Params;
	ax = axes('Units','centimeters');
	ax.Position = [0.6,figh-2.1,2.3,1.75];
		plotComparison(parsNav15,@nav15_NB,clrs(1,:),'square'); 
		parsNav15(1) = parsNav15(1)/5;
		% parsNav15(2) = parsNav15(2)/5;
		plotComparison(parsNav15,@nav15_NB,clrs(2,:),'v'); 
		xlim([-140,30]); ylim([0,1]); ylabel(''); box off;
		set(get(gca,'yaxis'),'visible','off');
		ch = get(gca,'Children');
		for i = 1:length(ch)
			ch(i).MarkerSize=3;
			ch(i).LineWidth=0.75;
		end
		set(gca,'xtick',[-100:50:0]);
		set(gca,'xticklabel',{'-100','-50','0 mV'});
		set(gca,'tickdir','out')
		set(gca,'fontsize',7);
		title('High CSI');
		set(gca,'LineWidth',0.75)

	% Low CSI Model
	load(fullfile(dataPath,'Nav15ParsNB'));
	parsNav15 = Params;
	parsNav15(1) = parsNav15(1)*0.01;
	parsNav15(2) = parsNav15(2)*100;
	ax = axes('Units','centimeters');
	ax.Position = [3.55,figh-2.1,2.3,1.75];
		plotComparison(parsNav15,@nav15_NB,clrs(1,:),'square'); 
		parsNav15(1) = parsNav15(1)/5;
		% parsNav15(2) = parsNav15(2)/100;
		% parsNav15(7) = parsNav15(7)/20;
		% parsNav15(8) = parsNav15(8)/20;
		plotComparison(parsNav15,@nav15_NB,clrs(2,:),'v');
		xlim([-140,30]); ylim([0,1]); ylabel(''); box off; set(gca,'tickdir','out')
		set(get(gca,'yaxis'),'visible','off');
		ch = get(gca,'Children');
		for i = 1:length(ch)
			ch(i).MarkerSize=3;
			ch(i).LineWidth=0.75;
		end
		set(gca,'xtick',[-100:50:0]);
		set(gca,'xticklabel',{'-100','-50','0 mV'});
		set(gca,'tickdir','out')
		set(gca,'fontsize',7);
		title('Low CSI');
		set(gca,'LineWidth',0.75)


	% L=labelpanel(0.125,0.24,'e');
	% 	annotation('textbox',[0.125+0.05,0.255,0.25,0.03], 'String',{'Nav1.5 model'}, 'LineStyle','none', 'FontWeight','normal', 'FontSize',7,'Margin',0);
	% L=labelpanel(0.55,0.24,'f');
	% 	annotation('textbox',[0.55+0.05,0.255,0.25,0.03], 'String',{'Low CSI model'}, 'LineStyle','none', 'FontWeight','normal', 'FontSize',7,'Margin',0);



function plotComparison(params,modelfunction,clr,markerType)

	basePath = fileparts(fileparts(mfilename('fullpath')));
	dataPath = fullfile(basePath,'figures','dependencies','data');
	load(fullfile(dataPath,'fittingTemplate.mat'));

	[Q,OpenPositions] = modelfunction(params);
	[~,~,~,~,~,I_A,I_I] = getchannelfitness(template,Q,OpenPositions,ones(1,8),[1,1]);


	FTA = FitBoltzman(template.Activation.Voltages,I_A',-10,-10,59,1);
	FTI = FitBoltzman2(template.Inactivation.Voltages,-I_I',-60,10,1);
	h = plot(template.Activation.Voltages,FTA.Gmx*I_A./(template.Activation.Voltages-FTA.ERev)','Marker',markerType,'Color',clr,'MarkerSize',5,'LineWidth',1.5); drawnow; 
	hold on;
	set(h.NodeChildren(1),'LineWidth',0.75);
	h1 = plot(template.Inactivation.Voltages,abs(I_I),'Marker',markerType','Color',clr,'MarkerSize',5,'LineWidth',1.5); drawnow; 
	set(h1.NodeChildren(1),'LineWidth',0.75);