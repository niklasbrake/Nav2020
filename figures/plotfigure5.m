function plotfigure5

	basePath = fileparts(fileparts(mfilename('fullpath')));
	dataPath = fullfile(basePath,'figures\dependencies\data');
	addpath(fullfile(basePath,'data_processing'));
	addpath(fullfile(basePath,'figures\dependencies\functions'));

	fig = figure('color','w','units','centimeters');
	fig.Position(3)= 8.5;
	fig.Position(4) = 9.5;
	figh = fig.Position(4);


	ax = axes('units','centimeters');
	ax.Position = [1.5,figh-3.5,6,0.28];
	ax.Position(4) = ax.Position(3)*range([0.5,3.2])/range([0.7,5.25]);
	ax.Units = 'normalized';
	buildmodeldiagram2(true); set(gca,'FontSize',7)
	% annotation('textbox',[0.1,0.96,0.05,0.03], 'String',{'a'}, 'LineStyle','none', 'FontWeight','bold', 'FontSize',8,'Margin',0.1,'VerticalAlignment','bottom');
	labelpanel(0.1,0.95,'a');

	simulateDIICN(dataPath);
	simulateDIVCN(dataPath);

function simulateDIICN(dataPath)
	fig = gcf;
	figh = fig.Position(4);

	temp = lines(3);
	clrs(1,:) = [0,0,0];
	clrs(2,:) = temp(3,:);

	% Nav1.5 model
	load(fullfile(dataPath,'Nav15ParsNB'));
	parsNav15 = Params;
	ax = axes('Units','centimeters');
	ax.Position = [-0.2+0.125+0.5,figh-6,2.25,1.75];
		plotComparison(parsNav15,@nav15_NB,clrs(1,:),'square'); 
		plotComparison(parsNav15,@nav15minusC1_NB,clrs(2,:),'v'); 
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

	% Nav1.4 model
	load(fullfile(dataPath,'Nav14ParsMGO'));
	parsNav14 = Params;
	ax = axes('Units','centimeters');
	ax.Position = [-0.2+0.125+3.25,figh-6,2.25,1.75];
		plotComparison(parsNav14,@nav14_MGO,clrs(1,:),'square'); 
		plotComparison(parsNav14,@nav14minusC1_MGO,clrs(2,:),'v'); 
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

	% Low CSI Model
	parsNav15(1) = parsNav15(1)*0.01;
	parsNav15(2) = parsNav15(2)*100;
	ax = axes('Units','centimeters');
	ax.Position = [-0.2+0.125+6,figh-6,2.25,1.75];
		plotComparison(parsNav15,@nav15_NB,clrs(1,:),'square'); 
		plotComparison(parsNav15,@nav15minusC1_NB,clrs(2,:),'v'); 
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

	L=labelpanel(0.04,0.555,'b');
		annotation('textbox',[0.04+0.05,0.57,0.25,0.03], 'String',{'Nav1.5 model'}, 'LineStyle','none', 'FontWeight','normal', 'FontSize',7,'Margin',0);
	L=labelpanel(0.35,0.555,'c');
		annotation('textbox',[0.35+0.05,0.57,0.25,0.03], 'String',{'Nav1.4 model'}, 'LineStyle','none', 'FontWeight','normal', 'FontSize',7,'Margin',0);
	L=labelpanel(0.665,0.555,'d');
		annotation('textbox',[0.67+0.05,0.57,0.25,0.03], 'String',{'Low CSI model'}, 'LineStyle','none', 'FontWeight','normal', 'FontSize',7,'Margin',0);


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
	ax.Position = [-0.2+4.4/3,figh-9,2.25,1.75];
		plotComparison(parsNav15,@nav15_NB,clrs(1,:),'square'); 
		parsNav15(1) = parsNav15(1)*20;
		parsNav15(2) = parsNav15(2)*20;
		parsNav15(7) = parsNav15(7)/20;
		parsNav15(8) = parsNav15(8)/20;
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

	% Low CSI Model
	load(fullfile(dataPath,'Nav15ParsNB'));
	parsNav15 = Params;
	parsNav15(1) = parsNav15(1)*0.01;
	parsNav15(2) = parsNav15(2)*100;
	ax = axes('Units','centimeters');
	ax.Position = [-0.2+8.9-4.4/3-2.25,figh-9,2.25,1.75];
		plotComparison(parsNav15,@nav15_NB,clrs(1,:),'square'); 
		parsNav15(1) = parsNav15(1)*20;
		parsNav15(2) = parsNav15(2)*20;
		parsNav15(7) = parsNav15(7)/20;
		parsNav15(8) = parsNav15(8)/20;
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


	L=labelpanel(0.125,0.24,'e');
		annotation('textbox',[0.125+0.05,0.255,0.25,0.03], 'String',{'Nav1.5 model'}, 'LineStyle','none', 'FontWeight','normal', 'FontSize',7,'Margin',0);
	L=labelpanel(0.55,0.24,'f');
		annotation('textbox',[0.55+0.05,0.255,0.25,0.03], 'String',{'Low CSI model'}, 'LineStyle','none', 'FontWeight','normal', 'FontSize',7,'Margin',0);



function plotComparison(params,modelfunction,clr,markerType)

	basePath = fileparts(fileparts(mfilename('fullpath')));
	dataPath = fullfile(basePath,'figures\dependencies\data');
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