function plotfigureS4

	basePath = fileparts(fileparts(mfilename('fullpath')));
	dataPath = fullfile(basePath,'figures','dependencies','data');
	addpath(fullfile(basePath,'data_processing'));
	addpath(fullfile(basePath,'figures','dependencies','functions'));

	fig=figure('color','w','units','centimeters');
	fig.Position(3)=8.5;
	% fig.Position(4) =11;
	fig.Position(4) = 9;

	delY = 1;


	axes('Position',[0.2,0.1,0.35,0.85]);
	annotation('textbox','String','High CSI','FontSize',8,'VerticalAlignment','bottom',...
				'HorizontalAlignment','center','Position',[0.2,0.94,0.35,0.1],'linestyle','none');
	default = load(fullfile(dataPath,'Nav15ParsNB.mat'));
	plotshifts(default,1,delY,dataPath)
	text(-160,1+7*delY/2,sprintf('Faster\nInactivation'),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',8);
	text(-160,1+5*delY/2,sprintf('Slower\nInactivation'),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',8);
	text(-160,1+3*delY/2,sprintf('Faster\nActivation'),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',8);
	text(-160,1+delY/2,sprintf('Slower\nActivation'),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',8);
	text(-160,0.5,'Baseline','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',8);


	axes('Position',[0.625,0.1,0.35,0.85]);
	annotation('textbox','String','Low CSI','FontSize',8,'VerticalAlignment','bottom',...
				'HorizontalAlignment','center','Position',[0.625,0.94,0.35,0.1],'linestyle','none');
	default.Params(1) = default.Params(1)*0.01;
	default.Params(2) = default.Params(2)*100;
	plotshifts(default,2,delY,dataPath)

function plotshifts(default,typ,delY,dataPath)

	load(fullfile(dataPath,'fittingTemplate.mat'));
	template.Activation.Voltages = linspace(-160,30,200);
	template.Inactivation.Voltages = linspace(-160,30,200);

	clrs(1,:) = [1,0,0];
	clrs(2,:) = [0,0.3,1];

	[Q,OpenPositions] = nav15_NB(default.Params);
	[~,~,~,~,~,I_A,I_I] = getchannelfitness(template,Q,OpenPositions,zeros(8,1),[1,1]);
	V = template.Activation.Voltages;
	FT = FitBoltzman(V,I_A',-20,-10,60,1);
	plot(V,I_A'*FT.Gmx./(V-FT.ERev),'LineWidth',2,'color',clrs(1,:));
	hold on;
	V = template.Inactivation.Voltages;
	plot(V,-I_I,'LineWidth',2,'color',clrs(2,:));
	FTI = FitBoltzman2(V,I_I',-60,10,-1);

	box off
	set(gca,'TickDir','out');
	set(get(gca,'yaxis'),'visible','off');
	xlim([-120,15]);

	for i = 1:4
		Params = default.Params;
		switch i
		case 1
			Params(1) = Params(1)*10;
			Params(2) = Params(2)*10;
		case 2
			Params(1) = Params(1)/10;
			Params(2) = Params(2)/10;
		case 3
			Params(23) = Params(23)*10;
		case 4
			Params(23) = Params(23)/10;
		end
		[Q,OpenPositions] = nav15_NB(Params);
		[~,~,~,~,~,I_A,I_I] = getchannelfitness(template,Q,OpenPositions,zeros(8,1),[1,1]);
		V = template.Activation.Voltages;
		V = template.Activation.Voltages;
		FTa = FitBoltzman(V,I_A',-20,-10,60,1);
		% plot(V,delY*i+I_A'*FTa.Gmx./(V-FTa.ERev),'LineWidth',2);
		V = template.Inactivation.Voltages;
		FTi = FitBoltzman2(V,I_I',-60,10,-1);
		% plot(V,delY*i+-I_I,'LineWidth',2);
		line([FTi.v50,FTa.v50],[delY*(4-i)+1.5,delY*(4-i)+1.5],'color','k','LineWidth',1)
		h=scatter(FTa.v50,delY*(4-i)+1.5,40,clrs(1,:),'filled');
		h=scatter(FTi.v50,delY*(4-i)+1.5,40,clrs(2,:),'filled');
	end
	ylim([0,delY*4+1]);
	yl = get(gca,'ylim');
	line(get(gca,'xlim'),[1,1]*yl(2),'color','k');
	set(gca,'FontSize',7);
	xlabel('Voltage (mV)')

	line([FT.v50,FT.v50],yl-[0,0.25],'linestyle','--','color',clrs(1,:),'LineWidth',1);
	line([FTI.v50,FTI.v50],yl-[0,0.25],'linestyle','--','color',clrs(2,:),'LineWidth',1);