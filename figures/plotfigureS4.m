function plotfigureS4

	basePath = fileparts(fileparts(mfilename('fullpath')));
	dataPath = fullfile(basePath,'figures','dependencies','data');
	addpath(fullfile(basePath,'data_processing'));
	addpath(fullfile(basePath,'figures','dependencies','functions'));

	fig=figure('color','w','units','centimeters');
	fig.Position = [0,0,8.5,13];


	delY = 1;


	default = load(fullfile(dataPath,'Nav15ParsNB.mat'));
	plotshifts(default,1,delY,dataPath)

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
	FT = FitBoltzman(template.Activation.Voltages,I_A',-20,-10,60,1);

	xactWT = template.Activation.Voltages;
	yactWT = I_A'*FT.Gmx./(xactWT-FT.ERev);
	yinactWT = -I_I';
	xinactWT = template.Inactivation.Voltages;
	ax = subplot(5,2,typ);
	ax.Position(1) = ax.Position(1)+0.1-0.02*(typ-1);
		plot(xactWT,yactWT,'k','LineWidth',1); hold on;
		plot(xinactWT,yinactWT,'k','LineWidth',1);
		xlim([-160,30]);
		ylim([0,1]);
		set(gca,'LineWidth',0.75)
		set(gca,'FontSize',7);
		if(typ==1)
			T = title('High CSI','FontSize',10);
			T.Position(2) = 1.2;
		else
			T = title('Low CSI','FontSize',10);
			T.Position(2) = 1.2;
		end
		T = text(10,-0.17,'mV','FontSize',6.6);

	for i = 1:4
		Params = default.Params;
		switch i
		case 1
			Params(1) = Params(1)*10;
			str = '10 \times \alpha_{4k}';
		case 2
			Params(1) = Params(1)/10;
			str = '\alpha_{4k} / 10';
		case 3
			Params(23) = Params(23)*10;
			str = '10 \times \gamma_k';
		case 4
			Params(23) = Params(23)/10;
			str = '\gamma_k / 10';
		end
		[Q,OpenPositions] = nav15_NB(Params);
		[~,~,~,~,~,I_A,I_I] = getchannelfitness(template,Q,OpenPositions,zeros(8,1),[1,1]);
		
		ax = subplot(5,2,2*i+typ);
		ax.Position(1) = ax.Position(1)+0.1-0.02*(typ-1);
			V = template.Activation.Voltages;
			FT = FitBoltzman(V,I_A',-20,-10,60,1);
			plot(xactWT,yactWT,'color','k','LineWidth',1,'LineStyle','--'); hold on;
			plot(V,I_A'*FT.Gmx./(V-FT.ERev),'LineWidth',1,'color','k');
			V = template.Inactivation.Voltages;
			FTi = FitBoltzman2(V,I_I',-60,10,-1);
			plot(xinactWT,yinactWT,'color','k','LineWidth',1,'LineStyle','--');
			plot(V,-I_I,'LineWidth',1,'color','k');
			xlim([-160,30]);
			ylim([0,1]);
			set(gca,'LineWidth',0.75)
			set(gca,'FontSize',7);
			if(typ == 1)
				text(-250,0.5,str,'FontSize',8,'HorizontalAlignment','center');
			end
			T = text(10,-0.17,'mV','FontSize',6.6);
	end
