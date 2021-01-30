function plotfigureS2

	basePath = fileparts(fileparts(mfilename('fullpath')));
	dataPath = fullfile(basePath,'figures','dependencies','data');
	addpath(fullfile(basePath,'data_processing'));
	addpath(fullfile(basePath,'figures','dependencies','functions'));

	load(fullfile(dataPath,'fittingTemplate.mat'));
	load(fullfile(dataPath,'Nav15ParsNB.mat'));
	[Q,OpenPositions,P] = nav15_NB(Params);



	%%% Numerically integrate system with voltage protocols %%%
	% Part 0. Initialize constants
		T_max = 500; % 500 frames with frame time of 0.1 microseconds, thus 5 miliseconds.
		N = length(Q(0));
	% Part 0. Initial baseline at -100 mV 
		dX_base = Q(-100*1e-3); % Get transition matrix for V = -100 mV
		temp = expsolver(dX_base,[1:100]*1e-3,[1 zeros(1,N-1)])'; % Integrate for 100 ms
		Xinit = temp(end,:)'; % Take final "steady-state" conformation of system
	% Part 1. Initilize constants for inactivaiton protocol
		V = template.Inactivation.Voltages; % Get voltage steps used in the inactivation experiments from template
		VSteps = length(V);	% Total number of voltage steps
		dX_Pulse = Q(-10*1e-3); % Test pulse is -10 mV
	% Part 1. Inactivation Protocol
		X1 = zeros(T_max,N,VSteps); % Allocate memory
		preX1 = zeros(100,N,VSteps); % Allocate memory
		for idx = 1:VSteps % For each voltage step
			V_temp = 1e-3*V(idx); % Scale voltage from mV to V
			dX = Q(V_temp); % Get transition matrix for this voltage
			preX1(:,:,idx) = expsolver(dX,[1:100]*1e-3,Xinit)'; % Integrate pre-pulse, lasting 100 ms
			X00 = preX1(end,:,idx)'; % Take "steady-state" conformation of system
			X1(:,:,idx) = expsolver(dX_Pulse,[1:T_max]*1e-5,X00)'; % Integrate test-pulse, for 500 ms.
		end
	% Part 2. Initilize constants for activaiton protocol
		V = template.Activation.Voltages; % Get voltage steps used in the activation experiments from template
		VSteps = length(V); % Total number of voltage steps
	% Part 2. Activation Protocol
		X2 = zeros(T_max,N,VSteps); % Allocates memory
		for idx = 1:VSteps % For each voltage step
			V_temp = 1e-3*V(idx); % Scale voltage from mV to V
			dX = Q(V_temp); % Get transition matrix
			X2(:,:,idx) = expsolver(dX,[1:T_max]*1e-5,Xinit)'; % Integrate voltage step for 500 ms
		end
	% Part 3. Recovery Initialization
		TSteps = template.Recovery.Delays;
		idx_temp = 21; %Index for V = -10 mV
		X00 = X2(end,:,idx_temp);
	% Part 3. Recovery Protocol
		X3 = zeros(T_max,N,length(TSteps));
		temp1 = expsolver(dX_Pulse,[1:80]*1e-3,X00')';
		for idx = 1:length(TSteps)	
			Xi3(:,:,idx) = expsolver(dX_base,linspace(0,TSteps(idx),100)*1e-3,temp1(end,:))';
			X3(:,:,idx) = expsolver(dX_Pulse,[1:T_max]*1e-5,Xi3(end,:,idx))';
		end

	%%% Handelling of numerical simulations %%%
	% Part 0. Sometimes matrix exp. converts to complex numbers with 0 imaginary part.
		X1 = real(X1); 
		X2 = real(X2);
		X3 = real(X3);
	% Part 1. Max Current calculations
		inActEst = squeeze(sum(X1(1:500,OpenPositions,:),2)).*(-10-template.ERev); % Take the computed conductance x densitiy in "open states") and scale w.r.t reversal potential
		actEst = squeeze(sum(X2(1:500,OpenPositions,:),2)).*(template.Activation.Voltages-template.ERev); % Repeat for activation protocol
		recovEst = squeeze(sum(X3(1:500,OpenPositions,:),2)); % Repeat for recovery protol (don't scale because we're looking at fraction recovery and always pulse to same voltage anyways)
		GMax = max(abs(actEst(:))); % Scale by peak conductance to get in range [0,1]
		actEst = actEst/GMax;
		GMax = max(abs(inActEst(:)));
		inActEst = inActEst/GMax;
		GMax = max(recovEst(:));
		recovEst = recovEst/GMax;
	% Part 2. Find peak magnitude of each current
		A = min(actEst',[],2); % Find max of each current
		B = max(actEst',[],2); % Find min of each current
		[a,b] = max([abs(A) B],[],2); % Take whichever has the larger amplitude
		X = [A B]; % This is the an n x 2 matrix of (x,y) coordinates for peak amplitude
		modelActivationGV = X(sub2ind([length(B) 2],[1:length(B)]',b)); % Get the IV curve
		I_A = modelActivationGV/max(abs(modelActivationGV)); % Scale it by the max
	% Part 2. Repeat for inactivation
		A = min(inActEst',[],2);
		B = max(inActEst',[],2);
		[a,b] = max([abs(A) B],[],2);
		X = [A B];
		modelInactivationGV = X(sub2ind([length(B) 2],[1:length(B)]',b));
		I_I = modelInactivationGV/max(abs(modelInactivationGV));
	% Part 2. Repeat for recovery
		I_R = max(recovEst)/max(max(recovEst));


	WT.A = load(fullfile(dataPath,'Nav1.5e','20170418c2','activation.mat'));
	WT.I = load(fullfile(dataPath,'Nav1.5e','20170418c3','inactivation.mat'));
	WT.R = load(fullfile(dataPath,'Nav1.5e','20170418c2','recovery1ms.mat'));
	WT.R2 = load(fullfile(dataPath,'Nav1.5e','20170418c2','recovery10ms.mat'));

	fig = figure('Position',[950.2000  229.0000  560.0000  420.0000],'color','w');
	fig.Units = 'centimeters';
	fig.Position(3) = 8.9;
	fig.Position(4) = 8.9;

	axes('Position',[0.03,3/4,0.3,0.23]);
	V = mean(WT.A.Voltage(:,WT.A.Epochs(4):WT.A.Epochs(4)+500)');
	WT.A.Current = activationleakcorrection(WT.A.Voltage,WT.A.Current,WT.A.Epochs);
	X = WT.A.Current(:,WT.A.Epochs(4)+33:WT.A.Epochs(4)+500+33)';
	[M,I] = max(abs(X));
	sig = sign(X(I+[0:length(I)-1]*501));
	FT = FitBoltzman(V,M.*sig/max(M(:)),-50,-10,60,1);
	X = X(:,find(V<FT.ERev-5));
	plot(X/max(abs(X(:))),'k'); ylim([-1,0]); xlim([0,500]);
	set(get(gca,'xaxis'),'visible','off'); set(get(gca,'yaxis'),'visible','off');

	axes('Position',[0.35,3/4,0.3,0.23]);
	V = mean(WT.I.Voltage(:,WT.I.Epochs(5)+30:WT.I.Epochs(5)+530)');
	WT.I.Current = inactivationleakcorrection(WT.I.Voltage,WT.I.Current,WT.I.Epochs);
	X = WT.I.Current(:,WT.I.Epochs(5)+30:WT.I.Epochs(5)+530)';
	[M,I] = max(abs(X));
	sig = sign(X(I+[0:length(I)-1]*501));
	FT = FitBoltzman(V,M.*sig/max(M(:)),-50,-10,60,1);
	X = X(:,find(V<FT.ERev));
	plot(X/max(abs(X(:))),'k'); ylim([-1,0]); xlim([0,500]);
	set(get(gca,'xaxis'),'visible','off'); set(get(gca,'yaxis'),'visible','off');

	axes('Position',[0.67,3/4,0.3,0.23]); hold on;
	X = WT.R.Current(:,WT.R.Epochs(5)-100:WT.R.Epochs(5)+100e2)';
	for i = 1:size(X,2)
		X(1:101+i*10*10+45,i) = nan;
	end
	pre = mean(WT.R.Current(:,WT.R.Epochs(4)+1:WT.R.Epochs(4)+500));
	mX = max(abs(pre));
	plot((1:500)-300,pre/mX,'k'); hold on;
	plot(X./mX,'k'); hold on;

	X = WT.R2.Current(:,WT.R.Epochs(5)-100:WT.R.Epochs(5)+100e2)';
	for i = 1:size(X,2)
		X(1:101+i*100*10+45,i) = nan;
	end
	xlim([-300,30e2]);
	ylim([-1,0]);
	set(get(gca,'xaxis'),'visible','off'); set(get(gca,'yaxis'),'visible','off');

	axes('Position',[0.03,0.48,0.3,0.23]);
	plot(actEst,'Color','b','LineWidth',0.5); xlim([0,500]);
	set(get(gca,'xaxis'),'visible','off'); set(get(gca,'yaxis'),'visible','off');
	axes('Position',[0.35,0.48,0.3,0.23]);
	plot(inActEst,'Color','b','LineWidth',0.5); xlim([0,500]);
	set(get(gca,'xaxis'),'visible','off'); set(get(gca,'yaxis'),'visible','off');
	axes('Position',[0.67,0.48,0.3,0.23]); 
	pre = sum(X2(:,OpenPositions,21)');
	mx = max(pre(:));
	plot(((1:300)-300)*1e-2,-pre(1:300)/mx); hold on;
	for idx = 1:size(X3,3)
		plot([linspace(0,TSteps(idx),100),TSteps(idx)+[1:T_max]*1e-2], ...
			-[sum(Xi3(:,OpenPositions,idx)')/mx, ...
			sum(X3(:,OpenPositions,idx)')/mx],'Color','b','LineWidth',0.5);
	end
	xlim([-3,30]);
	set(get(gca,'xaxis'),'visible','off'); set(get(gca,'yaxis'),'visible','off');


	load(fullfile(dataPath,'ephys_SummaryData.mat'));

	axes('Position',[0.065,0.1,0.4,0.3]);
	errorbar(Nav15e.activation(:,1),Nav15e.activation(:,2),Nav15e.activation(:,3),'Color','k','Marker','square','MarkerFaceColor','k','LineWidth',0.75,'LineStyle','none'); hold on;
	FTA = FitBoltzman(template.Activation.Voltages,I_A',-10,-10,59,1);
	h = plot(template.Activation.Voltages,FTA.Gmx*I_A./(template.Activation.Voltages-FTA.ERev)','Marker','square','Color','b','LineWidth',1); drawnow; 
	set(h.NodeChildren(1),'LineWidth',0.75);
	xlim([-130 30]); ylim([0,1]);

	errorbar(Nav15e.inactivation(:,1),Nav15e.inactivation(:,2),Nav15e.inactivation(:,3),'Color','k','Marker','v','MarkerFaceColor','k','LineWidth',0.75,'LineStyle','none'); hold on;
	FTI = FitBoltzman2(template.Inactivation.Voltages,-I_I',-60,10,1);
	h1 = plot(template.Inactivation.Voltages,abs(I_I),'Marker','v','Color','b','LineWidth',1); drawnow; 
	set(h1.NodeChildren(1),'LineWidth',0.75);
	xlim([-130 30]); ylim([0,1]);
	xlabel('Voltage(mV)'); set(gca,'FontSize',7); set(gca,'LineWidth',1);
	box off; set(gca,'TickDir','out');


	axes('Position',[0.57,0.1,0.4,0.3]);
	T_Steps = template.Recovery.Delays;
	errorbar(template.Recovery.Delays,template.Recovery.m,template.Recovery.s/sqrt(template.Recovery.N),'Color','k', ...
	'LineWidth',1,'Marker','v','MarkerFaceColor','k','LineStyle','none');
	hold on;
	h1 = plot(T_Steps,I_R,'Marker','v','Color','b'); drawnow;
	set(gca,'XScale','log');xlim([1 150]);ylim([0,1]);
	xlabel('Inter-pulse interval (ms)'); set(gca,'FontSize',7); set(gca,'LineWidth',1);
	box off; set(gca,'TickDir','out');
	set(gca,'xtick',[1,10,100]);

	
	labelpanel(0.0,0.95,'a');
	labelpanel(0.0,0.68,'b');
	labelpanel(0.0,0.4,'c');
	labelpanel(0.52,0.4,'d');