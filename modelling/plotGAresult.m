function [max_F mean_F FitnessTrack Pidx] = plotGAresult(Q,OpenPositions,WT,template)

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
		inactEst = squeeze(sum(X1(1:500,OpenPositions,:),2)).*(-10-template.ERev); % Take the computed conductance x densitiy in "open states") and scale w.r.t reversal potential
		actEst = squeeze(sum(X2(1:500,OpenPositions,:),2)).*(template.Activation.Voltages-template.ERev); % Repeat for activation protocol
		recovEst = squeeze(sum(X3(1:500,OpenPositions,:),2)); % Repeat for recovery protol (don't scale because we're looking at fraction recovery and always pulse to same voltage anyways)
		GMax = max(abs(actEst(:))); % Scale by peak conductance to get in range [0,1]
		actEst = actEst/GMax;
		GMax = max(abs(inactEst(:)));
		inactEst = inactEst/GMax;
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
		A = min(inactEst',[],2);
		B = max(inactEst',[],2);
		[a,b] = max([abs(A) B],[],2);
		X = [A B];
		modelInactivationGV = X(sub2ind([length(B) 2],[1:length(B)]',b));
		I_I = modelInactivationGV/max(abs(modelInactivationGV));
	% Part 2. Repeat for recovery
		I_R = max(recovEst)/max(max(recovEst));



figure(2);
axes('Position',[0,2/3,1/3,1/3]);
V = mean(WT.A.Voltage(:,WT.A.Epochs(4):WT.A.Epochs(4)+500)');
WT.A.Current = activationleakcorrection(WT.A.Voltage,WT.A.Current,WT.A.Epochs);
X = WT.A.Current(:,WT.A.Epochs(4):WT.A.Epochs(4)+500)';
[M,I] = max(abs(X));
sig = sign(X(I+[0:length(I)-1]*501));
FT = FitBoltzman(V,M.*sig/max(M(:)),-50,-10,60,1);
X = X(:,find(V<FT.ERev-5));
plot(X/max(abs(X(:)))); ylim([-1,0]); xlim([0,500]);
set(get(gca,'xaxis'),'visible','off'); set(get(gca,'yaxis'),'visible','off');

axes('Position',[0,1/3,1/3,1/3]);
V = mean(WT.I.Voltage(:,WT.I.Epochs(5)+30:WT.I.Epochs(5)+530)');
WT.I.Current = inactivationleakcorrection(WT.I.Voltage,WT.I.Current,WT.I.Epochs);
X = WT.I.Current(:,WT.I.Epochs(5)+30:WT.I.Epochs(5)+530)';
[M,I] = max(abs(X));
sig = sign(X(I+[0:length(I)-1]*501));
FT = FitBoltzman(V,M.*sig/max(M(:)),-50,-10,60,1);
X = X(:,find(V<FT.ERev));
plot(X/max(abs(X(:)))); ylim([-1,0]); xlim([0,500]);
set(get(gca,'xaxis'),'visible','off'); set(get(gca,'yaxis'),'visible','off');

axes('Position',[0,0,1/3,1/3]);
X = WT.R.Current(:,WT.R.Epochs(5)-100:WT.R.Epochs(5)+3000)';
for i = 1:size(X,2)
	X(1:101+i*10*10+45,i) = nan;
end
X = X(102:end,:);
mX = max(abs(WT.R.Current(:,WT.R.Epochs(4):WT.R.Epochs(4)+500)'));
plot(X./mX); hold on;
xlim([0,150e2]);
set(get(gca,'xaxis'),'visible','off'); set(get(gca,'yaxis'),'visible','off');

axes('Position',[1/3,2/3,1/3,1/3]);
plot(actEst);
set(get(gca,'xaxis'),'visible','off'); set(get(gca,'yaxis'),'visible','off');
axes('Position',[1/3,1/3,1/3,1/3]);
plot(inactEst);
set(get(gca,'xaxis'),'visible','off'); set(get(gca,'yaxis'),'visible','off');
axes('Position',[1/3,0,1/3,1/3]); hold on;
for idx = 1:size(X3,3)
	plot([linspace(0,TSteps(idx),100),TSteps(idx)+[1:T_max]*1e-2], ...
		-[sum(Xi3(:,OpenPositions,idx),2)/max(sum(X2(:,OpenPositions,21),2)); ...
		sum(X3(:,OpenPositions,idx),2)/max(sum(X2(:,OpenPositions,21),2))]);
end
xlim([0,150]);
set(get(gca,'xaxis'),'visible','off'); set(get(gca,'yaxis'),'visible','off');


basePath = fileparts(fileparts(mfilename('fullpath')));
dataPath = fullfile(basePath,'figures\dependencies\data');

load(fullfile(dataPath,'AdamoData.mat'));

axes('Position',[2/3,2/3,1/3,1/3]);
errorbar(Nav15e.activation(:,1),Nav15e.activation(:,2),Nav15e.activation(:,3),'Color','k','Marker','square','MarkerFaceColor','k','LineWidth',0.75,'LineStyle','none'); hold on;
FTA = FitBoltzman(template.Activation.Voltages,I_A',-10,-10,59,1);
h = plot(template.Activation.Voltages,FTA.Gmx*I_A./(template.Activation.Voltages-FTA.ERev)','Marker','square','Color','b','LineWidth',1); drawnow; 
set(h.NodeChildren(1),'LineWidth',0.75);
xlim([-110 30]); ylim([0,1]);
set(get(gca,'xaxis'),'visible','off'); set(get(gca,'yaxis'),'visible','off');


axes('Position',[2/3,1/3,1/3,1/3]);
errorbar(Nav15e.inactivation(:,1),Nav15e.inactivation(:,2),Nav15e.inactivation(:,3),'Color','k','Marker','v','MarkerFaceColor','k','LineWidth',0.75,'LineStyle','none'); hold on;
FTI = FitBoltzman2(template.Inactivation.Voltages,-I_I',-60,10,1);
h1 = plot(template.Inactivation.Voltages,abs(I_I),'Marker','v','Color','b','LineWidth',1); drawnow; 
set(h1.NodeChildren(1),'LineWidth',0.75);
xlim([-160 30]); ylim([0,1]);
set(get(gca,'xaxis'),'visible','off'); set(get(gca,'yaxis'),'visible','off');


axes('Position',[2/3,0,1/3,1/3]);
T_Steps = template.Recovery.Delays;
errorbar(template.Recovery.Delays,template.Recovery.m,template.Recovery.s/sqrt(template.Recovery.N),'Color','k', ...
'LineWidth',1,'Marker','v','MarkerFaceColor','k','LineStyle','none');
hold on;
h1 = plot(T_Steps,I_R,'Marker','v','Color','b'); drawnow;
set(gca,'XScale','log');xlim([1 150]);ylim([0,1]);
set(get(gca,'xaxis'),'visible','off'); set(get(gca,'yaxis'),'visible','off');


annotation('line',[0,1],[1/3,1/3]);
annotation('line',[0,1],[2/3,2/3]);
annotation('line',[1/3,1/3],[0,1]);
annotation('line',[2/3,2/3],[0,1]);

drawnow;