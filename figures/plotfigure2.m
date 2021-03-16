function plotfigure2

basePath = fileparts(fileparts(mfilename('fullpath')));
dataPath = fullfile(basePath,'figures','dependencies','data');
addpath(fullfile(basePath,'data_processing'));
addpath(fullfile(basePath,'figures','dependencies','functions'));
addpath(fullfile(basePath,'modelling'));

fig = figure('Units','centimeters','color','w');
fig.Position(3) = 12;
fig.Position(4) = 8.5;
% fig.Position(3) = 11.4;
% fig.Position(4) = 8;


ax = axes('Position',[0.03,0.62,0.45,0.36]);
ax.Units = 'centimeters';
ax.Position(4) = ax.Position(3)*range([0.5,3.2])/range([0.7,5.25]);
ax.Units = 'normalized';
buildmodeldiagram(false); set(gca,'FontSize',7)

plotmodelpredictions(dataPath);
plotdata(dataPath);

labelpanel(0.01,0.95,'a');
labelpanel(0.01,0.54,'b');
labelpanel(0.01,0.34,'c');
labelpanel(0.265,0.34,'d');
labelpanel(0.535,0.95,'e');
labelpanel(0.535,0.7,'f');

function plotmodelpredictions(dataPath)
	clrs = lines(6);
	clrs=clrs([2:4,6],:);
	clrs(4,2)=0.3;

	load(fullfile(dataPath,'Nav15ParsNB.mat'));
	Fs = 1e-6;
	% Fs = 1e-6;
	[actEst1,V,I_A] = simActivation(Params,Fs);
	ax1 = axes('position',[0.34,0.1,0.16,0.27],'LineWidth',0.75,'fontsize',7,'box','off','tickdir','out','yscale','log','color','none'); hold on;
		[a1,b] = max(abs(actEst1));
		b(b==5e4) = nan;
		h=plot(V,b*1e-4,'LineWidth',1,'color','k');
		xlim([-70 35]);
		ylim([0.25 4]*1e-7/Fs); 
		XL = xlabel('Voltage (mV)');XL.Units = 'normalized'; XL.Position(2) = -0.2;
		ylabel('Time (ms)');
		set(gca,'ytick',[0.25,0.5,1,2,4]*1e-7/Fs)
		set(gca,'yticklabel',[0.25,0.5,1,2,4])
		set(gca,'xtick',[-70:35:50]);
	ax2 = axes('position',[0.08,0.1,0.16,0.27],'LineWidth',0.75,'fontsize',7,'box','off','tickdir','out','color','none'); hold on;
		FTA = FitBoltzman(V(1:20:200),I_A(1:20:200)',0,-10,59,1);
		plot(V,FTA.Gmx*I_A(:)./(V(:)-FTA.ERev),'LineWidth',1,'color',h.Color);
		xlim([-70 35]); ylim([0,1.1]);
		XL = xlabel('Voltage (mV)');XL.Units = 'normalized'; XL.Position(2) = -0.2;
		YL = ylabel('G/Gmax');
		YL.Position(1) = (YL.Position(1)+70)*0.9-70;
		set(gca,'xtick',[-70:35:50]);
	ax3 = axes('position',[0.205,0.4,0.125,0.15],'color','none'); 
		plot(actEst1(:,[1:10:200]),'Color',h.Color); hold on; 
		xlim([0,5e-3/Fs]); set(get(gca,'xaxis'),'visible','off'); set(get(gca,'yaxis'),'visible','off');
		TL=title('(\gamma*, \alpha*)','FontWeight','normal','fontsize',6); TL.Units = 'normalized'; TL.Position(2)=1.05;
		ylim([min(actEst1(:)),max(actEst1(:))]);


	Params(23) = Params(23)*5;
	Params(5) = Params(5)*5;
	[actEst2,V,I_A] = simActivation(Params,Fs);
	axes(ax1);
		[a2,b] = max(abs(actEst2));
		b(b==5e4) = nan;
		h=plot(V,b*1e-4,'LineWidth',1,'color',clrs(1,:)); set(gca,'yscale','log');
	axes(ax2)
		FTA = FitBoltzman(V(1:20:200),I_A(1:20:200)',0,-10,59,1);
		plot(V,FTA.Gmx*I_A(:)./(V(:)-FTA.ERev),'LineWidth',1,'Color',h.Color);

	ax3 = axes('position',[0.05,0.4,0.125,0.15],'color','none'); 
		plot(actEst2(:,[1:10:200]),'Color',h.Color); 
		xlim([0,5e-3/Fs]); set(get(gca,'xaxis'),'visible','off'); set(get(gca,'yaxis'),'visible','off');
		TL=title('(5\gamma*, 5\alpha*)','FontWeight','normal','fontsize',6); TL.Units = 'normalized'; TL.Position(2)=1.05;
		mX = min(actEst2(:));
		ylim([min(actEst2(:)),max(actEst2(:))]);
		line([3,4]*1e-3/Fs,[0.6 0.6]*mX,'LineWidth',1,'Color','k')
		text(3.5*1e-3/Fs,0.7*mX,'1 ms','FontSize',6,'HorizontalAlignment','center','margin',0.1);

	Params(23) = Params(23)/20;
	Params(5) = Params(5)/20;
	[actEst3,V,I_A] = simActivation(Params,Fs);
	axes(ax1);
		[a3,b] = max(abs(actEst3));
		b(b==5e4) = nan;
		h=plot(V,b*1e-4,'LineWidth',1,'color',clrs(4,:)); set(gca,'yscale','log');
	axes(ax2); 
		FTA = FitBoltzman(V(1:20:200),I_A(1:20:200)',0,-10,59,1);
		plot(V,FTA.Gmx*I_A(:)./(V(:)-FTA.ERev),'LineWidth',1,'Color',h.Color);
	ax3 = axes('position',[0.36,0.4,0.125,0.15],'color','none'); 
		plot(actEst3(:,[1:10:200]),'Color',h.Color); 
		xlim([0,5e-3/Fs]); set(get(gca,'xaxis'),'visible','off'); set(get(gca,'yaxis'),'visible','off');
		TL=title('(0.2\gamma*, 0.2\alpha*)','FontWeight','normal','fontsize',6); TL.Units = 'normalized'; TL.Position(2)=1.05;
		ylim([min(actEst3(:)),max(actEst3(:))]);


function [actEst,V,I_A] = simActivation(Params,Fs)
	[Q,OpenPositions,P] = nav15_NB(Params);
	N = length(Q(0));
	dX_base = Q(-100*1e-3); % Get transition matrix for V = -100 mV
	temp = expsolver(dX_base,[1:100]*1e-3,[1 zeros(1,N-1)])'; % Integrate for 100 ms
	Xinit = temp(end,:)'; % Take final "steady-state" conformation of system
	VSteps = 200;
	V = linspace(-75,50,VSteps);
	T_max = 5e-3/Fs; % 500 frames with frame time of 0.1 microseconds, thus 5 miliseconds.

	X2 = zeros(T_max,N,VSteps); % Allocates memory
	for idx = 1:VSteps % For each voltage step
		V_temp = 1e-3*V(idx); % Scale voltage from mV to V
		dX = Q(V_temp); % Get transition matrix
		X2(:,:,idx) = expsolver(dX,[1:T_max]*Fs,Xinit)'; % Integrate voltage step for 500 ms
	end
	actEst = squeeze(sum(X2(:,OpenPositions,:).*[1,1],2)).*(V-62); % Repeat for activation protocol

	A = min(actEst',[],2); % Find max of each current
	B = max(actEst',[],2); % Find min of each current
	[a,b] = max([abs(A) B],[],2); % Take whichever has the larger amplitude
	X = [A B]; % This is the an n x 2 matrix of (x,y) coordinates for peak amplitude
	modelActivationGV = X(sub2ind([length(B) 2],[1:length(B)]',b)); % Get the IV curve
	I_A = modelActivationGV/max(abs(modelActivationGV)); % Scale it by the max


function plotdata(dataPath)

	load(fullfile(dataPath,'ephys_SummaryData.mat'));
	load(fullfile(dataPath,'TimeToPeak.mat'));

	D1 = load(fullfile(dataPath,'Nav1.5e-D1','20180307c1','activation.mat'));
	D1.Current = activationleakcorrection(D1.Voltage,D1.Current,D1.Epochs); % corrects for the leak current
	repD1 = D1.Current(:,D1.Epochs(4)-75:D1.Epochs(4)+600);
	M = max(abs(repD1(:,100:end)),[],2);
	repD1 = repD1./M; % Normalize to peak at 1
	repD1(repD1>0.1) = nan;

	D2 = load(fullfile(dataPath,'Nav1.5e-D2','20180316c1','activation.mat'));
	D2.Current = activationleakcorrection(D2.Voltage,D2.Current,D2.Epochs); % corrects for the leak current
	repD2 = D2.Current(:,D2.Epochs(4)-75:D2.Epochs(4)+600);
	M = max(abs(repD2(:,100:end)),[],2);
	repD2 = repD2./M; % Normalize to peak at 1
	repD2(repD2>0.1) = nan;

	D4 = load(fullfile(dataPath,'Nav1.5e-D4','20190318c3','activation.mat'));
	D4.Current = activationleakcorrection(D4.Voltage,D4.Current,D4.Epochs); % corrects for the leak current
	repD4 = D4.Current(:,D4.Epochs(4)-75:D4.Epochs(4)+600);
	M = max(abs(repD4(:,100:end)),[],2);
	repD4 = repD4./M; % Normalize to peak at 1
	repD4(repD4>0.1) = nan;


	WT = load(fullfile(dataPath,'Nav1.5e','20170418c2','activation.mat'));
	WT.Current = activationleakcorrection(WT.Voltage,WT.Current,WT.Epochs); % corrects for the leak current
	repWT = WT.Current(:,WT.Epochs(4)-75:WT.Epochs(4)+600);
	repWT = repWT./max(abs(repWT),[],2); % Normalize to peak at 1
	repWT(repWT>0.1) = nan;


	clrs = lines(6);
	clrs=clrs([2:4,6],:);

	V = [-110:5:60];
	vIdcs = [16:2:20];
	% vIdcs = [17,17,17];
	T = ((1:length(repWT))-101)*1e-2;
	for i = 1:length(vIdcs)
		axes('position',[0.56+0.15*(i-1),0.75,0.12,0.24],'FontSize',7);
		hold on;
		iV = vIdcs(i);
		cl = clrs(1,:);
		% if(i==3)
		% 	repD1(iV,:) = repD4(iV,:);
		% 	cl = clrs(4,:);
		% end
		% if(i==2)
		% 	repD1(iV,:) = repD2(iV,:);
		% 	cl = clrs(2,:);
		% end
		plot(T(1:120),repWT(iV,1:120),'.k','MarkerSize',4);	 
		plot(T(120:600),repWT(iV,120:600),'-k','LineWidth',1);
		h=plot(T(1:140),repD1(iV,1:140),'.','MarkerSize',4,'color',cl); 
		plot(T(140:600),repD1(iV,140:600),'-','LineWidth',1,'Color',cl);
		xlim([-1,5]);
		ylim([-1.1,0.2]);
		text(0.5,0.15,[int2str(V(iV)) ' mV'],'FontWeight','normal','FontSize',7,'HorizontalAlignment','left','VerticalAlignment','middle');
		set(get(gca,'xaxis'),'visible','off')
		set(get(gca,'yaxis'),'visible','off')
		scatter(0,0.16,14,[0,0,0],'v','filled');
	end
	temp = get(gcf,'children');
	axes(temp(1));
	line([3,5],[-1,-1],'LineWidth',1,'Color','k')
	text(4,-1,'2 ms','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',7)


	%% Find voltages where normalized conductance is greater than 0.05 %%
	% (WT has larger currents, so use 0.025)
	% (DIV-CN has smaller currents, so use 0.1)
	vThresh = Nav15e.activation(Nav15e.activation(:,2)>0.025,1);
	idcs.WT = find(V>=vThresh(1));

	vThresh = minusD1.activation(minusD1.activation(:,2)>0.05,1);
	idcs.D1 = find(V>=vThresh(1));

	vThresh = minusD2.activation(minusD2.activation(:,2)>0.05,1);
	idcs.D2 = find(V>=vThresh(1));
	
	vThresh = minusD3.activation(minusD3.activation(:,2)>0.05,1);
	idcs.D3 = find(V>=vThresh(1));
	
	vThresh = minusD4.activation(minusD4.activation(:,2)>0.1,1);
	idcs.D4 = find(V>=vThresh(1));

	fn = fieldnames(idcs);
	slot = [nan,2,3,5,6];
	rn = {'WT','DI','DII','DIII','DIV'};
	cn = [nan,1,3,4,2];
	ydist = [0.45,0.1];
	for k = 2:length(fn)
		[tempA,tempB] = ind2sub([2,2],k-1);
		ax = axes('position',[0.6+0.22*(tempA-1),ydist(tempB),0.16,0.27],'fontsize',7, ...
					'LineWidth',0.75,'TickDir','out');
		hold on;
		h(1) = plotwitherror(V(idcs.WT),time2Peak.WT(idcs.WT,:),false,'Color','k','LineWidth',1);
		h(2) = plotwitherror(V(idcs.(fn{k})),time2Peak.(fn{k})(idcs.(fn{k}),:),false,'Color',clrs(k-1,:),'LineWidth',1);
		L = legend(h,{'WT', [rn{k} '-CN']},'FontSize',6);
		L.ItemTokenSize = [10,5];
		L.Position(1) = ax.Position(1)+ax.Position(3)-L.Position(3)+0.02;
		L.Position(2) = ax.Position(2)+ax.Position(4)-L.Position(4);
		L.Box = 'off';
		xlim([-70,35]);
		set(gca,'yscale','log')
		ylim([0.5,4]);
		set(gca,'ytick',[0.5,1,2,4]);
		set(gca,'xtick',[-70:35:35]);
		if(tempA==1)
			ylabel('Time (ms)')
		end
		if(tempB==2)
			XL = xlabel('Voltage (mV)');XL.Units = 'normalized'; XL.Position(2) = -0.2;
		end
	end