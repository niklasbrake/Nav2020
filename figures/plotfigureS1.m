function plotfigureS1

basePath = fileparts(fileparts(mfilename('fullpath')));
dataPath = fullfile(basePath,'figures\dependencies\data');

pt2cm = 0.0352778;

fig = figure('units','centimeters','color','w');
fig.Position(3) = 8.5;
fig.Position(4) = 7;
fig.Position(2) = fig.Position(2)/2;

clrs = lines(6);
clrs=clrs([2,3,4,6],:);

setPanLab = @(x,y,str) labelpanel(x/fig.Position(3),y/fig.Position(4),str);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% PANEL A %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1 = 0.35;
y1 = fig.Position(4)+1.9-3.7;
axw = 1.45;
axh = 1.75;

setPanLab(0.1,fig.Position(4)-0.35,'a');
axes('units','centimeters','Position',[x1,y1,axw,axh]); hold on;
set(get(gca,'xaxis'),'visible','off');
set(get(gca,'yaxis'),'visible','off');
WT = load(fullfile(dataPath,'\mH1\20170311c3\activation.mat'));
	V = mean(WT.Voltage(:,WT.Epochs(4):WT.Epochs(4)+500)');
	WT.Current = activationleakcorrection(WT.Voltage,WT.Current,WT.Epochs);
	X = WT.Current(:,WT.Epochs(4):WT.Epochs(4)+500)';
	[M,I] = max(abs(X));
	sig = sign(X(I+[0:length(I)-1]*501));
	FT = FitBoltzman(V,M.*sig/max(M(:)),-50,-10,60,1);
	X = X(:,find(V<FT.ERev-5));
	mX = max(abs(X(:)));
plot(X,'LineWidth',0.5,'color','k');
xlim([1,500]);
ylim([-1.1,0.2]*mX);
h=title('Nav1.5','FontSize',7,'FontWeight','normal'); yl = get(gca,'ylim');
h.Position(2) =yl(2)-0.12*range(yl);
line([300,400],[-0.9 -0.9]*mX,'LineWidth',1,'Color','k')
text(300,-mX,'1 ms','FontSize',6);
fr = 2e3/mX;
line([300,300],[-0.9 -0.9+fr]*mX,'LineWidth',1,'Color','k')
text(250,-0.9*mX,['2 nA'],'FontSize',6,'Rotation',90);

axes('units','centimeters','Position',[x1+axw+0.2,y1,axw,axh]); hold on;
set(get(gca,'xaxis'),'visible','off');
set(get(gca,'yaxis'),'visible','off');
	D1 = load(fullfile(dataPath,'\mH1-D1\20180301c6\activation.mat'));
	V = mean(D1.Voltage(:,D1.Epochs(4):D1.Epochs(4)+750)');
	D1.Current = activationleakcorrection(D1.Voltage,D1.Current,D1.Epochs);
	X = D1.Current(:,D1.Epochs(4):D1.Epochs(4)+750)';
	[M,I] = max(abs(X));
	sig = sign(X(I+[0:length(I)-1]*751));
	FT = FitBoltzman(V,M.*sig/max(M(:)),-50,-10,60,1);
	X = X(:,find(V<FT.ERev-5));
	mX = max(abs(X(:)));
plot(X,'LineWidth',0.5,'color',clrs(1,:));
xlim([1,750]);
ylim([-1.1,0.2]*mX);
h=title('DI-CN','FontSize',7,'FontWeight','normal'); yl = get(gca,'ylim');
h.Position(2) =yl(2)-0.12*range(yl);
line([500,700],[-0.9 -0.9]*mX,'LineWidth',1,'Color','k')
text(500,-mX,'2 ms','FontSize',6);
fr = 2e3/mX;
line([500,500],[-0.9 -0.9+fr]*mX,'LineWidth',1,'Color','k')
text(425,-0.9*mX,['2 nA'],'FontSize',6,'Rotation',90);


axes('units','centimeters','Position',[x1+2*(axw+0.2),y1,axw,axh]); hold on;
set(get(gca,'xaxis'),'visible','off');
set(get(gca,'yaxis'),'visible','off');
	D2 = load(fullfile(dataPath,'\mH1-D2\20180711c3\activation.mat'));
	V = mean(D2.Voltage(:,D2.Epochs(4):D2.Epochs(4)+500)');
	D2.Current = activationleakcorrection(D2.Voltage,D2.Current,D2.Epochs);
	X = D2.Current(:,D2.Epochs(4):D2.Epochs(4)+500)';
	[M,I] = max(abs(X));
	sig = sign(X(I+[0:length(I)-1]*501));
	FT = FitBoltzman(V,M.*sig/max(M(:)),-50,-10,60,1);
	X = X(:,find(V<FT.ERev-5));
	mX = max(abs(X(:)));
plot(X,'LineWidth',0.5,'color',clrs(2,:));
xlim([1,500]);
ylim([-1.1,0.2]*mX);
h=title('DII-CN','FontSize',7,'FontWeight','normal'); yl = get(gca,'ylim');
h.Position(2) =yl(2)-0.12*range(yl);
line([300,400],[-0.9 -0.9]*mX,'LineWidth',1,'Color','k')
text(300,-mX,'1 ms','FontSize',6);
fr = 2e3/mX;
line([300,300],[-0.9 -0.9+fr]*mX,'LineWidth',1,'Color','k')
text(250,-0.9*mX,['2 nA'],'FontSize',6,'Rotation',90);

axes('units','centimeters','Position',[x1+3*(axw+0.2),y1,axw,axh]); hold on;
set(get(gca,'xaxis'),'visible','off');
set(get(gca,'yaxis'),'visible','off');
	D3 = load(fullfile(dataPath,'\mH1-D3\20180309c12\activation.mat'));
	V = mean(D3.Voltage(:,D3.Epochs(4):D3.Epochs(4)+500)');
	D3.Current = activationleakcorrection(D3.Voltage,D3.Current,D3.Epochs);
	X = D3.Current(:,D3.Epochs(4):D3.Epochs(4)+500)';
	[M,I] = max(abs(X));
	sig = sign(X(I+[0:length(I)-1]*501));
	FT = FitBoltzman(V,M.*sig/max(M(:)),-50,-10,60,1);
	X = X(:,find(V<FT.ERev-5));
	mX = max(abs(X(:)));
plot(X,'LineWidth',0.5,'color',clrs(3,:));
xlim([1,500]);
ylim([-1.1,0.2]*mX);
h=title('DIII-CN','FontSize',7,'FontWeight','normal'); yl = get(gca,'ylim');
h.Position(2) =yl(2)-0.12*range(yl);
line([300,400],[-0.9 -0.9]*mX,'LineWidth',1,'Color','k')
text(300,-mX,'1 ms','FontSize',6);
fr = 1e3/mX;
line([300,300],[-0.9 -0.9+fr]*mX,'LineWidth',1,'Color','k')
text(250,-0.9*mX,['1 nA'],'FontSize',6,'Rotation',90);


axes('units','centimeters','Position',[x1+4*(axw+0.2),y1,axw,axh]); hold on;
set(get(gca,'xaxis'),'visible','off');
set(get(gca,'yaxis'),'visible','off');
	D4 = load(fullfile(dataPath,'\mH1-D4\20181022c7\activation.mat'));
	V = mean(D4.Voltage(:,D4.Epochs(4):D4.Epochs(4)+500)');
	D4.Current = activationleakcorrection(D4.Voltage,D4.Current,D4.Epochs);
	X = D4.Current(:,D4.Epochs(4):D4.Epochs(4)+500)';
	[M,I] = max(abs(X));
	sig = sign(X(I+[0:length(I)-1]*501));
	FT = FitBoltzman(V,M.*sig/max(M(:)),-50,-10,60,1);
	X = X(:,find(V<50));
	mX = max(abs(X(:)));
plot(X,'LineWidth',0.5,'color',clrs(4,:));
xlim([1,500]+25);
ylim([-1.1,0.2]*mX);
h=title('DIV-CN','FontSize',7,'FontWeight','normal'); yl = get(gca,'ylim');
h.Position(2) =yl(2)-0.12*range(yl);
line([350,450],[-0.9 -0.9]*mX,'LineWidth',1,'Color','k')
text(350,-mX,'1 ms','FontSize',6);
fr = 0.5e3/mX;
line([350,350],[-0.9 -0.9+fr]*mX,'LineWidth',1,'Color','k')
text(300,-0.95*mX,['0.5 nA'],'FontSize',6,'Rotation',90);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% PANEL B %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setPanLab(0.1,fig.Position(4)+1.8-4.1-pt2cm*8/2,'b');
axes('Units','centimeters','Position',[1.75,fig.Position(4)+1.9-8,6,4],'color','w',...
	'box','off','tickdir','out','FontSize',7,'LineWidth',1); hold on;
xlim([-110,35]); 	xlabel('Test Pulse Potential (mV)')
ylim([0,1.1]); 		ylabel('Peak Conductance (Norm.)')


load(fullfile(dataPath,'\mH1\GV_Curves.mat'));
G = 0.95*A_I'./(0.5+A_V'-ActFit(:,1)')./ActFit(:,2)';
V = mean(A_V);
eb(1) = errorbar(mean(A_V),mean(G'),stderror(G'),'squarek', ...
	'Color','k','MarkerFaceColor','k','MarkerSize',5);
G = G(V<35,:);
V = V(V<35);
FT = FitBoltzmanCurve2(V,mean(G'),-50,-5);
xtemp = V;
ytemp = FT(xtemp);
plot(xtemp,ytemp,'Color','k');

load(fullfile(dataPath,'\mH1-D1\GV_Curves.mat'));
G = 0.9*A_I'./(5+A_V'-ActFit(:,1)')./ActFit(:,2)';
V = mean(A_V);
eb(2) = errorbar(mean(A_V),mean(G'),stderror(G'),'^', ...
	'Color',clrs(1,:),'MarkerFaceColor',clrs(1,:),'MarkerSize',5);
G = G(V<35,:);
V = V(V<35);
FT = FitBoltzmanCurve2(V,mean(G'),-50,-5);
xtemp = V;
ytemp = FT(xtemp);
plot(xtemp,ytemp,'Color',clrs(1,:));

load(fullfile(dataPath,'\mH1-D2\GV_Curves.mat'));
G = A_I'./(A_V'-ActFit(:,1)')./ActFit(:,2)';
V = mean(A_V);
eb(3) = errorbar(mean(A_V),mean(G'),stderror(G'),'v', ...
	'Color',clrs(2,:),'MarkerFaceColor',clrs(2,:),'MarkerSize',5);
G = G(V<35,:);
V = V(V<35);
FT = FitBoltzmanCurve2(V,mean(G'),-50,-5);
xtemp = V;
ytemp = FT(xtemp);
plot(xtemp,ytemp,'Color',clrs(2,:));

load(fullfile(dataPath,'\mH1-D3\GV_Curves.mat'));
ActFit(5,:) = []; A_I(5,:) = []; A_V(5,:) = [];
ActFit(10,:) = []; A_I(10,:) = []; A_V(10,:) = [];
G = A_I'./(A_V'-ActFit(:,1)')./ActFit(:,2)';
V = mean(A_V);
eb(4) = errorbar(mean(A_V),mean(G'),stderror(G'),'d', ...
	'Color',clrs(3,:),'MarkerFaceColor',clrs(3,:),'MarkerSize',5);
G = G(V<35,:);
V = V(V<35);
FT = FitBoltzmanCurve2(V,mean(G'),-50,-5);
xtemp = V;
ytemp = FT(xtemp);
plot(xtemp,ytemp,'Color',clrs(3,:));

load(fullfile(dataPath,'\mH1-D4\GV_Curves.mat'));
G = A_I'./(A_V'-ActFit(:,1)')./ActFit(:,2)';
V = mean(A_V);
eb(5) = errorbar(mean(A_V),mean(G'),stderror(G'),'o', ...
	'Color',clrs(4,:),'MarkerFaceColor',clrs(4,:),'MarkerSize',5);
G = G(V<35,:);
V = V(V<35);
FT = FitBoltzmanCurve2(V,mean(G'),-50,-15);
xtemp = V;
ytemp = FT(xtemp);
plot(xtemp,ytemp,'Color',clrs(4,:));

set(gca,'xtick',[-100:25:25])
L = legend(eb,{'Nav1.5','DI-CN','DII-CN','DIII-CN','DIV-CN'},'box','off','FontSize',7, ...
	'TextColor','k','Units','centimeters');
L.Position=[1.8,3,1.5,1.5];
L.ItemTokenSize = [10,7];

axes('Units','centimeters','Position',[5.8,0.9,2.5,1.8],'color','none');
set(get(gca,'xaxis'),'visible','off')
set(get(gca,'yaxis'),'visible','off')
