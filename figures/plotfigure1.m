function plotfigure1

basePath = fileparts(fileparts(mfilename('fullpath')));
dataPath = fullfile(basePath,'figures\dependencies\data');
addpath(fullfile(basePath,'data_processing'));
addpath(fullfile(basePath,'figures\dependencies\functions'));

fig = figure('units','centimeters','color','w');
fig.Position(3) = 11.4;
fig.Position(4) = 6.5;
fhy =fig.Position(4); 

clrs = lines(6);
clrs=clrs([2,3,4,6],:);

setPanLab = @(x,y,str) labelpanel(x/fig.Position(3),y/fig.Position(4),str);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% PANEL A %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
labelpanel(0.0,0.95,'a');
axes('units','centimeters','Position',[0,fhy-7+5.25,4.9,1.5]); hold on;
set(get(gca,'xaxis'),'visible','off'); xlim([0,1]);
set(get(gca,'yaxis'),'visible','off'); ylim([0,1]);

hlfHt = 0.165;
hlfHt_nm = hlfHt/1.5;
y1 = hlfHt_nm;
y4 = 1-hlfHt_nm;
y2 = (y4-y1)/3+y1;
y3 = 2*(y4-y1)/3+y1;

text(0,y1,'Nav1.5e DIV','HorizontalAlignment', ...
	'left','VerticalAlignment','middle','FontSize',6.7)
text(0,y2,'Nav1.5e DIII','HorizontalAlignment', ...
	'left','VerticalAlignment','middle','FontSize',6.7)
text(0,y3,'Nav1.5e DII','HorizontalAlignment', ...
	'left','VerticalAlignment','middle','FontSize',6.7)
text(0,y4,'Nav1.5e DI','HorizontalAlignment', ...
	'left','VerticalAlignment','middle','FontSize',6.7)

typesetAA = @(x,y,str) text(x,y,str,'HorizontalAlignment','center', ...
	'VerticalAlignment','middle', 'FontSize',7,'FontWeight','bold');

seqDI = {'L','S','A','L','R','T','F','R','V','L','R','A','L','K','T'};
seqDII = {'L','S','V','L','R','S','F','R','L','L','R','V','F','K','L'};
seqDIII = {'M','G','P','I','K','S','L','R','T','L','R','A','L','R','P'};
seqDIV = {'P','T','L','F','R','V','I','R','L','A','R','I','G','R','I'};
chargeLocations = [5,8,11];
x1 = 0.32;
wdth = (1-x1)/14;
for i = 1:length(seqDIV)
	tx(1)=typesetAA(x1+wdth*(i-1),y4,seqDI{i});
	tx(2)=typesetAA(x1+wdth*(i-1),y3,seqDII{i});
	tx(3)=typesetAA(x1+wdth*(i-1),y2,seqDIII{i});
	tx(4)=typesetAA(x1+wdth*(i-1),y1,seqDIV{i});
	if(sum(i==chargeLocations))
		for j = 1:4
			tx(j).Color = clrs(j,:);
		end
	end
end
for i = 1:length(chargeLocations)
	x1s = x1+wdth*(chargeLocations(i)-1-0.48);
	x2s = x1+wdth*(chargeLocations(i)-1+0.48);
	F=fill([x1s,x2s,x2s,x1s],[0,0,1,1],'w');
	F.FaceColor = 'none';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% PANEL B %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1 = 5.35;
y1 = fhy-7+5.18;
axw = 1.13;
axh = 1.75;

labelpanel(0.445,0.95,'b');
axes('units','centimeters','Position',[x1,y1,axw,axh]); hold on;
set(get(gca,'xaxis'),'visible','off');
set(get(gca,'yaxis'),'visible','off');
WT = load(fullfile(dataPath,'Nav1.5e\20170418c2\activation.mat'));
	V = mean(WT.Voltage(:,WT.Epochs(4):WT.Epochs(4)+600)');
	WT.Current = activationleakcorrection(WT.Voltage,WT.Current,WT.Epochs);
	X = WT.Current(:,WT.Epochs(4):WT.Epochs(4)+600)';
	[M,I] = max(abs(X));
	sig = sign(X(I+[0:length(I)-1]*501));
	FT = FitBoltzman(V,M.*sig/max(M(:)),-50,-10,60,1);
	X = X(:,find(V<FT.ERev-15));
	mX = max(abs(X(:)));
plot(X(:,1:2:end),'LineWidth',0.5,'color','k');
xlim([1,600]);
ylim([-1.1,0.2]*mX);
h=title('Nav1.5e','FontSize',7,'FontWeight','normal'); yl = get(gca,'ylim');
h.Position(2) =yl(2)-0.12*range(yl);
line([350,550],[-0.9 -0.9]*mX,'LineWidth',1,'Color','k')
text(350,-mX,'2 ms','FontSize',6);
fr = 5e3/mX;
line([350,350],[-0.9 -0.9+fr]*mX,'LineWidth',1,'Color','k')
text(250,-0.9*mX,['5 nA'],'FontSize',6,'Rotation',90);

axes('units','centimeters','Position',[x1+axw+0.1,y1,axw,axh]); hold on;
set(get(gca,'xaxis'),'visible','off');
set(get(gca,'yaxis'),'visible','off');
	D1 = load(fullfile(dataPath,'\Nav1.5e-D1\20180307c1\activation.mat'));
	V = mean(D1.Voltage(:,D1.Epochs(4):D1.Epochs(4)+750)');
	D1.Current = activationleakcorrection(D1.Voltage,D1.Current,D1.Epochs);
	X = D1.Current(:,D1.Epochs(4):D1.Epochs(4)+750)';
	[M,I] = max(abs(X));
	sig = sign(X(I+[0:length(I)-1]*751));
	FT = FitBoltzman(V,M.*sig/max(M(:)),-50,-10,60,1);
	X = X(:,find(V<FT.ERev-15));
	mX = max(abs(X(:)));
plot(X(:,1:2:end),'LineWidth',0.5,'color',clrs(1,:));
xlim([1,600]);
ylim([-1.1,0.2]*mX);
h=title('DI-CN','FontSize',7,'FontWeight','normal'); yl = get(gca,'ylim');
h.Position(2) =yl(2)-0.12*range(yl);
% line([500,700],[-0.9 -0.9]*mX,'LineWidth',1,'Color','k')
% text(500,-mX,'2 ms','FontSize',6);
fr = 2e3/mX;
line([450,450],[-0.9 -0.9+fr]*mX,'LineWidth',1,'Color','k')
text(350,-0.9*mX,['2 nA'],'FontSize',6,'Rotation',90);


axes('units','centimeters','Position',[x1+2*(axw+0.1),y1,axw,axh]); hold on;
set(get(gca,'xaxis'),'visible','off');
set(get(gca,'yaxis'),'visible','off');
	D2 = load(fullfile(dataPath,'Nav1.5e-D2\20180316c1\activation.mat'));
	V = mean(D2.Voltage(:,D2.Epochs(4):D2.Epochs(4)+600)');
	D2.Current = activationleakcorrection(D2.Voltage,D2.Current,D2.Epochs);
	X = D2.Current(:,D2.Epochs(4):D2.Epochs(4)+600)';
	[M,I] = max(abs(X));
	sig = sign(X(I+[0:length(I)-1]*501));
	FT = FitBoltzman(V,M.*sig/max(M(:)),-50,-10,60,1);
	X = X(:,find(V<FT.ERev-15));
	mX = max(abs(X(:)));
plot(X(:,1:2:end),'LineWidth',0.5,'color',clrs(2,:));
xlim([1,600]);
ylim([-1.1,0.2]*mX);
h=title('DII-CN','FontSize',7,'FontWeight','normal'); yl = get(gca,'ylim');
h.Position(2) =yl(2)-0.12*range(yl);
% line([300,500],[-0.9 -0.9]*mX,'LineWidth',1,'Color','k')
% text(300,-mX,'1 ms','FontSize',6);
fr = 2e3/mX;
line([450,450],[-0.9 -0.9+fr]*mX,'LineWidth',1,'Color','k')
text(350,-0.9*mX,['2 nA'],'FontSize',6,'Rotation',90);

axes('units','centimeters','Position',[x1+3*(axw+0.1),y1,axw,axh]); hold on;
set(get(gca,'xaxis'),'visible','off');
set(get(gca,'yaxis'),'visible','off');
	D3 = load(fullfile(dataPath,'Nav1.5e-D3\20180712c3\activation.mat'));
	V = mean(D3.Voltage(:,D3.Epochs(4):D3.Epochs(4)+600)');
	D3.Current = activationleakcorrection(D3.Voltage,D3.Current,D3.Epochs);
	X = D3.Current(:,D3.Epochs(4):D3.Epochs(4)+600)';
	[M,I] = max(abs(X));
	sig = sign(X(I+[0:length(I)-1]*501));
	FT = FitBoltzman(V,M.*sig/max(M(:)),-50,-10,60,1);
	X = X(:,find(V<FT.ERev-15));
	mX = max(abs(X(:)));
plot(X(:,1:2:end),'LineWidth',0.5,'color',clrs(3,:));
xlim([1,600]);
ylim([-1.1,0.2]*mX);
h=title('DIII-CN','FontSize',7,'FontWeight','normal'); yl = get(gca,'ylim');
h.Position(2) =yl(2)-0.12*range(yl);
% line([300,400],[-0.9 -0.9]*mX,'LineWidth',1,'Color','k')
% text(300,-mX,'1 ms','FontSize',6);
fr = 1e3/mX;
line([450,450],[-0.9 -0.9+fr]*mX,'LineWidth',1,'Color','k')
text(350,-0.9*mX,['1 nA'],'FontSize',6,'Rotation',90);


axes('units','centimeters','Position',[x1+4*(axw+0.1),y1,axw,axh]); hold on;
set(get(gca,'xaxis'),'visible','off');
set(get(gca,'yaxis'),'visible','off');
	D4 = load(fullfile(dataPath,'Nav1.5e-D4\20190318c3\activation.mat'));
	V = mean(D4.Voltage(:,D4.Epochs(4):D4.Epochs(4)+600)');
	D4.Current = activationleakcorrection(D4.Voltage,D4.Current,D4.Epochs);
	X = D4.Current(:,D4.Epochs(4):D4.Epochs(4)+600)';
	% X(1:45,:) = nan;
	[M,I] = max(abs(X));
	sig = sign(X(I+[0:length(I)-1]*501));
	FT = FitBoltzman(V,M.*sig/max(M(:)),-50,-10,60,1);
	X = X(:,find(V<50));
	mX = max(abs(X(:)));
plot(X(:,1:2:end),'LineWidth',0.5,'color',clrs(4,:));
xlim([1,600]+30);
ylim([-1.1,0.2]*mX);
h=title('DIV-CN','FontSize',7,'FontWeight','normal'); yl = get(gca,'ylim');
h.Position(2) =yl(2)-0.12*range(yl);
% line([350,450],[-0.9 -0.9]*mX,'LineWidth',1,'Color','k')
% text(350,-mX,'1 ms','FontSize',6);
fr = 0.3e3/mX;
line([450,450],[-0.9 -0.9+fr]*mX,'LineWidth',1,'Color','k')
text(350,-0.95*mX,['0.3 nA'],'FontSize',6,'Rotation',90);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% PANEL C %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

labelpanel(0.5,0.65,'c');
axes('Units','centimeters','Position',[7,0.9,4.25,3.6],'color','w',...
	'box','off','tickdir','out','FontSize',7,'LineWidth',1); hold on;
xlim([-110,35]); 	xlabel('Test Pulse Potential (mV)')
ylim([0,1.1]); 		ylabel('Peak Conductance (Norm.)')


load(fullfile(dataPath,'ephys_SummaryData.mat'));
eb(1) = errorbar(Nav15e.activation(:,1),Nav15e.activation(:,2),Nav15e.activation(:,3),'squarek', ...
	'Color','k','MarkerFaceColor','k','MarkerSize',4);
FT = FitBoltzmanCurve2(Nav15e.activation(:,1),Nav15e.activation(:,2),-50,-5);
xtemp = Nav15e.activation(:,1);
ytemp = FT(xtemp);
plot(xtemp,ytemp,'Color','k');

eb(2) = errorbar(minusD1.activation(:,1),minusD1.activation(:,2),minusD1.activation(:,3),'^', ...
	'Color',clrs(1,:),'MarkerFaceColor',clrs(1,:),'MarkerSize',4);
FT = FitBoltzmanCurve2(minusD1.activation(:,1),minusD1.activation(:,2),-50,-5);
xtemp = minusD1.activation(:,1);
ytemp = FT(xtemp);
plot(xtemp,ytemp,'Color',clrs(1,:));

eb(3) = errorbar(minusD2.activation(:,1),minusD2.activation(:,2),minusD2.activation(:,3),'v', ...
	'Color',clrs(2,:),'MarkerFaceColor',clrs(2,:),'MarkerSize',4);
FT = FitBoltzmanCurve2(minusD2.activation(:,1),minusD2.activation(:,2),-50,-5);
xtemp = minusD2.activation(:,1);
ytemp = FT(xtemp);
plot(xtemp,ytemp,'Color',clrs(2,:));

eb(4) = errorbar(minusD3.activation(:,1),minusD3.activation(:,2),minusD3.activation(:,3),'d', ...
	'Color',clrs(3,:),'MarkerFaceColor',clrs(3,:),'MarkerSize',4);
FT = FitBoltzmanCurve2(minusD3.activation(:,1),minusD3.activation(:,2),-50,-5);
xtemp = minusD3.activation(:,1);
ytemp = FT(xtemp);
plot(xtemp,ytemp,'Color',clrs(3,:));

eb(5) = errorbar(minusD4.activation(:,1),minusD4.activation(:,2),minusD4.activation(:,3),'o', ...
	'Color',clrs(4,:),'MarkerFaceColor',clrs(4,:),'MarkerSize',4);
FT = FitBoltzmanCurve2(minusD4.activation(:,1),minusD4.activation(:,2),-50,-5);
xtemp = minusD4.activation(:,1);
ytemp = FT(xtemp);
plot(xtemp,ytemp,'Color',clrs(4,:));

xlim([-100,35]);
set(gca,'xtick',[-100:25:25])
L = legend(eb,{'Nav1.5e','DI-CN','DII-CN','DIII-CN','DIV-CN'},'box','off','FontSize',7, ...
	'TextColor','k','Units','centimeters');
L.Position=[7,2.8,1.5,1.5];
L.ItemTokenSize = [10,7];

axes('Units','centimeters','Position',[9.9,0.9,1.4,1.2],'color','none');
set(get(gca,'xaxis'),'visible','off')
set(get(gca,'yaxis'),'visible','off')

line([0,1/3],[-130 -130],'Color','k','LineWidth',0.5);
for i = [-100,35]
	line([1/3 1],[i i],'Color','k','LineWidth',0.5);
end
line([1/3 1/3],[-130 35],'Color','k','LineWidth',0.5);
text(2/3,75,['35 mV'],'FontSize',6,'HorizontalAlignment','center')
text(2/3,-150,['-100 mV'],'FontSize',6,'HorizontalAlignment','center')
text(1/4,-120,['-130 mV'],'FontSize',6,'HorizontalAlignment','right','VerticalAlignment','bottom')
text(2/3,-31,'\vdots','FontSize',7,'HorizontalAlignment','center','VerticalAlignment','middle','interpreter','latex')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% PANEL A, Structure %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ax = axes;
ax.Units = 'centimeters';
ax.Position(1)=0;
ax.Position(2) = 0.15;
ax.Position(3) = 5.5;
ax.Position(4)=4.42;
I = imread(fullfile(dataPath,'NavStructure.png'));
imshow(I);
A = annotation('textbox',[0.3,0.12,0.02,0.02],...
	'String','DI','LineStyle','none','FontSize',7,'color',clrs(1,:), ...
	'VerticalAlignment','middle','HorizontalAlignment','center');
A = annotation('textbox',[0.09,0.09,0.02,0.02],...
	'String','DII','LineStyle','none','FontSize',7,'color',clrs(2,:), ...
	'VerticalAlignment','middle','HorizontalAlignment','center');
A = annotation('textbox',[0.13,0.6,0.02,0.02],...
	'String','DIII','LineStyle','none','FontSize',7,'color',clrs(3,:), ...
	'VerticalAlignment','middle','HorizontalAlignment','center')
A = annotation('textbox',[0.28,0.64,0.02,0.02],...
	'String','DIV','LineStyle','none','FontSize',7,'color',clrs(4,:), ...
	'VerticalAlignment','middle','HorizontalAlignment','center')

annotation('arrow',[0.27,0.33],[0.36,0.36],'HeadWidth',4,'HeadLength',4)

annotation('ellipse',[0.29,0.335,0.01,0.05])
annotation('rectangle',[0.295,0.365,0.005,0.03],'LineStyle','none','FaceColor','w')
annotation('arrow',[0.3,0.3],[0.375,0.38],'HeadWidth',2,'HeadLength',2,'HeadStyle','plain')

hwr = fig.Position(3)/fig.Position(4);


co = [0.41,0.21]; ang = pi/2;
A = annotation('arrow',[co(1),co(1)+0.005*cos(ang)],[co(2),co(2)+0.005*hwr*sin(ang)],'HeadWidth',4, ...
'HeadLength',4,'HeadStyle','plain');
annotation('textbox',[co - 0.025*[cos(ang),hwr*sin(ang)] - 0.025,0.05,0.05],...
'String','S1','LineStyle','none','FontSize',6, ...
'VerticalAlignment','middle','HorizontalAlignment','center');

co = [0.46,0.285]; ang = -pi;
A = annotation('arrow',[co(1),co(1)+0.005*cos(ang)],[co(2),co(2)+0.005*hwr*sin(ang)],'HeadWidth',4, ...
'HeadLength',4,'HeadStyle','plain');
annotation('textbox',[co - 0.025*[cos(ang),hwr*sin(ang)] - 0.025,0.05,0.05],...
'String','S2','LineStyle','none','FontSize',6, ...
'VerticalAlignment','middle','HorizontalAlignment','center');

co = [0.48,0.37]; ang = -pi;
A = annotation('arrow',[co(1),co(1)+0.005*cos(ang)],[co(2),co(2)+0.005*hwr*sin(ang)],'HeadWidth',4, ...
'HeadLength',4,'HeadStyle','plain');
annotation('textbox',[co - 0.025*[cos(ang),hwr*sin(ang)] - 0.025,0.05,0.05],...
'String','S3','LineStyle','none','FontSize',6, ...
'VerticalAlignment','middle','HorizontalAlignment','center');

co = [0.445,0.425]; ang = -pi*3/4;
A = annotation('arrow',[co(1),co(1)+0.005*cos(ang)],[co(2),co(2)+0.005*hwr*sin(ang)],'HeadWidth',4, ...
'HeadLength',4,'HeadStyle','plain');
annotation('textbox',[co - 0.025*[cos(ang),hwr*sin(ang)] - 0.025,0.05,0.05],...
'String','S4','LineStyle','none','FontSize',6, ...
'VerticalAlignment','middle','HorizontalAlignment','center');

co = [0.39,0.46]; ang = -pi*3/4;
A = annotation('arrow',[co(1),co(1)+0.005*cos(ang)],[co(2),co(2)+0.005*hwr*sin(ang)],'HeadWidth',4, ...
'HeadLength',4,'HeadStyle','plain');
annotation('textbox',[co - 0.025*[cos(ang),hwr*sin(ang)] - 0.025,0.05,0.05],...
'String','S5','LineStyle','none','FontSize',6, ...
'VerticalAlignment','middle','HorizontalAlignment','center');

co = [0.35,0.29]; ang = 0;
annotation('arrow',[co(1),co(1)+0.005*hwr*cos(ang)],[co(2),co(2)+0.005*sin(ang)],'HeadWidth',4, ...
'HeadLength',4,'HeadStyle','plain');
annotation('textbox',[co - 0.015*[hwr*cos(ang),sin(ang)] - 0.025,0.05,0.05],...
'String','S6','LineStyle','none','FontSize',6, ...
'VerticalAlignment','middle','HorizontalAlignment','center');