function plotfigure2b

basePath = fileparts(fileparts(mfilename('fullpath')));
dataPath = fullfile(basePath,'figures\dependencies\data');
load(fullfile(dataPath,'Nav15ParsNB.mat'));
load(fullfile(dataPath,'ParameterSensitivity(gamma_alpha).mat'));

offset = [-0.6,0.2,0,0];
fig = figure('Units','centimeters','color','w');
fig.Position(3) = 8.9;
fig.Position(4) = 4;

ax = axes('Position',[0.09,0.2,0.29,0.65]);
pos = [ax.Position(1) ax.Position(2)+ax.Position(4)-0.08 0.1 0.1];
n1 = length(Metrics_2par.var1);
n2 = length(Metrics_2par.var2);
x1 = Metrics_2par.var1;
y1 = Metrics_2par.var2;
x2 = interp1(1:n1,x1,linspace(1,n1,100));
y2 = interp1(1:n2,y1,linspace(1,n2,100));
[X1,Y1] = meshgrid(x1,y1);
[X2,Y2] = meshgrid(x2,y2);
Z = interp2(X1,Y1,Metrics_2par.act.v50',X2,Y2);
imagesc(x2,y2,Z);
set(gca,'LineWidth',0.75)
set(gca,'FontSize',7);
xlim([0,max(X2(:))]);
ylim([0,1e6]);
set(gca,'CLim',[-50,0]);
xlabel([char(945) '*']);
ylabel([char(947) '*']); axis xy;
ylab = get(get(gca,'yaxis'),'label');
ylab.Position(1) = ylab.Position(1)*1.1;
B = annotation(fig,'textbox',pos, 'String',{'E'}, 'LineStyle','none', 'FontWeight','bold', 'FontSize',8,'Margin',0);
B.Units = 'centimeter';
B.Position = B.Position + offset;
C = colorbar;
C.Position(1) = ax.Position(1)+ax.Position(3)+0.005;
C.Position(2) = ax.Position(2)+2*ax.Position(4)/3;
C.Position(4) = ax.Position(4)/3;
C.Limits = [-50 0];
C.Ticks = [-50,-25,0];
C.TickLabels = {'-50','-25','0 mV'};
C.Title.String = 'V_{1/2} (GV)';
C.Title.FontSize = 7;
C.Title.HorizontalAlignment = 'left';
C.Title.Margin = 0.1;
C.Title.Units = 'Normalized';
C.Title.Position = [-1,1.1,1];
set(gca,'xtick',[0:5:10]*1e4)
set(gca,'xticklabel',[0,0.5,1]);
xlab = get(get(gca,'xaxis'),'label');
xlab.Position(2) = xlab.Position(2)*0.85;
set(gca,'ytick',[0:5:10]*1e5); set(get(gca,'yaxis'),'exponent',6);
pos = [C.Position(1)-0.015,ax.Position(2),0.1,0.1];
annotation(fig,'textbox',pos, 'String',{'\times10^5'}, 'LineStyle','none', 'FontSize',7);


ax = axes('Position',[0.59,0.2,0.29,0.65]);
pos = [ax.Position(1) ax.Position(2)+ax.Position(4)-0.08 0.1 0.1];
Z = interp2(X1,Y1,Metrics_2par.inact.v50',X2,Y2);
imagesc(x2,y2,Z);
set(gca,'LineWidth',0.75)
set(gca,'FontSize',7);
xlim([0,max(X2(:))]);
ylim([0,1e6]);
set(gca,'CLim',[-110,-30]);
xlabel([char(945) '*']);
ylabel([char(947) '*']); axis xy;
ylab = get(get(gca,'yaxis'),'label');
ylab.Position(1) = ylab.Position(1)*1.1;
A = annotation(fig,'textbox',pos, 'String',{'F'}, 'LineStyle','none', 'FontWeight','bold', 'FontSize',8,'Margin',0);
A.Units = 'centimeter';
A.Position = A.Position + offset;
C = colorbar;
C.Position(1) = ax.Position(1)+ax.Position(3)+0.005;
C.Position(2) = ax.Position(2)+2*ax.Position(4)/3;
C.Position(4) = ax.Position(4)/3;
C.Ticks = [-110,-70,-30];
C.TickLabels = {'-110','-70','-30 mV'};
C.Title.String = 'V_{1/2} (SSI)';
C.Title.FontSize = 7;
C.Title.HorizontalAlignment = 'left';
C.Title.Margin = 0.1;
C.Title.Units = 'Normalized';
C.Title.Position = [-1,1.1,1];
set(gca,'xtick',[0:5:10]*1e4)
set(gca,'xticklabel',[0,0.5,1]);
xlab = get(get(gca,'xaxis'),'label');
xlab.Position(2) = xlab.Position(2)*0.85;
set(gca,'ytick',[0:5:10]*1e5); set(get(gca,'yaxis'),'exponent',6);
pos = [C.Position(1)-0.015,ax.Position(2),0.1,0.1];
annotation(fig,'textbox',pos, 'String',{'\times10^5'}, 'LineStyle','none', 'FontSize',7);

CM = parula;
colormap(0.2+0.8*CM)