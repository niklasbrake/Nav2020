function main(minusC1)

% fig = figure;
% fig.Position = [-1599 -204.2000 1600 750.4000];

temp = lines(3);
clrs(2,:) = temp(3,:);

X1 = text(1,3,'I_{4}','FontSize',7,'HorizontalAlignment','center','VerticalAlignment','top');
X2 = text(2,3,'I_{14}','FontSize',7,'HorizontalAlignment','center','VerticalAlignment','top');
X3 = text(3,3,'I_{24}','FontSize',7,'HorizontalAlignment','center','VerticalAlignment','top');
X4 = text(4,3,'I_{34}','FontSize',7,'HorizontalAlignment','center','VerticalAlignment','top');
X5 = text(5,3,'I_{}','FontSize',7,'HorizontalAlignment','center','VerticalAlignment','top');
X6 = text(1,2,'C_{4}','FontSize',7,'HorizontalAlignment','center','VerticalAlignment','top');
X7 = text(2,2,'C_{14}','FontSize',7,'HorizontalAlignment','center','VerticalAlignment','top');
X8 = text(3,2,'C_{24}','FontSize',7,'HorizontalAlignment','center','VerticalAlignment','top');
X9 = text(4,2,'C_{34}','FontSize',7,'HorizontalAlignment','center','VerticalAlignment','top');
X10 = text(5,2,'O_{4}','FontSize',7,'HorizontalAlignment','center','VerticalAlignment','top');
X11 = text(1,1,'C_{0}','FontSize',7,'HorizontalAlignment','center','VerticalAlignment','top');
X12 = text(2,1,'C_{1}','FontSize',7,'HorizontalAlignment','center','VerticalAlignment','top');
X13 = text(3,1,'C_{2}','FontSize',7,'HorizontalAlignment','center','VerticalAlignment','top');
X14 = text(4,1,'C_{3}','FontSize',7,'HorizontalAlignment','center','VerticalAlignment','top');
X15 = text(5,1,'O_{}','FontSize',7,'HorizontalAlignment','center','VerticalAlignment','top');
set(get(gca,'XAxis'),'Visible','off')
set(get(gca,'YAxis'),'Visible','off')
xlim([0.7,5.25]);
ylim([0.5 3.2]);

hold on;

horizontalLines(X1,X2,'\alpha x_\alpha y_\alpha','\beta x_\beta y_\beta');
horizontalLines(X2,X3,'\alpha x_\alpha y_\alpha','\beta x_\beta y_\beta');
horizontalLines(X3,X4,'\alpha x_\alpha y_\alpha','\beta x_\beta y_\beta');
horizontalLines(X4,X5,'\gamma_i','\delta_i');

horizontalLines(X6,X7,'\alpha x_\alpha ','\beta x_\beta ');
horizontalLines(X7,X8,'\alpha x_\alpha ','\beta x_\beta ');
horizontalLines(X8,X9,'\alpha x_\alpha ','\beta x_\beta ');
horizontalLines(X9,X10,'\gamma_4','\delta_4');

horizontalLines(X11,X12,'\alpha','\beta');
horizontalLines(X12,X13,'\alpha','\beta');
horizontalLines(X13,X14,'\alpha','\beta');
horizontalLines(X14,X15,'\gamma','\delta');

verticalLines(X1,X6,'r','i');
verticalLines(X2,X7,'r y_\beta','i y_\alpha');
verticalLines(X3,X8,'r y_\beta^2','i y_\alpha^2');
verticalLines(X4,X9,'r y_\beta^3','i y_\alpha^3');
verticalLines(X5,X10,'r_O','i_O');

verticalLines(X6,X11,'\beta_4','\alpha_4');
verticalLines(X7,X12,'\beta_4 x_\beta','\alpha_4 x_\alpha');
verticalLines(X8,X13,'\beta_4 x_\beta^2','\alpha_4 x_\alpha^2');
verticalLines(X9,X14,'\beta_4 x_\beta^3','\alpha_4 x_\alpha^3');
verticalLines(X10,X15,'\beta_{4O}','\alpha_{4O}');

if(~minusC1)
	return;
end

pos = get(gca,'Position');
pos = get(gca,'Position');
xl = get(gca,'xlim');
yl = get(gca,'ylim');
figCox = @(x) pos(1)+pos(3)*(x-xl(1))/range(xl);
figCoy = @(y) pos(2)+pos(4)*(y-yl(1))/range(yl);

% for i = 1:3
% 	pos2 = [1.5,0.55+(i-1),0.33,0.65];
% 	if(i~=2)
% 		pos2(3) = pos2(3)+0.035;
% 	end
% 	f2 = annotation('rectangle',[figCox(pos2(1)),figCoy(pos2(2)),pos(3)*pos2(3)/range(xl),pos(4)*pos2(4)/range(yl)],'FaceColor','w','LineStyle','none');
% 	f2.FaceAlpha = 0.9;
% end
x = [0.75 1.5 1.5 0.75];
y = [0.6 0.6 3.08 3.08];
pos2 = [0.75,0.6,0.75,2.48];
f = annotation('rectangle',[figCox(pos2(1)),figCoy(pos2(2)),pos(3)*pos2(3)/range(xl),pos(4)*pos2(4)/range(yl)],'FaceColor','w','LineStyle','--');
f.Color = clrs(2,:);
f.LineWidth = 1.5;
% f.FaceAlpha = 0.9;
f.FaceAlpha = 0.8;

verticalLine(X6,X11);
verticalLine(X7,X12);
verticalLine(X8,X13);
verticalLine(X9,X14);
verticalLine(X10,X15);

function horizontalLines(X1,X2,string1,string2)

pos = get(gca,'Position');
pos = get(gca,'Position');
xl = get(gca,'xlim');
yl = get(gca,'ylim');
figCox = @(x) pos(1)+pos(3)*(x-xl(1))/range(xl);
figCoy = @(y) pos(2)+pos(4)*(y-yl(1))/range(yl);

y1 = X1.Extent(2)+X1.Extent(4)/2+0.05;
y2 = X2.Extent(2)+X2.Extent(4)/2;
x1 = X1.Extent(1)+X1.Extent(3);
x2 = X2.Extent(1);
% ws=0.1;
offset = 0.04;

annotation('arrow',figCox([x1+offset,x2-offset]),figCoy([y1+0.025 y1+0.025]),'HeadWidth',3,'HeadLength',3,'HeadStyle','vback2');
annotation('arrow',figCox([x2-offset,x1+offset]),figCoy([y1-0.025 y1-0.025]),'HeadWidth',3,'HeadLength',3,'HeadStyle','vback2');


% quiver(x1+ws+offset,y1+0.025,x2-x1-1.5*ws,0,'LineWidth',0.75,'Color','k','MaxHeadSize',0.5);
% quiver(x2-ws-offset,y1-0.025,x1-x2+1.5*ws,0,'LineWidth',0.75,'Color','k','MaxHeadSize',0.5);
text((x2+x1)/2,y1+0.035,string1,'FontSize',5,'HorizontalAlignment','center','VerticalAlignment','bottom');
text((x2+x1)/2,y1-0.025,string2,'FontSize',5,'HorizontalAlignment','center','VerticalAlignment','top');


function verticalLines(X1,X2,string1,string2)

pos = get(gca,'Position');
pos = get(gca,'Position');
xl = get(gca,'xlim');
yl = get(gca,'ylim');
figCox = @(x) pos(1)+pos(3)*(x-xl(1))/range(xl);
figCoy = @(y) pos(2)+pos(4)*(y-yl(1))/range(yl);

y1 = X1.Extent(2);
y2 = X2.Extent(2)+X2.Extent(4);
x1 = X1.Position(1);
x2 = X2.Position(1);
ws=0.1;
offset = 0.01;


annotation('arrow',figCox([x1+0.025,x2+0.025]),figCoy([y1-offset y2+offset]),'HeadWidth',3,'HeadLength',3,'HeadStyle','vback2');
annotation('arrow',figCox([x1-0.025,x1-0.025]),figCoy([y2+offset y1-offset]),'HeadWidth',3,'HeadLength',3,'HeadStyle','vback2');

% quiver(x1+0.025,y1-ws-offset,0,y2-y1+1.5*ws,'LineWidth',0.75,'Color','k','MaxHeadSize',0.5);
% quiver(x1-0.025,y2+ws+offset,0,y1-y2-1.5*ws,'LineWidth',0.75,'Color','k','MaxHeadSize',0.5);
text(x1+0.075,(y2+y1)/2,string1,'FontSize',5,'HorizontalAlignment','left','VerticalAlignment','middle');
text(x1-0.075,(y2+y1)/2,string2,'FontSize',5,'HorizontalAlignment','right','VerticalAlignment','middle');


function verticalLine(X1,X2)

clrs = lines(6);
clrs(1,:) = [0,0,0];
clrs(2,:) = clrs(6,:);

pos = get(gca,'Position');
pos = get(gca,'Position');
xl = get(gca,'xlim');
yl = get(gca,'ylim');
figCox = @(x) pos(1)+pos(3)*(x-xl(1))/range(xl);
figCoy = @(y) pos(2)+pos(4)*(y-yl(1))/range(yl);

y1 = X1.Extent(2);
y2 = X2.Extent(2)+X2.Extent(4);
x1 = X1.Position(1);
x2 = X2.Position(1);
ws=0.1;
offset = 0.01;


% A = annotation('arrow',figCox([x1+0.025,x2+0.025]),figCoy([y1-offset y2+offset]),'HeadWidth',6,'HeadLength',6,'HeadStyle','vback2');
A = annotation('arrow',figCox([x1,x1]),figCoy([y2+offset y1-offset]),'HeadWidth',10,'HeadLength',10,'HeadStyle','vback2');
A.Color = clrs(2,:);
A.LineWidth = 2;
