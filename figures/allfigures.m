function allfigures

warning('off','MATLAB:print:FigureTooLargeForPage');
figs = string(int2str((1:7)'));
for i = 1:length(figs)
	eval(['plotfigure' figs{i}]);
	fig = gcf;
	fig.Units = 'inches';
	if(i==1 || i == 2)
		set(gcf,'PaperPositionMode','Auto','PaperUnits','inches','PaperSize',[5,fig.Position(4)],'Renderer','Painters');
	else
		set(gcf,'PaperPositionMode','Auto','PaperUnits','inches','PaperSize',[3.5,fig.Position(4)],'Renderer','Painters');
	end
	% print(['E:\Documents\Work\PhD\Nav15\Manuscript\figures\figure' figs{i} '.png'],'-dpng');
	% print(['E:\Documents\Work\PhD\Nav15\Manuscript\figures\figure' figs{i} '.pdf'],'-dpdf');
	print(['E:\Documents\Work\PhD\Nav15\Manuscript\figures\figure' figs{i} '.svg'],'-dsvg');
	close all;
end
warning('on','MATLAB:print:FigureTooLargeForPage');
