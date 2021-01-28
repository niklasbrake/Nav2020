function allfigures

warning('off','MATLAB:print:FigureTooLargeForPage');
figs = string(int2str((1:7)'));
figs = [figs;"S1";"S2";"S3";"S4";"S5"];
figs = ["S1","S3"];
for i = 1:length(figs)
	eval(['plotfigure' figs{i}]);
	fig = gcf;
	set(gcf,'PaperPositionMode','Auto','PaperUnits','Centimeters','PaperSize',[fig.Position(3),fig.Position(4)],'Renderer','Painters');
	% print(['E:\Documents\Work\PhD\Nav15\Manuscript\figures\figure' figs{i} '.pdf'],'-dpdf');
	print(['E:\Documents\Work\PhD\Nav15\Manuscript\figures\figure' figs{i} '.svg'],'-dsvg');
	% if(~isempty(strfind(figs{i},'S')))
	% 	print(['E:\Documents\Work\PhD\Nav15\Manuscript\figures\Supplemental Figure ' figs{i} '.png'],'-dpng','-r600');
	% else
	% 	print(['E:\Documents\Work\PhD\Nav15\Manuscript\figures\Figure ' figs{i} '.png'],'-dpng','-r600');
	% end
	close all;
end
warning('on','MATLAB:print:FigureTooLargeForPage');
