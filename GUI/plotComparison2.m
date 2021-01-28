function F_Breakdown = plotComparison(varargin)
ip = inputParser;

paramName = 'Parameters';
defaultval = [];
errorMsg = 'Input must be numeric.';
validationFcn = @(x) assert(isnumeric(x),errorMsg);
addParameter(ip,paramName,defaultval,validationFcn);

paramName = 'TemplateFolder';
defaultval = 'E:\Documents\Work\PhD\Nav15\Experiments\Data\Nav1.5e';
errorMsg = 'Input must be nummeric.';
validationFcn = @(x) assert(exist(x)==7,errorMsg);
addParameter(ip,paramName,defaultval,validationFcn);

paramName = 'ModelFunction';
defaultval = @CGO_Model;
errorMsg = 'error with model function.';
validationFcn = @(x) assert(true,errorMsg);
addParameter(ip,paramName,defaultval,validationFcn);

paramName = 'axeshandle';
defaultval = [];
errorMsg = '';
validationFcn = @(x) assert(true,errorMsg);
addParameter(ip,paramName,defaultval,validationFcn);

parse(ip,varargin{:});
Parameters = ip.Results.Parameters;
TemplateFolder = ip.Results.TemplateFolder;
ModelFunction = ip.Results.ModelFunction;
axeshandle = ip.Results.axeshandle;
if(isempty(axeshandle))
	ax{1} = subplot(2,2,1);
	ax{2} = subplot(2,2,2); 
	ax{3} = subplot(2,2,3); 
	ax{4} = subplot(2,2,4); 
else
	ax = axeshandle;
end

load(fullfile(TemplateFolder,'AverageResponseCurves.mat'))
load(fullfile(TemplateFolder,'Recovery_Curves.mat'));
load(fullfile(TemplateFolder,'GV_Curves.mat'));

[Q,OpenPositions] = ModelFunction(Parameters);
Template = constructtemplate(TemplateFolder);
[~,F_Breakdown,actEst,inactEst,recovEst,I_A,I_I,I_R] = getchannelfitness(Template,Q,OpenPositions,ones(1,7),ones([1,length(OpenPositions)]));

axes(ax{1})
plot([1:500]/100,activationCurrents(1001:1500,1:34)/max(max(abs(activationCurrents))),'r')
hold on;
plot([41:500]/100,actEst(1:460,:),'-k')
xlim([1 500]/100); ylim([-1.1 0.1]); axis square;

axes(ax{2})
plot([41:500]/100,inactEst(1:460,:),'-k')
hold on;
plot([1:500]/100,inactivationCurrents(11001:11500,:)/max(max(abs(inactivationCurrents))),'r')
xlim([1 500]/100); ylim([-1.1 0.1]); axis square;

axes(ax{3})
errorbar(Template.Activation.Voltages,Template.Activation.m,Template.Activation.s/sqrt(Template.Activation.N),'Color','r'); hold on;
errorbar(Template.Inactivation.Voltages,Template.Inactivation.m,Template.Inactivation.s/sqrt(Template.Inactivation.N),'Color','b');
plot(Template.Activation.Voltages,I_A,'*r');
plot(Template.Inactivation.Voltages,I_I,'*b');
xlim([-150 35]); ylim([-1 0]); axis square;

axes(ax{4})
T_Steps = Template.Recovery.Delays;
% plot(a,b,'.r','MarkerSize',15);
errorbar(Template.Recovery.Delays,Template.Recovery.m,Template.Recovery.s/sqrt(Template.Recovery.N),'Color','r');
hold on;
plot(T_Steps,I_R,'*k','MarkerSize',5); axis square;
set(gca,'XScale','log');xlim([1 150]);