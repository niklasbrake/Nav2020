function F_Breakdown = plotComparison(varargin)
ip = inputParser;

paramName = 'Parameters';
defaultval = [];
errorMsg = 'Input must be numeric.';
validationFcn = @(x) assert(isnumeric(x),errorMsg);
addParameter(ip,paramName,defaultval,validationFcn);

paramName = 'TemplateFolder';
defaultval = 'Nav1.5e';
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

paramName = 'GWeights';
defaultval = [1,1];
errorMsg = 'Input must be nummeric.';
validationFcn = @(x) assert(true,errorMsg);
addParameter(ip,paramName,defaultval,validationFcn);

paramName = 'ModelColor';
defaultval = 'r';
errorMsg = 'Input must be valid plot color.';
validationFcn = @(x) assert(true,errorMsg);
addParameter(ip,paramName,defaultval,validationFcn);

paramName = 'DataColor';
defaultval = 'k';
errorMsg = 'Input must be plot color.';
validationFcn = @(x) assert(true,errorMsg);
addParameter(ip,paramName,defaultval,validationFcn);

paramName = 'ActivationMarker';
defaultval = '*';
errorMsg = 'Input must be valid plot color.';
validationFcn = @(x) assert(true,errorMsg);
addParameter(ip,paramName,defaultval,validationFcn);

paramName = 'InactivationMarker';
defaultval = 'x';
errorMsg = 'Input must be plot color.';
validationFcn = @(x) assert(true,errorMsg);
addParameter(ip,paramName,defaultval,validationFcn);


parse(ip,varargin{:});
Parameters = ip.Results.Parameters;
TemplateFolder = ip.Results.TemplateFolder;
ModelFunction = ip.Results.ModelFunction;
axeshandle = ip.Results.axeshandle;
dataColor = ip.Results.DataColor;
modelColor = ip.Results.ModelColor;
ActivationMarker = ip.Results.ActivationMarker;
InactivationMarker = ip.Results.InactivationMarker;


load(fullfile(TemplateFolder,'AverageResponseCurves.mat'))
load(fullfile(TemplateFolder,'Recovery_Curves.mat'));
load(fullfile(TemplateFolder,'GV_Curves.mat'));

[Q,OpenPositions] = ModelFunction(Parameters);
Template = constructtemplate(TemplateFolder);
[~,F_Breakdown,actEst,inactEst,recovEst,I_A,I_I,I_R] = getchannelfitness(Template,Q,OpenPositions,ones(1,8),ip.Results.GWeights);


offset = 40;
MarkerSize = 5;

T_Steps = Template.Recovery.Delays;
errorbar(Template.Recovery.Delays,Template.Recovery.m,Template.Recovery.s/sqrt(Template.Recovery.N),'Color',dataColor, ...
	'MarkerSize',MarkerSize,'LineWidth',1,'Marker',ActivationMarker,'MarkerFaceColor',dataColor,'LineStyle','none');
hold on;
h1 = plot(T_Steps,I_R,'Marker',ActivationMarker,'Color',modelColor,'MarkerSize',MarkerSize); drawnow;
set(h1.NodeChildren(1),'LineWidth',1);
set(gca,'XScale','log');xlim([1 150]);
set(gca,'FontSize',12); set(get(gca,'XAxis'),'LineWidth',1); set(get(gca,'YAxis'),'LineWidth',1);
xlabel('Inter-pulse interval (ms)');
ylabel('Normalized Peak Conductance')
