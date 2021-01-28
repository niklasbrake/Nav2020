function [F_Breakdown] = plotComparison(varargin)
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
defaultval = @nav14_MGO;
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

% I_A = I_A./(Template.Activation.Voltages-Template.ERev)';
% I_A = I_A/max(I_A);

offset = 40;
MarkerSize = 5;

% errorbar(Template.Activation.Voltages,abs(Template.Activation.m),Template.Activation.s/sqrt(Template.Activation.N),'Color',dataColor,'Marker',ActivationMarker,'MarkerFaceColor',dataColor,'MarkerSize',MarkerSize,'LineWidth',0.75); hold on;
% errorbar(Template.Inactivation.Voltages,abs(Template.Inactivation.m),Template.Inactivation.s/sqrt(Template.Inactivation.N),'Color',dataColor,'Marker',InactivationMarker,'MarkerFaceColor',dataColor,'MarkerSize',MarkerSize,'LineWidth',0.75);

load('ephys_SummaryData.mat');
errorbar(Nav15e.activation(:,1),Nav15e.activation(:,2),Nav15e.activation(:,3),'Color',dataColor,'Marker',ActivationMarker,'MarkerFaceColor',dataColor,'MarkerSize',MarkerSize,'LineWidth',0.75,'LineStyle','none'); hold on;
errorbar(Nav15e.inactivation(:,1),Nav15e.inactivation(:,2),Nav15e.inactivation(:,3),'Color',dataColor,'Marker',InactivationMarker,'MarkerFaceColor',dataColor,'MarkerSize',MarkerSize,'LineWidth',0.75,'LineStyle','none');


if(~strcmp(modelColor,'none'))
	% h = plot(Template.Activation.Voltages,abs(I_A),[ActivationMarker modelColor],'MarkerSize',MarkerSize,'LineWidth',0.75); drawnow; set(h.NodeChildren(1),'LineWidth',0.75);
	% hold on;
	% h1 = plot(Template.Inactivation.Voltages,abs(I_I),[InactivationMarker modelColor],'MarkerSize',MarkerSize,'LineWidth',0.75); drawnow; set(h1.NodeChildren(1),'LineWidth',0.75);
	FTA = FitBoltzman(Template.Activation.Voltages,I_A',-10,-10,59,1);
	FTI = FitBoltzman2(Template.Inactivation.Voltages,-I_I',-60,10,1);
	h = plot(Template.Activation.Voltages,FTA.Gmx*I_A./(Template.Activation.Voltages-FTA.ERev)','Marker',ActivationMarker,'Color',modelColor,'MarkerSize',MarkerSize,'LineWidth',1.5); drawnow; 
	set(h.NodeChildren(1),'LineWidth',0.75);
	% plot(Template.Activation.Voltages,1./(1+exp((Template.Activation.Voltages-FTA.v50)/FTA.k)),'Color',modelColor);
	h1 = plot(Template.Inactivation.Voltages,abs(I_I),'Marker',InactivationMarker','Color',modelColor,'MarkerSize',MarkerSize,'LineWidth',1.5); drawnow; 
	set(h1.NodeChildren(1),'LineWidth',0.75);
	% h1 = plot(Template.Inactivation.Voltages,FTI.Gmx*I_I./(Template.Inactivation.Voltages-FTI.ERev)',[InactivationMarker modelColor],'MarkerSize',MarkerSize,'LineWidth',0.75); drawnow; set(h1.NodeChildren(1),'LineWidth',0.75);
	% plot(Template.Inactivation.Voltages,1./(1+exp((Template.Inactivation.Voltages-FTI.v50)/FTI.k)),'Color',modelColor);
end

xlim([-150 30]); ylim([0,1]);
set(get(gca,'XAxis'),'LineWidth',0.75); set(get(gca,'YAxis'),'LineWidth',0.75);
set(get(gca,'XAxis'),'LineWidth',0.75); set(get(gca,'YAxis'),'LineWidth',0.75);
set(gca,'XTick',[-150,-100,-50 0]);
set(gca,'XTickLabels',{'-150','-100','-50','0 mV'});
ylabel('I/I_{max}');
set(gca,'FontSize',10);