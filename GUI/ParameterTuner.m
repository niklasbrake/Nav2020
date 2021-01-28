function varargout = ParameterTuner(varargin)
% PARAMETERTUNER MATLAB code for ParameterTuner.fig
%      PARAMETERTUNER, by itself, creates a new PARAMETERTUNER or raises the existing
%      singleton*.
%
%      H = PARAMETERTUNER returns the handle to a new PARAMETERTUNER or the handle to
%      the existing singleton*.
%
%      PARAMETERTUNER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PARAMETERTUNER.M with the given input arguments.
%
%      PARAMETERTUNER('Property','Value',...) creates a new PARAMETERTUNER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ParameterTuner_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ParameterTuner_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ParameterTuner

% Last Modified by GUIDE v2.5 05-Mar-2019 14:48:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ParameterTuner_OpeningFcn, ...
                   'gui_OutputFcn',  @ParameterTuner_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before ParameterTuner is made visible.
function ParameterTuner_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ParameterTuner (see VARARGIN)

% Choose default command line output for ParameterTuner
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ParameterTuner wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ParameterTuner_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fileName,pathName] = uigetfile('*.MModel');
[~,fileName] = fileparts(fileName);
copyfile(fullfile(pathName,[fileName '.MModel']),fullfile(pathName,[fileName '.mat']))
load(fullfile(pathName,[fileName '.mat']));
delete(fullfile(pathName,[fileName '.mat']));

handles.text4.String = [fileName '.MModel'];

% handles.Im = Im;

params = constructModelCode(Edges,Parameters,'modelfile.m',pathName,OpenPositions);

handles.params = params;
handles.uitable1.Data = {};
for i = 1:length(handles.params)
	handles.uitable1.Data{i,1} = handles.params{i};
end
guidata(hObject,handles);

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = PickWorkspaceVariable;
waitfor(h);
handles.template = evalin('base','privatevar');
evalin('base','clear privatevar');
[~,s1,s2] = fileparts(handles.template.path);
handles.text3.String = [s1 s2];
guidata(hObject,handles);


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.axes1); hold off; cla;
axes(handles.axes2); hold off; cla;
axes(handles.axes3); hold off; cla;
% axes(handles.axes4); hold off; cla;
axes(handles.axes5); hold off; cla;

MB = plotComparison2('axeshandle',{handles.axes1,handles.axes2,handles.axes3,handles.axes5},'TemplateFolder',handles.template.path,'Parameters',handles.params,'ModelFunction',@modelfile);
handles.text5.String = num2str(MB);

% playActivation(hObject,eventdata,handles)
guidata(hObject,handles);

% --- Executes when entered data in editable cell(s) in uitable1.
function uitable1_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

if(~isnumeric(eventdata.NewData))
	hObject.Data{eventdata.Indices(1),eventdata.Indices(2)} = eventdata.PreviousData;
	return;
end
handles.params(eventdata.Indices(1)) = eventdata.NewData;
guidata(hObject,handles)


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


h = PickWorkspaceVariable;
waitfor(h);
handles.params = evalin('base','privatevar');
evalin('base','clear privatevar');
guidata(hObject,handles);
if(length(handles.params) ~= size(handles.uitable1.Data,1))
	disp('eee');
	return;
end

for i = 1:length(handles.params)
	handles.uitable1.Data{i,2} = handles.params(i);
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fileName,pathName] = uigetfile({'*.png';'*.jpg'});
[cdata,colormap] = imread(fullfile(pathName,fileName));
axes(handles.axes4); imshow(cdata,colormap);
guidata(hObject,handles);


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.axes6);
playActivation(hObject,eventdata,handles)


function playActivation(hObject,eventdata,handles)

[Q,OpenPositions] = modelfile(handles.params);
T_max = 500;
N = length(Q(0));
Template = constructTemplate(handles.template.path);

% Initial baseline at -100 mV 
dX_base = Q(-100*1e-3);
[T,temp] = ode15s(@(t,X) dX_base*X,[1:100]*1e-3,[1 zeros(1,N-1)]);
Xinit = temp(end,:)';

% Activation Initialization
VIMAX = 30;
V = Template.Activation.Voltages(1:VIMAX);
VSteps = length(V);
V_init = 1e-3*(-100);
dX = Q(V_init);

X2 = zeros(T_max,N,VSteps);

% Activation Protocol
for idx = 1:VSteps-1
	V_temp = 1e-3*V(idx);
	dX = Q(V_temp);
    [T,temp] = ode15s(@(t,X) dX*X,[1:T_max]*1e-5,Xinit);
	X2(:,:,idx) = temp;
end
actEst = squeeze(sum(X2(1:end,OpenPositions,:),2));
GMax = max(max(abs(actEst)));
actEst = actEst/GMax;

GV = actEst(1:500,:).*(Template.Activation.Voltages(1:VIMAX)-62);
% D1 = squeeze(sum(sum(X2(:,[2,3,5,6,8,9],:))));
% D4 =  squeeze(sum(sum(X2(:,[1,2,3,4,5,6],:))));

% figure;
% yyaxis left;
% plot(D1); hold on;
% plot(D4);
% yyaxis right;
% plot(max(abs(actEst)));

XX = zeros(3,5,500,VIMAX);
for i = 1:VIMAX
	XX(:,:,:,i) = permute(reshape(X2(:,:,i)',[5 3 500]),[2 1 3]);
end

% axes(handles.axes6);
for m = 7:6:VIMAX
% m = VIMAX;
for t = [1:5:500]
% for t = 80:400
	imagesc(XX(:,:,t,m)); drawnow;
end
end

% fprintf('%.4e\n',handles.params);