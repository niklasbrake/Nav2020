function varargout = PickWorkspaceVariable(varargin)
% PICKWORKSPACEVARIABLE MATLAB code for PickWorkspaceVariable.fig
%      PICKWORKSPACEVARIABLE, by itself, creates a new PICKWORKSPACEVARIABLE or raises the existing
%      singleton*.
%
%      H = PICKWORKSPACEVARIABLE returns the handle to a new PICKWORKSPACEVARIABLE or the handle to
%      the existing singleton*.
%
%      PICKWORKSPACEVARIABLE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PICKWORKSPACEVARIABLE.M with the given input arguments.
%
%      PICKWORKSPACEVARIABLE('Property','Value',...) creates a new PICKWORKSPACEVARIABLE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PickWorkspaceVariable_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PickWorkspaceVariable_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PickWorkspaceVariable

% Last Modified by GUIDE v2.5 27-Feb-2019 10:14:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PickWorkspaceVariable_OpeningFcn, ...
                   'gui_OutputFcn',  @PickWorkspaceVariable_OutputFcn, ...
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


% --- Executes just before PickWorkspaceVariable is made visible.
function PickWorkspaceVariable_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PickWorkspaceVariable (see VARARGIN)

% Choose default command line output for PickWorkspaceVariable
handles.output = hObject;

T = evalin('base','whos');
for i = 1:length(T)
  handles.listbox2.String{i} = [T(i).name '     ' int2str(T(i).size(1)) 'x' int2str(T(i).size(2)) '    ' T(i).class];
end
handles.basevars = T;


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PickWorkspaceVariable wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PickWorkspaceVariable_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2
S = get(hObject,'Value');
varname = handles.basevars(S).name;
var = evalin('base',varname);
assignin('base','privatevar',var);
close;

% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
