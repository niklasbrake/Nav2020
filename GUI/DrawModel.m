function varargout = DrawModel(varargin)
	%DRAWMODEL MATLAB code file for DrawModel.fig
	%      DRAWMODEL, by itself, creates a new DRAWMODEL or raises the existing
	%      singleton*.
	%
	%      H = DRAWMODEL returns the handle to a new DRAWMODEL or the handle to
	%      the existing singleton*.
	%
	%      DRAWMODEL('Property','Value',...) creates a new DRAWMODEL using the
	%      given property value pairs. Unrecognized properties are passed via
	%      varargin to DrawModel_OpeningFcn.  This calling syntax produces a
	%      warning when there is an existing singleton*.
	%
	%      DRAWMODEL('CALLBACK') and DRAWMODEL('CALLBACK',hObject,...) call the
	%      local function named CALLBACK in DRAWMODEL.M with the given input
	%      arguments.
	%
	%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
	%      instance to run (singleton)".
	%
	% See also: GUIDE, GUIDATA, GUIHANDLES

	% Edit the above text to modify the response to help DrawModel

	% Last Modified by GUIDE v2.5 26-Feb-2019 10:56:06

	% Begin initialization code - DO NOT EDIT
	gui_Singleton = 1;
	gui_State = struct('gui_Name',       mfilename, ...
	                   'gui_Singleton',  gui_Singleton, ...
	                   'gui_OpeningFcn', @DrawModel_OpeningFcn, ...
	                   'gui_OutputFcn',  @DrawModel_OutputFcn, ...
	                   'gui_LayoutFcn',  [], ...
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

% --- Executes just before DrawModel is made visible.
function DrawModel_OpeningFcn(hObject, eventdata, handles, varargin)
	% This function has no output args, see OutputFcn.
	% hObject    handle to figure
	% eventdata  reserved - to be defined in a future version of MATLAB
	% handles    structure with handles and user data (see GUIDATA)
	% varargin   unrecognized PropertyName/PropertyValue pairs from the
	%            command line (see VARARGIN)

	% Choose default command line output for DrawModel
	handles.output = hObject;
	handles.gridScale = 2.5e-2;
	set(gca,'DataAspectRatio',[1,1,1])
	% set( gca, 'DataAspectRatioMode', 'auto' )
	set(get(gca,'XAxis'),'Visible','off'); xlim([-0.25 1.25]);
	set(get(gca,'YAxis'),'Visible','off'); ylim([0 1]);
	temp = xlim; x1 = temp(1); x2 = temp(2);
	temp = ylim; y1 = temp(1); y2 = temp(2);
	set(gca,'XTick',[x1:handles.gridScale:x2]);
	set(gca,'YTick',[y1:handles.gridScale:y2]);
	handles.Vertex = {};
	handles.Edges = {};
	handles.deletePoints = [];
	hold on;
	grid on;

	% set(hObject,'toolbar','figure');
	% set(hObject,'menubar','figure');

	% Update handles structure
	guidata(hObject, handles);

	% UIWAIT makes DrawModel wait for user response (see UIRESUME)
	% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = DrawModel_OutputFcn(hObject, eventdata, handles)
	% varargout  cell array for returning output args (see VARARGOUT);
	% hObject    handle to figure
	% eventdata  reserved - to be defined in a future version of MATLAB
	% handles    structure with handles and user data (see GUIDATA)

	% Get default command line output from handles structure
	varargout{1} = handles.output;

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

	idcs = eventdata.Indices;
	try
		str2sym(eventdata.NewData);
		params = symvar(eventdata.NewData);
		for i = 1:length(params)
			idxtemp = find(strcmp(params{i},handles.uitable2.Data(:,1)));
			if(isempty(idxtemp))
				handles.uitable2.Data{end+1,1} = params{i};
				handles.uitable2.Data{end,2} = params{i};
			end
		end

		source =  handles.Edges{idcs(1)}.source;
		target =  handles.Edges{idcs(1)}.target;
		S = get(source,'Children'); S = S(3);
		T = get(target,'Children'); T = T(3);
		sourcePos = [S.XData S.YData];
		targetPos = [T.XData T.YData];
		Tv = targetPos'- sourcePos';
		Tv_norm = Tv/norm(Tv);
		Nv = [-Tv(2);Tv(1)];
		Nv_norm = Nv/norm(Nv);
		midPoint = mean([sourcePos;targetPos]);
		eps = 1e-2;

		T = handles.Edges{idcs(1)}.text;
		thet = atan(Tv(2)/Tv(1))*180/pi;
		T.Position(1:2) = midPoint;
		T.Rotation = thet;
		T.String = parseRateString(eventdata.NewData);
		posHat = [T.Extent(1)+0.5*T.Extent(3) T.Extent(2)+0.5*T.Extent(4)];
		difPos = posHat - midPoint;
		T.Position = midPoint-difPos+2.5*eps*Nv_norm';

		handles.uitable2.Data = updateUITable2(handles.uitable1,handles.uitable2);
	catch
		warning('Parse error. Input must have correct mathematical syntax.');
		handles.uitable1.Data{idcs(1),idcs(2)} = eventdata.PreviousData;
	end
	guidata(hObject,handles);

function Table = updateUITable2(uitable1,uitable2)
	allParams = {};
	for i = 1:size(uitable1.Data,1)
		vars = symvar(uitable1.Data{i,3});
		for j = 1:length(vars)
			if(isempty(find(strcmp(vars{j},allParams))))
				allParams{end+1} = vars{j};
			end
		end
	end
	removeParams = [];
	Ntemp = size(uitable2.Data,1);
	for i = 4:Ntemp
		if(isempty(find(strcmp(uitable2.Data{i,1},allParams))))
			removeParams(end+1) = i;
		end
	end
	keepIdcs = setdiff([1:Ntemp],removeParams);
	Table = uitable2.Data(keepIdcs,:);

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
	% hObject    handle to pushbutton2 (see GCBO)
	% eventdata  reserved - to be defined in a future version of MATLAB
	% handles    structure with handles and user data (see GUIDATA)

	[fileName,pathName] = uiputfile('*.m');
	if(fileName)
		constructModelCode(handles.uitable1.Data,handles.uitable2.Data,fileName,pathName,handles.OpenPositions);
	end
	open(fullfile(pathName,fileName));

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
	% hObject    handle to pushbutton3 (see GCBO)
	% eventdata  reserved - to be defined in a future version of MATLAB
	% handles    structure with handles and user data (see GUIDATA)

	[fileName,pathName] = uiputfile('*.MModel');
	if(fileName)
		N = length(handles.Vertex);
		aliveVertices = setdiff([1:N],handles.deletePoints);
		Edges = handles.uitable1.Data;
		for i = 1:length(aliveVertices)
			imap = aliveVertices(i);
			HPC = get(handles.Vertex{imap},'Children'); HPC = HPC(3);
			VertexPositions(i,:) = [HPC.XData HPC.YData];
			VertexTags{i} = get(handles.Vertex{imap},'Tag');
		end
		Parameters = handles.uitable2.Data;
		[Im,map] = frame2im(getframe(gca));
		OpenPositions = handles.OpenPositions;
		save(fullfile(pathName,fileName),'Edges','VertexPositions','Parameters','VertexTags','Im','OpenPositions');
	end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
	% hObject    handle to pushbutton4 (see GCBO)
	% eventdata  reserved - to be defined in a future version of MATLAB
	% handles    structure with handles and user data (see GUIDATA)

	[fileName,pathName] = uigetfile('*.MModel');
	if(fileName)
		[~,fileName] = fileparts(fileName);
		copyfile(fullfile(pathName,[fileName '.MModel']),fullfile(pathName,[fileName '.mat']))
		load(fullfile(pathName,[fileName '.mat']));
		delete(fullfile(pathName,[fileName '.mat']));

		cla
		handles.Vertex = {};
		handles.Edges = {};
		handles.deletePoints = [];
		handles.uitable1.Data = Edges;
		handles.uitable2.Data = Parameters;
		handles.OpenPositions = {};
		if(exist('OpenPositions'))
			handles.OpenPositions = OpenPositions;
		end
		for i = 1:size(VertexPositions,1)
			N = i-1;
			hPoint = impoint(handles.axes1,VertexPositions(i,:));
			setString(hPoint,VertexTags{i});
			hPointChildren = get(hPoint,'Children');
			hPointChildren(1).FontWeight = 'bold';
			hPointChildren(1).BackgroundColor = 'none';
			hPointChildren(1).FontSize = 6;
			hPointChildren(2).Visible = 'off';
			hPointChildren(3).Marker = 'square';
			hPointChildren(3).MarkerEdgeColor = 'k';
			hPointChildren(3).MarkerFaceColor = 'k';
			set(hPoint,'Tag',VertexTags{i});
			handles.Vertex{i} = hPoint;
			menu = get(hPoint,'UIContextMenu');
			temp = get(menu,'Children');
			set(temp(1),'MenuSelectedFcn',@(source,eventdata) deletePoint(source,eventdata,hPoint));
			temp = find(strcmp(handles.OpenPositions,VertexTags{i}));
			submenu = uimenu(menu,'Label','Add Connection');
			for j = 1:N
				uimenu(submenu,'Label',get(handles.Vertex{j},'Tag'),'Callback',@(source,eventdata) connectEdge(source,eventdata,i,j));
				menu2 = get(handles.Vertex{j},'UIContextMenu');
				submenu2 = get(menu2,'Children');
				submenu2 = submenu2(2);
				uimenu(submenu2,'Label',['x' int2str(i)],'Callback',@(source,eventdata) connectEdge(source,eventdata,j,i));
			end
			if(isempty(temp))
				uimenu(menu,'Label','Toggel Conducting State (N)','Callback',@(source,eventdata) toggleConducting(source,eventdata,hPoint));
			else
				uimenu(menu,'Label','Toggel Conducting State (Y)','Callback',@(source,eventdata) toggleConducting(source,eventdata,hPoint));
				hPointChildren(3).MarkerEdgeColor = 'b';
				hPointChildren(3).MarkerFaceColor = 'b';
			end
			addNewPositionCallback(hPoint,@(pos) fcn2(pos,hPoint,handles.gridScale));
		end
		for i = 1:size(Edges,1)
			sourceIdx = find(strcmp(Edges{i,1},VertexTags));
			targetIdx = find(strcmp(Edges{i,2},VertexTags));
			source = handles.Vertex{sourceIdx};
			target = handles.Vertex{targetIdx};

			S = get(source,'Children'); S = S(3);
			T = get(target,'Children'); T = T(3);
			sourcePos = [S.XData S.YData];
			targetPos = [T.XData T.YData];
			Tv = targetPos'- sourcePos';
			Tv_norm = Tv/norm(Tv);
			Nv = [-Tv(2);Tv(1)];
			Nv_norm = Nv/norm(Nv);
			midPoint = mean([sourcePos;targetPos]);
			eps = 1e-2;
			arrow = quiver(sourcePos(1)+0.05*Tv(1)+0.25*eps*Nv_norm(1),sourcePos(2)+0.05*Tv(2)+0.25*eps*Nv_norm(2),Tv(1)+0.25*eps*Nv_norm(1),Tv(2)+0.25*eps*Nv_norm(2),'k');
			label = text(midPoint(1),midPoint(2),'','FontSize',6);
			handles.Edges{i}.text = label;
			handles.Edges{i}.text.String = parseRateString(Edges{i,3});
			handles.Edges{i}.source = source;
			handles.Edges{i}.target = target;
			handles.Edges{i}.arrow = arrow;
			thet = atan(Tv(2)/Tv(1))*180/pi;
			handles.Edges{i}.text.Position(1:2) = midPoint;
			handles.Edges{i}.text.Rotation = thet;
			posHat = [handles.Edges{i}.text.Extent(1)+0.5*handles.Edges{i}.text.Extent(3) handles.Edges{i}.text.Extent(2)+0.5*handles.Edges{i}.text.Extent(4)];
			difPos = posHat - midPoint;
			handles.Edges{i}.text.Position = midPoint-difPos+2.5*eps*Nv_norm';
			addNewPositionCallback(source,@(pos) fcn(pos,arrow,source,target,label,1));
			addNewPositionCallback(target,@(pos) fcn(pos,arrow,source,target,label,-1));
		end
	end

	guidata(hObject,handles);

function toggleConducting(source,eventdata,hPoint)
	ax = get(hPoint,'Parent');
	handles = guidata(ax);
	tag = get(hPoint,'Tag');
	N = find(strcmp(handles.OpenPositions,tag));
	hPointChildren = get(hPoint,'Children');
	if(isempty(N))
		source.Text = 'Toggel Conducting State (Y)';
		handles.OpenPositions{end+1} = tag;
		hPointChildren(3).MarkerEdgeColor = 'b';
		hPointChildren(3).MarkerFaceColor = 'b';
	else
		source.Text = 'Toggel Conducting State (N)';
		handles.OpenPositions = handles.OpenPositions(setdiff(1:length(handles.OpenPositions),N));
		hPointChildren(3).MarkerEdgeColor = 'k';
		hPointChildren(3).MarkerFaceColor = 'k';
	end
	
	guidata(get(ax,'Parent'),handles);


% --------------------------------------------------------------------
function drawpoint_Callback(hObject, eventdata, handles)
	% hObject    handle to drawpoint (see GCBO)
	% eventdata  reserved - to be defined in a future version of MATLAB
	% handles    structure with handles and user data (see GUIDATA)

	N = length(handles.Vertex);
	aliveVertices = setdiff([1:N],handles.deletePoints);
	if(~isempty(handles.deletePoints))
		temp = sort(handles.deletePoints);
		M = temp(1);
		handles.deletePoints = setdiff(handles.deletePoints,M);
	else
		M = N+1;
	end

	hPoint = impoint();
	setString(hPoint,['x' int2str(M)]);
	hPointChildren = get(hPoint,'Children');
	hPointChildren(1).FontWeight = 'bold';
	hPointChildren(1).BackgroundColor = 'none';
	hPointChildren(1).FontSize = 6;
	hPointChildren(2).Visible = 'off';
	hPointChildren(3).Marker = 'square';
	hPointChildren(3).MarkerEdgeColor = 'k';
	hPointChildren(3).MarkerFaceColor = 'k';
	set(hPoint,'Tag',['x' int2str(M)]);

	
	% label = text(hPointChildren(2).XData,hPointChildren(2).YData,['x' int2str(M)],'FontSize',6);
	handles.Vertex{M} = hPoint;
	handles.OpenPositions{M} = 0;

	menu = get(hPoint,'UIContextMenu');
	temp = get(menu,'Children');
	set(temp(1),'MenuSelectedFcn',@(source,eventdata) deletePoint(source,eventdata,hPoint));
	submenu = uimenu(menu,'Label','Add Connection');
	for i = aliveVertices
		uimenu(submenu,'Label',get(handles.Vertex{i},'Tag'),'Callback',@(source,eventdata) connectEdge(source,eventdata,M,i));
		menu2 = get(handles.Vertex{i},'UIContextMenu');
		submenu2 = get(menu2,'Children');
		submenu2 = submenu2(2);
		uimenu(submenu2,'Label',['x' int2str(M)],'Callback',@(source,eventdata) connectEdge(source,eventdata,i,M));
	end

	uimenu(menu,'Label','Toggel Conducting State (N)','Callback',@(source,eventdata) toggleConducting(source,eventdata,hPoint));
	addNewPositionCallback(hPoint,@(pos) fcn2(pos,hPoint,handles.gridScale));
	guidata(hObject,handles);


function fcn2(pos,source,gridScale)
	S = get(source,'Children'); L = S(1); S = S(3);
	tempx = round(S.XData/gridScale)*gridScale;
	tempy = round(S.YData/gridScale)*gridScale;
	difx = S.XData-tempx;
	dify = S.YData-tempy;
	L.Position = L.Position - [difx dify 0];
	S.XData = tempx;
	S.YData = tempy;

function deletePoint(source,eventdata,hPoint)
	ax = get(hPoint,'Parent');
	handles = guidata(ax);
	delTag = get(hPoint,'Tag');
	edgeLedger = handles.uitable1.Data(:,[1,2]);
	[edgesToRemove,~] = find(strcmp(edgeLedger,delTag));
	if(~isempty(edgesToRemove))
		edgesToKeep = setdiff([1:size(edgeLedger,1)],edgesToRemove);
		handles.uitable1.Data = handles.uitable1.Data(edgesToKeep,:);
		for i = 1:length(edgesToRemove)
			delete(handles.Edges{edgesToRemove(i)}.arrow);
			delete(handles.Edges{edgesToRemove(i)}.text);
		end
		handles.Edges = handles.Edges(edgesToKeep);
	end
	aliveVertices = setdiff([1:length(handles.Vertex)],handles.deletePoints);
	for i = aliveVertices
		p = handles.Vertex{i};
		menu = get(get(p,'UIContextMenu'),'Children');
		submenu = get(menu(1),'Children');
		for j = 1:length(submenu)
			if(strcmp(get(submenu(j),'Label'),delTag))
				delete(submenu(j));
			end
		end
		if(strcmp(get(p,'Tag'),delTag))
			handles.deletePoints(end+1) = i;
		end
	end
	N = find(strcmp(handles.OpenPositions,delTag));
	if(~isempty(N))
		handles.OpenPositions = handles.OpenPositions(setdiff(1:length(handles.OpenPositions),N));
	end

	delete(hPoint);
	handles.uitable2.Data = updateUITable2(handles.uitable1,handles.uitable2);
	guidata(get(ax,'Parent'),handles);


function connectEdge(menuSource,eventdata,sourceIdx,targetIdx)
	handles = guidata(menuSource);
	source = handles.Vertex{sourceIdx};
	target = handles.Vertex{targetIdx};
	S = get(source,'Children'); S = S(3);
	T = get(target,'Children'); T = T(3);

	sourcePos = [S.XData S.YData];
	targetPos = [T.XData T.YData];
	Tv = targetPos'- sourcePos';
	Tv_norm = Tv/norm(Tv);
	Nv = [-Tv(2);Tv(1)];
	Nv_norm = Nv/norm(Nv);
	midPoint = mean([sourcePos;targetPos]);
	eps = 1e-2;

	arrow = quiver(sourcePos(1)+0.05*Tv(1)+0.25*eps*Nv_norm(1),sourcePos(2)+0.05*Tv(2)+0.25*eps*Nv_norm(2),Tv(1)+0.25*eps*Nv_norm(1),Tv(2)+0.25*eps*Nv_norm(2),'k');
	
	n = length(handles.Edges)+1;
	handles.Edges{n}.source = source;
	handles.Edges{n}.target = target;
	label = text(midPoint(1),midPoint(2),'','FontSize',6);
	handles.Edges{n}.text = label;
	handles.Edges{n}.arrow = arrow;
	handles.uitable1.Data{n,1} = get(source,'Tag');
	handles.uitable1.Data{n,2} = get(target,'Tag');
	addNewPositionCallback(source,@(pos) fcn(pos,arrow,source,target,label,1));
	addNewPositionCallback(target,@(pos) fcn(pos,arrow,source,target,label,-1));
	guidata(menuSource,handles)
	delete(menuSource);

function fcn(pos,arrow,source,target,label,direction)
	try
		handles = guidata(arrow);
	catch
		return;
	end

	S = get(source,'Children'); S = S(3);
	T = get(target,'Children'); T = T(3);

	if(direction == 1)
		S.XData = round(S.XData/handles.gridScale)*handles.gridScale;
		S.YData = round(S.YData/handles.gridScale)*handles.gridScale;
	else
		T.XData = round(T.XData/handles.gridScale)*handles.gridScale;
		T.YData = round(T.YData/handles.gridScale)*handles.gridScale;
	end

	sourcePos = [S.XData S.YData];
	targetPos = [T.XData T.YData];

	Tv = targetPos'- sourcePos';
	Tv_norm = Tv/norm(Tv);
	Nv = [-Tv(2);Tv(1)];
	Nv_norm = Nv/norm(Nv);

	midPoint = mean([sourcePos;targetPos]);
	eps = 1e-2;

	arrow.XData = sourcePos(1)+0.05*Tv(1)+0.25*eps*Nv_norm(1);
	arrow.YData = sourcePos(2)+0.05*Tv(2)+0.25*eps*Nv_norm(2);
	arrow.UData = Tv(1)+0.25*eps*Nv_norm(1);
	arrow.VData = Tv(2)+0.25*eps*Nv_norm(2);

	
	T = label;
	thet = atan(Tv(2)/Tv(1))*180/pi;
	T.Position(1:2) = midPoint;
	T.Rotation = thet;
	posHat = [T.Extent(1)+0.5*T.Extent(3) T.Extent(2)+0.5*T.Extent(4)];
	difPos = posHat - midPoint;
	T.Position = midPoint-difPos+2.5*eps*Nv_norm';
	guidata(arrow,handles);

% --------------------------------------------------------------------
function cm_axes_Callback(hObject, eventdata, handles)
	% hObject    handle to cm_axes (see GCBO)
	% eventdata  reserved - to be defined in a future version of MATLAB
	% handles    structure with handles and user data (see GUIDATA)


function text = parseRateString(string)
	greekAlphabet = {'alpha';'theta';'tau';'beta';'vartheta';'pi'; ...
	'gamma';'varpi';'phi';'delta';'kappa';'rho';'varphi';'epsilon'; ...
	'lambda';'varrho';'chi';'varepsilon';'mu';'sigma';'psi';'zeta';'nu'; ...
	'varsigma';'omega';'eta';'xi'};
	
	vars = symvar(string);
	repIdcs = [];
	for i = 1:length(vars)
		vars2 = strsplit(vars{i},'_');
		for j = 1:length(vars2)
			idxtemp = find(strcmp(vars2{j},greekAlphabet));
			if(~isempty(idxtemp) && ~sum(idxtemp==repIdcs))
				repIdcs(end+1) = idxtemp;
				string = strrep(string,greekAlphabet{idxtemp},['\' greekAlphabet{idxtemp}]);
			end
		end
	end


	text = strrep(string,'*','');


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2

contents = cellstr(get(hObject,'String'));
fontsize = str2num(contents{get(hObject,'Value')});
for i = 1:length(handles.Vertex)
	child = get(handles.Vertex{i},'Children');
	set(child(1),'FontSize',fontsize);
end
for i = 1:length(handles.Edges)
	handles.Edges{i}.text.FontSize = fontsize;
end
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
