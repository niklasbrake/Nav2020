function params = constructModelCode(EdgeTable,ParameterTable,fileName,filePath,OpenPositions)

	Vertices = unique(EdgeTable(:,1));
	temp = cellfun(@(x) str2num(x(2:end)),Vertices,'UniformOutput',false);
	[~,I] = sort([temp{:}]);
	Vertices = Vertices(I);

	% Construct Edges table
	V = length(Vertices);
	E = size(EdgeTable,1);
	for i = 1:E
		Edges(i,1) = find(strcmp(EdgeTable{i,1},Vertices));
		Edges(i,2) = find(strcmp(EdgeTable{i,2},Vertices));
	end

	for i = 1:length(OpenPositions)
		O{i} = int2str(sort(find(strcmp(Vertices,OpenPositions{i}))));
	end
	O = O(cellfun(@(x) ~isempty(x),O));

	% Deal with microscopic reversibility
	G = repmat({[]},[V 1]);
	for i = 1:V
		Idcs = find(strcmp(EdgeTable(:,1),Vertices{i}));
		for j = 1:length(Idcs)
			idcs2 = find(strcmp(EdgeTable{Idcs(j),2},Vertices));
			G{i} = [G{i} idcs2];
		end
	end
	circuits = findCircuits(G);

	domainIdcs = [];
	domainVars = {};
	MRIdcs = [];
	for j = 1:size(ParameterTable,1)	
		if(strcmp(ParameterTable{j,2},'independent'))
			domainIdcs(end+1) = j;
			domainVars{end+1} = ParameterTable{j,1};
		elseif(strcmp(ParameterTable{j,2}(1:2),'MR'))
			MRIdcs(end+1) = j;
		end
	end

	funCallStr = '(';
	for k = 1:length(domainVars)
		if(k == 1)
			funCallStr = [funCallStr domainVars{k}];
		else
			funCallStr = [funCallStr ',' domainVars{k}];
		end
	end
	funCallStr = [funCallStr ')'];

	Q = cell(V);
	MR_Edges = [];
	MROrder = [];
	for i = 1:size(EdgeTable,1)
		S = str2sym(EdgeTable{i,3});
		vars = symvar(S);
		for j = 1:length(vars)
			idxtemp = find(strcmp(ParameterTable(:,1),char(vars(j))));
			if(~isempty(idxtemp) && isempty(find(domainIdcs==idxtemp)) && isempty(find(MRIdcs==idxtemp)))
				S = subs(S,vars(j),str2sym(ParameterTable{idxtemp,2}));
			end
			if(~isempty(find(MRIdcs==idxtemp)))
				MR_Edges(end+1,:) = Edges(i,[1:2]);
				MROrder(end+1) = str2num(ParameterTable{idxtemp,2}(3:end));
			end
		end
		vars = symvar(S);
		for j = 1:length(vars)
			idxtemp = find(strcmp(ParameterTable(:,1),char(vars(j))));
			if(~isempty(idxtemp))
				repValue = str2num(ParameterTable{idxtemp,2});
				if(~isempty(repValue))
					S = subs(S,vars(j),repValue);
				end
			end
		end
		theta{i} = char(vpa(S,8));
		Q{Edges(i,2),Edges(i,1)} = ['@' funCallStr ' ' theta{i}];
	end
	[~,SortI] = sort(MROrder,'descend');
	MR_Edges = MR_Edges(SortI,:);
    if(length(MR_Edges)>0)
        MR_Idx = getMRCircuits(MR_Edges,circuits);
    else
        MR_Idx = [];
    end
	for i = length(MR_Idx):-1:1
		swithIJ = 0;
		C = circuits{MR_Idx(i)};
		I = [];
		for j = 1:length(C)-1
			if(C(j:j+1)==MR_Edges(i,:))
				i_rep = find(sum(Edges == C(j:j+1),2)==2);
				continue;
			end
			I(end+1) = find(sum(Edges == C(j:j+1),2)==2);
		end
		CR = fliplr(C);
		J = [];
		for j = 1:length(CR)-1
			if(CR(j:j+1)==MR_Edges(i,:))
				i_rep = find(sum(Edges == CR(j:j+1),2)==2);
				swithIJ = 1;
				continue;
			end
			J(end+1) = find(sum(Edges == CR(j:j+1),2)==2);
		end
		if(swithIJ)
			MRString = ['(' strjoin(theta(I),'*') '/ (' strjoin(theta(J),'*') '))'];
		else
			MRString = ['(' strjoin(theta(J),'*') '/ (' strjoin(theta(I),'*') '))'];
		end
		theta{i_rep} = MRString;
		Q{Edges(i_rep,2),Edges(i_rep,1)} = ['@' funCallStr ' ' theta{i_rep}];
	end


	funPartN = length(funCallStr) + 2;
	for i = 1:V
		outStr = [];
		for j = 1:V
			if(~isempty(Q{j,i}))
				outStr = [outStr '-'  Q{j,i}(funPartN+1:end)];
			end
		end
		Q{i,i} = ['@' funCallStr ' ' outStr];
	end

	params = symvar(strjoin(theta,'+'));
	idxRemove = [];
	for i = 1:length(domainVars)
		temp = find(strcmp(params,domainVars{i}));
		if(~isempty(temp))
			idxRemove(end+1) = temp;
		end
	end
	idcs = setdiff([1:length(params)],idxRemove);
	params = params(idcs);


	% Write to file
	fid = fopen(fullfile(filePath,fileName),'w');
	[~,funName] = fileparts(fileName);


	addpath(fileparts(mfilename('fullpath')));

	fprintf(fid,['function [Q,OpenPositions,P] = ' funName '(Params) \n']);
	fprintf(fid,['%% ' funName ' Defined with the DrawModel GUI and programmatically generated \n%% through constructModelCode. \n']);
	fprintf(fid,['%% \t [Q,OpenPositions,P] = ' funName '(Params) Generate transition matrix Q parameterized by\n%%\t input Params (length=' int2str(length(params)) '). \n']);

	N = [1];
	paramText = params{N};
	while(length(paramText)<40 && N(end) < length(params));
		N = [N length(N)+1];
		paramText = strjoin(params(N),', ');
	end
	fprintf(fid,['%% \t Parameter order: ' paramText]);
	while(N(end) < length(params))
		N = [N(end)+1];
		paramText = params{N};
		while(length(paramText)<35 && N(end) < length(params));
			N = [N N(end)+1];
			paramText = strjoin(params(N),', ');
		end
		fprintf(fid,['\n%% \t \t ' paramText]);
	end
	fprintf(fid,['.\n']);

	
	fprintf(fid,'%% \n');
	fprintf(fid,'%% See also constructModelCode, DrawModel. \n');

	
	for i = 1:length(params)
		Str = [params{i} ' = Params(' int2str(i) ');\n'];
		fprintf(fid,Str);
	end
	
	Str = [];
	fprintf(fid,['preQ = repmat({@(V) 0},[' int2str(V) ' ' int2str(V) ']);\n']);
	for i = 1:V
		for j = 1:V
			if(~isempty(Q{j,i}))
				Str = ['preQ{' int2str(j) ',' int2str(i) '} = ' Q{j,i} ';\n'];
				fprintf(fid,Str);
			end
		end
	end

	Str = ['Q = @(v) cellfun(@(f)f(v),preQ);\n'];
	fprintf(fid,Str);

	Str = ['OpenPositions = [' strjoin(O,',') ']; \n'];
	fprintf(fid,Str);


	for i = 1:length(params)
		Str = ['P.' params{i} ' = Params(' int2str(i) ');\n'];
		fprintf(fid,Str);
	end

	fclose(fid);