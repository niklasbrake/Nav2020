function MR_Idx = getMRCircuits(MR_Edges,circuits)
	for e = 1:size(MR_Edges,1)
		E = [',' int2str(MR_Edges(e,1)) ',' int2str(MR_Edges(e,2))  ','];
		ER = [',' int2str(MR_Edges(e,2)) ',' int2str(MR_Edges(e,1))  ','];
		n = []; I = [];
		for i = 1:length(circuits)
			C = [',' char(strjoin(string(circuits{i}),',')) ','];
			idcs = strfind(C,E);
			idcs2 = strfind(C,ER);
			if(~isempty(idcs) || ~isempty(idcs2))
				n(end+1) = length(circuits{i});
				I(end+1) = i;
			end
		end
		[n,idcs] = sort(n);
		I = I(idcs);
		N = unique(n);
		smallLoopN = N(2);
		idcs2 = find(n == smallLoopN);
		MR_Idx(e) = 0;
		noOtherMREdges = ones(length(idcs2),1);
		for i = 1:length(idcs2)
			if(e == 1)
				continue;
			end
			for e2 = 1:e-1
				E2 = [',' int2str(MR_Edges(e2,1)) ',' int2str(MR_Edges(e2,2))  ','];
				E2R = [',' int2str(MR_Edges(e2,2)) ',' int2str(MR_Edges(e2,1))  ','];
				C = [',' char(strjoin(string(circuits{I(idcs2(i))}),',')) ','];
				if(~isempty(strfind(C,E2)) || ~isempty(strfind(C,E2R)))
					noOtherMREdges(i) = 0;
				end
			end
		end
		j = find(noOtherMREdges);
		MR_Idx(e) = I(idcs2(j(1)));
		circuits{MR_Idx(e)};
	end
end


