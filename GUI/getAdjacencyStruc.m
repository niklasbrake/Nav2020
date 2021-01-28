function G = getAdjacencyStruc(E)
	V = unique(E);

	G = repmat({[]},[length(V),1]);
	for i = 1:length(V)
		N = find(E(:,1) == i);
		for j = 1:length(N)
			G{i} = [G{i} E(N(j),2)];
		end
	end