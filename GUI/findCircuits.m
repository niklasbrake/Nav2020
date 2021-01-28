function circuits = findCircuits(G)
	n = length(G);
	s = 1;
	B = repmat({[]},[n,1]);
	circuits = {[]};
	while s < n
		F = getInducedGraph(G,[s:n]);
		connnectedV = find(cellfun(@(v) ~isempty(v),F));
		if(~isempty(connnectedV))
			s = connnectedV(1);
			for i = connnectedV';
				blocked(i) = 0;
				B{i} = [];
			end
			[f,A,B,blocked,vStack,circuits] = CIRCUIT(s,F,B,blocked,[],s,circuits);
			s = s + 1;
		else
			s = n;
		end
	end
end


function [f,A,B,blocked,vStack,circuits] = CIRCUIT(v,A,B,blocked,vStack,s,circuits)
	f = 0;
	vStack = [vStack,v];
	blocked(v) = 1;
	Av = A{v};
	for i = 1:length(Av)
		w = Av(i);
		if(w == s)
			circuits{end+1} = [vStack s];
			f = 1;
		elseif(~blocked(w))
			[fw,A,B,blocked,vStack,circuits] = CIRCUIT(w,A,B,blocked,vStack,s,circuits);
			if(fw)
				f = 1;
			end
		end
	end
	if(f)
		[B,blocked] = UNBLOCK(v,B,blocked);
	else
		for i = 1:length(Av)
			w = Av(i);
			if(isempty(find(v==B{w})))
				B{w}(end+1) = v;
			end
		end
	end
	vStack = vStack([1:end-1]);
end




function [B,blocked] = UNBLOCK(u,B,blocked)
	Bu = B{u};
	blocked(u) = 0;
	while(length(Bu)>0)
		w = Bu(1);
		Bu = Bu(2:end);
		if(blocked(w))
			[B,blocked] = UNBLOCK(w,B,blocked);
		end
	end
	B{u} = Bu;
end


function F = getInducedGraph(G,W)
	notW = setdiff([1:length(G)],W);
	F = repmat({[]},[length(G) 1]);
	for i = 1:length(W)
		F{W(i)} = setdiff(G{W(i)},notW);
	end
end



function getStrongComponent(G)
	N = length(G);
	preOrderNum = zeros(N,1);
	for v = 1:N
		if(preOrderNum(v)==0)
			[stack1,stack2,preOrderNum,C,CC] = findCC(G,v,[],[],preOrderNum,1,zeros(N,1));
		end
	end
end

function [stack1,stack2,preOrderNum,C,CC] = findCC(G,v,stack1,stack2,preOrderNum,C,CC)
	preOrderNum(s) = C;
	C = C + 1;
	stack1 = [stack1 v];
	stack2 = [stack2 v];
	Ns = G{v};
	for i= 1:length(Ns)
		w = Ns(i);
		if(preOrderNum(w) == 0)
			[stack1,stack2,preOrderNum,C,CC] = findCC(G,v,stack1,stack2,preOrderNum,C,CC)
		else
			if(CC(w) == 0)
				while(stack2(end) >= preOrderNum(w))
					stack2 = stack2(1:end-1);
				end
			end
		end
	end
	if(stack2(end) == v)
		newComp = max(CC)+1;
		while(stack1(end)~=v)
			CC(stack1(end)) = newComp;
			stack1 = stack1(1:end-1);
		end
		stack2 = stack2(1:end-1);
	end
end