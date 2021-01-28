function responses = inactivationleakcorrection(commands,responses,Epochs)
	% Elementwise rather than mean...
	for i = 1:size(responses,1)
		y1 = median(responses(i,1:Epochs(3)));
		x1 = median(commands(i,1:Epochs(3)));
		y2 = median(responses(i,Epochs(3):Epochs(4)));
		x2 = median(commands(i,Epochs(3):Epochs(4)));
		m = (y2-y1)/(x2-x1);
		C = @(V) (V-x1)*m+y1;
		responses(i,:) = responses(i,:) - C(commands(i,:));
		responses(i,Epochs(4):Epochs(5)) = responses(i,Epochs(4):Epochs(5)) - responses(i,Epochs(5)-100);
		responses(i,Epochs(6):Epochs(6)) = responses(i,Epochs(6):Epochs(6)) - responses(i,Epochs(6)-100);
	end