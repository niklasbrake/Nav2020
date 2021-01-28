function responses = recoveryleakcorrection(commands,responses,Epochs,lag)
	% Elementwise rather than mean...
	C1 = zeros(size(responses));
	for i = 1:size(responses,1)
		y1 = mean(responses(i,1:Epochs(3)));
		x1 = mean(commands(i,1:Epochs(3)));
		y2 = mean(responses(i,Epochs(3):Epochs(4)));
		x2 = mean(commands(i,Epochs(3):Epochs(4)));
		m = (y2-y1)/(x2-x1);
		C = @(V) (V-x1)*m+y1;
		responses(i,:) = responses(i,:) - C(commands(i,:));
		C1(i,Epochs(4):Epochs(5)) = responses(i,Epochs(5)-100);
		C1(i,Epochs(6)+100*lag*(i-1):Epochs(7)+100*lag*(i-1)) = responses(i,Epochs(7)+100*lag*(i-1)-100);
		responses(i,:) = responses(i,:) - C1(i,:); 
	end