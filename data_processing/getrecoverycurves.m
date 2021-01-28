function [a,b,RecovFit] = getrecoverycurves(Folder,OverideLoad)
% GETRECOVERYCURVES creates the file Recovery_Curves.mat used in constructtemplate.
% 	[a,b,RecovFit] = GETRECOVERYCURVES(Folder) takes the Folder containing subfolders for 
% 	all the experiments of a single condition. (a,b) are the lag time and fraction recovered,
% 	recpetively. RecovFit is an n x 3 matrix, where each row is a different replicate and
% 	each column is A, gamma1, and gamma2, respectively, fitted parameters of the function 
% 		i = (1-exp(-t*gamma1))*A + (1-exp(-t*gamma2))*(1-A)
% 
% 	See also constructtemplate, FitBiExponential.

% If Recovery_Curves already exists, and override == 0, that is, you don't want to recompute it,
% Then load the file and return the contents.
	if(nargin == 1)
		OverideLoad = 0;
	end
	if(exist(fullfile(Folder,'Recovery_Curves.mat')) && ~OverideLoad)
		load(fullfile(Folder,'Recovery_Curves.mat'));
		return;
	end

	F = dir(Folder); F = F(3:end);

	a = [];
	b = [];
	for i = 1:length(F)
		% Since activation protocol is run first, all experiments that have been processed
		% must have the activation.mat file
		if(~exist(fullfile(Folder,F(i).name,'activation.mat')))
			try
				output = python('pyloadABF.py',fullfile(Folder,F(i).name));
				temp = split(output,'/');
				flag = double(temp{2});
				if(flag == 48)
					continue;
				end
			catch
				continue;
			end
		end

		try
			load(fullfile(Folder,F(i).name,'recovery1ms.mat'))
			Current = recoveryleakcorrection(Voltage,Current,Epochs,1);
			Current2 = Current;
			pre = min(Current2(:,41775:4.5e4),[],2);
			for j = 1:size(Current,1)
				% Current2(j,Epochs(6)+100*(j-1):Epochs(6)+100*(j-1)+36) = 0;
				b(i,j) = min(Current2(j,Epochs(6)+100*(j-1)+35:Epochs(6)+100*(j-1)+150))/pre(j);
				a(i,j) = j;
			end
			% [temp1 temp2] = findpeaks(min(Current2(:,Epochs(6)+ 50:Epochs(6)+100*size(Current,1)+50),[],1)/pre, ...
			% 	'MinPeakDistance',80);
			% b(2*i-1,1:25) = temp1;
			% a(2*i-1,1:25) = temp2 + double(Epochs(6)) +  50 - double(Epochs(5));
		catch
			disp([F(i).name ': no recovery1ms data.'])
			continue
		end	

		try
			load(fullfile(Folder,F(i).name,'recovery10ms.mat'))
			Current = recoveryleakcorrection(Voltage,Current,Epochs,10);
			Current2 = Current;
			pre = min(Current2(:,41925:4.5e4),[],2);
			for j = 1:size(Current,1)
				% Current2(j,Epochs(6)+100*10*(j-1):Epochs(6)+100*10*(j-1)+36) = 0;
				% Current2(j,Epochs(7)+100*10*(j-1):Epochs(7)+100*10*(j-1)+36) = 0;
				b(i,25+j) = min(Current2(j,Epochs(6)+1000*(j-1)+35:Epochs(6)+1000*(j-1)+150))/pre(j);
				a(i,25+j) = 10*j;
			end
			% [temp1 temp2] = findpeaks(min(Current2(:,Epochs(6)+ 50:Epochs(6)+100*10*size(Current,1)+50),[],1)/pre, ...
			% 	'MinPeakDistance',980);
			% b(2*i-1,26:40) = temp1;
			% a(2*i-1,26:40) = temp2 + double(Epochs(6)) + 50 - double(Epochs(5));
		catch
			disp([F(i).name ': no recovery10ms data.'])
			continue
		end
	end

	a = a(b(:,end)~=0,:);
	b = b(b(:,end)~=0,:);

	[~,idcs] = sort(a(1,:));
	a = a(:,idcs);
	b = b(:,idcs);

	% a = repmat(sort([1:25 10:10:150]),[size(b,1) 1]);

	for i = 1:size(a,1)
		[fitresult1, gof] = FitBiExponential(a(i,:),b(i,:));
		RecovFit(i,:) = coeffvalues(fitresult1);
	end

	save(fullfile(Folder,'Recovery_Curves.mat'),'a','b','RecovFit');