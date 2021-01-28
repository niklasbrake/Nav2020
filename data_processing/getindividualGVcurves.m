function [InactFit ActFit I_V I_I A_V A_I] = getindividualGVcurves(Folder,override)
% GETINDIVIDUALGVCURVES extracts the I-V curves for activation and inactivation protocols 
% and fits Boltzmans to get the V_1/2 and k of activation and steady-state inactvation.
% 	[InactFit ActFit I_V I_I A_V A_I] = GETINDIVIDUALGVCURVES(Folder) takes the Folder
% 	containing subfolders for all the experiments of a single condition. Saves the output
% 	variables in the file GV_Curves.mat. 
% 	InactFit and ActFit are n x 4 matrices, where each row is another replicate and the columns 
% 	are the fitted reversal potential, peak conductance	and the k and V_1/2 parameters of a 
% 	Boltzman fit to the steady-state inactivation and activaiton data, respectively. I_V, I_I, 
% 	A_V, A_I are the inactivation and activation step voltages and peak currents. That is, 
% 	plot(A_V',A_I') will plot the I-V curve of every experiment.
% 
% 	For the assumed file structure, see the <a href="matlab:web('https://github.com/niklasbrake/Brake-et-al-2019/wiki')">documentation</a>.
% 
% 	** The function calls a python function which requires the packages: pyabf, numpy, mat4py.



% If GV_Curves already exists, and override == 0, that is, you don't want to recompute it,
% Then load the file and return the contents.
	if(nargin == 1)
		override = 0;
	end
	if(exist(fullfile(Folder,'GV_Curves.mat')) && ~override)
		load(fullfile(Folder,'GV_Curves.mat'));
		return;
	end

% Initialize stuff
	F = dir(Folder); F = F([F(:).isdir]); F = F(3:end);
	VClamp = [];
	inactivationResponse = [];
	exResponse = [];
	count1 = 0;
	count2 = 0;

	ActFit = [];
	InactFit = [];

	I_V = [];
	I_I = [];

	A_V = [];
	A_I = [];

% For each subfolder in Folder
	for i = 1:length(F)
	% We start each experiment with the activation protocol, so if there is not activation protocol
	% file, this either isn't an experiment folder, or the data hasn't been extracted into a *.mat file
		if(~exist(fullfile(Folder,F(i).name,'activation.mat')))
			try
				% There was a *.abf extraction function on MATLAB file exchange but it didn't work. Oh well, I found
				% a Python implementation that works, so I wrote a python file (pyloadABF) and call it from matlab
				% This requires the packages: os, pyabf, numpy, mat4py, sys.
				output = python('pyloadABF.py',fullfile(Folder,F(i).name));
			catch
				disp('error')
				% continue;
			end
		end

	% Now there should defintely be an activation.mat file if this is an experiment folder. If there are any
	% issues, skip this folder
		try
			load(fullfile(Folder,F(i).name,'activation.mat'));
			Current = activationleakcorrection(Voltage,Current,Epochs); % corrects for the leak current
			V = mean(Voltage(:,Epochs(4):Epochs(5)),2)'; % The voltage step occured bewteen Epochs 4 and 5-60.
			responses = Current(:,Epochs(4)+40:Epochs(5)-1); % See documentation for note on the offset
			[~,iMax] = min(responses');
			idcs = sub2ind(size(responses),1:size(responses,1),iMax);
			I = responses(idcs); % Get the IV curve for each voltage protocol
			idcs = find(V>30);
			FT = fitlm(V(idcs),I(idcs));
			ERev = -FT.Coefficients{1,1}./FT.Coefficients{2,1};
			G = I./(V-ERev);
			FT = FitBoltzmanCurve(V(V<40),G(V<40),median(G(30:end)),-10,-15);
			[FT2,gof] = FitBoltzmanCurve2(V(V<40),G(V<40)/FT.Gmx,-10,-15);
			ActFit(i-count2,1) = ERev;
			ActFit(i-count2,2) = FT.Gmx;
			ActFit(i-count2,3) = FT2.k;
			ActFit(i-count2,4) = FT2.v50;
			if(length(V)<40)
				V(end+1:40) = nan;
				I(end+1:40) = nan;
			end
			A_V(i-count2,:) = V(1:40);
			A_I(i-count2,:) = I(1:40);
		catch
			disp([F(i).name ': different activation protocol']);
			count2 = count2 +1;
		end

		% Repeat with the inactivation file. The voltages are always below the reversal potential of sodium, so
		% no shananigans when calculating the IV curve.
		try
			load(fullfile(Folder,F(i).name,'inactivation.mat'))
			Current = inactivationleakcorrection(Voltage,Current,Epochs);
			% Current = capactianceCorrection(Voltage,Current,Epochs);
			I_V(i-count1,:) = mean(Voltage(:,Epochs(4):Epochs(5)),2);
			responses = Current(:,Epochs(5)+40:Epochs(6));
			I_I(i-count1,:) = min(responses,[],2);
			v50Pre = interp1(I_I(i-count1,:),I_V(i-count1,:),0.5*min(I_I(i-count1,:)));
			FT = FitBoltzman2(I_V(i-count1,:),I_I(i-count1,:),v50Pre,10,-1);
			InactFit(i-count1,1) = 60;
			InactFit(i-count1,2) = FT.Gmx;
			InactFit(i-count1,3) = FT.k;
			InactFit(i-count1,4) = FT.v50;
		catch
			disp([F(i).name ': no inactivation data.'])
			count1 = count1 + 1;
		end

	end

	save(fullfile(Folder,'GV_Curves.mat'),'InactFit','ActFit','I_V','I_I','A_V','A_I');