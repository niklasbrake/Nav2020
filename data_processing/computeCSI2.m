function [OSI,CSI,Inact,O] = computeCSI(Folder)
% COMPUTECSI calculates the fraction of steady-state inactivation
% for closed-states and open-states using the inactivation protocol.
% 	[OSI,CSI,Inact] = computeCSI(Folder) takes the Folder	containing 
% 	subfolders for all the experiments of a single condition and saves the
% 	output in the file CSI.mat. Outputs are the fraction of open-state
% 	inactivation (OSI), closed-state inactivation (CSI), total steady-
% 	state inactivation (Inact) for each voltage step. These are calculated 
% 	using the formulas presented in Armtrong C. (2006) "Na channel 
% 	inactivation from open and closed states". PNAS. This function
% 	also plots the output. 


F = dir(Folder); F = F(3:end);
idcsE = 0;

for i = 1:length(F)
	if(exist(fullfile(Folder,F(i).name,'inactivation.mat'))~=2)
		idcsE=idcsE+1;
		continue;
	end
	if(exist(fullfile(Folder,F(i).name,'activation.mat'))~=2)
		idcsE=idcsE+1;
		continue;
	end
	load(fullfile(Folder,F(i).name,'activation.mat')); % Load inactivation data
	V = mean(Voltage(:,Epochs(4):Epochs(4)+500)');
	Current = activationleakcorrection(Voltage,Current,Epochs);
	X = Current(:,Epochs(4):Epochs(4)+500)';
	[~,idxMax] = min(X);
	for j = 1:40
		id = idxMax(j);
		while id>2 && X(id,j)<X(idxMax(j),j)*0.05
			id = id-1;
		end
		idxStart(j) = id;
	end
	idxStart = ceil(median(idxStart(find(and(V>-60,V<55)))));
	X = X(idxStart:end,:);
	[Imax,idxMax] = min(X);

	[~,idx] = min(Imax);
	Vlo = V(idx)+10;
	if(Vlo>50)
		idcsE=idcsE+1;
		continue;
	end
	Vhi = 55;
	idcsLinear = and(V>Vlo,V<Vhi);
	FT = fitlm(V(idcsLinear),Imax(idcsLinear));
	ERev = -FT.Coefficients.Estimate(1)/FT.Coefficients.Estimate(2);
	if(ERev <= 40 || ERev >= 80)
		disp(['Error: Erev for ' F(i).name  ' out of range.']);
		idcsE = idcsE + 1;
		continue;
	end

	load(fullfile(Folder,F(i).name,'inactivation.mat')); % Load inactivation data
	Current = inactivationleakcorrection(Voltage,Current,Epochs); % Correct for leak current
	Current = Current(7:end,:);
	I_V = mean(Voltage(7:end,Epochs(4):Epochs(4)+500),2); % Get the voltage of pre pulse
	% preCurrent = Current(:,Epochs(4)+20:Epochs(4)+20+16e2); % Get the current in response to the pre-pulse (20 used to get rid of capactive transient)
	preCurrent = Current(:,Epochs(4)+idxStart:Epochs(5)-50); % Get the current in response to the pre-pulse (20 used to get rid of capactive transient)
	preCurrent = preCurrent-median(preCurrent,2); % Responses are baselined to return to 0
	preG = preCurrent ./(I_V-ERev); % Correct for reversal potential to estimate open probability of channels
	postCurrent = Current(:,Epochs(5)+idxStart:Epochs(6)); % Get currents in response to test-pulse
	I_I = min(postCurrent,[],2) - min(Current(:,Epochs(5)+4000:Epochs(6)),[],2); % get peak current minimum (they should all be negative)
	FT = FitBoltzman2(I_V,I_I,-60,10,min(I_I));
	if(max(abs(I_I./FT.Gmx))>2)
		disp(['Error: SSI scale for ' F(i).name  ' out of range.']);
		idcsE = idcsE + 1;
		continue;
	end
	Inact(i-idcsE,:) = 1-I_I./FT.Gmx; % Get the fraction of unavailable channels.

	OSI(i-idcsE,:) = sum(preG,2);
	FT = FitBoltzmanCurve(I_V,OSI(i-idcsE,:)',-10,-5,max(OSI(i-idcsE,:)));
	OSI(i-idcsE,:)=OSI(i-idcsE,:)/FT.Gmx;

	CSI(i-idcsE,:) = Inact(i-idcsE,:) - OSI(i-idcsE,:); % CSI are channels is amount of unavailble channels that can't be accounted 
														% for by the pre-pulse current
end
V = floor(I_V);
save(fullfile(Folder,'CSI.mat'),'CSI','OSI','Inact','V');

figure
errorbar(V,mean(CSI),std(CSI)/sqrt(size(CSI,1)),'LineWidth',0.75)
hold on;
errorbar(V,mean(OSI),std(OSI)/sqrt(size(CSI,1)),'LineWidth',0.75)
errorbar(V,mean(Inact),std(Inact)/sqrt(size(CSI,1)),'LineWidth',0.75)
% legend({'CSI','OSI','Inactivation'},'Location','northwest');
box off;
xlim([-120,-5])
set(gca,'TickDir','out')
set(gca,'FontSize',12);
% set(gca,'LineWidth',0.75);
xlabel('Holding Potential');
