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
% for i = 1:10
	try
		load(fullfile(Folder,F(i).name,'activation.mat')); % Load inactivation data
		V = mean(Voltage(:,Epochs(4):Epochs(4)+500)');
		Current = activationleakcorrection(Voltage,Current,Epochs);
		X = Current(:,Epochs(4)+30:Epochs(4)+530)';
		[~,idxMax] = max(abs(X));
		for j = 1:size(X,2)
			Imax(j) = X(idxMax(j),j);
		end
		[~,idx] = min(Imax);
		Vlo = V(idx)+10;
		Vhi = 65;
		idcsLinear = and(V>Vlo,V<Vhi);
		FT = fitlm(V(idcsLinear),Imax(idcsLinear));
		ERev = -FT.Coefficients.Estimate(1)/FT.Coefficients.Estimate(2);
		if(ERev <= 50 || ERev >= 70)
			idcsE = idcsE + 1;
			continue;
		end

		load(fullfile(Folder,F(i).name,'inactivation.mat')); % Load inactivation data
		Current = inactivationleakcorrection(Voltage,Current,Epochs); % Correct for leak current
		V = mean(Voltage(:,Epochs(4):Epochs(4)+500),2); % Get the voltage of pre pulse
		% preCurrent = Current(:,Epochs(4)+42:Epochs(4)+42+16e2); % Get the current in response to the pre-pulse (42 used to get rid of capactive transient)
		preCurrent = Current(:,Epochs(4)+42:Epochs(5)-50); % Get the current in response to the pre-pulse (42 used to get rid of capactive transient)
		preCurrent = preCurrent-median(preCurrent,2); % Responses are baselined to return to 0
		preG = preCurrent ./(V-ERev); % Correct for reversal potential to estimate open probability of channels
		% peakG = max(abs(preG(end,:))); % Get peak conductance
		% O = preCurrent./peakG; % Normalize to get a probablility
		% OSI(i-idcsE,:) = sum(preG,2)/max(sum(preG,2)); % OSI is proportional to the total occupancy of the open state
		OSI(i-idcsE,:) = sum(preG,2);
		FT = FitBoltzmanCurve(V,OSI(i-idcsE,:)',-10,-5,max(OSI(i-idcsE,:)));
		OSI(i-idcsE,:)=OSI(i-idcsE,:)/FT.Gmx;

		I_V = mean(Voltage(:,Epochs(4):Epochs(5)),2); % Get the voltages of the test-pulses
		I_V = round(I_V/5)*5; % Round 
		postCurrent = Current(:,Epochs(5)+42:Epochs(6)); % Get currents in response to test-pulse
		I_I = min(postCurrent,[],2); % get peak current minimum (they should all be negative)
		Inact(i-idcsE,:) = 1 - I_I/min(I_I); % Get the fraction of unavailable channels.

		CSI(i-idcsE,:) = Inact(i-idcsE,:) - OSI(i-idcsE,:); % CSI are channels is amount of unavailble channels that can't be accounted 
															% for by the pre-pulse current
	catch
		idcsE = idcsE + 1;
	end
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
set(gca,'LineWidth',0.75);
xlabel('Holding Potential');
