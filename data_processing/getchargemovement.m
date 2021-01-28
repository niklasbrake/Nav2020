function [Charge_Activation,Charge_Inactivation] = getchargemovement(Folder)
% GETCHARGEMOVEMENT creates the file ChargeMovement.mat used in constructtemplate.
% 	Charge = GETCHARGEMOVEMENT(Folder) takes the Folder containing subfolders for 
% 	all the experiments of a single condition. Charge_Activation and Charge_Inactivation
% 	are vectors containing the integral of the current corresponding the different
% 	voltage steps in the activation and inactivation protocols, respectively.
% 
% 	See also constructtemplate.

DIR = dir(Folder);
DIR = DIR(3:end);
Idcs = find([DIR.isdir]);

for j = 1:length(Idcs)
	try
		SubFolder = DIR(Idcs(j)).name;
		load(fullfile(Folder,SubFolder,'activation.mat'))
		Current = activationleakcorrection(Voltage,Current,Epochs);
		temp = Current(:,Epochs(4)+41:Epochs(4)+500)'; % Offset of 41 to avoid integrating capacitive transients.
		temp = temp/max(max(abs(temp)));
		Charge_Activation(j,:) = sum(temp);
	catch
	end
	try
		load(fullfile(Folder,SubFolder,'inactivation.mat'))
		Current = inactivationleakcorrection(Voltage,Current,Epochs);
		temp = Current(:,Epochs(5)+41:Epochs(5)+500)'; % Offset of 41 to avoid integrating capacitive transients.
		temp = temp/max(max(abs(temp)));
		Charge_Inactivation(j,:) = sum(temp);
	catch
	end
end

save(fullfile(Folder,'ChargeMovement.mat'),'Charge_Activation','Charge_Inactivation');