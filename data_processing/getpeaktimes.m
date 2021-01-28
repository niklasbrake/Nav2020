function [ActPeak,InactPeak] = getpeaktimes(Folder);
% GETPEAKTIMES creates the file PeakTimes.mat used in constructtemplate.
% 	[ActPeak,InactPeak] = GETPEAKTIMES(Folder) takes the Folder
% 	containing subfolders for all the experiments of a single condition. Saves the output
% 	in the file PeakTimes.mat. ActPeak and InactPeak are n x m  matrices where each 
% 	row corresponds to a replicate and each column correponds to the time between
% 	voltage step and peak conductance following the activation and inactivation protocol,
%  	respectively. The first column is thereform the time to peak with the most negative
% 	voltage step. If the detected peak is less than 5% of the max current, time to peak
% 	is set to 500 (5 ms).

folderList = dir(Folder);
folderList = folderList(3:end);
Idcs = find([folderList.isdir]);

missing_Counter = 0;
for j = 1:length(Idcs)
	folder2 = folderList(Idcs(j)).name;
	load(fullfile(Folder,folder2,'activation.mat'))
	Current = activationleakcorrection(Voltage,Current,Epochs);
	temp = Current(1:34,Epochs(4)+51:Epochs(4)+550)';
	temp = temp/max(max(abs(temp)));
	[x I] = max(abs(temp));
	I(x<0.05) = 450;
	ActPeak(j,:) = I+50;
	try
		load(fullfile(Folder,folder2,'inactivation.mat'))
	catch
		missing_Counter = missing_Counter + 1;
		continue;
	end
	Current = inactivationleakcorrection(Voltage,Current,Epochs);
	temp = Current(1:32,Epochs(5)+51:Epochs(5)+550)';
	temp = temp/max(max(abs(temp)));
	[x I] = max(abs(temp));
	I(x<0.05) = 550;
	InactPeak(j-missing_Counter,:) = I+50;
end

save(fullfile(Folder,'PeakTimes.mat'),'ActPeak','InactPeak')