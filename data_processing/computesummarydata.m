function computesummarydata(Folder)
% COMPUTESUMMARYDATA generates the files necessary for constructtemplate.
% 	COMPUTESUMMARYDATA(Folder) produces the following files from a folde
% 	containing subfolders for each experiment that contains the files
% 	activation.abf, inactivation.abf, recovery1ms.abf, and recovery10ms.abf.
% 		GV_Curves.mat
% 		Recovery_Curves.mat
% 		AverageResponseCurves.mat
% 		PeakTimes.mat
% 		ChargeMovement.mat
% 
% 	For more information, see <a href="matlab:web('https://github.com/niklasbrake/Brake-et-al-2019/wiki')">the reference page</a>.
% 
% 	See also getindividualGVcurves, getrecoverycurves, getresponsecurves, getpeaktimes, getchargemovement, constructtemplate.

% Build GV_Curves.mat
	getindividualGVcurves(Folder);
% Build Recovery_Curves.mat
	getrecoverycurves(Folder,1);
% Build AverageResponseCurves.mat
	getresponsecurves(Folder);
% Build PeakTimes.mat
	getpeaktimes(Folder);
% Build ChargeMovement.mat
	getchargemovement(Folder);
% Build template for fitting algorithm
	Template = constructtemplate(Folder)
	save(fullfile(Folder,'Template.mat'),'Template');