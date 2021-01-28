function Template = constructtemplate(Folder)
% CONSTRUCTTEMPLATE returns a struct of template data used for fitting.
% 	Template = CONSTRUCTTEMPLATE(folder). folder is a path to a folder that
% 	must contain the following files:
% 		GV_Curves.mat
% 		ChargeMovement.mat
% 		Recovery_Curves.mat
% 		AverageResponseCurves.mat
% 		PeakTimes.mat
% 	Template is a struct that is used for fitting.
% 
% 	For more information, see <a href="matlab:web('https://github.com/niklasbrake/Brake-et-al-2019/wiki')">the reference page</a>.
% 		
% See also computesummarydata.


%%% All features have 3 values: mean (m), standard deviation (s), number of samples (N) %%%

%%% Summary data quantifications %%%
	load(fullfile(Folder,'GV_Curves.mat'));
% Activaiton. Boltzman fit parameters
	Template.v50A.m = mean(ActFit(:,4));
	Template.v50A.s = std(ActFit(:,4));
	Template.v50A.N = length(ActFit(:,4));
	Template.kA.m = mean(ActFit(:,3));
	Template.kA.s = std(ActFit(:,3));
	Template.kA.N = length(ActFit(:,3));
	Template.AV = mean(A_V(:,1:34));
	Template.ERev = mean(ActFit(:,1));
% Activation. I-V curve
	A_Curve = A_I'./max(abs(A_I')); A_Curve = A_Curve(1:34,:);
	Template.Activation.Voltages =  mean(A_V(:,1:34));
	Template.Activation.m = mean(A_Curve,2);
	Template.Activation.s = std(A_Curve,[],2);
	Template.Activation.N = size(A_Curve,2);
% Activation. G-V curve
	Template.Activation.G.m = mean(ActFit(:,2).*A_I(:,1:34)./(A_V(:,1:34)-ActFit(:,1)));
	Template.Activation.G.s = std(ActFit(:,2).*A_I(:,1:34)./(A_V(:,1:34)-ActFit(:,1)));
	Template.Activation.G.N =  size(A_Curve,2);
% Inactivation. Boltzman fit parameters.
	Template.v50I.m = mean(InactFit(:,4));
	Template.v50I.s = std(InactFit(:,4));
	Template.v50I.N = length(InactFit(:,4));
	Template.kI.m = mean(InactFit(:,3));
	Template.kI.s = std(InactFit(:,3));
	Template.kI.N = length(InactFit(:,3));
	Template.IV = mean(I_V);
% Inactivation. I-V curve
	I_Curve = I_I'./max(abs(I_I')); I_Curve = I_Curve(1:32,:);
	Template.Inactivation.Voltages = mean(I_V);
	Template.Inactivation.m = mean(I_Curve,2);
	Template.Inactivation.s = std(I_Curve,[],2);
	Template.Inactivation.N = size(I_Curve,2);
% Recovery. Biexponential fit parameters (gamma = 1/tau)
	load(fullfile(Folder,'Recovery_Curves.mat'));
	Template.gamma1.m = mean(RecovFit(:,2));
	Template.gamma1.s = std(RecovFit(:,2));
	Template.gamma1.N = length(RecovFit(:,2));
	Template.gamma2.m = mean(RecovFit(:,3));
	Template.gamma2.s = std(RecovFit(:,3));
	Template.gamma2.N = length(RecovFit(:,3));
	weightedGamma = RecovFit(:,1).*RecovFit(:,2)+(1-RecovFit(:,1)).*RecovFit(:,3);
	Template.weightedGamma.m = mean(weightedGamma);
	Template.weightedGamma.s = std(weightedGamma);
	Template.weightedGamma.N = length(weightedGamma(:,1));
% Recovery. Fraction recovered curve
	[Template.Recovery.Delays,I] = sort(a(1,:));
	Template.Recovery.m = mean(b(:,I));
	Template.Recovery.N = size(b(:,I),1);
	Template.Recovery.s = std(b(:,I));

%%% Total charge movement (integral of current) %%%
	load(fullfile(Folder,'ChargeMovement.mat'));
% Activation.
	Template.charge_activation.m = mean(Charge_Activation(:,1:34));
	Template.charge_activation.s = std(Charge_Activation(:,1:34));
	Template.charge_activation.N = size(Charge_Activation(:,1:34),1);
% Inactivation.
	Template.charge_inactivation.m = mean(Charge_Inactivation(:,1:32));
	Template.charge_inactivation.s = std(Charge_Inactivation(:,1:32));
	Template.charge_inactivation.N = size(Charge_Inactivation(:,1:32),1);

%%% Average activation current waveforms %%%
	load(fullfile(Folder,'AverageResponseCurves.mat'))
	Template.Y = activationCurrents(1040:1539,1:34)/max(max(abs(activationCurrents(1040:1539,1:34))));

%%% Time to peak current %%%
	load(fullfile(Folder,'PeakTimes.mat'))
% Activation
	Template.ActPeak.m = mean(ActPeak(:,1:34));
	Template.ActPeak.s = std(ActPeak(:,1:34));
	Template.ActPeak.N = length(ActPeak(:,1));
% Inactivation
	Template.InactPeak.m = mean(InactPeak(:,1:32));
	Template.InactPeak.s = std(InactPeak(:,1:32));
	Template.InactPeak.N = length(InactPeak(:,1));

%%% Default model parameters. These are taken from Capes et al. 2013 and relate to the Nav1.4 model presented in that paper %%%
	Template.defaultP0 = [2100.00,1410.00,0.21,0.21,14910.00,0.32,1890000.00,170.00,-2.39,-2.39,800.00,-0.91,1000.00,-0.50,1000.00,-0.50,700.00,-0.50,14640.00,1.99,24880.00,22500.00,0.01,0.01,600.00,1.60,-0.30,-0.30,1.00,0.12,1.00,0.86];
	Template.defaultOrder = {'alpha_k0';'alpha_q';'x_alpha';'y_alpha';'beta_k0';'beta_q';'x_beta';'y_beta';'alpha4_k0';'alpha4_q';'beta4_k0';'beta4_q';'alpha4O_k0';'alpha4O_q';'beta4O_k0';'beta4O_q';'gamma_k0';'gamma_q';'delta_k0';'delta_q';'delta4_k0';'delta4_q';'deltai_k0';'deltai_q';'i_k0';'i_q';'r_k0';'r_q';'iO_k0';'iO_q';'rO_k0';'rO_q'};
	Template.idcsSort = [9,13,14,10,1,2,11,15,16,12,5,6,21,22,23,24,19,20,17,18,25,29,30,26,27,31,32,28,3,7,4,8]; % I don't think this is used for anything anymore, but I'll keep it just in case...
	Template.defaultNonNegative =  [1:2:32 4 8]; % Which of the above parameters must be non negative

%%% Path to the folder containing the experimental data %%%
	Template.path = Folder;



