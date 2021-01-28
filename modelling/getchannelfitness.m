function [F,F_Breakdown,actEst,inActEst,recovEst,I_A,I_I,I_R] = getchannelfitness(Template,Q,OpenPositions,FittingWeights,GWeights)
% GETCHANNELFITNESS computes a number representing how well a channel model fits ephys data.
% 	[F,F_Breakdown,actEst,inActEst,recovEst,I_A,I_I,I_R] = GETCHANNELFITNESS(Template,Q,OpenPositions,FittingWeights,GWeights) fits
% 	the model defined by the voltage-dependent transition matrix Q, the states corresponding to open conformations, OpenPositions,
% 	the conductance weighted for each open state, GWeights, to the data in the struct Template. The importance of each feature
% 	is defined in FittingWeights. 
% 	F is the overall fitness of the model whereas F_Breakdown contains the error computed against each feature. actEst, inActEst,
% 	and recovEst are the simulated currents in response to the activation, inactivaiton, and recovery protocols, respectively.
% 	I_A, I_I, and I_R are the simulated GV, SS-inactvation, and recovery curves respectively.
% 	
% 	For more information, see <a href="matlab:web('https://github.com/niklasbrake/Brake-et-al-2019/wiki')">the reference page</a>.
% 
% 	See also expsolver, constructModelCode, constructtemplate.

%%% Numerically integrate system with voltage protocols %%%
% Part 0. Initialize constants
	T_max = 500; % 500 frames with frame time of 0.1 microseconds, thus 5 miliseconds.
	N = length(Q(0));
% Part 0. Initial baseline at -100 mV 
	dX_base = Q(-100*1e-3); % Get transition matrix for V = -100 mV
	temp = expsolver(dX_base,[1:100]*1e-3,[1 zeros(1,N-1)])'; % Integrate for 100 ms
	Xinit = temp(end,:)'; % Take final "steady-state" conformation of system
% Part 1. Initilize constants for inactivaiton protocol
	V = Template.Inactivation.Voltages; % Get voltage steps used in the inactivation experiments from template
	VSteps = length(V);	% Total number of voltage steps
	dX_Pulse = Q(-10*1e-3); % Test pulse is -10 mV
% Part 1. Inactivation Protocol
	X1 = zeros(T_max,N,VSteps); % Allocate memory
	preX1 = zeros(100,N,VSteps); % Allocate memory
	for idx = 1:VSteps % For each voltage step
		V_temp = 1e-3*V(idx); % Scale voltage from mV to V
		dX = Q(V_temp); % Get transition matrix for this voltage
		preX1(:,:,idx) = expsolver(dX,[1:100]*1e-3,Xinit)'; % Integrate pre-pulse, lasting 100 ms
		X00 = preX1(end,:,idx)'; % Take "steady-state" conformation of system
		X1(:,:,idx) = expsolver(dX_Pulse,[1:T_max]*1e-5,X00)'; % Integrate test-pulse, for 500 ms.
	end
% Part 2. Initilize constants for activaiton protocol
	V = Template.Activation.Voltages; % Get voltage steps used in the activation experiments from template
	VSteps = length(V); % Total number of voltage steps
% Part 2. Activation Protocol
	X2 = zeros(T_max,N,VSteps); % Allocates memory
	for idx = 1:VSteps % For each voltage step
		V_temp = 1e-3*V(idx); % Scale voltage from mV to V
		dX = Q(V_temp); % Get transition matrix
		X2(:,:,idx) = expsolver(dX,[1:T_max]*1e-5,Xinit)'; % Integrate voltage step for 500 ms
	end
% Part 3. Recovery Initialization
	TSteps = Template.Recovery.Delays;
	idx_temp = 21; %Index for V = -10 mV
	X00 = X2(end,:,idx_temp);
% Part 3. Recovery Protocol
	X3 = zeros(T_max,N,length(TSteps));
	for idx = 1:length(TSteps)	
		temp1 = expsolver(dX_Pulse,[1:80]*1e-3,X00')';
		temp2 = expsolver(dX_base,linspace(0,TSteps(idx),100)*1e-3,temp1(end,:))';
		X3(:,:,idx) = expsolver(dX_Pulse,[1:T_max]*1e-5,temp2(end,:))';
	end

%%% Handelling of numerical simulations %%%
% Part 0. Sometimes matrix exp. converts to complex numbers with 0 imaginary part.
	X1 = real(X1); 
	X2 = real(X2);
	X3 = real(X3);
% Part 1. Max Current calculations
	inActEst = squeeze(sum(X1(1:500,OpenPositions,:).*GWeights,2)).*(-10-Template.ERev); % Take the computed conductance (GWeights x densitiy in "open states") and scale w.r.t reversal potential
	actEst = squeeze(sum(X2(1:500,OpenPositions,:).*GWeights,2)).*(Template.Activation.Voltages-Template.ERev); % Repeat for activation protocol
	recovEst = squeeze(sum(X3(1:500,OpenPositions,:).*GWeights,2)); % Repeat for recovery protol (don't scale because we're looking at fraction recovery and always pulse to same voltage anyways)
	GMax = max(abs(actEst(:))); % Scale by peak conductance to get in range [0,1]
	actEst = actEst/GMax;
	GMax = max(abs(inActEst(:)));
	inActEst = inActEst/GMax;
	GMax = max(recovEst(:));
	recovEst = recovEst/GMax;
% Part 2. Find peak magnitude of each current
	A = min(actEst',[],2); % Find max of each current
	B = max(actEst',[],2); % Find min of each current
	[a,b] = max([abs(A) B],[],2); % Take whichever has the larger amplitude
	X = [A B]; % This is the an n x 2 matrix of (x,y) coordinates for peak amplitude
	modelActivationGV = X(sub2ind([length(B) 2],[1:length(B)]',b)); % Get the IV curve
	I_A = modelActivationGV/max(abs(modelActivationGV)); % Scale it by the max
% Part 2. Repeat for inactivation
	A = min(inActEst',[],2);
	B = max(inActEst',[],2);
	[a,b] = max([abs(A) B],[],2);
	X = [A B];
	modelInactivationGV = X(sub2ind([length(B) 2],[1:length(B)]',b));
	I_I = modelInactivationGV/max(abs(modelInactivationGV));
% Part 2. Repeat for recovery
	I_R = max(recovEst)/max(max(recovEst));
% Part 3. Calculate time to peak conductance
	offset = 22; % This offset is used to align the model and the experimental data.
	temp = actEst;
	[x I] = max(abs(temp)); % Get the peak conductances
	I(x<0.05) = 500-offset; % Set an conductances are smaller than 0.05 to "infinite lag" since they wouldn't stand above the noise in the experimental recordings
	ActPeak = I+offset; % Add this lag to allow direct comparison with experimentally determined time-to-peak
	idcsAct = find(Template.ActPeak.s > 0); % Get the indicies where experimentally-determined time-to-peak is meaningful (variance is non-zero).
% Part 3. Repeat for inactivation
	temp = inActEst;
	[x I] = max(abs(temp));
	I(x<0.05) = 500-offset;
	InactPeak = I+offset;
	idcsInact = find(Template.InactPeak.s > 0);

	if(sum(abs(FittingWeights))==0)
		F=0;
		F_Breakdown = zeros(8,1);
		return
	end

%%% Calculate Error %%%
% Activation
	F_Breakdown(1) = mean(abs(Template.Activation.m - I_A)./sqrt(Template.Activation.s.^2/Template.Activation.N)); % Difference between IV curves, normalized by standard error of template
	F_Breakdown(2) = mean(abs(Template.charge_activation.m - sum(actEst))./sqrt(Template.charge_activation.s.^2/Template.charge_activation.N)); % Integral of current
	F_Breakdown(3) = mean(abs(Template.ActPeak.m(idcsAct) - ActPeak(idcsAct))./sqrt(Template.ActPeak.s(idcsAct).^2/length(idcsAct))); % Time to peak 
% Inactivation
	F_Breakdown(4) = mean(abs(Template.Inactivation.m - I_I)./sqrt(Template.Inactivation.s.^2/Template.Inactivation.N)); % IV curves
	F_Breakdown(5) = mean(abs(Template.charge_inactivation.m - sum(inActEst))./sqrt(Template.charge_inactivation.s.^2/Template.charge_inactivation.N)); % Integral of current
	F_Breakdown(6) =  mean(abs(Template.InactPeak.m(idcsInact) - InactPeak(idcsInact))./sqrt(Template.InactPeak.s(idcsInact).^2/length(idcsInact))); % Time to peak
% Recovery
	F_Breakdown(7) = mean(abs(Template.Recovery.m - I_R)./sqrt(Template.Recovery.s.^2/Template.Recovery.N)); % Difference in recovery curves
	F_Breakdown(8) = 100*sum(sum(temp2(:,OpenPositions)));
% Weight each component
	F_Breakdown = F_Breakdown.*FittingWeights; % Scale each error by weight
	F = norm(F_Breakdown,2); % Take the L2 norm
% Set infinite error to a model that breaks anything in this function
	if(isnan(F))
		F = Inf;
	end
% Fitness is the reciprocal of the error
	F = F^(-1);