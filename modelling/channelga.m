function [Generation Fitness FitnessTrack]= channelga(pc,pm1,pm2,varargin)
% CHANNELGA uses an evolutationary-type algorithm to fit Markov model parameters.
% 	[Generation Fitness FitnessTrack] = CHANNELGA(pc,pm1,pm2) where pc is the probability 
% 	that at least one recomination event will occur for each parents pair, pm1 is the probablility 
% 	of mutatation 1 and pm2 is the probability of mutation 2.
% 		e.g. channelga(0.8,0.3,0.4)
% 	Generation is an n x m matrix of parameters values assicated with the final 
% 	generation. Columns are parameters and rows are individuals, sorted by fitness.
% 	Fitness is the fitness value associated with each individual (n x 1 vector). 
% 	FitnessTrack is the maximum fitness of every past generation.
% 
% 	[Generation Fitness FitnessTrack] = CHANNELGA(pc,pm1,pm2,Name,Value)
% 	 
% 	See also getchannelfitness, gendist, constructtemplate, nav15_NB.


	ip = inputParser;

	paramName = 'OutputFolder';
	defaultval = '.';
	errorMsg = 'Input must be a valid path.';
	validationFcn = @(x) assert(exist(x)==7,errorMsg);
	addParameter(ip,paramName,defaultval,validationFcn);

	paramName = 'TemplateFolder';
	defaultval = '.';
	errorMsg = 'Input must be nummeric.';
	validationFcn = @(x) assert(exist(x)==7,errorMsg);
	addParameter(ip,paramName,defaultval,validationFcn);

	paramName = 'StartPoint';
	defaultval = [];
	errorMsg = 'Input must be nummeric.';
	validationFcn = @(x) assert(isnumeric(x),errorMsg);
	addParameter(ip,paramName,defaultval,validationFcn);

	paramName = 'FixedParams';
	defaultval = [];
	errorMsg = 'Input must be nummeric.';
	validationFcn = @(x) assert(isnumeric(x),errorMsg);
	addParameter(ip,paramName,defaultval,validationFcn);

	paramName = 'FittingWeights';
	defaultval = ones(1,8);
	errorMsg = 'Input must be nummeric.';
	validationFcn = @(x) assert(isnumeric(x),errorMsg);
	addParameter(ip,paramName,defaultval,validationFcn);

	paramName = 'ModelFunction';
	defaultval = @nav14_MGO;
	errorMsg = 'Input must be nummeric.';
	validationFcn = @(x) assert(true,errorMsg);
	addParameter(ip,paramName,defaultval,validationFcn);

	paramName = 'NonNegative';
	defaultval = [1,2,5,7,8,11:2:21,22,25,26,29:32];
	errorMsg = 'Input must be nummeric.';
	validationFcn = @(x) assert(true,errorMsg);
	addParameter(ip,paramName,defaultval,validationFcn);

	paramName = 'GWeights';
	defaultval = [];
	errorMsg = 'Input must be nummeric.';
	validationFcn = @(x) assert(true,errorMsg);
	addParameter(ip,paramName,defaultval,validationFcn);

	parse(ip,varargin{:});
	OutputFolder = ip.Results.OutputFolder;
	TemplateFolder = ip.Results.TemplateFolder;
	P0 = ip.Results.StartPoint;
	FittingWeights = ip.Results.FittingWeights;
	FixedParams = ip.Results.FixedParams;
	ModelFunction = ip.Results.ModelFunction;
	NonNegative = ip.Results.NonNegative;
	GWeights = ip.Results.GWeights;

% Use TemplateFolder to create template struct from experimental data
	template = constructtemplate(TemplateFolder);
% If not startpoint was given, use the defaultP0 from constructtemplate (parameters from Capes et al. 2013 for Nav1.4. See function nav14_MGO.m)
	if(isempty(P0))
		P0 = template.defaultP0;
	end
% If no conductance weights are input, assume all open states have equal conductance
	if(isempty(GWeights))
		[Q,OpenPositions] = ModelFunction(P0);
		GWeights = ones(1,length(OpenPositions));
	end

	if(size(P0,1)>1 && size(P0,2)>1) % If entire first generation is specified
		[N1,N2] = size(P0);
		Generation = P0; % Set first generation to StartPoint
		PP0 = P0(1,:);
		FloatingParams = setdiff([1:size(PP0,2)],FixedParams); % If there are any parameters that are fixed
		P0 = P0(1,FloatingParams); % reflect this in the variable P0
	else % If there is only one initial condition
		PP0 = P0;
		FloatingParams = setdiff([1:size(PP0,2)],FixedParams);
		P0 = P0(1,FloatingParams);

		N1 = 75; 			% Number of individuals
		N2 = length(P0); 	% Number of parameters
		Generation = P0.*(1+0.3*randn(N1,N2)); % Intialize first generation with initial condition plus random Gaussian noise
		Generation(1,:) = P0;
	end
	N3 = 10; % Parents to next generation (the top n individuals are not mutated or killed)

	WT.A = load('E:\Documents\Work\PhD\Nav15\Experiments\Data\Nav1.5e\ASM_20170418_cell2\activation.mat');
	WT.I = load('E:\Documents\Work\PhD\Nav15\experiments\Data\Nav1.5e\ASM_20170418_cell3\inactivation.mat');
	WT.R = load('E:\Documents\Work\PhD\Nav15\experiments\Data\Nav1.5e\ASM_20170418_cell2\recovery1ms.mat');


% Turn off warnings on each worker
	% spmd
	% 	warning('off','MATLAB:ode15s:IntegrationTolNotMet');
	% 	warning('off','signal:findpeaks:largeMinPeakHeight');
	% 	warning('off','MATLAB:illConditionedMatrix')
	% 	warning('off','curvefit:prepareFittingData:removingNaNAndInf');
	% 	warning('off','curvefit:fit:invalidStartPoint');
	% end
	% warning('off','curvefit:prepareFittingData:removingNaNAndInf');

% Initialize stuff
	TotalTime = tic();
	max_F = [0];
	mean_F = [];
	FitnessTrack = [];
	X_Running = [];
	Y_Running = [];
	Ptemp2 = zeros(1,N2);

% Begin the selection process
	n = 1;
	while(1)
	% Weird stuff with parallel processing said these variables couldn't be found, but this fixes the problem?
		F = zeros(1,N1);
		Generation = Generation;
		fw = rand(8,N1);
		fw = (fw./sum(fw)*8);
		for i = 1:N1
			try
				Ptemp = PP0; % Generation one
				Ptemp(FloatingParams) = Generation(i,:); % update all the non-fixed parameters to the current individual
				[Q,OpenPositions] = ModelFunction(Ptemp); % Get the corresponding transition matrix
				F(i) = getchannelfitness(template,Q,OpenPositions,fw(:,i),GWeights); % Calculate fitness of individual
			catch ME
				F(i) = 0; % If anything breaks, give the model a fitness of 0
			end
		end

	% Quicksave (in case I want to stop early)
		save(fullfile(OutputFolder,'temp_results.mat'),'F','Generation');

	% Quick analysis of generation
		max_F(n) = max(F); % Max fitness
		mean_F(n) = mean(F); % Average fitness
		[F, Pidx] = sort(F,'descend'); % Sort fitnesses
		FitnessTrack(n) = F(1); % Update max fitness tracker
		Generation = Generation(Pidx,:); % Sort generation by fitness
	% Get the transition matrix from the best individual

		figure(1);
		subplot(1,2,1); hold off;
		histogram(F,'NumBins',50); title('Fitness Distribution'); box off;
		subplot(1,2,2); hold off;
		plot(max_F,'b*-'); hold on; plot(mean_F,'r*-'); 
		% legend('Max Fitness','Mean Fitness','location','northwest'); box off;
		xlabel('Generation')

		if(max_F(n)>max(max_F) || n == 1)
			Ptemp2 = PP0;
			Ptemp2(FloatingParams) = Generation(1,:);
			[Q,OpenPositions] = ModelFunction(Ptemp2);
		% Update progress tracker
			plotGAresult(Q,OpenPositions,WT,template);
		end
	% End fitting when the max fitness does not increase by at least 10% over 50 generations
		if(n > 50)
			if(max_F(end) < max_F(end-40)*1.1)
				break;
			end
		end

	% Get parents based on fitness
		F2 = max(1+(F-mean(F))/(2*std(F)),0);
		Parents = gendist(F2/sum(F2),2,N1/2);
		Children = zeros(N1,N2); % Initialize children
	% Send top unique parents directly to next generation...
		temp = unique(Generation,'rows','stable');
		Children(1:N3,:) = temp(1:N3,:);

	% Estimate the weighted variance of each parameter
	if(n<15)
		wSig = eye(N2,N2).*abs(Generation(1,:))/60;
	else
		wMu = sum(Generation .* F2',1)'/sum(F2); % Weighted average (samples of parameter space that led to higher fitness are weighted proportionally more)
		wSig = zeros(N2,N2); % Compute weighted variance estimate
		for i = 1:N1
			wSig = wSig + F2(i)*(Generation(i,:)'-wMu)*(Generation(i,:)-wMu');
		end
		wSig = wSig/(sum(F2));
	end

	% Perform mutations
		Prob_M1 = rand(N1,1);
		Prob_M2 = rand(N1,1);
		for i = N3/2+1:N1/2 % Loop through each parent pair (exluding parents that were sent directly as children to the next generation)
			idx1 = 2*i-1;
			idx2 = 2*i;

		% Copy parent's parameters directly to children
			Children(idx1,:) = Generation(Parents(1,i),:);
			Children(idx2,:) = Generation(Parents(2,i),:);

		% MUTATION 1 (convex recombination). Each parent is an N2-dimensional point in Euclidean space. Consider the orthogonal projection
		% of these points onto a subspace, defined by the parameters that will be mixed. In this subspace, draw a line
		% between both parents and select at random a point along this line. Then embed this point back in the original
		% space. This is how recombination of the parent's parameters will occur.
			temp = rand(1,N2)<pc/N2; % Each parameter in this pair will be swapped with probability pc/N2. This defines the subspace
			idx_swap = find(temp); % Find the indicies of these parameters
		% Recombination will occur in each child with probability pm1.
			if(Prob_M1(idx1) < pm1) 
				Children(idx1,idx_swap) = Generation(Parents(1,i),idx_swap) + rand()*(Generation(Parents(2,i),idx_swap) - Generation(Parents(1,i),idx_swap));
			end
			if(Prob_M1(idx2) < pm1)
				Children(idx2,idx_swap) = Generation(Parents(2,i),idx_swap) + rand()*(Generation(Parents(1,i),idx_swap) - Generation(Parents(2,i),idx_swap));
			end

		% MUTATION 2 (random perturbation). With probability pm2, each of the child's parameters will be perturbed randomly. The magnitude
		% of the perturbation is scaled by the estimated variance of that parameter, i.e. a random variable sampled from a
		% multivariate Gaussian with mean zero and variance 3*wSig, the estimated variance of the locally well-behaved manifold 
		% of parameter space
			if(Prob_M2(idx1) < pm2)
				Children(idx1,:) = mvnrnd(Children(idx1,:),3*wSig,1);
			end
			if(Prob_M2(idx2) < pm2)
				Children(idx2,:) = mvnrnd(Children(idx2,:),3*wSig,1);
			end
		end

	% Ensure constraints on parameters. i.e. non-negativity of rate constants
		Children = constrainParameters(Children,FloatingParams,NonNegative,length(PP0));
	
	% Send these children to the next generation of parents
		Generation = Children;

	% Keep a running tab of unique individuals and their fitness
		[temp,I] = unique(Generation,'rows');
		X_Running = [X_Running; temp];
		Y_Running = [Y_Running; F(I)'];
		n = n+1;
	end

% Display how long this took
	disp(['Total Elapsed Time: ' num2str(toc(TotalTime)) 's']);

% Get final generation parmaeter values
	temp = repmat(PP0,[N1,1]);
	temp(:,FloatingParams) = Generation;
	Generation = temp;

% Save algorithm hyperparameters
	inputParameters = ip.Results;

% Save results
	if(~exist(fullfile(OutputFolder,'GA_Results.mat')))
		save(fullfile(OutputFolder,'GA_Results.mat'),'Generation');
	else
		n = 2;
		while(exist(fullfile(OutputFolder,['GA_Results(' int2str(n) ').mat'])))
			n = n + 1;
		end
		save(fullfile(OutputFolder,['GA_Results(' int2str(n) ').mat']),'Generation');
	end

function Children = constrainParameters(Children,FloatingParams,NonNegative,N)
% This is a function that makes 0 a sticky boundary for parameters that should be
% non negative as determined by NonNegative.
	temp = zeros(N,1); % Initialize
	temp(NonNegative) = 1; % Get true/false encoding of non-negativity
	temp = temp(FloatingParams); % Get the subset of these are non-fixed parameters
	NonNegative = find(temp); % Find the indicies that are supposed to be non-negative
	temp = Children(:,NonNegative); % temporary variable
	temp(temp < 0) = 0; % Set any of these that are negative to be 0
	Children(:,NonNegative) = temp; % Update the children with constrained parameters.
