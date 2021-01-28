function getresponsecurves(Folder)
% GETRESPONSECURVES calculates average waveforms for activation,
% inactivation, and recovery protocols.
% 	GETRESPONSECURVES(Folder) uses the experiments in Folder to 
% 	calculate average current waveforms and saves them in the file
% 	AverageResponseCurves.mat.


	if(nargin == 0)
		Folder = pwd;
	end

	F = dir(Folder); F = F(3:end);

	% Get Inactivation currents
	VClamp_A = [];
	VClamp_I = [];
	activationResponse = [];
	inactivationResponse = [];
	exResponse = [];
	for i = 1:length(F)
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
			load(fullfile(Folder,F(i).name,'activation.mat'));
			Current = activationLeakCorrection(Voltage,Current,Epochs);
			% Get average voltage during clamp period
			temp = mean(Voltage(:,Epochs(4):Epochs(5)),2);
			% Append to running vector
			VClamp_A(end+1:end+length(temp)) = temp;
			% Get current response during clamp period (+ extra because of inactivaiton)
			temp2 = Current(:,Epochs(4)-1000:Epochs(5)-1);
			% activationResponse(end+1:end+length(temp),:) = temp2;
			activationResponse(end+1:end+length(temp),:) = temp2/max(max(abs(temp2)));
		catch
			disp([F(i).name ': no activation data.'])
		end

		try
			load(fullfile(Folder,F(i).name,'inactivation.mat'))
			Current = inactivationLeakCorrection(Voltage,Current,Epochs);
			% Get average voltage during clamp period
			temp = mean(Voltage(:,Epochs(4):Epochs(5)),2);
			% Append to running vector
			VClamp_I(end+1:end+length(temp)) = temp;
			% Get current response during clamp period (+ extra because of inactivaiton)
			temp2 = Current(:,Epochs(4)-1000:Epochs(6)-1);
			% inactivationResponse(end+1:end+length(temp),:) = temp2;
			inactivationResponse(end+1:end+length(temp),:) = temp2/max(max(abs(temp2)));
		catch
			disp([F(i).name ': no inactivation data.'])
		end
	end
	% Save the activation responses of inactivation protocol for activation data
	% exVClamp = VClamp;
	% exResponse = inactivationResponse(:,1:Epochs(5)-Epochs(4)+1);
	% Get all the unique voltage protocol magnitudes
	inactivationPotentials = unique(round(VClamp_I/5)*5);
	for i = 1:length(inactivationPotentials)
		% Find all the currents corresponding to each voltage command
		idx = find(round(VClamp_I/5)*5==inactivationPotentials(i));
		% Get average current response
		inactivationCurrents(i,:) = mean(inactivationResponse(idx,:),1);
	end
	inactivationCommands = VClamp_I;

	% Get all the unique voltage protocol magnitudes
	activationPotentials = unique(round(VClamp_A/5)*5);
	for i = 1:length(activationPotentials)
		% Find all the currents corresponding to each voltage command
		idx = find(round(VClamp_A/5)*5==activationPotentials(i));
		% idx2 = find(round(exVClamp_A/5)*5==activationPotentials(i));
		% activationCurrents(i,:) = mean([exResponse(idx2,:);activationResponse(idx,:)],1);
		% Get average current response
		activationCurrents(i,:) = mean(activationResponse(idx,:),1);
	end
	activationCommands = VClamp_A;

	% Get delay currents (both 1ms and 10ms)
	upTimes = [];
	downTimes = [];
	recoveryResponse = [];
	Delay = [];
	for i = 1:length(F)
		try
			load(fullfile(Folder,F(i).name,'recovery1ms.mat'))
			Current = recoveryLeakCorrection(Voltage,Current,Epochs,1);
		catch
			disp([F(i).name ': no recovery1ms data.'])
			continue
		end	
		k = length(upTimes);
		% Loop over each voltage command
		for j = 1:size(Voltage,1)
			Delay(i,end+1) = 100*(i-1);
			% Find end of first pulse
			temp = find(Voltage(j,Epochs(4)+300:Epochs(6)+300)'<-60);
			downTimes(k+j) = temp(1)+300;
			% Find beginning of second pulse
			temp = find(Voltage(j,Epochs(4)+300+downTimes(k+j):Epochs(8))'>-60);
			% Append beginning of second pulse (time starting from end of first pule) to Time
			upTimes(k+j) = temp(1)+300+downTimes(k+j);
			% Save current response to first pulse
			temp2 = [Current(j,Epochs(4)-1000:Epochs(4)+upTimes(k+j)+4999) Current(j,Epochs(4)+upTimes(k+j)+4999)*ones(1,32000-1-(upTimes(k+j)+4999))];
			% recoveryResponse(k+j,:) = temp2;
			recoveryResponse(k+j,:) = temp2/max(max(abs(temp2)));
		end
		
		% Do the same thing for the 10ms data
		try
			load(fullfile(Folder,F(i).name,'recovery10ms.mat'))
			Current = recoveryLeakCorrection(Voltage,Current,Epochs,10);
		catch
			disp([F(i).name ': no recovery10ms data.'])
			continue
		end
		k = length(upTimes);
		for j = 1:size(Voltage,1)
			% Find end of first pulse
			temp = find(Voltage(j,Epochs(4)+300:Epochs(6)+300)'<-60);
			downTimes(k+j) = temp(1)+300;
			% Find beginning of second pulse
			temp = find(Voltage(j,Epochs(4)+300+downTimes(k+j):Epochs(8))'>-60);
			% Append beginning of second pulse (time starting from end of first pule) to Time
			upTimes(k+j) = temp(1)+300+downTimes(k+j);
			% Save current response to first pulse
			temp2 = [Current(j,Epochs(4)-1000:Epochs(4)+upTimes(k+j)+4999) Current(j,Epochs(4)+upTimes(k+j)+4999)*ones(1,32000-1-(upTimes(k+j)+4999))];
			% recoveryResponse(k+j,:) = temp2;
			recoveryResponse(k+j,:) = temp2/max(max(abs(temp2)));
		end
	end
	% Unique lags between first and second pulse (in frames)
	delayTimeFrames = sort(unique(round(upTimes)));
	for i = 1:length(delayTimeFrames)
		% Find recordings corresponding to this lag time
		idx = find(round((upTimes))==delayTimeFrames(i));
		% Average responses to second pulse
		recoveryCurrents(i,:) = mean(recoveryResponse(idx,:),1);
	end
	% Convert from frames to ms
	upTimes = upTimes/100;
	downTimes = downTimes/100;
	upTime = sort(unique(round(upTimes)));
	downTime = sort(unique(round(downTimes)));
	

	% Transpose data to make plotting easier
	inactivationCurrents = inactivationCurrents';
	activationCurrents = activationCurrents';
	recoveryCurrents = recoveryCurrents';
	inactivationResponse = inactivationResponse';
	activationResponse = activationResponse';
	recoveryResponse = recoveryResponse';

	save(fullfile(Folder,'ResponseBundle.mat'),'inactivationResponse','inactivationCommands','activationResponse','activationCommands','recoveryResponse','upTimes','downTimes');
	save(fullfile(Folder,'AverageResponseCurves.mat'),'inactivationCurrents','inactivationPotentials','activationCurrents','activationPotentials','recoveryCurrents','upTime','downTime');

function responses = activationLeakCorrection(commands,responses,Epochs)
	% Elementwise rather than mean...
	for i = 1:size(responses,1)
		y1 = mean(responses(i,1:Epochs(3)));
		x1 = mean(commands(i,1:Epochs(3)));
		y2 = mean(responses(i,Epochs(3):Epochs(4)));
		x2 = mean(commands(i,Epochs(3):Epochs(4)));
		m = (y2-y1)/(x2-x1);
		C = @(V) (V-x1)*m+y1;
		responses(i,:) = responses(i,:) - C(commands(i,:));
		responses(i,Epochs(4):Epochs(5)) = responses(i,Epochs(4):Epochs(5)) - responses(i,Epochs(5)-100);
	end

function responses = inactivationLeakCorrection(commands,responses,Epochs)
	% Elementwise rather than mean...
	for i = 1:size(responses,1)
		y1 = mean(responses(i,1:Epochs(3)));
		x1 = mean(commands(i,1:Epochs(3)));
		y2 = mean(responses(i,Epochs(3):Epochs(4)));
		x2 = mean(commands(i,Epochs(3):Epochs(4)));
		m = (y2-y1)/(x2-x1);
		C = @(V) (V-x1)*m+y1;
		responses(i,:) = responses(i,:) - C(commands(i,:));
		responses(i,Epochs(4):Epochs(5)) = responses(i,Epochs(4):Epochs(5)) - responses(i,Epochs(5)-100);
		responses(i,Epochs(5):Epochs(6)) = responses(i,Epochs(5):Epochs(6)) - responses(i,Epochs(6)-100);
	end

function responses = recoveryLeakCorrection(commands,responses,Epochs,lag)
	% Elementwise rather than mean...
	C1 = zeros(size(responses));
	for i = 1:size(responses,1)
		y1 = mean(responses(i,1:Epochs(3)));
		x1 = mean(commands(i,1:Epochs(3)));
		y2 = mean(responses(i,Epochs(3):Epochs(4)));
		x2 = mean(commands(i,Epochs(3):Epochs(4)));
		m = (y2-y1)/(x2-x1);
		C = @(V) (V-x1)*m+y1;
		responses(i,:) = responses(i,:) - C(commands(i,:));
		C1(i,Epochs(4):Epochs(5)) = responses(i,Epochs(5)-100);
		C1(i,Epochs(6)+100*lag*(i-1):Epochs(7)+100*lag*(i-1)) = responses(i,Epochs(7)+100*lag*(i-2));
		responses(i,:) = responses(i,:) - C1(i,:); 
	end

function responses2 = capactianceCorrection(commands,responses,Epochs)
	Pre = mean(mean(Voltage(:,Epochs(2):Epochs(3)),2));
	Post = mean(mean(Voltage(:,Epochs(3):Epochs(4)),2));
	C_Correct = zeros(size(responses,2),1);
	C_Correct(Epochs(4):Epochs(4)+1000) = mean(responses(:,Epochs(3):Epochs(3)+1000),1)/(Post-Pre);
	DeltaV = Command-Post;
	responses2 = responses'- C_Correct*DeltaV';

