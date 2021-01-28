function [OSI,Inact,CSI] = computemodelcsi(Q,OpenPositions)


T_max = 10e-3;
N = length(Q(0));

% Initial baseline at -100 mV 
dX_base = Q(-100*1e-3);
[~,temp] = ode15s(@(t,X) dX_base*X(:),[1,100]*1e-3,[1 zeros(1,N-1)]);
Xinit = temp(end,:)';
% Activation Initialization
V = [-120:0];
% V = [-120:035];
VSteps = length(V);

% Activation Protocol
dt = 1e-5; % 10,000 Hz
X2 = zeros(ceil(T_max/dt),N,VSteps);
OSI = zeros(VSteps,1);
Inact = zeros(VSteps,1);
CSI = zeros(VSteps,1);
for idx = 1:VSteps	
	V_temp = 1e-3*V(idx);
	dX = Q(V_temp);
	[~,X2(:,:,idx)] = ode15s(@(t,X) dX*X(:),dt:dt:T_max,Xinit);
	OSI(idx) = sum(dX(N/3,OpenPositions(1))*X2(:,OpenPositions(1),idx)*dt);
	for i = 1:4
		runCSI(i,:) = dX(i,i+5)*X2(:,i+5,idx)*dt - dX(i+5,i)*X2(:,i,idx)*dt;
	end
	CSI(idx) = sum(runCSI(:));
	% Inact(idx) = squeeze(sum(X2(end,[1:N/3],idx),2));
end
% CSI = Inact-OSI;
Inact=OSI+CSI;


plot(-120:0,Inact,'LineWidth',1)
hold on;
plot(-120:0,OSI,'LineWidth',1)
plot(-120:0,CSI,'LineWidth',1)
xlabel('Voltage (mV)');
