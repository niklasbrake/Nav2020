function computetimetopeak

fittingResults = load('Par_nav15_NB.mat');
Params = fittingResults.Generation(1,:);

figure; 
[actEst1,V] = simActivation(Params);
[a1,b] = max(abs(actEst1));
b(b==5e4) = nan;
subplot(1,3,2); h=plot(V,b*1e-4,'LineWidth',1); hold on; set(gca,'yscale','log');
subplot(3,3,4); plot(actEst1(:,[1:10:200]),'Color',h.Color); hold on; 
xlim([0,50000]); set(get(gca,'xaxis'),'visible','off'); set(get(gca,'yaxis'),'visible','off');
subplot(1,3,3); plot(V,max(abs(actEst1./(V-62)))/max(max(abs(actEst1./(V-62)))),'LineWidth',1); hold on;

Params(21) = Params(21)*5;
Params(5) = Params(5)*5;
[actEst2,V] = simActivation(Params);
[a2,b] = max(abs(actEst2));
b(b==5e4) = nan;
subplot(1,3,2); h=plot(V,b*1e-4,'LineWidth',1); set(gca,'yscale','log');
subplot(3,3,7); plot(actEst2(:,[1:10:200]),'Color',h.Color); 
xlim([0,50000]); set(get(gca,'xaxis'),'visible','off'); set(get(gca,'yaxis'),'visible','off');
subplot(1,3,3); plot(V,max(abs(actEst2./(V-62)))/max(max(abs(actEst2./(V-62)))),'LineWidth',1);


Params(21) = Params(21)/25;
Params(5) = Params(5)/25;
[actEst3,V] = simActivation(Params);
[a3,b] = max(abs(actEst3));
b(b==5e4) = nan;
subplot(1,3,2); h=plot(V,b*1e-4,'LineWidth',1); set(gca,'yscale','log');
subplot(3,3,1); plot(actEst3(:,[1:10:200]),'Color',h.Color); 
xlim([0,50000]); set(get(gca,'xaxis'),'visible','off'); set(get(gca,'yaxis'),'visible','off');
subplot(1,3,3); plot(V,max(abs(actEst3./(V-62)))/max(max(abs(actEst3./(V-62)))),'LineWidth',1);

subplot(1,3,2);
set(gca,'LineWidth',1);
legend({'(\gamma_k^0, \alpha_k^0)','(5\gamma_k^0, 5\alpha_k^0)','(0.2\gamma_k^0, 0.2\alpha_k^0)'});
xlabel('Voltage (mV)');
ylabel('Time (ms)');

% subplot()


function [actEst,V] = simActivation(Params)
	[Q,OpenPositions,P] = nav15_NB(Params);
	N = length(Q(0));
	dX_base = Q(-100*1e-3); % Get transition matrix for V = -100 mV
	temp = expsolver(dX_base,[1:100]*1e-3,[1 zeros(1,N-1)])'; % Integrate for 100 ms
	Xinit = temp(end,:)'; % Take final "steady-state" conformation of system
	VSteps = 200;
	V = linspace(-75,50,VSteps);
	T_max = 50000; % 500 frames with frame time of 0.1 microseconds, thus 5 miliseconds.

	X2 = zeros(T_max,N,VSteps); % Allocates memory
	for idx = 1:VSteps % For each voltage step
		V_temp = 1e-3*V(idx); % Scale voltage from mV to V
		dX = Q(V_temp); % Get transition matrix
		X2(:,:,idx) = expsolver(dX,[1:T_max]*1e-7,Xinit)'; % Integrate voltage step for 500 ms
	end
	actEst = squeeze(sum(X2(:,OpenPositions,:).*[1,1],2)).*(V-62); % Repeat for activation protocol
