function [y,T] = expsolver(A,T,y0)
% EXPSOLVER solves a system of linear ODEs using the matrix exponential.
% 	[y,T] = EXPSOLVER(A,T,y0) where y is the solution, starting from initial
% 	conditions y0, T is a uniformly spaced time vector (used for calculating dt
% 	and length of integration), and A is the transition	matrix corresponding 
% 	to the ODE system. 

	dt = T(2)-T(1);
	y = zeros(length(y0),length(T));
	y(:,1) = y0;
	E = expmGA(A * dt);	% Compute matrix exponential e^(A*dt)
	for i = 1:length(T)-1
	    y(:,i+1) = E*y(:,i); 
	end