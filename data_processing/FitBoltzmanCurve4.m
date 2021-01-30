function [fitresult, gof] = FitBoltzman(V, I)
% [fitresult, gof] = FitBoltzman(V,I,Vhalf,K0) fits (V,I) with the function
% 	i = 1/(1+exp((v-v50)/k))
% fitresult is ordered as [k v50].

[xData, yData] = prepareCurveData( V, I );

% Set up fittype and options.
ft = fittype( '1/(1+exp(-(v-v50A)/kA)) - 1/(1+exp(-(v-v50B)/kB))', 'independent', 'v', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Algorithm = 'Trust-Region';
opts.DiffMaxChange = 100;
opts.Display = 'Off';
opts.MaxFunEvals = 600;
opts.MaxIter = 400;
opts.StartPoint = [20,20,-80,-20];
opts.TolFun = 1e-6;
opts.Lower = [0,0,-200,-80];
opts.Upper = [100,100,0,50];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
