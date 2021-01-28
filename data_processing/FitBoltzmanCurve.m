function [fitresult, gof] = FitBoltzman(V, I,Vhalf0,K0,Gmx)
% [fitresult, gof] = FitBoltzman(V,I,Vhalf,K0) fits (V,I) with the function
% 	i = 1/(1+exp((v-v50)/k))
% fitresult is ordered as [k v50].

[xData, yData] = prepareCurveData( V, I );

% Set up fittype and options.
ft = fittype( 'Gmx/(1+exp((v-v50)/k))', 'independent', 'v', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Algorithm = 'Trust-Region';
opts.DiffMaxChange = 100;
opts.Display = 'Off';
opts.MaxFunEvals = 600;
opts.MaxIter = 400;
opts.StartPoint = [Gmx,K0,Vhalf0];
opts.TolFun = 1e-6;
opts.Lower = [0 -20 -150];
opts.Upper = [Inf 20 150];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
