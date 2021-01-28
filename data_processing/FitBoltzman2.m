function [fitresult, gof] = FitBoltzman2(V,I,Vhalf,K0,Gmx)
% [fitresult, gof] = FitBoltzman2(V,I,Vhalf,K0,Gmx) fits (V,I) with the function
% 	i = 1/(Gmx*(1+exp((v-v50)/k)))
% fitresult is ordered as [Gmx k v50]


[xData, yData] = prepareCurveData( V, I );

% Set up fittype and options.
ft = fittype( '1/(Gmx*(1+exp((v-v50)/k)))', 'independent', 'v', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Algorithm = 'Trust-Region';
opts.DiffMaxChange = 100;
opts.Display = 'Off';
opts.MaxFunEvals = 600;
opts.MaxIter = 400;
opts.StartPoint = [Gmx K0 Vhalf];
opts.TolFun = 1e-6;
opts.Lower = [-1e5 -20 -200];
opts.Upper = [0 20 40];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'M vs. V', 'untitled fit 1', 'Location', 'NorthEast' );
% % Label axes
% xlabel V
% ylabel M
% grid on


