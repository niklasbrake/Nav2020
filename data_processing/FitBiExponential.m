function [fitresult, gof] = FitExponential(T, I)
% [fitresult, gof] = FitExponential(T,I) fits the function
% 	i = (1-exp(-t*gamma1))*A + (1-exp(-t*gamma2))*(1-A)
% to the data (T,I). fitresult is ordered as [A,gamma1,gamm2].


[xData, yData] = prepareCurveData( T, I );

% Set up fittype and options.
ft = fittype( '(1-exp(-(t)*gamma1))*A + (1-exp(-(t)*gamma2))*(1-A)', 'independent', 't', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Algorithm = 'Trust-Region';
opts.DiffMaxChange = 0.1;
opts.Display = 'Off';
opts.MaxFunEvals = 1000;
opts.MaxIter = 1000;
opts.StartPoint = [1 1 1];
opts.Lower = [0 0 0];
opts.Upper = [1 Inf Inf];
opts.TolFun = 1e-6;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
