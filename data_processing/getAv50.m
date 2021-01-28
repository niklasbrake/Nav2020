function [WT,B1,B3] = getAv50

figure;
subplot(1,3,1);
WT = load('E:\Documents\Work\PhD\Nav15\Experiments\Data\Nav1.5e\GV_Curves.mat');
	x = WT.A_V(1,:)';
	y = transpose(WT.A_I);
	ys = -y./min(y);
	outliers = [6,16,45];
	ys(:,outliers) = [];
	WT.pars = fitpars(x,ys);

subplot(1,3,2);
B1 = load('E:\Documents\Work\PhD\Nav15\Experiments\Data\Nav1.5e+B1\GV_Curves.mat');
	x = B1.A_V(1,:)';
	y = transpose(B1.A_I);
	ys = -y./min(y);
	B1.pars = fitpars(x,ys);

subplot(1,3,3);
B3 = load('E:\Documents\Work\PhD\Nav15\Experiments\Data\Nav1.5e+B3\GV_Curves.mat');
	x = B3.A_V(1,:)';
	y = transpose(B3.A_I);
	ys = -y./min(y);
	outliers = 10;
	ys(:,outliers) = [];
	B3.pars = fitpars(x,ys);


function pars = fitpars(x,ys)

	preG = zeros(size(ys));
	for i = 1:size(ys,2)
		FT = FitBoltzman(x,ys(:,i),-15,-10,59,1);
		preG(:,i) = ys(:,i)*FT.Gmx./(x-FT.ERev);
	end
	V = x(1:32);
	G = preG(1:32,:);
	G = G./median(G(28:end,:));

	plot(V,G);

	for i = 1:size(ys,2)
		FT = FitBoltzmanCurve(V,G(:,i),-15,-10);
		pars(i,1) = FT.v50;
		pars(i,2) = FT.k;
	end