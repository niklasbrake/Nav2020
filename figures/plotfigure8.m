	fig = figure('color','w','units','centimeters');
	fig.Position(3) = 8.5;
	fig.Position(4) = 8;

	T = [1,-55.5,1.3,-54.2,0.8,-18.3,1.5,-19.0,1.8,1;
	1,-35.0,1.0,-43.0,1.0,-17.0,1.0,-18.0,1.0,1;
	2,-51.7,2.6,-58.1,1.3,-18.9,1.4,-20.1,3.1,1;
	2,-42.0,2.0,-52.0,1,-18.0,2,-21.0,2.0,1;
	3,-64.9,1.5,-59.9,0,-25.5,1.6,-25.5,0,0;
	4,-54.0,1.0,-60.0,1.0,-27.0,1.0,-26.0,1.0,1;
	4,-54.0,0.4,-62.0,0.1,-25.0,0.5,-33.0,0.3,1;
	4,-66.1,0.6,-66.5,0.5,-22.8,0.9,-23.9,1.2,0;
	4,-74.2,1.9,-65.3,1.6,-24.8,1.2,-24.7,2.3,0;
	4,-67.7,1.3,-66.0,2.2,-19.2,2.0,-19.1,4.7,0;
	4,-64.0,0.7,-59.9,0.5,-8.0,0.7,-6.1,0.56,0;
	5,-76.0,1.0,-78.0,1.0,-36.0,1.0,-37.0,1.0,1;
	5,-52.0,4.2,-56.1,3.5,-27.8,2.3,-28.7,1.6,1; % Qu two-microelectrode
	5,-65.3,0.9,-65.9,0.8,-18.6,4.2,-23.3,2.0,1;
	5,-84.8,2.5,-74.0,2.4,-35.8,1.4,-34.6,1.9,1;
	5,-77.1,0.5,-72.9,1.0,-25.6,1.0,-24.9,1.0,0;
	5,-82.0,0.5,-73.5,1.0,-16.8,0.5,-14.8,0.7,0;
	6,-74.3,2.3,-72.2,0.6,-36.7,1.1,-34.8,1.7,0;
	6,-51.5,0.4,-50.8,0.3,-13.4,0.8,-14.2,0.3,1;
	6,-54.7,0.7,-51.4,1.3,-10.4,0.65,-9.9,0.7,1;
	7,-70.9,0.5,-65.7,0.5,-18.6,0.4,-17.4,1.8,0;
	7,-68.2,0.4,-69.2,0.4,-22.0,2.7,-27.7,1.3,1;
	% 8,-54.8,1.7,-62.6,2.3,4.7,0.7,-3.3,1.0,1; % Outlier (see Methods)
	8,-43.2,2.0,-47.8,1.5,-12.5,1.7,-16.5,1.4,0];

	mrks = {'o','v','square','square'};
	mrks = {'o','v','o'};

	clrs = lines(8);
	clrs(8,:) = clrs(4,:);
	clrs(4,:) = [0,1,1];
	clrs(5,:) = [0,0,1];
	clrs(6,:) = [0.5,1,0];

	axes('Position',[0.6083 0.6660 0.3537 0.2871]);
		for i = [1,2,3,7,8,6,4,5]
			tempT = T(T(:,1)==i,:);
			dSSI = tempT(:,4)-tempT(:,2);
			sdSSI = sqrt(tempT(:,5).^2+tempT(:,3).^2);
			L = size(tempT,1);
			locs = ([1:L]/2)-mean(1:L)/2;
			if(L>1)
				locs = locs*0.2;
			end
			for k = 1:L
				x = i+locs(k);
				y = dSSI(k);
				sy = sdSSI(k);
				h = plot([x,x],[-sy sy]+y,'color',0*clrs(i,:)); hold on;
				h.Color = [h.Color 0.6];
				scatter(x,y,10,0*clrs(i,:),'filled','Marker',mrks{1+tempT(k,end)},'MarkerFaceAlpha',0.6);
			end
		end
		xlim([0.5,8.5]);
		xticks([1:8]);
		xlabel('Nav1.x');
		ylabel(['\Delta SSI (mV)']);
		box off
		set(gca,'TickDir','out')
		set(gca,'LineWidth',1)
		line([0.5,8.5],[0,0],'color','k','LineStyle','--','LineWidth',1)
		set(gca,'FontSize',7);
		ylim([-15,15])
		h(1) = plot(nan,nan,'ok','MarkerSize',4,'MarkerFaceColor','k');
		h(2) = plot(nan,nan,'vk','MarkerSize',4,'MarkerFaceColor','k');

	axes('Position',[0.1177 0.6621 0.3470 0.2903]);
		for i = [1,2,3,7,8,6,4,5]
			tempT = T(T(:,1)==i,:);
			dSSI = tempT(:,8)-tempT(:,6);
			sdSSI = sqrt(tempT(:,7).^2+tempT(:,5).^2);
			L = size(tempT,1);
			locs = ([1:L]/2)-mean(1:L)/2;
			if(L>1)
				locs = locs*0.2;
			end
			for k = 1:L
				x = i+locs(k);
				y = dSSI(k);
				sy = sdSSI(k);
				h = plot([x,x],[-sy sy]+y,'color',0*clrs(i,:)); hold on;
				h.Color = [h.Color 0.6];
				scatter(x,y,10,0*clrs(i,:),'filled','Marker',mrks{1+tempT(k,end)},'MarkerFaceAlpha',0.6);
			end
		end
		xlim([0.5,8.5]);
		xticks([1:8]);
		xlabel('Nav1.x');
		ylabel(['\Delta GV (mV)']);
		box off
		set(gca,'TickDir','out')
		set(gca,'LineWidth',1)
		line([0.5,8.5],[0,0],'color','k','LineStyle','--','LineWidth',1)
		set(gca,'FontSize',7);
		ylim([-15,15])
		h(1) = plot(nan,nan,'ok','MarkerSize',3,'MarkerFaceColor','k');
		h(2) = plot(nan,nan,'vk','MarkerSize',3,'MarkerFaceColor','k');
		L = legend(h,'HEK293','X. Oocyte');
		L.ItemTokenSize = [10,10];
		L.Position = [0.1138 0.8662 0.1916 0.0906];
		L.Box = 'off'

	axes('position',[0.2248 0.1105 0.4820 0.4106]);
		for i = [1,2,3,7,8,6,4,5]
			tempT = T(T(:,1)==i,:);
			sCSI = sqrt(tempT(:,7).^2+tempT(:,3).^2);
			CSI = tempT(:,6)-tempT(:,2);
			dSSI = tempT(:,4)-tempT(:,2);
			sdSSI = sqrt(tempT(:,5).^2+tempT(:,3).^2);
			L = size(tempT,1);
			for k = 1:L
				x = CSI(k); sx = sCSI(k); y = dSSI(k); sy = sdSSI(k);
				h = plot([-sx sx]+x,[y, y],'color',clrs(i,:)); hold on;
				h.Color = [h.Color,0.6];
				h = plot([x x],[-sy, sy]+y,'color',clrs(i,:));
				h.Color = [h.Color,0.6];

				scatter(x,y,20,clrs(i,:),'filled','Marker',mrks{1+tempT(k,end)},'MarkerFaceAlpha',0.6);
			end
		end
		xlabel(['CSI* (mV)']);
		ylabel('\Delta SSI (mV)');
		box off
		set(gca,'TickDir','out')
		set(gca,'LineWidth',1)
		xlim([10,80])
		line(get(gca,'xlim'),[0,0],'color','k','LineStyle','--','LineWidth',1)
		ylim([-15,15])
		dSSI = T(:,4)-T(:,2);
		sdSSI = sqrt(T(:,5).^2+T(:,3).^2);
		FT = fitlm((T(:,6)-T(:,2)),dSSI,'Weights',1./sdSSI)
		plot(10:80,FT.predict((10:80)'),'color','k')
		set(gca,'FontSize',7);
		R2 = FT.Rsquared.Ordinary;
		pV = fix(log(FT.anova{1,end})/log(10));
		text(60,-8,['R^2 = ' num2str(R2,2)],'FontSize',6);
		text(60,-12,sprintf(['p < 10^{' int2str(pV) '}']),'FontSize',6);

		for i = 1:8
			text(90,10-20*i/8,['Nav1.' int2str(i)],'color',clrs(i,:),'FontSize',7);
		end
			% h(i) = plot(nan,nan,'ok','MarkerSize',4,'MarkerFaceColor','k');


		L = labelpanel(0.013,0.93,'a');
		L = labelpanel(0.5,0.93,'b');
		L = labelpanel(0.12,0.5,'c');