function A = labelpanel(x,y,str,forcetype)

if(nargin<4)
	% str = lower(str); % NPG
	str = upper(str); % CellPress
else
	if(strcmp(forcetype,'NPG'))
		str = lower(str); % NPG
	else
		str = upper(str); % CellPress
	end
end

A = annotation('textbox', [x,y,0.05,0.05],'String',str, 'LineStyle', ...
	'none', 'FontWeight','bold', 'FontSize',8,'Margin',0);