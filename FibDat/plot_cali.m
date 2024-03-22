cal_data = {};
cal_data_aero = {};
cal_out = {};
cal_out_aero = {};
fib_path = dir("Proc/fib");
aero_path = dir("Proc/aero");

contains = @(str, pattern) ~cellfun('isempty', strfind(str, pattern));

function fig_handle = read_n_plot(fib_raw, fib_cal, aero_raw, aero_cal)

	contains = @(str, pattern) ~cellfun('isempty', strfind(str, pattern));
	name = strsplit(aero_raw, '/');
	name = name(3);

	fig = figure();
	set(fig,'position',[10,10,225,225]);
	ax1 = axes(fig);
	hold on
	ax2 = axes(fig);
	hold on
	set(ax2, 'XAxisLocation', 'top');
	set(ax2, 'YAxisLocation', 'right');
	set(ax2, 'Color', 'none');
	set(ax2, 'YTick', []);
	axp = get(ax1, 'Position');
	axp = [0.17, 0.13, axp(3), 0.7];
	set(ax1, 'Position', axp);
	set(ax2, 'Position', axp);

	cal_data = dlmread(fib_raw, ',', 1, 0);
	cal_out = dlmread(fib_cal, '', 0, 0);
	cal_data_aero = dlmread(aero_raw, ',', 1, 0);
	cal_out_aero = dlmread(aero_cal, '', 0, 0);
	m_data = dlmread('MaterialData_water_1p31.csv', ',', 1, 2);
	cal_out = cal_out(:, end);
	cal_out_aero = cal_out_aero(:, end);

	theory = interp1(cal_data(:, end), cal_data(:, 5), ...
		linspace(0, (4095-cal_out(2))/cal_out(1), 1000), ...
		'linear', 'extrap');

	p1 = plot(ax1, linspace(0, (4095-cal_out(2))/cal_out(1), 1000)*1e6, theory, ...
		'--k','LineWidth',0.5,'DisplayName','Theoretical Response');
	p2 = plot(ax1, linspace((0-cal_out(2))/cal_out(1), ...
		(4095-cal_out(2))/cal_out(1), 1000)*1e6, linspace(0, 4095, 1000), ...
		'-.k','LineWidth',0.5, "DisplayName", "Fibre Fit");
	p3 = plot(ax1, cal_data(:, end-1)*1e6, cal_data(:, 5), ...
		'xk','LineWidth',0.5,'DisplayName','Fibre Data');
	p4 = plot(ax1, linspace((0-cal_out_aero(2))/cal_out_aero(1), ...
		(4095-cal_out_aero(2))/cal_out_aero(1), 1000)*1e6,...
		linspace(0, 4095, 1000), ...
		'-k','LineWidth',0.5, "DisplayName", "Aerosol Fit");
	p5 = plot(ax1, cal_data_aero(:, end-1)*1e6, cal_data_aero(:, 5), ...
		'*k','LineWidth',0.5,'DisplayName','Aerosol Data');

	handles = {p1, p2, p3, p4, p5};
	order = [2, 3, 1];

	for i=1:size(cal_data, 1)
		p1 = plot(ax1, [cal_data(i, end-1)*1e6, cal_data(i, end-1)*1e6],...
			[0, 4095], ':k', "HandleVisibility", "off", ...
			"LineWidth", 0.3);
		txt = strcat('f', num2str(order(i)));
		p2 = text(ax1, cal_data(i, end-1)*1e6, 4000, txt,...
			'HorizontalAlignment', 'right',...
			'Rotation',90,'fontname','times','VerticalAlignment', ...
			'top',"fontsize", 8);
		handles = [handles, {p1, p2}];
	end
	for i=1:size(cal_data_aero, 1)
		p1 = plot(ax1, [cal_data_aero(i, end-1)*1e6, ...
			cal_data_aero(i, end-1)*1e6], [0, 4095], ...
			':k', "HandleVisibility", "off");
		txt = strcat('a', num2str(order(i)));
		p2 = text(ax1, cal_data_aero(i, end-1)*1e6, 3600, txt,...
			'HorizontalAlignment','right', 'Rotation',90,'fontname',...
			'times','VerticalAlignment','top',"fontsize", 8);
		handles = [handles, {p1, p2}];
	end

	p1 = text(0.35e-4, 300, name, ...
		'fontname','times',"fontsize", 8);
	handles = [handles, {p1}];

	set(ax1,'fontname','times');
	set(ax1, "fontsize", 6)
	set(ax2,'fontname','times');
	set(ax2, "fontsize", 6)
	set(ax1, 'Box', 'off');
	set(ax2, 'Box', 'off');

	ylim(ax1, "manual");
	ylim(ax1, [0, 4095])
	xlim(ax1, "manual");
	%xlim(ax1, [0, ((4095-cal_out(2))/cal_out(1))*1e6])
	if contains(name, '-AA-')
		x2 = 0.9e-4;
	else
		x2 = 1.9e-4;
	end
	xlim(ax1, [0, x2])
	sca_samples = linspace(0, x2/1e6, 5);
	rad_samples = interp1(m_data(:,1), m_data(:,2), sca_samples, 'linear');
	rad_samples(1) = 0;
	rad_samples = round(rad_samples);
	xlim(ax2, "manual");
	xlim(ax2, [0, rad_samples(end)])
	set(ax2, 'XTickMode', 'manual')
	set(ax2,'XTick',rad_samples(1:end-1))

	set(fig, "resize", "off")
	set(fig, 'PaperPositionMode', 'manual');

	fig_handle = fig;
end

function p_out = del_dots(p_in)
	row_del = [];
	for i=1:size(p_in, 2)
		if (strcmp(p_in{i}, '.') == 1) || ...
			(strcmp(p_in{i}, '..') == 1)
			row_del = [row_del, i];
		end
	end
	p_in(row_del) = [];
	p_out = p_in;
end

function p_out = path_finder(p_in)
	contains = @(str, pattern) ~cellfun('isempty', strfind(str, pattern));
	p_out = cell(2, size(p_in, 2));
	for i=1:size(p_in, 2)
		file_list = dir(p_in{i});
		file_list = {file_list.name};
		row_del = [];
		for j=1:size(file_list, 2)
			if (strcmp(file_list{j}, '.') == 1) || ...
				(strcmp(file_list{j}, '..') == 1)
				row_del = [row_del, j];
				continue
			elseif ~contains({file_list{j}}, 'CalData')
				row_del = [row_del, j];
			end
		end
		file_list(row_del) = [];
		i
		p_in{i}
		file_list
		p_out{1, i} = strcat(p_in{i}, '/', file_list{1});
		p_out{2, i} = strcat(p_in{i}, '/', file_list{2});
	end
end

fib_path = {fib_path.name};
aero_path = {aero_path.name};

fib_path = strcat("Proc/fib/", del_dots(fib_path))
aero_path = strcat("Proc/aero/", del_dots(aero_path))
fib_path = path_finder(fib_path);
aero_path = path_finder(aero_path);

fig_store = cell(1, size(fib_path, 2));
for i=1:size(fib_path, 2)
	fig_store{i} = read_n_plot(fib_path{2,i}, fib_path{1,i}, ...
		aero_path{2,i}, aero_path{1,i});
end

lgd_fig = figure();
set(lgd_fig,'position',[10,10,255,255]);
lgd_ax = axes(lgd_fig);
p1 = plot(lgd_ax, 1, 1, ...
	'--k','LineWidth',1,'DisplayName','Theoretical Response'); hold on;
p2 = plot(lgd_ax, 1, 1, ...
	'-.k','LineWidth',1, "DisplayName", "Fibre Fit");
p3 = plot(lgd_ax, 100, 100, ...
	'xk','LineWidth',1,'DisplayName','Fibre Data','Tag','OFF');
p4 = plot(lgd_ax, 1, 1, ...
	'-k','LineWidth',1, "DisplayName", "Aerosol Fit");
p5 = plot(lgd_ax, 100, 100, ...
	'*k','LineWidth',1,'DisplayName','Aerosol Data','Tag','OFF');
lgd = legend(lgd_ax, 'show');
set(lgd,'fontname','times');
set(lgd, 'Box', 'off')
set(lgd, "fontsize", 8)
ylim(lgd_ax, "manual");
ylim(lgd_ax, [0, 1])
xlim(lgd_ax, "manual");
xlim(lgd_ax, [0, 1])
axis off;

fig_store = [fig_store, {lgd_fig}];

for i=1:size(fig_store, 2)
	figure(fig_store{i})
	print("-dsvg", strcat("./plots/", num2str(i), ".svg"))
end

close all

%fig = figure();
%subs = {subplot(3,2,1), subplot(3,2,2), subplot(3,2,3), ...
%	subplot(3,2,4), subplot(3,2,5)};
%hold on
%for i=1:size(fig_store, 2)
%	for j=1:size(fig_store{i}, 2)
%		fig_objs = get(fig_store{i}{j}, 'children');
%		for k=1:size(fig_objs)
%			copyobj(fig_objs(k), subs{i});
%		end
%	end
%	ax = subs{i};	
%end























