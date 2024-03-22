glmt_data = {};
mie_data = {};
path = dir(".");

contains = @(str, pattern) ~cellfun('isempty', strfind(str, pattern));

path = {path.name};
row_del = [];
for i=1:size(path, 2)
    if ~isfolder(path{i}) == true
        row_del = [row_del, i];
    elseif ~contains({path{i}}, 'job')
        row_del = [row_del, i];
    end
end
path(row_del) = [];

p_out = cell(2, size(path, 2));
for i=1:size(path, 2)
    file_list = dir(path{i});
    file_list = {file_list.name};
    row_del = [];
    for j=1:size(file_list, 2)
        if ~contains({file_list{j}}, '.csv')
            row_del = [row_del, j];
        end
    end
    file_list(row_del) = [];
    p_out{1, i} = strcat(path{i}, '/', file_list{1});
    p_out{2, i} = strcat(path{i}, '/', file_list{2});
end
path = p_out

fig = figure(1)
set(gcf,'position',[10,10,454,300]);
%title('Effect of Beam Width on {\sigma}', "fontsize", 10)
ax1 = axes(fig);
set(ax1, 'Position', [0.12, 0.3, 0.38, 0.59]);
ax2 = axes(fig);
set(ax2, 'Position', [0.6, 0.3, 0.38, 0.59]);
lgd_ax = axes(fig);
set(lgd_ax, 'Position', [0.1, 0, 0.8, 0.15]);
plot(lgd_ax, 1, 1, '-k', 'LineWidth', 1, 'DisplayName', 'LMT'), hold on;
plot(lgd_ax, 1, 1, ':k', 'LineWidth', 1, 'DisplayName', ...
	'GLMT {b}_{w}={20}{\mu}m'), hold on;
plot(lgd_ax, 1, 1, '--k', 'LineWidth', 1, 'DisplayName', ...
	'GLMT {b}_{w}={40}{\mu}m'), hold on;
plot(lgd_ax, 1, 1, '-.k', 'LineWidth', 1, 'DisplayName', ...
	'GLMT {b}_{w}={100}{\mu}m'), hold on;

for i=1:size(path, 2)
    glmt_data{i} = dlmread(path{1, i}, ',');
    mie_data{i} = dlmread(path{2, i}, ',');
end
glmt_rad_40 = glmt_data{1}(:, 1);
glmt_scs_40 = glmt_data{1}(:, 2)/pi;
mie_rad = mie_data{1}(:, 1);
mie_scs = mie_data{1}(:, 2);
glmt_rad_100 = glmt_data{5}(:, 1);
glmt_scs_100 = glmt_data{5}(:, 2)/pi;
glmt_rad_20 = glmt_data{6}(:, 1);
glmt_scs_20 = glmt_data{6}(:, 2)/pi;

plot(ax1, pi*glmt_rad_40(10:end/2).^2/1e6, glmt_scs_40(10:end/2)/1e6, ...
	'--k', 'LineWidth', 1, 'DisplayName', 'GLMT {b}_{w}={40}_{\mu}m')
hold on
plot(ax1, pi*glmt_rad_100(10:end/2).^2/1e6, glmt_scs_100(10:end/2)/1e6, ...
	'-.k', 'LineWidth', 1, 'DisplayName', 'GLMT {b}_{w}={100}_{\mu}m')
hold on
plot(ax1, pi*glmt_rad_20(10:end/2).^2/1e6, glmt_scs_20(10:end/2)/1e6, ...
	':k', 'LineWidth', 1, 'DisplayName', 'GLMT {b}_{w}={20}_{\mu}m')
hold on
plot(ax1, pi*mie_rad(10:end/2).^2/1e6, mie_scs(10:end/2)/1e6, ...
	'-k', 'LineWidth', 1, 'DisplayName', 'LMT')
hold on
loglog(ax2, glmt_rad_40(10:end/2), glmt_scs_40(10:end/2)/1e6, ...
	'--k', 'LineWidth', 1, 'DisplayName', 'GLMT {b}_{w}={40}_{\mu}m')
hold on
loglog(ax2, glmt_rad_100(10:end/2), glmt_scs_100(10:end/2)/1e6, ...
	'-.k', 'LineWidth', 1, 'DisplayName', 'GLMT {b}_{w}={100}_{\mu}m')
hold on
loglog(ax2, glmt_rad_20(10:end/2), glmt_scs_20(10:end/2)/1e6, ...
	':k', 'LineWidth', 1, 'DisplayName', 'GLMT {b}_{w}={20}_{\mu}m')
hold on
loglog(ax2, mie_rad(10:end/2), mie_scs(10:end/2)/1e6, ...
	'-k', 'LineWidth', 1, 'DisplayName', 'LMT')
hold on

set(fig, "resize", "off")
set(fig, 'PaperPositionMode', 'manual');

set(ax1,'fontname','times');
set(ax1, "fontsize", 6)
set(ax2,'fontname','times');
set(ax2, "fontsize", 6)

ylabel(ax1, '{\sigma}_{sca} (mm^{2})', "fontsize", 8)
xlabel(ax1, '{\sigma}_{geo} (mm^{2})', "fontsize", 8)
ylabel(ax2, '{\sigma}_{sca} (mm^{2})', "fontsize", 8)
xlabel(ax2, 'Droplet Radius ({\mu}m)', "fontsize", 8)

lgd = legend(lgd_ax, 'show');
set(lgd, 'location', 'south');
set(lgd,'fontname','times');
set(lgd, 'Box', 'off')
set(lgd, "fontsize", 8)
set(lgd, "numcolumns", 2)
set(lgd_ax,'visible','off');

limx1 = get(ax1, 'xlim');
limy1 = get(ax1, 'ylim');

annotation('textbox',[.15 .82 .1 .1], ...
    'String','(a)','EdgeColor','none','fontname','times', "fontsize", 8)
annotation('textbox',[.65 .82 .1 .1], ...
    'String','(b)','EdgeColor','none','fontname','times', "fontsize", 8)
annotation('textbox',[.33 .92 .1 .1], ...
    'String','Effect of Beam Width on {\sigma}',...
	'EdgeColor','none','fontname','times', "fontsize", 10, ...
	'FontWeight', 'bold')

xlim(ax1, [limx1(1), pi*mie_rad(end/2).^2/1e6])
ylim(ax1, [limy1(1), mie_scs(end/2)/1e6])
xlim(ax2, [2e0, mie_rad(end/2)])
ylim(ax2, [3e-6, mie_scs(end/2)/1e6])

xticks(ax2, [3E0, 1E1, 3E1])

