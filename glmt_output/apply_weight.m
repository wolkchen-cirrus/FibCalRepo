function scs_tab = apply_weight(data_dir)

	warning('off', 'all');
	inputs = data_dir;
	contains = @(str, pattern) ~cellfun('isempty', strfind(str, pattern));
	folder = fileparts(mfilename('fullpath'));
	addpath(genpath(folder));
	data_dir = strcat("./", inputs);
	path = dir(data_dir);
	path = {path.name};
	row_del = [];
	for i=1:size(path, 2)
		if isfolder(path{i}) == true
			row_del = [row_del, i];
		elseif ~contains({path{i}}, '.xls')
			row_del = [row_del, i];
		end
	end
	path(row_del) = [];
	[rad, scs_tot_glmt] = textread(strcat(data_dir, "CEXT_GLMT"));
	scs_tot_glmt = [rad, scs_tot_glmt];
	[rad, scs_tot_mie] = textread(strcat(data_dir, "CEXT_MIE"));
	scs_tot_mie = [rad, scs_tot_mie];
	c_sca_mie = zeros(size(path, 2), 1);
	c_sca_glmt = zeros(size(path, 2), 1);
	rad = zeros(size(path, 2), 1);
	for i=1:size(path, 2)
		data = dlmread(path{i}, '\t');
		if isempty(data)
			continue
		end
		len = floor(size(data, 1)/2);
		ang = flip(abs(data(1:len, 1)));
		p11 = flip(data(1:len, 3));
		for j=2:size(p11, 1)
			step = ang(j)*pi/180-ang(j-1)*pi/180;
			p11(j) = step*0.5*(p11(j)+p11(j-1));
		end
		p11 = p11./sum(p11);
		rad(i) = str2double(strjoin(regexp(path{i}, '\d*', 'Match'), '.'));
		disp(rad(i))
		[~, idx] = min(abs(scs_tot_mie(:, 1) - rad(i)), [], 1);
		sig_mie = scs_tot_mie(idx, 2);
		sig_dot_mie = sig_mie.*p11;
		[~, idx] = min(abs(scs_tot_glmt(:, 1) - rad(i)), [], 1);
		sig_glmt = scs_tot_glmt(idx, 2);
		sig_dot_glmt = sig_glmt.*p11;
		residuals = trapz(ang.*(pi/180), ...
			Weight(ang.*(pi/180)).*sin(ang.*(pi/180)));
		mie_out = sum(sig_dot_mie.*residuals);
		glmt_out = sum(sig_dot_glmt.*residuals);
		c_sca_mie(i) = mie_out;
		c_sca_glmt(i) = glmt_out;
	end
	scs_tab = [rad, c_sca_mie];
	scs_tab = scs_tab(all(scs_tab,2),:);
	scs_tab = sortrows(scs_tab);
	csvwrite(strcat(data_dir, "weighted_scs_mie.csv"), scs_tab)
	scs_tab = [rad, c_sca_glmt];
	scs_tab = scs_tab(all(scs_tab,2),:);
	scs_tab = sortrows(scs_tab);
	csvwrite(strcat(data_dir, "weighted_scs_glmt.csv"), scs_tab)

endfunction

