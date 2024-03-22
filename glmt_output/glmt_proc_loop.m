%folders = {"job1/", "job2/", "job3/", "job4/", "job5/"};
folders = {"job6/"};

data_store = cell(1, size(folders, 2));
i = 1;
for folder=folders
    sprintf("Processing data from %s", folder{1})
    weighted_data = apply_weight(folder{1});
    data_store{i} = weighted_data;
	i ++;
end

