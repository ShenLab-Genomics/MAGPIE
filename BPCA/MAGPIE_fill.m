data = csvread(filename, 1, 7);% skip header and info columns
data = fillmissing(data, 'constant', 999.0);
data_filled = BPCAfill(data);
[~, name, ~] = fileparts(filename);
rootPath = '../data/temp/';
currentFolder = pwd;
fullPath = fullfile(currentFolder, rootPath, name+"_bpca.csv");
writematrix(data_filled, fullPath)
