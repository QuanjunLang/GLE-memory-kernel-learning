

temp = matlab.desktop.editor.getActiveFilename;
parent_dir_temp = strsplit(temp, '/');
parent_dir = strjoin(parent_dir_temp(1:end-1), filesep);

fig_dir = [parent_dir, '/figure/'];
data_dir = [parent_dir, '/data/'];

addpath(genpath(fileparts(which(mfilename))));

