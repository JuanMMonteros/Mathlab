clear; close all; clc;

%% Output directory

out_dir = mfilename('fullpath');
out_dir = out_dir(1:end-length(mfilename));
out_dir = [out_dir, 'out/'];

if ~exist(out_dir,'dir')
    mkdir(out_dir);
end
ebno_ber2e2 = 6.7;
lw=500e3;
file_name = strcat('o_data_', num2str(lw), '.mat');
file = [out_dir, file_name];
load(file);
snr_loss=zeros(length(filter_lengths),1);
for idx=1:length(filter_lengths)
   snr_loss(idx)= snr_loss_v(idx,4);
end
filename = strcat(out_dir, 'bsp_16.mat');
save(filename, 'snr_loss');