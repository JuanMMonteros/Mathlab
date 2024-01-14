clear
close all 
clc
%% Output directory

out_dir = mfilename('fullpath');
out_dir = out_dir(1:end-length(mfilename));
out_dir = [out_dir, 'out/'];

if ~exist(out_dir,'dir')
    mkdir(out_dir);
end

%% General configuration
    % -- Tx --
    config_s.tx_s.n_sym                  = 1e5     ;   % Number of symbols
    config_s.tx_s.n_pil                  = 100     ;   % Number of pilots symbols
    config_s.tx_s.n_payload              = 127     ;   % Symbols between pilots
    
    config_s.tx_s.M                      = 16      ;   % Modulation levels
    config_s.tx_s.BR                     = 32e9   ;   % Baud rate [Bd]
    
    % -- Laser --
    config_s.lw                          = 100e3       ;   % LW [Hz]
    
    % -- Noise --
    config_s.ebno_db                     = 10       ;   % EbNo [dB] 
    
    % -- BPS --
    % -- BPS --
    cfg_s.bps_s.n_phases            = 16      ;   % BPS phases
    cfg_s.bps_s.filter_lengt        = 64      ;   % BPS filter len
    
    % -- CS corrector --
    config_s.en_cs_corrector             = 1       ;   % CS corrector for BPS path
    config_s.cs_cor_s.window_length      = 64      ;   % CS corrector window len
   %%
% Valores para barrer
phases = [4, 8, 12, 16];
filter_lengths = [5, 15, 39, 87, 131, 201];
lw=[0, 10e3, 100e3, 500e3, 1e6];
ebno=3;
n_lw = length(lw);
n_filter = length(filter_lengths);
n_ph = length(phases);

out_c = cell(n_filter,n_ph);
ber_sim_v = zeros(n_filter,n_ph);
ebno_sim_v = zeros(n_filter,n_ph);

for idx_lw = 1:n_lw
    cfg_s = config_s;
    cfg_s.lw = lw(idx_lw);
    for idx_ph = 1:n_ph
        cfg_s.bps_s.n_phases = phases(idx_ph);
        parfor idx_filter = 1:n_filter
            while ber_sim_v(idx_filter, idx_ph)<=1e-3 && ber_sim_v(idx_filter, idx_ph)>=5e-2
            cfg_s.bps_s.filter_length = filter_lengths(idx_filter);
            cfg_s.bps_s.ebno_db = ebno_db;
            o_data = bps_sim(cfg_s);
            if o_data.ber_est_bps >5e-2
                ebno=ebno-1
            else if o_data.ber_est_bps < 1e-3
                    ebno=ebno+1
                end
            end
            ber_sim_v(idx_filter, idx_ph) = o_data.ber_est_bps;
            ebno_sim_v(idx_filter, idx_ph) = o_data.ebno_db;
        end
    end
    file_name = strcat(out_dir, 'o_data_', num2str(lw(idx_lw)), '.mat');
    save(file_name, 'ber_sim_v','ebno_sim_v','filter_lengths','phases','cfg_s');
 end






