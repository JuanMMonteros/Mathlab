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

%-- plots --$
config_s.plots=1;
% -- Tx --
config_s.tx_s.n_sym                  = 1e5     ;   % Number of symbols
config_s.tx_s.n_pil                  = 100     ;   % Number of pilots symbols
config_s.tx_s.n_payload              = 127     ;   % Symbols between pilots

config_s.tx_s.M                      = 16      ;   % Modulation levels
config_s.tx_s.BR                     = 32e9   ;   % Baud rate [Bd]

% -- Laser --
config_s.lw                          = 100e3     ;   % Tx LW [Hz]

% -- Noise --
config_s.ebno_db                     = 10      ;   % EbNo [dB]

% -- CPR --
config_s.n_taps_cpr                  = 65      ;   % Taps to filter CPR error

% Valores para barrer
payload = [10,20,50,100,150];
filter_lengths = [5,9,13,15,21];
lw=[0, 10e3, 100e3, 500e3, 1e6];
target_ber = 2e-2;  % Target BER
n_lw = length(lw);
n_filter = length(filter_lengths);
n_pl = length(payload);

out_c = cell(n_filter,n_pl);
ber_sim_v = zeros(n_filter,n_pl);
snr_loss_v = zeros(n_filter,n_pl);
for idx_lw = 1:n_lw
    for idx_pl = 1:n_pl
        parfor idx_filter = 1:n_filter
            cfg_s = config_s;
            cfg_s.lw = lw(idx_lw);
            cfg_s.tx_s.n_payload = payload(idx_pl);
            if filter_lengths(idx_filter) > 10
                ebno = 6.7;  % Es el valor teorico 
            else 
                ebno = 11;
            end
            while true
                cfg_s.n_taps_cpr = filter_lengths(idx_filter);
                cfg_s.ebno_db = ebno;
                o_data = cpr_sim(cfg_s);
                
                if o_data.ber_est > target_ber+0.001
                    ebno = ebno + 0.1;
                elseif o_data.ber_est < target_ber - 0.001
                    ebno = ebno - 0.1;
                else
                    break;  % Exit the loop when target BER is achieved
                end
            end
            
            ber_sim_v(idx_filter, idx_pl) = o_data.ber_est;
            snr_loss_v(idx_filter, idx_pl) = o_data.ebno_db + 10*log10(log2(cfg_s.tx_s.M));
        end
    end
    file_name = strcat(out_dir, 'o_data_', num2str(lw(idx_lw)), '.mat');
    save(file_name, 'ber_sim_v','snr_loss_v','filter_lengths','payload','config_s');
end


