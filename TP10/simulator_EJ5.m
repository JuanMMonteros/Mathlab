function [odata_s] = simulator_EJ5(config_s)

    if exist('config_s','var')
        
        % General
        
        fs_ch = config_s.fs_ch;
        fs_dsp = config_s.fs_dsp;
        
        n_fires = config_s.n_fires;
        
        en_noise = config_s.en_noise;
        en_plots = config_s.en_plots;
        
        % TX
        chirp_bw = config_s.chirp_bw;
        chirp_T = config_s.chirp_T;
        chirp_P = config_s.chirp_P;
        
        f0 = config_s.f0;
        pw_tx_dbm = config_s.pw_tx_dbm;
        
        % Target and channel
        range = config_s.range;
        speed = config_s.speed;
        
        range_max = config_s.range_max;
        speed_max = config_s.speed_max;
        
        % RX
        
        snr_db = config_s.snr_db;
        fft_zp = config_s.fft_zp;
        
        n_thr = config_s.n_thr;
        
    else
        
        % General
        
        fs_ch = 400e6; % Frecuancia del canal
        fs_dsp = 100e6;% Frecuencia del dsp
        n_fires = 1e2; % Disparos al target 
        en_noise = 0; % activador de  ruido 
        en_plots = 1; % activador de  graficos 
        
        % TX
        chirp_bw = 100e6; %
        chirp_T = 12e-6;
        chirp_P = 100;
        
        f0 = 24e9;
        pw_tx_dbm = 13;
        
        % Target and channel
        
        range = 200;
        speed = 70;
        range_1=50;
        speed_1 =30;
        range_max = 300;
        speed_max = 100;
        
        % RX
        
        snr_db = 16;
        fft_zp = 8;
        
        n_thr = 20;
        
    end        

    %% Parametros generales y variables dependientes
    
    if en_noise == 0
        n_fires = 1;
    end
    
    ts = 1/fs_ch;
    
    % Tx
    time_total = chirp_P*chirp_T;
    t_chirp_v = 0:ts:chirp_T - ts;
    t_v = 0:ts:time_total-ts;
    
    pw_tx = 1e-3*10^(pw_tx_dbm/10);
    
    chirp_slope = chirp_bw / chirp_T;
    chirp_len = length(t_chirp_v);
    
    % LO
    c = 3e8;
    lambda = c/f0;
    
    % Target
    tau_max = 2*range_max/c;
    fd_max = 2*speed_max/lambda;
    fr_max = tau_max * chirp_slope;
    fb_max = fr_max + fd_max;
    tau_max_samples = fix(tau_max*fs_dsp);
    
    t_meas = chirp_T - tau_max;
    
    tau = 2*range/c;
    fd = 2*speed/lambda;
    fr = tau * chirp_slope;
    fb = fr + fd;
    
    col_tar_loc = fix(fb*t_meas) + 1;
    row_tar_loc = fix(fd*chirp_T*chirp_P)+1;
    
    % Modifico un poco el target para que caiga en el centro de la celda
    fb = (col_tar_loc-1) / t_meas;
    fd = (row_tar_loc - 1) / (chirp_T*chirp_P);
    tau = (fb-fd) / chirp_slope;
    range = tau * c / 2;
    speed = fd * lambda / 2;
    %%
    delta_R=c*chirp_T/(2*chirp_bw*t_meas);
    delta_V=lambda/(2*chirp_P*chirp_T);
    range1=range + range_1*delta_R;
    speed1=speed + speed_1*delta_V;
    tau1 = 2*range1/c; 
    %Target Aux
    fd1 = 2*speed1/lambda;
    fr1 = tau1 * chirp_slope;
    fb1 = fr1 + fd1;
    
    col_tar_loc = fix(fb1*t_meas) + 1;
    row_tar_loc = fix(fd1*chirp_T*chirp_P)+1;
    
    % Modifico un poco el target para que caiga en el centro de la celda
    fb1 = (col_tar_loc-1) / t_meas;
    fd1 = (row_tar_loc - 1) / (chirp_T*chirp_P);
    tau1 = (fb1-fd1) / chirp_slope;
    range1 = tau1 * c / 2;
    speed1 = fd1 * lambda / 2;
    
    % Rx
    snr = 10^(snr_db/10);
    
    h_aaf_v = fir1(100, 1.5*fb_max/(fs_ch/2));
    h_aaf_delay = round((length(h_aaf_v)-1)/2);
    
    OVR = fix(fs_ch/fs_dsp);
    
    tp_m = zeros(n_fires, n_thr);
    fn_m = zeros(n_fires, n_thr);
    fp_m = zeros(n_fires, n_thr);
    tn_m = zeros(n_fires, n_thr);

    %% Warnings
    
    % Nyquist in channel
    if fs_ch/2 <= chirp_bw
        error('Ch Nyquist');
    end
    
    % Nyquist in DSP
    if fs_dsp/2 <= fb_max
        error('Ch Nyquist');
    end

    % Chirp length enough
    if t_meas <= 0
        error('T<tau_max');
    end
    
    %% Process
    
    % Tx
    freq_ramp_v = chirp_slope * t_chirp_v ;
    freq_ramps_v = repmat(freq_ramp_v, 1, chirp_P);
    
    phase_parabs_v = 2*pi*cumsum(freq_ramps_v)*ts;
    
    chirp_frame_v = exp(1j*phase_parabs_v);
    
    tx_v = sqrt(pw_tx) * chirp_frame_v;
    
    % Ch & Tg
    tau_samples = fix(tau*fs_ch);
    tau_samples1 = fix(tau1*fs_ch);
    rx_v0 = [zeros(1,tau_samples), tx_v(1:end-tau_samples)] .* ...
                                                    exp(-1j*2*pi*fd*t_v);
    rx_v1 = [zeros(1,tau_samples1), tx_v(1:end-tau_samples1)] .* ...
                                                    exp(-1j*2*pi*fd1*t_v);
    rx_v=rx_v1+rx_v0;
                                                
    % Rx
    mix_v = tx_v .* conj(rx_v);
    pw_mix = var(mix_v);
    
    if en_noise
        sn = pw_mix / snr * chirp_P * t_meas;
    else
        sn = 0;
    end
    
    pn = sn * fs_ch;
    
    psd_m = 0;
   
    est_range=zeros(n_fires,1);
    est_speed=zeros(n_fires,1);
    for idx = 1:n_fires
        
        if ~mod(idx, 50)
            fprintf('Running fire %d of %d ... \n', idx, n_fires);
        end
        
        n_v = sqrt(pn/2) * (randn(size(mix_v))+1j*randn(size(mix_v)));
        mix_noisy_v = mix_v + n_v;
        
        mix_fil_v = filter(h_aaf_v, 1, [mix_noisy_v, zeros(1,h_aaf_delay)]);
        mix_fil_v = mix_fil_v(h_aaf_delay+1:end);
        
        dsp_v = mix_fil_v(1:OVR:end);
        dsp_len = length(dsp_v);
       
        dsp_m = reshape(dsp_v, round(dsp_len/chirp_P), chirp_P).';
        dsp_cut_m = dsp_m(:, tau_max_samples+1:end);
        
        NFFT_zp_v = fft_zp * size(dsp_cut_m);
        NFFT_v = size(dsp_cut_m);
        
        fft_zp_m = fft2(dsp_cut_m, NFFT_zp_v(1), NFFT_zp_v(2)) / ...
                                                    (NFFT_v(1)*NFFT_v(2));
                                                
        
        N_beat_cells = round(NFFT_zp_v(2)*fb_max/fs_dsp);
        
        i = 1;
        for j = 1:fft_zp:NFFT_zp_v(1)
            fft_m(i,:) = fft_zp_m(j, 1:fft_zp:N_beat_cells);
            i = i+1;
        end
           % Encuentra el valor máximo en toda la matriz fft_m
           maxValue = max(abs(fft_zp_m), [], 'all');

           % Encuentra los índices de los valores máximos
            [rowIndex, colIndex, pageIndex] = ind2sub(size(fft_zp_m), find(abs(fft_zp_m) == maxValue));
            est_range(idx) = (colIndex-1)*fs_dsp / size(fft_zp_m, 2) /chirp_slope* c / 2;
            est_speed(idx) = (rowIndex-1) * fs_dsp / (size(fft_zp_m, 1))*(t_meas*chirp_P)*lambda/2;
        
        NFFT_v = size(fft_m);
        
        psd_m = psd_m + (abs(fft_m).^2/n_fires);
        
        if idx==1 && en_plots
            figure;
             X = linspace(0,  fb_max, size(fft_m, 2));
             Y = linspace(0, fd_max, size(fft_m, 1));


            % Graficar los datos
            surf(X,Y,abs(fft_m), 'EdgeColor', 'none');

            % Etiquetar los ejes
            xlabel('Beta freq');
            ylabel('Doppler freq');
            zlabel('Magnitud');
            xlim([0,fb_max]);
            ylim([0,fd_max]);

            % Agregar título
            title('Gráfico 3D de la transformada de Fourier 2D');
            
        end
        target_comb = sub2ind(NFFT_v, row_tar_loc, col_tar_loc);
        
        if idx == 1
            max_th = 1.2 * max(abs(fft_m),[],'all');
            min_th = max_th/n_thr;
            delta_th = (max_th-min_th) / n_thr;
            thr_v = min_th:delta_th:max_th;
        end
        
        for th_idx = 1:n_thr
            
            assert_comb_m = abs(fft_m) > thr_v(th_idx);
            
            detec_comb_m = find(assert_comb_m>0);
            
            tp_m(idx, th_idx) = ismember(target_comb, detec_comb_m);
            fp_m(idx, th_idx) = length(detec_comb_m) - tp_m(idx, th_idx);
        end
        
        fn_m(idx, :) = 1 - tp_m(idx, :);
        tn_m(idx, :) = (NFFT_v(1)*NFFT_v(2)) - fp_m(idx, :); 
        
    end
    
    % SNR estimation
    max_psd = max(psd_m, [], 'all');
    max_idx = find(psd_m == max_psd);
    [idx_r, idx_c] = ind2sub(size(psd_m), max_idx);
    if idx_c < NFFT_v(2)/2
        pn_est = mean(psd_m(idx_r, fix(end/2):end));
    else
        pn_est = mean(psd_m(idx_r, 1:fix(end/2)));
    end
    ps_est = max_psd - pn_est;
    snr_est = ps_est/pn_est;
    snr_est_db = 10*log10(snr_est);
    
    % PD & PFA estimation
    pd_est_v = sum(tp_m, 1) ./ (sum(tp_m, 1) + sum(fn_m, 1)); 
    pfa_est_v = sum(fp_m, 1) ./ (sum(fp_m, 1) + sum(tn_m, 1)); 
    
    % PD & PFA theo
    [pd_th_v, pfa_th_v] = myrocsnr(snr_db);
    %preccion 
    range_sim_prec = std(est_range);
    speed_sim_prec = std(est_speed);
    range_theo_prec = 1/(2*pi*t_meas)*sqrt(6/snr_est)*chirp_T*c/(2*chirp_bw);
    speed_theo_prec = 1/(2*pi*chirp_T*chirp_P)*sqrt(6/snr_est)*lambda/2;
   
    %% Output
    odata_s.snr_est_db = snr_est_db;
    odata_s.pd_est_v = pd_est_v;
    odata_s.pfa_est_v = pfa_est_v;
    odata_s.pd_th_v = pd_th_v;
    odata_s.pfa_th_v = pfa_th_v;
    odata_s.range_sim_prec = range_sim_prec;
    odata_s.speed_sim_prec = speed_sim_prec;
    odata_s.range_theo_prec = range_theo_prec;
    odata_s.speed_theo_prec = speed_theo_prec;
    
end
