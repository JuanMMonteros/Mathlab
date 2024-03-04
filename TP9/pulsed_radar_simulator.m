
function [odata] = pulsed_radar_simulator(config_s)


    if exist('config_s','var')
        tau = config_s.tau;
        Po = config_s.Po;
        range = config_s.range;
        max_range = config_s.max_range;
        No = config_s.No;
        Niters = config_s.Niters;
        NOS = config_s.NOS;           % Factor de sobremuestreo del canal
    else
        tau = 5e-9;       % Ancho del pulso [s]
        Po = 5e3 ;      % Potencia instantanea [W]
        range=700;       % Distancia al target [m]
        max_range = 2.5e3;  % Distancia maxima alcanzable [m]
        No = 1*(.4e-9)^2; % PSD del ruido one-side [W/Hz]
        Niters = 2500;    % Iteraciones del simulador
        NOS = 16;           % Factor de sobremuestreo del canal
    end        

    %% Parametros generales
    
    c = 3e8;            % Velocidad de la luz

    GT = 35;            % Ganancia direccional de la antena TX [dBi]
    GR = 35;            % Ganancia direccional de la antena RX [dBi]
    
    f0 = 77e9;          % Frequencia de portadora [Hz]
    lambda = c/f0;      % Longitud de onda [m]
    RCS = 1;            % Radar cross section [m^2]

    gain_detector = 5;  % ganancia [-]

    
    fs = (1/tau)*NOS;   % Frecuencia de muestreo de la simulacion

    %% Transmisor
    
    Lpulso = tau*fs;
    t_tx = (0:Lpulso-1)*1/fs;
    x_t = ones(Lpulso,1);   % Pulso en el TX
    %plot(t_tx,x_t)
    s_t = sqrt(Po)*x_t;     % Señal transmitida
    %plot(t_tx,abs(s_t).^2)

    %% Canal
    
    % Atenuacion (alpha)
    % GT y GR estan en dBi (los tengo que pasar a escala lineal)
    gt_lin = 10^(GT/10);
    gr_lin = 10^(GR/10);
    alpha = sqrt(gt_lin*gr_lin*lambda^2*RCS/ ((4*pi)^3 * range^4));
    
    if alpha > 1
        alpha = 1;
    end
    
    fprintf("\n");
    fprintf("Potencia instantanea TX: %2.2f W \n", Po);
    fprintf("Channel Gain: %2.2f dB \n", 20*log10(alpha));
    fprintf("Channel Attenuation: %2.2f dB \n", -20*log10(alpha));

    % Delay
    delay0 = 2*range/c; % segundos
    delay_samples0 = round(delay0*fs);
    real_delay0 = delay_samples0/fs; % Es el que voy a usar para comparar la performancce

    % Cambio de fase
    phase_change0 = exp(1j*2*pi*f0*real_delay0);

    ch_output_aux0 = phase_change0.*alpha.*[zeros(delay_samples0,1); s_t];
    zeros_que_faltan0 = round(2*max_range/c*fs)-length(ch_output_aux0);
    ch_output0 = [ch_output_aux0; zeros(zeros_que_faltan0,1)];
    %plot(ch_output0);

%     % Aux para hacer otro target
%     delta_range = -300;
%     alpha1 = sqrt(gt_lin*gr_lin*lambda^2*RCS/ ((4*pi)^3 * (range+delta_range)^4));
%     delay1 = 2*(range+delta_range)/c; % segundos
%     delay_samples1 = round(delay1*fs);
%     real_delay1 = delay_samples1/fs; % Es el que voy a usar para comparar la performancce
%     phase_change1 = exp(1j*2*pi*f0*real_delay1);
%     ch_output_aux1 = phase_change1.*alpha1.*[zeros(delay_samples1,1); s_t];
%     zeros_que_faltan1 = round(2*max_range/c*fs)-length(ch_output_aux1);
%     ch_output1 = [ch_output_aux1; zeros(zeros_que_faltan1,1)];

    y_t = ch_output0;

%     figure
%     subplot(211);
%     plot(abs(s_t).^2)
%     subplot(212);
%     plot(abs(y_t).^2)

    Ptx_dBm= 10*log10(Po/1e-3);
    Prx = max(abs(y_t).^2);
    Prx_dBm= 10*log10(Prx/1e-3);
    channel_att_measured = Ptx_dBm - Prx_dBm;
    
    fprintf("Measured Channel Attenuation: %2.2f dB \n", channel_att_measured);
    
    snr_teo = gain_detector.^2 * Prx * tau/No;
    fprintf("SNR: %2.2fdB \n", 10*log10(snr_teo));
    fprintf("\n")
    
    max_delay = 2*max_range/c;
    max_delay_samples = round(max_delay*fs);
    real_max_delay = max_delay_samples/fs; 
    real_max_range = real_max_delay/2 * c; 

    %% Front End
    
    deltaR = tau * 3e8 / 2;
    Ncells = ceil(max_range/deltaR)-1;
    bw_noise = fs;
    noise_power = No*bw_noise;
    PRX_peak = max(abs(gain_detector*y_t).^2);
    
    % Aproximacion de threshold maximo para pulso rectangular
    max_threshold = (sqrt(PRX_peak) + 3.5*sqrt(noise_power/NOS)).^2; 
    Nthrs = 100;
    thresholds = (max_threshold/Nthrs:max_threshold/Nthrs:max_threshold);
    
    fn_vector = zeros(Nthrs,1);
    fp_vector = zeros(Nthrs,1);
    tp_vector = zeros(Nthrs,1);
    tn_vector = zeros(Nthrs,1);

    est_range = zeros(Niters,1);
    
    y_mf_esperanza = 0;
    for niter = 1:Niters

        noise = sqrt(noise_power/2)*(randn(size(y_t)) + 1j.*randn(size(y_t)));
        z_t = gain_detector*y_t + noise;

        % Matched filter
        h_mf = conj(x_t(end:-1:1));
        h_mf = h_mf/sum(h_mf); % Normalizo para energia unitaria
        y_mf = filter(h_mf, 1, z_t);

        [~, imax] = max(abs(y_mf));
        est_range(niter) = (imax-length(h_mf))*real_max_range/(length(y_mf));

    %     plot(abs(y_mf).^2);
    %     hold all

        phase_decim = mod(delay_samples0-1,NOS);
        y_mf_decim_sq = abs(y_mf(1+phase_decim:NOS:end)).^2;
        
        y_mf_esperanza = y_mf_esperanza + y_mf_decim_sq/Niters;

        cell_of_interest = 1+ceil(range/deltaR); % El +1 es por matlab
        % Cuando cuente PFA, descarto la primer celda y la ultima

        %
%         figure
%         timeline = 1/fs*(0:length(y_mf)-1);
%         plot(timeline, abs(y_mf).^2); hold all
%         timeline2 = 1/(fs)*(phase_decim:NOS:length(y_mf)-1);
%         plot(timeline2, y_mf_decim_sq,'o'); hold all

        % Contar TP, FN, TN, FP

        X = (y_mf_decim_sq>thresholds).';
        Z = X(:,cell_of_interest);
        tps = Z; 
        fns = 1-Z;
        fps = (sum(X(:,2:end-1),2)-Z);
        tns = Ncells-1-fps;

        tp_vector=tp_vector+tps;
        tn_vector=tn_vector+tns;
        fp_vector=fp_vector+fps;
        fn_vector=fn_vector+fns;

    end
    
    % Estimacion de SNR usando la esperanza de la salida del MF
    [value_max, idx_max] = max(y_mf_esperanza);
    if idx_max > length(y_mf_esperanza)/2
        Pn_est = mean(y_mf_esperanza(1:100));
    else
        Pn_est = mean(y_mf_esperanza(end-100:end));
    end
    Ps_est = value_max - Pn_est;
    snr_est = Ps_est/Pn_est;

    range_sim_prec = std(est_range);
    range_theo_res = tau/2*c;
    range_theo_prec = range_theo_res/snr_teo;
    
    pfa_vector = fp_vector./(fp_vector+tn_vector);
    pfa_vector( fp_vector<20)=0;
    pd_vector = tp_vector./(tp_vector+fn_vector);

    %plot(y_mf_decim_sq);  
    %plot(ones(size(y_mf)).*max_threshold);

    odata.thresholds = thresholds;
    odata.pd_vector=pd_vector;
    odata.pfa_vector=pfa_vector;
    odata.tp_vector=tp_vector;
    odata.tn_vector=tn_vector;
    odata.fp_vector=fp_vector;
    odata.fn_vector=fn_vector;
    odata.est_range = est_range;
    odata.snr_teo = snr_teo;
    odata.snr_est = snr_est;
    odata.range_sim_prec = range_sim_prec;
    odata.range_theo_prec = range_theo_prec;
    odata.range_theo_res = range_theo_res;

end
