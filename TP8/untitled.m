% Parametros
BR = 32e9;      %Bd
L = 50000;      % Simulation Length
N = 4;          % Oversampling rate
fs = N*BR;      % Sampling rate to emulate analog domain
T = 1/BR;       % Time interval between two consecutive symbols
Ts = 1/fs;      % Time between 2 conseutive samples at Tx output

ps_taps = 10;   % Pulse shaping taps
M=4;            % M-QAM

rolloff_values = [0.1, 0.5, 0.8];
num_rolloffs = length(rolloff_values);

figure;

for i = 1:num_rolloffs
    rolloff = rolloff_values(i);
    
    % Two symbols generation
    decsymbs = randi([0 M-1], L, 1);
    x = qammod(decsymbs, M);

    % Upsampling to change sampling rate
    xup = upsample(x, N);

    % Filter to interpolate the signal
    h = rcosine(BR, fs, 'sqrt', rolloff, ps_taps);
    yup = filter(h, 1, xup);

    % Error de timing
    Lyup = length(yup);
    line = (0:Lyup-1).';

    % Tiempo ideal
    time_ideal = line .* Ts;
    
    % PARAMETROS DEL CLOCK
    clk_phase = 0.3;
    clk_ppm = 100;

    % Modelo ppm
    fs_real = fs * (1+clk_ppm*1e-6);
    Ts_real = 1/fs_real;

    % Modelo jitter random
    time_real = clk_phase*1/BR + line.*Ts_real;

    % Aplicar los fenómenos resampleando
    yrs = interp1(time_ideal, yup, time_real, 'spline', 0);

    % Receptor
    yrx = filter(h, 1, yup);
    yrx_rs = filter(h, 1, yrs);

    % Downsampling
    ybd = downsample(yrx, N);
    ybd_rs = downsample(yrx_rs, N);

    % Corrección de errores de timing usando Gardner
    Ntaps = 63;
    buffer = zeros(Ntaps, 1);
    Ti = Ts;

    interp_out = zeros(L * N, 1);
    uk = 0;
    uk_d = 0;
    signo_d = 0;
    uk_log = zeros(L, 1);
    extra_mem = 0;
    init_timer = Ntaps * 2;
    int_err = 0;
    Wm = 0;
    timing_error = zeros(L, 1);
    ek_log = 0;

    % TR PLL
    Kp = 0 * 8e-3;
    Ki = Kp / 1000;

    for m = 1:L * N - 300
        index = m + extra_mem;
        buffer = yrx_rs(index:index + Ntaps - 1);
        taps = interp_filter(Ntaps, uk);

        yf = conv(taps, buffer);
        yf = yf(Ntaps);

        interp_out(m) = yf;

        if mod(m, N) == 1 && m > init_timer
            
            %TED
            y_n_actual = interp_out(m);
            y_n_delay_1_muestra = interp_out(m - N);
            y_n_delay_media_muestra = interp_out(m - N / 2);

            ek = real(conj(y_n_delay_media_muestra) * ...
                (y_n_delay_1_muestra - y_n_actual));
            
            ek_log(end + 1) = ek;

            % PLL
            prop_err = Kp * ek;
            int_err = int_err + Ki * ek;
            Wm = prop_err + int_err;
            timing_error(fix(m / N)) = Wm;

            uk = uk + Wm;
            uk_log(fix(m / N)) = uk;
            overflow = 0;
            underflow = 0;

            if uk > 0.5
                uk = uk - 1;
                overflow = 1;
            elseif uk < -0.5
                uk = uk + 1;
                underflow = 1;
            end

            if underflow == 1
                extra_mem = extra_mem - 1;
            elseif overflow == 1
                extra_mem = extra_mem + 1;
            end
        end
    end

    subplot(num_rolloffs, 1, i);

    % Calcular la curva S
    differential_phase_error = ek_log;
    performance_curve = ek_log/T;

    % Graficar la curva S
    plot(performance_curve, differential_phase_error, 'LineWidth', 2);

    % Lengendas
    xlabel('Differential Phase Error (DPE)');
    ylabel('Performance del TED');
    title(['Curva S para rolloff = ' num2str(rolloff)]);
    grid on;  
end






