%%

clear 
close all
clc;

% Parametros 
BR = 32e9;      %Bd
L = 10000;      % Simulation Length
N = 4;          % Oversampling rate
fs = N*BR;      % Sampling rate to emulate analog domain
T = 1/BR;       % Time interval between two consecutive symbols
Ts = 1/fs;      % Time between 2 conseutive samples at Tx output
rolloff = 0.3;  % Pulse shaping rolloff 
ps_taps = 10;   % Pulse shaping taps
M=4;            % M-QAM
k=log2(M);

% PARAMETROS DEL CLOCK
clk_phase = 0;
clk_ppm = 100;

%% Two symbols generation 
decsymbs = randi([0 M-1], L,1);
x=qammod(decsymbs, M);

% Upsampling to change sampling rate
xup = upsample(x,N);

% Filter to interpolate the signal
h = rcosine(BR, fs, 'sqrt', rolloff, ps_taps);
yup = filter(h, 1, xup); 

%% Errores de timing
Lyup = length(yup);
line = (0:Lyup-1).';

% Tiempo ideal
time_ideal =  line .* Ts;

% Modelo ppm
fs_real = fs * (1+clk_ppm*1e-6);
Ts_real = 1/fs_real;

% Modelo jitter random
time_real = clk_phase*1/BR + line.*Ts_real; 

% Aplico los fenomenos resampleando 
yrs = interp1(time_ideal,yup,time_real, 'spline', 0);

%% RECEPTOR
% MF
yrx = filter(h,1,yup);
yrx_rs = filter(h,1,yrs);

% Downsampling
ybd = downsample(yrx,N);
ybd_rs = downsample(yrx_rs,N);

%% CORRECCION DE ERRORES DE TIMING USANDO GARDNER

Ntaps = 63; % taps del filtro interpolador en simbolos
buffer = zeros(Ntaps,1);
Ti = Ts; % Implica TED ideal

interp_out = zeros(L*N,1); % Salida del interpolador
uk=0; % mu
uk_d=0; % mu anterior
signo_d = 0; % signo anterior
uk_log = zeros(L,1);
extra_mem=0; % Implica descartar o repetir una muestra
init_timer = Ntaps*2;
int_err=0; Wm=0;
timing_error = zeros(L,1);
ek_log = 0;

% TR PLL
Kp=8e-3; Ki=Kp/1000;

for m=1:L*N-300
    index=m+extra_mem;
    buffer = yrx_rs(index:index+Ntaps);
    taps = interp_filter(Ntaps,uk); %uk entre -0.5 y 0.5
    
    yf = conv(taps,buffer);
    yf=yf(Ntaps+1); 
    interp_out(m) = yf;
    
    if mod(m,N)==1 && m>init_timer
        
        % TED
        y_n_actual = interp_out(m);
        y_n_delay_1_muestra = interp_out(m-N);
        y_n_delay_media_muestra = interp_out(m-N/2);
        
        ek = real(conj(y_n_delay_media_muestra) * ...
                            (y_n_delay_1_muestra-y_n_actual));
        ek_log(end+1) = ek;
        
        % PLL
        prop_err = Kp*ek;
        int_err = int_err+Ki*ek;
        Wm = prop_err+int_err;
        timing_error(fix(m/N)) = Wm;
        
        uk=uk+Wm;
        uk_log(fix(m/N))=uk;
        overflow=0; underflow=0;
        if uk>0.5
            uk=uk-1;
            overflow=1;
        elseif uk<-0.5
            uk=uk+1;
            underflow=1;
        end
            
        if underflow==1
              extra_mem=extra_mem-1;
        elseif overflow==1
              extra_mem= extra_mem+1;  
        end
    end
end

%% PLOTS

eyediagram(yrx(L*N/2:end-1000),2*N);
title('Diagrama de ojo sin ppm')
eyediagram(yrx_rs(L*N/2:end-1000),2*N);
title('Diagrama de ojo con ppm')

figure
subplot(1,2,1)
plot(real(ybd(1:L/2)),'.');
title('Simbolos recibidos sin ppm')

subplot(1,2,2)
plot(real(ybd_rs(1:L/2)),'.');
title('Simbolos recibidos con errores')

figure; 
subplot(1,2,1)
plot(time_ideal,'r','Linewidth',2);
hold on;
plot(time_real, '--b','Linewidth',2);
legend('Tiempo ideal','Tiempo real');
title('Tiempo Real vs Ideal')
subplot(1,2,2)
plot(time_ideal-time_real, '-k','Linewidth',2);
title('Diferencia')

ybd2 = downsample(interp_out,N,0);
figure
s2plot=ybd2; Lsp = length(s2plot);
t=(0:Lsp-1).*Ts;
plot(real(s2plot),'.');
title('Simbolos recuperados luego de la correccion')

figure
plot((uk_log(1:end-1000)))
title('Entrada del interpolador');
