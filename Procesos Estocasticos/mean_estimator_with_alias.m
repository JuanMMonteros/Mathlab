%-----------------------------------------------------------------------------%
%				                FUNDACION FULGOR
%-----------------------------------------------------------------------------%

% This script computes the mean of a signal sampled with alias

clear;      % Clear all variables
close all;  % Close all figures
clc;        % Clear command window

%% Parameters

fs = 1e6;               % Sampling frequency [Hz]
Ts = 1/fs;              % Sampling time [s]

OVS = 4;                % Oversampling factor (fs_ch/fs)

var_r = 1;              % Noise variance
mean_r = 0.5;           % Real mean to estimate
 
n_samples = 1e5;        % Sim length in digital domain

NFFT = 4*n_samples;     % FFT size

fz = 15;                % Plots Font Size

%% Processing

% Signal: random values (mean 0) + mean
x_v = sqrt(var_r)*randn(n_samples, 1);
x_v = x_v - mean(x_v) + mean_r; 
t_v = (0:n_samples-1)*Ts;
f_v = (0:NFFT-1)*fs/NFFT;

% Downsampling factors
OVS_v = [1,10:10:50];

% Alloc memory for output
mean_est_v = zeros(1, length(OVS_v));
leg_c = cell(1,length(OVS_v));

for idx = 1:length(OVS_v)
    
    % Select factor
    OVS = OVS_v(idx);
    
    % Downsampling
    xd_v = x_v(1:OVS:end);
    xd_len = length(xd_v);
    fsd = fs/OVS;
    td_v = (0:xd_len-1)*1/fsd;
    fd_v = (-NFFT/2:NFFT/2-1)*fsd/NFFT;

    % Mean estimator
    mean_est_v(idx) = mean(xd_v);
    
    % Plot input in time
    figure(1); 
    leg_c{idx} = sprintf('OVS:%d',OVS);
    plot(td_v,xd_v,'Linewidth',1.5);
    hold on;
    
    % Plot signal in freq
    figure(2); 
    plot(fd_v,fftshift(abs(fft(xd_v,NFFT)))/(n_samples/OVS),'Linewidth',1.5);
    hold on;
    
    % Plot signal in freq (ZOOM)
    figure(3); 
    plot(fd_v,fftshift(abs(fft(xd_v,NFFT)))/(n_samples/OVS),'Linewidth',1.5);
    hold on;
    
end

%% Plots 

% Complete time plot
figure(1)
title('xd[n]','Interpreter','latex','FontSize', fz);
xlabel('Time [s]', 'Interpreter','latex','FontSize', fz);
ylabel('Amplitude', 'Interpreter','latex','FontSize', fz);
legend(leg_c,'Interpreter','latex','Location','sw','FontSize',fz-2);
grid on;
set(gcf, 'Position', [50 50 500 500],'Color', 'w');

% Complete freq plot
figure(2)
title('Xd[w]','Interpreter','latex','FontSize', fz);
xlabel('Freq [Hz]', 'Interpreter','latex','FontSize', fz);
ylabel('Amplitude', 'Interpreter','latex','FontSize', fz);
legend(leg_c,'Interpreter','latex','Location','nw','FontSize',fz-2);
grid on; xlim([-6e4,6e4])
set(gcf, 'Position', [50 50 500 500],'Color', 'w');

% Complete freq plot [ZOOM]
figure(3)
title('Xd[w]','Interpreter','latex','FontSize', fz);
xlabel('Freq [Hz]', 'Interpreter','latex','FontSize', fz);
ylabel('Amplitude', 'Interpreter','latex','FontSize', fz);
legend(leg_c,'Interpreter','latex','Location','nw','FontSize',fz-2);
grid on; xlim([-10,10])
set(gcf, 'Position', [50 50 500 500],'Color', 'w');

% Estimation compare
figure;
plot(OVS_v, mean_r*ones(1,length(OVS_v)) ,'--k','Linewidth',1);
hold on
plot(OVS_v, mean_est_v ,'-xr','Linewidth',2);
tit = sprintf('Mean estimator with alias');
title(tit,'Interpreter','latex','FontSize', fz);
xlabel('OVS', 'Interpreter','latex','FontSize', fz);
ylabel('Mean', 'Interpreter','latex','FontSize', fz);
legend({'Original','Measured'},'Interpreter','latex','Location','ne','FontSize',fz-2);
grid on; ylim([0,mean_r*2]); xticks(OVS_v)
set(gcf, 'Position', [50 50 500 500],'Color', 'w');