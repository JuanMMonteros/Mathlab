clc; close all; clear
%%
fz = 15;
n_exp=1;
n_sam = 100e3; % numero de muestras
mu_x = 0; %media 
var_x = 2; %Variansa
x_v = sqrt(var_x) * randn(n_sam,n_exp) + mu_x; % WGN
% Filter
n_taps = 25;
h_v = ones(1, n_taps) / n_taps; %Filtro de promedios moviles
y_v = filter(h_v, 1, x_v);
% Plot
figure
plot(x_v ,'-r','Linewidth',2);
hold on;
plot(y_v ,'--b','Linewidth',2);
tit = sprintf('Signals in time');
title(tit,'Interpreter','latex','FontSize', fz);
xlabel('Discrete time', 'Interpreter','latex','FontSize', fz);
ylabel('Amplitude', 'Interpreter','latex','FontSize', fz);
legend({'x(n)','y(n)'}, 'Interpreter','latex','FontSize', fz-2);
grid on;
set(gcf, 'Position', [550 50 500 500],'Color', 'w');
%% PDS
NFFT = 1024;
win_v = hanning(128);
noverlap = 0;

Sx_theo_v = var_x/(2*pi) * ones(NFFT, 1);
[Sx_est_v, wd_v]=pwelch(x_v-mean(x_v), win_v, noverlap, NFFT, 'twoside');

Sy_theo_v = Sx_theo_v .* abs(fft(h_v,NFFT)).^2;
[Sy_est_v, ~]=pwelch(y_v-mean(y_v), win_v, noverlap, NFFT, 'twoside');

% Plot
figure
plot(wd_v,Sx_theo_v,'-k','Linewidth',1.5);
hold on;
plot(wd_v,Sy_theo_v,'-m','Linewidth',1.5);
plot(wd_v,Sx_est_v,'--r','Linewidth',2);
plot(wd_v,Sy_est_v,'--b','Linewidth',2);
tit = sprintf('PSDs');
title(tit,'Interpreter','latex','FontSize', fz);
xlabel('Discrete frequency [rad]', 'Interpreter','latex','FontSize', fz);
ylabel('PSD amplitude', 'Interpreter','latex','FontSize', fz);
legend({'Sx','Sy','\^{Sx}','\^{Sy}'},'Interpreter','latex','location','s',...
                                                        'FontSize', fz-2);
grid on; xlim([0,2*pi])
set(gcf, 'Position', [550 50 500 500],'Color', 'w');

% Compute variance
x_var = var(x_v);
x_psd_area = sum(Sx_est_v)*(2*pi/NFFT);

y_var = var(y_v);
y_psd_area = sum(Sy_est_v)*(2*pi/NFFT);

fprintf('\n - x(n) variance: %.2f', x_var)
fprintf('\n - Sx area: %.2f\n\n', x_psd_area)

fprintf('\n - y(n) variance: %.2f', y_var)
fprintf('\n - Sy area: %.2f\n\n', y_psd_area)