%-----------------------------------------------------------------------------%
%				                FUNDACION FULGOR
%-----------------------------------------------------------------------------%

% This script shows how to generate different random distributions

clear;      % Clear all variables
close all;  % Close all figures
clc;        % Clear command window

% Fix the seed to simulate always the same values
rng(1)

%% Discrete uniform random distribution

% Parameters
n_du = 1e4;
min_du = 1;
max_du = 6;
n_levels = max_du - min_du + 1;
level_v = 0:n_levels-1;

% Generate values
x_du_v = randi(n_levels, 1, n_du) + min_du - 1;

% PMF and CDF
n_plot = 10000;
x_v = linspace(min_du, max_du, n_plot); % Vector used to plot theo pdf and cdf
x_step = n_plot/(n_levels-1);
th_pmf_du_v = zeros(n_plot,1);
for idx = 1:n_levels
    idx_impulse = fix((level_v(idx))*x_step) + 1;
    idx_impulse(idx_impulse>n_plot) = n_plot;
    th_pmf_du_v(idx_impulse) = 1/n_levels;
end
th_cdf_du_v = cumsum(th_pmf_du_v); 

% Plot
figure('Name','Discrete uniform random distribution'); 
fz = 15;

% Histogram
subplot(1,3,1);
histogram(x_du_v, n_levels);
xlabel("x", 'Interpreter','latex','FontSize', fz);
ylabel('Repetitions', 'Interpreter','latex','FontSize', fz);
title('Histogram', 'Interpreter','latex','FontSize', fz);
grid on; 

% PDF/PMF
subplot(1,3,2);
histogram(x_du_v, n_levels,'Normalization','pdf','BinMethod','integers');
hold on; grid on;
plot(x_v,th_pmf_du_v,'-or')
xlabel("x", 'Interpreter','latex','FontSize', fz);
ylabel('Px', 'Interpreter','latex','FontSize', fz);
title('PDF/PMF', 'Interpreter','latex','FontSize', fz);
legend({'Sim','Th'},'Location','se','Interpreter','latex','FontSize', fz-2);
 
% CDF
subplot(1,3,3);
histogram(x_du_v, n_levels,'Normalization','cdf','BinMethod','integers');
hold on; grid on;
plot(x_v,th_cdf_du_v,'r')
xlabel("x", 'Interpreter','latex','FontSize', fz);
ylabel('Fx', 'Interpreter','latex','FontSize', fz);
title('CDF', 'Interpreter','latex','FontSize', fz);
legend({'Sim','Th'},'Location','se','Interpreter','latex','FontSize', fz-2);
 
set(gcf, 'Position', [50 50 1200 400],'Color', 'w');

%% Continuous uniform random distribution

clear;

% Parameters
n_cu = 1e4;
min_cu = 1;
max_cu = 6;
n_levels = max_cu - min_cu + 1;

% Generate values
x_cu_v = (max_cu-min_cu)*rand(n_cu,1) + min_cu;

% PMF and CDF
n_plot = 1000;
x_v = linspace(min_cu, max_cu, n_plot); % Vector used to plot theo pdf and cdf
dx = n_levels / n_plot;
th_pdf_cu_v = 1/n_levels*ones(n_plot,1);
th_cdf_cu_v = cumsum(th_pdf_cu_v) * dx; 

% Plot
figure('Name','Continuous uniform random distribution'); 
fz = 15;
hist_size = 50;

% Histogram
subplot(1,3,1);
histogram(x_cu_v, hist_size);
xlabel("x", 'Interpreter','latex','FontSize', fz);
ylabel('Repetitions', 'Interpreter','latex','FontSize', fz);
title('Histogram', 'Interpreter','latex','FontSize', fz);
grid on; 

% PDF/PMF
subplot(1,3,2);
histogram(x_cu_v, hist_size,'Normalization','pdf');
hold on; grid on;
plot(x_v,th_pdf_cu_v,'-r')
xlabel("x", 'Interpreter','latex','FontSize', fz);
ylabel('Px', 'Interpreter','latex','FontSize', fz);
title('PDF/PMF', 'Interpreter','latex','FontSize', fz);
legend({'Sim','Th'},'Location','se','Interpreter','latex','FontSize', fz-2);
 
% CDF
subplot(1,3,3);
histogram(x_cu_v, hist_size,'Normalization','cdf');
hold on; grid on;
plot(x_v,th_cdf_cu_v,'r')
xlabel("x", 'Interpreter','latex','FontSize', fz);
ylabel('Fx', 'Interpreter','latex','FontSize', fz);
title('CDF', 'Interpreter','latex','FontSize', fz);
legend({'Sim','Th'},'Location','se','Interpreter','latex','FontSize', fz-2);
 
set(gcf, 'Position', [50 50 1200 400],'Color', 'w');

%% Gaussian Random Distribution

clear;

% Parameters
n_ga = 1e4;
mean_ga = 5;
var_ga = 6;

% Generate values
x_ga_v = sqrt(var_ga).*randn(n_ga,1) + mean_ga;

% PMF and CDF
n_plot = 1000;
x_v = linspace(min(x_ga_v),max(x_ga_v), n_plot); % Vector used to plot 
dx = (max(x_ga_v) - min(x_ga_v)) / n_plot;
th_pdf_ga_v = 1/(sqrt(2*pi*var_ga)).*exp(-((x_v-mean_ga).^2)./(2*var_ga));
th_cdf_ga_v = cumsum(th_pdf_ga_v) * dx; 

% Plot
figure('Name','Gaussian Distribution'); 
fz = 15;
hist_size = 50;

% Histogram
subplot(1,3,1);
histogram(x_ga_v, hist_size);
xlabel("x", 'Interpreter','latex','FontSize', fz);
ylabel('Repetitions', 'Interpreter','latex','FontSize', fz);
title('Histogram', 'Interpreter','latex','FontSize', fz);
grid on; 

% PDF/PMF
subplot(1,3,2);
histogram(x_ga_v, hist_size,'Normalization','pdf');
hold on; grid on;
plot(x_v,th_pdf_ga_v,'-r')
xlabel("x", 'Interpreter','latex','FontSize', fz);
ylabel('Px', 'Interpreter','latex','FontSize', fz);
title('PDF/PMF', 'Interpreter','latex','FontSize', fz);
legend({'Sim','Th'},'Location','se','Interpreter','latex','FontSize', fz-2);
 
% CDF
subplot(1,3,3);
histogram(x_ga_v, hist_size,'Normalization','cdf');
hold on; grid on;
plot(x_v,th_cdf_ga_v,'r')
xlabel("x", 'Interpreter','latex','FontSize', fz);
ylabel('Fx', 'Interpreter','latex','FontSize', fz);
title('CDF', 'Interpreter','latex','FontSize', fz);
legend({'Sim','Th'},'Location','se','Interpreter','latex','FontSize', fz-2);
 
set(gcf, 'Position', [50 50 1200 400],'Color', 'w');

%% Rician distribution

% Rician distribution is the absolute value of a complex normal
% distribution where each component R & I have variance "var_ri"
% The complex mean is "mean_ri"

clear;

% Parameters
n_ri = 1e4;
mean_ri = 2 + 1j;
var_ri = 1;

% Generate values
x_nr_v = sqrt(var_ri) .* (randn(n_ri,1)) + real(mean_ri);
x_ni_v = sqrt(var_ri) .* (randn(n_ri,1)) + imag(mean_ri);
x_ri_v = abs(x_nr_v + 1j*x_ni_v);

% PMF and CDF
n_plot = 1000;
x_v = linspace(min(x_ri_v),max(x_ri_v), n_plot); % Vector used to plot
dx = (max(x_ri_v) - min(x_ri_v)) / n_plot;
mean_abs = abs(mean_ri);
th_pdf_ri_v = x_v/var_ri .* exp(-(x_v.^2+mean_abs.^2)/(2*var_ri)) .* ...
                                                besseli(0,x_v.*mean_abs/var_ri);
th_cdf_ri_v = cumsum(th_pdf_ri_v) * dx; 

% Plot
figure('Name','Rician distribution'); 
fz = 15;
hist_size = 50;

% Histogram
subplot(1,3,1);
histogram(x_ri_v, hist_size);
xlabel("x", 'Interpreter','latex','FontSize', fz);
ylabel('Repetitions', 'Interpreter','latex','FontSize', fz);
title('Histogram', 'Interpreter','latex','FontSize', fz);
grid on; 

% PDF/PMF
subplot(1,3,2);
histogram(x_ri_v, hist_size,'Normalization','pdf');
hold on; grid on;
plot(x_v,th_pdf_ri_v,'-r')
xlabel("x", 'Interpreter','latex','FontSize', fz);
ylabel('Px', 'Interpreter','latex','FontSize', fz);
title('PDF/PMF', 'Interpreter','latex','FontSize', fz);
legend({'Sim','Th'},'Location','se','Interpreter','latex','FontSize', fz-2);
 
% CDF
subplot(1,3,3);
histogram(x_ri_v, hist_size,'Normalization','cdf');
hold on; grid on;
plot(x_v,th_cdf_ri_v,'r')
xlabel("x", 'Interpreter','latex','FontSize', fz);
ylabel('Fx', 'Interpreter','latex','FontSize', fz);
title('CDF', 'Interpreter','latex','FontSize', fz);
legend({'Sim','Th'},'Location','se','Interpreter','latex','FontSize', fz-2);
 
set(gcf, 'Position', [50 50 1200 400],'Color', 'w');

%% Rayleigh distribution

% Rayleigh distribution is the absolute value of a complex normal
% distribution where each component R & I have variance "var_ray"
% The complex mean is 0. 
% Rician with mean 0 -> Rayleigh

clear;

% Parameters
n_ray = 1e4;
var_ray = 1;

% Generate values
x_nr_v = sqrt(var_ray) .* (randn(n_ray,1));
x_ni_v = sqrt(var_ray) .* (randn(n_ray,1));
x_ray_v = abs(x_nr_v + 1j*x_ni_v);

% PMF and CDF
n_plot = 1000;
x_v = linspace(min(x_ray_v),max(x_ray_v), n_plot); % Vector used to plot
dx = (max(x_ray_v) - min(x_ray_v)) / n_plot;
th_pdf_ray_v = x_v/var_ray .* exp(-x_v.^2/(2*var_ray));
th_cdf_ray_v = cumsum(th_pdf_ray_v) * dx; 

% Plot
figure('Name','Rayleigh distribution'); 
fz = 15;
hist_size = 50;

% Histogram
subplot(1,3,1);
histogram(x_ray_v, hist_size);
xlabel("x", 'Interpreter','latex','FontSize', fz);
ylabel('Repetitions', 'Interpreter','latex','FontSize', fz);
title('Histogram', 'Interpreter','latex','FontSize', fz);
grid on; 

% PDF/PMF
subplot(1,3,2);
histogram(x_ray_v, hist_size,'Normalization','pdf');
hold on; grid on;
plot(x_v,th_pdf_ray_v,'-r')
xlabel("x", 'Interpreter','latex','FontSize', fz);
ylabel('Px', 'Interpreter','latex','FontSize', fz);
title('PDF/PMF', 'Interpreter','latex','FontSize', fz);
legend({'Sim','Th'},'Location','se','Interpreter','latex','FontSize', fz-2);
 
% CDF
subplot(1,3,3);
histogram(x_ray_v, hist_size,'Normalization','cdf');
hold on; grid on;
plot(x_v,th_cdf_ray_v,'r')
xlabel("x", 'Interpreter','latex','FontSize', fz);
ylabel('Fx', 'Interpreter','latex','FontSize', fz);
title('CDF', 'Interpreter','latex','FontSize', fz);
legend({'Sim','Th'},'Location','se','Interpreter','latex','FontSize', fz-2);
 
set(gcf, 'Position', [50 50 1200 400],'Color', 'w');

%% Reorder figures
figure(5);
figure(4);
figure(3);
figure(2);
figure(1);