%-----------------------------------------------------------------------------%
%				               FUNDACION FULGOR
%-----------------------------------------------------------------------------%

clc; close all; clear

%% Parameters and vectors creation 

% Fontsize
fz = 15;

% WGN
n_exp = 1e3;
n_sam = 1000;
mu_x = 0.5;
var_x = 5;

x_m = sqrt(var_x) * randn(n_sam, n_exp) + mu_x;

% Filter
h_v = rcosine(0.25,1,'normal',0.5,2);
n_taps = length(h_v);
h_mf_v = conj(h_v(end:-1:1));

y_m = filter(h_v, 1, x_m);

%% Autocovariance of x
% Cx(tau) = E{(x(t))*(x(t-tau))}
% g(tau) = x(t)*x(t-tau)
Cx_theo_v = zeros(2*n_sam-1, 1);
Cx_theo_v((length(Cx_theo_v)+1)/2) = Cx_theo_v((length(Cx_theo_v)+1)/2) + var_x;

% Estimation
g_tau_m = zeros(2*n_sam-1 , n_exp);

for idx_exp = 1:n_exp
    
    x_v = x_m(:,idx_exp);
    g_tau_m(:,idx_exp) = xcorr(x_v - mu_x, 'unbiased');
    
end

Cx_est_v = 1/n_exp * sum(g_tau_m, 2);
tau_v = (-n_sam+1:n_sam-1);

% Plot
figure
plot(tau_v,Cx_theo_v ,'-r','Linewidth',2);
hold on;
plot(tau_v,Cx_est_v ,'--b','Linewidth',2);
tit = sprintf('Autocovariance of x');
title(tit,'Interpreter','latex','FontSize', fz);
xlabel('Tau', 'Interpreter','latex','FontSize', fz);
ylabel('Amplitude', 'Interpreter','latex','FontSize', fz);
legend({'Theo','Estimated'}, 'Interpreter','latex','FontSize', fz-2);
grid on;
set(gcf, 'Position', [550 50 500 500],'Color', 'w');

%% Autocorrelation of x
% Rx(tau) = E{(x(t)-mu)*(x(t-tau)-mu)}
% g(tau) = x(t)*x(t-tau)
Rx_theo_v = mu_x^2 * ones(2*n_sam-1, 1);
Rx_theo_v((length(Rx_theo_v)+1)/2) = Rx_theo_v((length(Cx_theo_v)+1)/2) + var_x;

% Estimation
g_tau_m = zeros(2*n_sam-1 , n_exp);

for idx_exp = 1:n_exp
    
    x_v = x_m(:,idx_exp);
    g_tau_m(:,idx_exp) = xcorr(x_v, 'unbiased');
    
end

Rx_est_v = 1/n_exp * sum(g_tau_m, 2);
tau_v = (-n_sam+1:n_sam-1);

% Plot
figure
plot(tau_v, Rx_theo_v ,'-r','Linewidth',2);
hold on;
plot(tau_v, Rx_est_v ,'--b','Linewidth',2);
tit = sprintf('Autocorrelation of x');
title(tit,'Interpreter','latex','FontSize', fz);
xlabel('Tau', 'Interpreter','latex','FontSize', fz);
ylabel('Amplitude', 'Interpreter','latex','FontSize', fz);
legend({'Theo','Estimated'}, 'Interpreter','latex','FontSize', fz-2);
grid on;
set(gcf, 'Position', [550 50 500 500],'Color', 'w');

%% Mean of y
% mu_y = mu_x * sum(h(t)) = mu_x * H(f=0);
mu_y_theo_v = mu_x * sum(h_v) * ones(n_sam, 1);
mu_y_est_v = 1/n_exp * sum(y_m, 2);

% Plot
figure
plot(mu_y_theo_v ,'-r','Linewidth',2);
hold on;
plot(mu_y_est_v ,'--b','Linewidth',2);
tit = sprintf('Mean of y');
title(tit,'Interpreter','latex','FontSize', fz);
xlabel('Time', 'Interpreter','latex','FontSize', fz);
ylabel('Amplitude', 'Interpreter','latex','FontSize', fz);
legend({'Theo','Estimated'}, 'Interpreter','latex','FontSize', fz-2);
grid on; 
set(gcf, 'Position', [550 50 500 500],'Color', 'w');

%% Autocovariance of y
% Cy(tau) = E{(y(t))*(y(t-tau))}
% g(tau) = y(t)*y(t-tau)
h_conv_h_mf_v = conv(h_v,h_mf_v);
Cy_theo_v = conv(Cx_theo_v, h_conv_h_mf_v,'same');

% Estimation
g_tau_m = zeros(2*n_sam-1 , n_exp);

for idx_exp = 1:n_exp
    
    y_v = y_m(:,idx_exp);
    g_tau_m(:,idx_exp) = xcorr(y_v - mu_y_theo_v, 'unbiased');
    
end

Cy_est_v = 1/n_exp * sum(g_tau_m, 2);
tau_v = (-n_sam+1:n_sam-1);

% Plot
figure
plot(tau_v,Cy_theo_v ,'-r','Linewidth',2);
hold on;
plot(tau_v,Cy_est_v ,'--b','Linewidth',2);
tit = sprintf('Autocovariance of y');
title(tit,'Interpreter','latex','FontSize', fz);
xlabel('Tau', 'Interpreter','latex','FontSize', fz);
ylabel('Amplitude', 'Interpreter','latex','FontSize', fz);
legend({'Theo','Estimated'}, 'Interpreter','latex','FontSize', fz-2);
grid on; xlim([-n_taps,n_taps]);
set(gcf, 'Position', [550 50 500 500],'Color', 'w');

%% Autocorrelation of y
% Ry(tau) = E{(y(t)-mu)*(y(t-tau)-mu)}
% g(tau) = y(t)*y(t-tau)
h_conv_h_mf_v = conv(h_v,h_mf_v);
Ry_theo_v = conv(Rx_theo_v, h_conv_h_mf_v,'same');

% Estimation
g_tau_m = zeros(2*n_sam-1 , n_exp);

for idx_exp = 1:n_exp
    
    y_v = y_m(:,idx_exp);
    g_tau_m(:,idx_exp) = xcorr(y_v, 'unbiased');
    
end

Ry_est_v = 1/n_exp * sum(g_tau_m, 2);
tau_v = (-n_sam+1:n_sam-1);

% Plot
figure
plot(tau_v,Ry_theo_v ,'-r','Linewidth',2);
hold on;
plot(tau_v,Ry_est_v ,'--b','Linewidth',2);
tit = sprintf('Autocorrelation of y');
title(tit,'Interpreter','latex','FontSize', fz);
xlabel('Tau', 'Interpreter','latex','FontSize', fz);
ylabel('Amplitude', 'Interpreter','latex','FontSize', fz);
legend({'Theo','Estimated'}, 'Interpreter','latex','FontSize', fz-2);
grid on; xlim([-n_taps,n_taps]);
set(gcf, 'Position', [550 50 500 500],'Color', 'w');

%% Variance of y
% var_y = E{|y-mu_y|^2} = E{y^2} - mu_y^2 = Cy(0) =
% = [conv(Cx,h(t),h*(-t)) valued in 0]
% = var_x * [conv(h(t),h*(-t)) valued in 0] 
% = var_x * sum(abs(h(t))^2) 

% g(tau) =|y-mu_y|^2

var_y_theo_v = zeros(n_sam,1);
var_y_theo_v(:) = var_x * sum(abs(h_v).^2);

% Estimation
g_tau_m = zeros(n_sam , n_exp);

for idx_exp = 1:n_exp
    
    y_v = y_m(:,idx_exp);
    g_tau_m(:,idx_exp) = abs(y_v-mu_y_theo_v).^2;
    
end

var_y_est_v = 1/n_exp * sum(g_tau_m, 2);

% Plot
figure
plot(var_y_theo_v ,'-r','Linewidth',2);
hold on;
plot(var_y_est_v ,'--b','Linewidth',2);
tit = sprintf('Variance of y');
title(tit,'Interpreter','latex','FontSize', fz);
xlabel('Time', 'Interpreter','latex','FontSize', fz);
ylabel('Amplitude', 'Interpreter','latex','FontSize', fz);
legend({'Theo','Estimated'}, 'Interpreter','latex','FontSize', fz-2);
grid on; ylim([0,2*var_y_theo_v(end)])
set(gcf, 'Position', [550 50 500 500],'Color', 'w');

