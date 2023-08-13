fz = 15;

n_exp = 10e3; % cantidad de experimentos 
n_sam = 1500; % cantidad de muestras
mu_x = 1; % media 
var_x = 2; % varianza

x_m = sqrt(var_x) * randn(n_sam, n_exp) + mu_x; % WGN

% Filter
n_taps = 25;
h_v = ones(1, n_taps) / n_taps; % Filtro de promedios moviles
h_mf_v = conj(h_v(end:-1:1));
y_m = filter(h_v, 1, x_m);

%% 

%Autocovariance of X
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
tit = sprintf('Autocovariance of X');
title(tit,'Interpreter','latex','FontSize', fz);
xlabel('Tau', 'Interpreter','latex','FontSize', fz);
ylabel('Amplitude', 'Interpreter','latex','FontSize', fz);
legend({'Theo','Estimated'}, 'Interpreter','latex','FontSize', fz-2);
grid on;
set(gcf, 'Position', [550 50 500 500],'Color', 'w');

%% Autocorrelation of X
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
tit = sprintf('Autocorrelation of X');
title(tit,'Interpreter','latex','FontSize', fz);
xlabel('Tau', 'Interpreter','latex','FontSize', fz);
ylabel('Amplitude', 'Interpreter','latex','FontSize', fz);
legend({'Theo','Estimated'}, 'Interpreter','latex','FontSize', fz-2);
grid on;
set(gcf, 'Position', [550 50 500 500],'Color', 'w');

%% Autocorrelation of Y
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
tit = sprintf('Autocorrelation of Y');
title(tit,'Interpreter','latex','FontSize', fz);
xlabel('Tau', 'Interpreter','latex','FontSize', fz);
ylabel('Amplitude', 'Interpreter','latex','FontSize', fz);
legend({'Theo','Estimated'}, 'Interpreter','latex','FontSize', fz-2);
grid on; xlim([-n_taps,n_taps]);
set(gcf, 'Position', [550 50 500 500],'Color', 'w');


%% Hope of Y
% mu_y = mu_x * sum(h(t)) = mu_x * H(f=0);
mu_y_theo_v = mu_x * sum(h_v) * ones(n_sam, 1);
mu_y_est_v = 1/n_exp * sum(y_m, 2);

% Plot
figure
plot(mu_y_theo_v ,'-r','Linewidth',2);
hold on;
plot(mu_y_est_v ,'--b','Linewidth',2);
tit = sprintf('Hope of Y');
title(tit,'Interpreter','latex','FontSize', fz);
xlabel('Time', 'Interpreter','latex','FontSize', fz);
ylabel('Amplitude', 'Interpreter','latex','FontSize', fz);
legend({'Theo','Estimated'}, 'Interpreter','latex','FontSize', fz-2);
grid on; 
set(gcf, 'Position', [550 50 500 500],'Color', 'w');

%% Variance of y
var_y_theo_v = zeros(n_sam,1);
var_y_theo_v(:)=var_x*sum(abs(h_v).^2);

% Estimation
g_tau_m=zeros(n_sam,n_exp);

for idx_exp =1:n_exp
    y_v= y_m(:,idx_exp);
    g_tau_m(:,idx_exp) = abs(y_v-mu_y_theo_v).^2;
end
var_y_est_v= 1/n_exp *sum(g_tau_m,1);

%Plot
figure
plot(var_y_theo_v,'-r','Linewidth',2);
hold on;
plot(var_y_est_v ,'--b','Linewidth',2);
tit = sprintf('Varianza of Y');
title(tit,'Interpreter','latex','FontSize', fz);
xlabel('Time', 'Interpreter','latex','FontSize', fz);
ylabel('Amplitude', 'Interpreter','latex','FontSize', fz);
legend({'Theo','Estimated'}, 'Interpreter','latex','FontSize', fz-2);
grid on; 
set(gcf, 'Position', [550 50 500 500],'Color', 'w');

