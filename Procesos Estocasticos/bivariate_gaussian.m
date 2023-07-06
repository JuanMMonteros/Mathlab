%-----------------------------------------------------------------------------%
%				                FUNDACION FULGOR
%-----------------------------------------------------------------------------%

% This script shows how to generate a bivariate gaussian function

clear;      % Clear all variables
close all;  % Close all figures
clc;        % Clear command window

% Fix the seed to simulate always the same values
rng(1)

%%

% Mean
mu_x = 0;
mu_y = 0;
mu_v = [mu_x, mu_y]; % Mean vector [meanX, meanY]

% Var
var_x = 0.2;
cov = 0.1; 
var_y = 1;
var_m = [var_x, cov; cov, var_y]; % Variance matrix [varX, cov; cov, varY]
sigma_m = sqrt(var_m);

% Variable vector
x_v = (-10:0.2:10)*sqrt(var_x) + mu_x;
x_len = length(x_v);
y_v = (-10:0.2:10)*sqrt(var_y) + mu_y;
y_len = length(y_v);

[X_m, Y_m] = meshgrid(x_v, y_v);
Z_m = [X_m(:), Y_m(:)];

% Generate pdf
F_m = mvnpdf(Z_m, mu_v, sigma_m);
F_m = reshape(F_m, y_len, x_len);

% 3D Plot
fz = 15;
figure; 
surf(x_v, y_v, F_m);
xlabel("x", 'Interpreter','latex','FontSize', fz);
ylabel('y', 'Interpreter','latex','FontSize', fz);
zlabel('pdf', 'Interpreter','latex','FontSize', fz);
title('Bivariate gaussian PDF', 'Interpreter','latex','FontSize', fz);
set(gcf, 'Position', [50 50 600 600],'Color', 'w');