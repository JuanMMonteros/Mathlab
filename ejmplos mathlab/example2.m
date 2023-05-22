%-----------------------------------------------------------------------------%
%				                FUNDACION FULGOR
%-----------------------------------------------------------------------------%

% The nomenclature "_v" at the end of a variable indicates that is a vector
% The nomenclature "_m" at the end of a variable indicates that is a matrix
% The nomenclature "_s" at the end of a variable indicates that is a struct
% The nomenclature "_c" at the end of a variable indicates that is a cell

clear;      % Clear all variables
close all;  % Close all figures
clc;        % Clear command window

%% Parameters

fs = 2e6;           % Sampling frequency
ts = 1/fs;          % Sampling time

n_samples = 1e6;

%% Processing

% Time vector
t_v = (0:n_samples-1)*ts;

% Signals creation
x1_v = 0;
x2_v = 0;
for idx = 1:5
    if mod(idx,2)
        x1_v = x1_v + exp(1j*2*pi*(10^idx)*t_v);
    else
        x2_v = x2_v + exp(1j*2*pi*(10^idx)*t_v);
    end
end

% Spectrum
NFFT = n_samples;
X1_v = fft(x1_v,NFFT)./n_samples;
X2_v = fft(x2_v,NFFT)./n_samples;

% Positive frequencies vector
f_pos_v = (0:NFFT/2-1)*(fs/2)/(NFFT/2);

%% Plots
fz = 15;

% Signal spectrum
figure;
semilogx(f_pos_v,abs(X1_v(1:end/2)),'-r','Linewidth',2)
hold on;
semilogx(f_pos_v,abs(X2_v(1:end/2)),'--b','Linewidth',2)

title('Signals in frequency','Interpreter','latex','FontSize', fz);
leg_c = {'X1','X2'};
legend(leg_c,'Location','nw','Interpreter','latex','FontSize', fz-2);
xlabel('Frequency [Hz]', 'Interpreter','latex','FontSize', fz);
ylabel('Module', 'Interpreter','latex','FontSize', fz);
grid on

set(gcf, 'Position', [50 50 700 500],'Color', 'w');