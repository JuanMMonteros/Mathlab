clear all
%close all

fs = 16e9; %GHz
Ts = 1/fs;
z = tf('z', Ts);

%%
L = 100000; % Simulation Length
t = [0:L-1].*Ts;

% System modelling
Kp = .005;
Ki = Kp/100;

H0 = Kp + Ki*z/(z-1);
NCO = z/(z-1);
H = feedback(H0*NCO,z^-0);

%% Step response
figure
step(H)

%% Ramp response
figure
step(H*z/(z-1))
hold all
step(z/(z-1))


% %%
figure
pzmap(H)

%%


