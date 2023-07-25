clear all
%close all

nexp=10e3;
nsamples=1000;
Ts=1e-2;
tline=Ts.*(0:nsamples-1);
omega0=2*pi*20;

phi = 2*pi*rand(nexp,1)-pi;
x = cos(tline+phi);
save('for_tp/signal_exercise_1.mat')

% phi=rand(nexp,1);
% x = phi.*[0:nsamples-1];
% save('for_tp/signal_exercise_2.mat')


%%
t2=500;
term1=x(:,t2);
z = x.*conj(term1);
correlation_for_t1 = sum(z,1)/nexp;

% figure
% plot(mean(x,1))

figure
plot(correlation_for_t1(1+t2:end));
hold all
