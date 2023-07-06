clear all
%close all

%% Genero una VA gaussiana de media 1, varianza 1.2
Nexp=1000;

mean_log=zeros(Nexp,1);
var_log=zeros(Nexp,1);


for n=1:Nexp
    %randn('seed',1);
    Ngaus=20000;
    my_var = 1.2;
    sigma = sqrt(my_var);
    my_mean = 1;
    x = sigma*randn(Ngaus,1)+my_mean;
    mean_log(n) = mean(x);
    var_log(n) = var(x); % std(x)=sqrt(var(x))
end

var(mean_log)

%%
figure
hist(mean_log,100)
xlim([0.95,1.05])

%%
% figure
% hist(var_log,100)

%%
% x = randn(1000,1);
% y = 0.25*x.^2;
% 
% figure
% plot(x,y,'.')
% 
% corr(x,y)
