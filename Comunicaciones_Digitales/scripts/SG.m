%clear all
%close all
rand('seed',1)
h_ref = fir1(4,0.3)';
Ncoeff=length(h_ref);
Lsim=100e3;

x=randn(Lsim, 1);
y_ref=filter(h_ref,1,x)+0.0025*randn(Lsim,1);

h=zeros(Ncoeff,1);
step=1e-2;
buffer=zeros(size(h));

y=zeros(size(x));
error=zeros(size(x));
mse_log=zeros(Lsim,1);
h_log=zeros(Lsim, Ncoeff);

for n=1:Lsim
    
    buffer(2:end)=buffer(1:end-1);
    buffer(1)=x(n);
    y(n)=sum(h.*buffer);
    
    error(n)=y(n)-y_ref(n);
    h=h-step*error(n)*buffer;
        
    h_log(n,:)=h;
    
end

plot(10*log10(filter(ones(2560,1)/2560,1,abs(error).^2)))
hold all
