%clear all
%close all

h_ref = fir1(4,0.3)';
Ncoeff=length(h_ref);
Lbatch=1e3;

Niters=10000;
c=zeros(Ncoeff,1);
Dtap=0.0001;
gradient=zeros(Ncoeff,1);


mse_log=zeros(Niters,1);
gradient_log = zeros(Niters, Ncoeff);
h_log=zeros(Niters, Ncoeff);

for iter=1:Niters
    x=randn(Lbatch, 1);
    y_ref=filter(h_ref,1,x);

    y=filter(c,1,x);
    mse0=mean(abs(y-y_ref).^2);
    mse_log(iter)=mse0;
    
    for nc=1:Ncoeff
        caux=c; caux(nc)=caux(nc)+Dtap;
        y=filter(caux,1,x);
        mse1=mean(abs(y-y_ref).^2);
        gradient(nc)=(mse1-mse0)/Dtap;
    end
    c=c-0.001*gradient;
    
    gradient_log(iter,:)=gradient;
    h_log(iter,:)=c;
    
end

figure
stem(h_ref);
hold all
stem(c);

figure
plot(h_ref-c,'-o')

figure
plot(h_log)

figure
plot(10*log10(mse_log))

figure
plot(gradient_log)

