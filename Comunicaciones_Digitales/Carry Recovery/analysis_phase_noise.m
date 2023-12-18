clear all
% close all

fs=16e9;
LW=50e3;


Lsim=1000e3;
Niters=10;


for niter=1:Niters
    fnoise=sqrt(LW*2*pi/fs).*randn(Lsim,1);
    pnoise=cumsum(fnoise);
     plot(pnoise)
     hold all
    
    if niter==1
        Pxx_mean = 1/Niters*pwelch(exp(1j*pnoise), hanning(1024),0,4096,fs);
    else
        Pxx_mean = Pxx_mean + 1/Niters*pwelch(exp(1j*pnoise), hanning(1024),0,4096,fs);
    end
end
% grid on

plot(fftshift(Pxx_mean))
