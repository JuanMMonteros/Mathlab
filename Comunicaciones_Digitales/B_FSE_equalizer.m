clear all
close all

%Basic TX QPSK 
BR = 32e9;  %Bd
L = 3e4;    % Simulati on Length
N = 2;      % Oversampling rate
fs = N*BR;  % Sampling rate to emulate analog domain
T = 1/BR;   % Time interval between two consecutive symbols
Ts = 1/fs;  % Time between 2 conseutive samples at Tx output
EbNo = 15; % dB
k = 2;      % k=log2(M)

% Two symbols generation (+1,-1)
xi = 2*randi([0,1],L,1)-1;
xq = 2*randi([0,1],L,1)-1;
x=xi+1j.*xq;

% Upsampling to change sampling rate
xup = upsample(x,N);

% Filter to interpolate the signal
g = rcosine(BR,fs,'sqrt', 0.2,25);
g = g/sum(g);
yup = filter(g,1,xup);

%% Agrego ISI y ruido

b = fir1(10,0.4); % Si el orden es impar hay que re-ajustar la fase de muestreo
%b = 1;
ych = filter(b,1,yup);

SNR_dB = EbNo - 10*log10(N) + 10*log10(k);
SNR = 10^(SNR_dB/10);
Pt = var(yup);
No = Pt/SNR;
sigma = sqrt(No/2);
ruido = sigma.*randn(length(yup),1)+1j.*sigma.*randn(length(yup),1);

rx = ych+ruido;

%%

NTAPS=101;
step=2e-3;
leak = 1e-3;
Xbuffer = zeros(NTAPS,1);
W = zeros(1, NTAPS);
W((NTAPS+1)/2)=1.0;
LY = length(rx);
out_eq = zeros(LY,1);
out_eq_down = zeros(LY/2,1);

% Plot
figure(100);
NFFT=256;
f = N*(-BR/2:BR/NFFT:BR/2-BR/NFFT);
plot(f/1e9,20.*log10(abs(fftshift(fft(conv(g,b), NFFT)))), 'LineWidth',2)
hold on;

for m=1:LY-NTAPS-1
    
    Xbuffer(2:end)=Xbuffer(1:end-1); Xbuffer(1)=rx(m);
    yeq = W*Xbuffer;
    out_eq(m)=yeq;
    
    if mod(m,N)==0
        error=yeq-slicer_qam(yeq,4);
        grad = error.*Xbuffer';
        W = W*(1-step*leak) - step*grad;
        out_eq_down(m/N)=yeq;
    end
    
    % Plot 
    if ~mod(m, 1000)
        figure(100); clf;
        HCH = abs(fftshift(fft(conv(g,b), NFFT)));
        plot(f/1e9,20.*log10(HCH),'LineWidth',2)
        hold on; grid on;
        HEQ = abs(fftshift(fft(W, NFFT)));
        plot(f/1e9,20.*log10(HEQ/(max(HEQ))),'LineWidth',2)
        plot(f/1e9,20.*log10(abs(HCH.*HEQ/(max(HEQ)))),'LineWidth',1.5)
        plot([BR BR]/2/1E9,ylim,'-.m','LineWidth',2)
        plot(-[BR BR]/2/1E9,ylim,'-.m','LineWidth',2)
       
        legend({'Channel','Equalizer','Folded Esp','+/- BR/2'},'Location','s')
        xlabel('Frequency [GHz]')
        ylabel('Amplitude [dB]')
        
        pause(.1);
        
    end
    
end

%%
scatterplot(out_eq_down(end-L/2:end-L/20))

figure; 
plot(real(out_eq_down),'.')

figure
plot(real(W), 'LineWidth', 2);
hold all
plot(imag(W)); % El canal es real, la parte imaginaria es siempre 0
grid on
xlabel('Tap Number')
ylabel('Tap Amplitude')
title('Impulse Response of equalizer')
