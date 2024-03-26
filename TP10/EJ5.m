
clear; close all; clc;

config_s.fs_ch =400e6;
config_s.fs_dsp=100e6;

config_s.n_fires=10;

config_s.en_noise=0;
config_s.en_plots=0;

% TX
config_s.chirp_bw=100e6;
config_s.chirp_T=12e-6;
config_s.chirp_P=100;

config_s.f0=23e9;
config_s.pw_tx_dbm=13;

% Target and channel

config_s.range=200;
config_s.speed=70;

config_s.range_max=300;
config_s.speed_max=100;

% RX

config_s.snr_db=10;
config_s.fft_zp=16;

config_s.n_thr=20;

%Segundo target 
 config_s.deltaR=10;
 config_s.deltaV =10;

deltaR_v = [3 1 0.5];
n_deltaR = length(deltaR_v);

deltaV_v = [3 1 0.5];


for idx = 1:n_deltaR
  

    config_s.deltaR = deltaR_v(idx);
    config_s.deltaV = deltaV_v(idx);
    
    o_data = simulator_EJ5(config_s);
    figure;
	subplot(2,1,1);
	plot(o_data.X, mean(abs(o_data.fft_zp_m),1),'Linewidth',2);
	grid on;
	title('FFT')
	xlabel('Beat frequency')
	ylabel('Amplitude')
	xlim([1.8e6 ,1.95e6]);

        
	subplot(2,1,2);
	plot(o_data.Y, mean(abs(o_data.fft_zp_m),2),'Linewidth',2);
	grid on;
	title('FFT')
	xlabel('Doppler frequency')
	ylabel('Amplitude ')
     xlim([1500 ,2400]);

end
