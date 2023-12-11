%-----------------------------------------------------------------------------%
%                                   FULGOR
%
% Programmer(s): Francisco G. Rainero
%-----------------------------------------------------------------------------%

clear 
close all
clc;

%% BER vs EbNo
EbNo_db_v = 0:15;

ber_pam2_v = berawgn(EbNo_db_v,'PAM',2);
ber_pam4_v = berawgn(EbNo_db_v,'PAM',4);
ber_qpsk_v = berawgn(EbNo_db_v,'QAM',4);
ber_qam16_v = berawgn(EbNo_db_v,'QAM',16);
ber_qam64_v = berawgn(EbNo_db_v,'QAM',64);

% Plot
figure
semilogy(EbNo_db_v,ber_pam2_v, 'x-k');
hold on
semilogy(EbNo_db_v,ber_pam4_v, 'x-b');
semilogy(EbNo_db_v,ber_qpsk_v, '--r');
semilogy(EbNo_db_v,ber_qam16_v, '--m');
semilogy(EbNo_db_v,ber_qam64_v, '--c');
title('BER vs EbNo')
xlabel('EbNo [dB]')
ylabel('BER')
legend({'PAM2','PAM4','QPSK','QAM16','QAM64'},'Location','sw'); 
grid on;
set(gcf, 'Position', [50 50 500 500],'Color', 'w');

%% BER vs EbNo

SNR_db_v = 0:15;

% PAM2
M = 2;
SNR_v = 10.^(SNR_db_v/10);
EbNo_v = SNR_v / log2(M);
EbNo_db_v = 10*log10(EbNo_v);
ber_pam2_v = berawgn(EbNo_db_v,'PAM',M);

% PAM4
M = 4;
SNR_v = 10.^(SNR_db_v/10);
EbNo_v = SNR_v / log2(M);
EbNo_db_v = 10*log10(EbNo_v);
ber_pam4_v = berawgn(EbNo_db_v,'PAM',M);

% QPSK
M = 4;
SNR_v = 10.^(SNR_db_v/10);
EbNo_v = SNR_v / log2(M);
EbNo_db_v = 10*log10(EbNo_v);
ber_qpsk_v = berawgn(EbNo_db_v,'QAM',M);

% QAM16
M = 16;
SNR_v = 10.^(SNR_db_v/10);
EbNo_v = SNR_v / log2(M);
EbNo_db_v = 10*log10(EbNo_v);
ber_qam16_v = berawgn(EbNo_db_v,'QAM',M);

% QAM64
M = 64;
SNR_v = 10.^(SNR_db_v/10);
EbNo_v = SNR_v / log2(M);
EbNo_db_v = 10*log10(EbNo_v);
ber_qam64_v = berawgn(EbNo_db_v,'QAM',M);

% Plot
figure
semilogy(SNR_db_v,ber_pam2_v, 'x-k');
hold on
semilogy(SNR_db_v,ber_pam4_v, 'x-b');
semilogy(SNR_db_v,ber_qpsk_v, '--r');
semilogy(SNR_db_v,ber_qam16_v, '--m');
semilogy(SNR_db_v,ber_qam64_v, '--c');
title('BER vs SNR')
xlabel('SNR [dB]')
ylabel('BER')
legend({'PAM2','PAM4','QPSK','QAM16','QAM64'},'Location','sw'); 
grid on;
set(gcf, 'Position', [50 50 500 500],'Color', 'w');
