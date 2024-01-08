clear
close all

tx_config.M = 64; % Cantidad de niveles de la modulacion
tx_config.NOS = 2; % Tasa de sobremuestreo
tx_config.Lsymbs = 100e3; % Cantidad de simbolos
tx_config.rolloff = 0.5; % Rolloff del filtro conformador
tx_config.pulse_shaping_ntaps = 200;
tx_config.pulse_shaping_type = 1; % 0: RRC, 1: RC

odata = transmisor_MQAM(tx_config);

tx_symbs = odata.tx_symbs;
tx_data = odata.oversampled_output;

scatterplot(tx_symbs)
eyediagram(tx_data(1:10e3), tx_config.NOS )

figure
pwelch(tx_data)
