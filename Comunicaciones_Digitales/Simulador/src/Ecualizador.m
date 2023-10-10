
function [o_data_s] = Ecualizador(i_config_s,rx_norm,RCMA)
%--------------------------%
%     DEFAULT SETTINGS
%--------------------------%
%% ECUALIZADOR FSE 
% Vamos a empezar solo con DD (slicer) sin FCR (carrier recovery)
config.NTAPS = 63; % Por comodidad impar
config.step_cma = 2^-11;%2e-3;
config.step_dd = 2^(-9);
config.tap_leak_gain = 1e-4;%1e-3;
config.force_cma_enable = 0;
config.N=2;
config.M=4;

 %--------------------------%
 %       REASSIGNMENT
 %--------------------------%

    fn = fieldnames(i_config_s);
    for k = 1:numel(fn)
        if isfield(config,(fn{k}))==1
            config.(fn{k})= i_config_s.(fn{k});
        else
            %error("%s: Parametro del simulador no valido", fn{k})
        end
    end


    %--------------------------%
    %         VARIABLES
    %--------------------------%
NTAPS=config.NTAPS;
step_cma=config.step_cma; 
step_dd=config.step_dd;
tap_leak_gain=config.tap_leak_gain; 
force_cma_enable=config.force_cma_enable; 
N=config.N;
M=config.M;
    
    %--------------------------%
    %         PROCESS
    %--------------------------%
   

Lrx = length(rx_norm);
buffer = zeros(NTAPS,1); %Columna
i_equalizer = rx_norm; % Columna
o_filter = zeros(Lrx,1); % salida del FSE (tasa 2)
yk = zeros(fix(Lrx/N),1); % salida del FSE (tasa 1)
ak = zeros(fix(Lrx/N),1); % salida del slicer (tasa 1)
W = zeros(NTAPS,1); W((NTAPS+1)/2)=8.0;


error = zeros(size(ak));
coeffs = zeros( [length(ak), length(W)]);

% Arranco con CMA, luego de 1/3 de la simulacion, paso a DD
% El ultimo tercio de la simulacion la uso para procesar BER

for m=1:Lrx
    
    if (m < fix(Lrx/3))
        enable_cma=1;
    else
        enable_cma=0;
    end
    
    % Actualizo buffer de entrada
    buffer(2:end)=buffer(1:end-1);
    buffer(1) = i_equalizer(m);
    
    % Filtro
    yfilt = W.'*buffer; %OJO: ' transpone y conjuga! usar .'
   
  
    if mod(m,N)==0 %Decimo por N
        m2 = ceil(m/N); % Con m cuento muestras y con m2 cuento simbolos
        % Este slicer es lentoooo, cambiarlo por uno mas eficiente
        yk(m2) = yfilt;
        o_filter(m) = yfilt;
       slicer_out = slicer(yfilt,M);
       ak(m2)=slicer_out;
        
        % Calculo del error
        if enable_cma || force_cma_enable
            Ek = yfilt*(abs(yfilt)-RCMA);
            step=step_cma;
        else
            Ek = yfilt-slicer_out;
            step=step_dd;
        end
        error(m2) = Ek;
        coeffs(m2,:) = W(:);
        
        % Emulo el sobre-muestreo del error, ejecutando el LMS
        % dentro de este if
        W = W*(1-step*tap_leak_gain) - step.*Ek.*conj(buffer);
    end
end
o_data_s.yk=yk;
