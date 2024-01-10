%-----------------------------------------------------------------------------%
%                                   FULGOR
%
% Programmer(s): Santiago F. Leguizamon
% Created on   : July 2023
% Description  : Blind phase search algorithm
%-----------------------------------------------------------------------------%

function [o_data_s, o_cfg_s] = m_bps(i_data_s, i_cfg_s)

    %--------------------------%
    %     DEFAULT SETTINGS
    %--------------------------%

    cfg_s.n_phases          = 16    ;
    cfg_s.filter_length     = 128   ;
    cfg_s.mapper_levels     = 4     ;

    %--------------------------%
    %       REASSIGNMENT
    %--------------------------%

    if nargin > 0
        cfg_s = overwrite_parameters(i_cfg_s, cfg_s);
    end

    %--------------------------%
    %           ERRORS
    %--------------------------%
    
    %--------------------------%
    %     LOCAL VARIABLES
    %--------------------------%

    n_phases = cfg_s.n_phases;
    filter_length = cfg_s.filter_length;
    mapper_levels = cfg_s.mapper_levels;
    group_delay = floor((filter_length - 1) / 2);

    tested_phases = ((0 : n_phases - 1) - n_phases / 2) * (pi / 2) / n_phases;
    bps_filter = ones(filter_length, 1) ./ filter_length;
    
    %--------------------------%
    %          PROCESS
    %--------------------------%

    % Rotation by test phases
    rot_data = i_data_s.signal_v.' * exp(1j * tested_phases);
    
    % Slicer
    decided_v = my_qam_slicer(rot_data, mapper_levels);
    
    % Error calculation
    squared_error_v = abs( (decided_v - rot_data) .^ 2);
    
    % Filtering error
    filtered_error_vector = filter(bps_filter, 1, squared_error_v);

    % Calculating the phase idx that results in the smallest error
    [~, min_indexes] = min(filtered_error_vector.'); % Traspuesta para sacar minimo por cada fila
    
    % Selecting the optimal phase
    min_phases = tested_phases(min_indexes).';
    bps_phase = mod(unwrap(min_phases * 4) / 4, 2 * pi);

    % Alingning data with the phase values
    delayed_data = [zeros(group_delay,1); i_data_s.signal_v(1:end-group_delay).'];
    
    % Final correction
    bps_output = delayed_data.* exp(1j.*bps_phase);

    %--------------------------%
    %          OUTPUT
    %--------------------------%

    o_data_s.bps_output = bps_output;
    o_data_s.bps_phase = bps_phase;
    o_data_s.bps_delayed_input = delayed_data;

    o_cfg_s = cfg_s;

    %--------------------------%
    %           PLOTS
    %--------------------------%

end
