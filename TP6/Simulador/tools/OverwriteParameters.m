%------------------------------------------------------------------------%
%				               FUNDACION FULGOR
%------------------------------------------------------------------------%

function s_output = OverwriteParameters (overwrite_s, default_config_s)
    
    fn = fieldnames(overwrite_s);
    for k = 1:numel(fn)
        if isfield(default_config_s,(fn{k}))
            if (~isstruct(overwrite_s.(fn{k})))
                default_config_s.(fn{k}) = overwrite_s.(fn{k});
            else
                default_config_s.(fn{k}) = ...
                    OverwriteParameters(overwrite_s.(fn{k}), default_config_s.(fn{k}));
            end
        else
            warning("%s: Invalid parameter.", fn{k});
        end
    end
    s_output = default_config_s;
end