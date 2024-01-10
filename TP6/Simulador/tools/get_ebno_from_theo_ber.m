% Get the OSNR to obtain the theoretical BER

function [ebno_db_v] = get_ebno_from_theo_ber(ber_v, M)

    if M>2
        mod = 'QAM';
    else
        mod = 'PAM';
    end
    
    ebno_db_v = zeros(size(ber_v));
    
    for idx = 1:length(ber_v)
    
        ber = ber_v(idx);
        done = 0;
        ebno_int_db_v = 0:0.1:15;

        % Iteration to obtain ber vector

        while done == 0

           ber_int_v = berawgn(ebno_int_db_v, mod, M);

           if (min(ber_int_v) < ber) && (max(ber_int_v) > ber)

               done = 1;

           elseif min(ber_int_v) > ber

               ebno_int_db_v = ebno_int_db_v + 0.1;

           elseif max(ber_int_v) < ber    

               ebno_int_db_v = ebno_int_db_v - 0.1;

           end

        end

        % Interpolation
        ebno_db = interp1(log10(ber_int_v), ebno_int_db_v, log10(ber));
        ebno_db_v(idx) = ebno_db;
        
    end
end
