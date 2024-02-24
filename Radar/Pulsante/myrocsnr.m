function [pd,pfa] = myrocsnr(snr_dB, varargin)

    p=inputParser;
    
    addRequired(p, 'snr_dB');
    addOptional(p, 'MinPfa',1e-10);
    addOptional(p, 'MaxPfa',1);
    addOptional(p, 'NumPoints',101);
    
    parse(p,snr_dB,varargin{:});
    
    MinPfa = log10(p.Results.MinPfa);
    MaxPfa = log10(p.Results.MaxPfa);
    NumPoints = p.Results.NumPoints;
    
    snr = 10.^(snr_dB./10);
    
    pfa = logspace(MinPfa,MaxPfa,NumPoints).';
    
    t0 = sqrt(2*snr);
    t1 = sqrt(-2.*log(pfa));
    pd = marcumq(t0,t1);
    
end

