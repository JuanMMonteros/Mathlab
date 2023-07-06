function [y] = get_limited_band_signal(L,max_bw)
% Returns a reals signal that is the combination of several tones with random amplitudes and phases
% that does not exceed the max_bw. max_bw is in the range [0, pi]

Nfreqs=5000;
nline=(0:L-1)';
%freqs=linspace(0,max_bw,Nfreqs);
freqs=max_bw*rand(1,Nfreqs);
amps = rand(1,Nfreqs)+0.5;
phases=2*pi*rand(1,Nfreqs);
y=sum(amps.*cos(freqs.*nline+phases),2);

end

