%TODO: Add frequency axis
function [] = plotPSD(nsigs,sampPerSym,symRate,Nc,Nsym,guard)

fs_desired = 1;
ts = 1/sampPerSym/symRate;
Nmin = ceil(1/(fs_desired*ts)); %minimum length DFT for desired frequency granularity
%for efficient computation, choose FFT size to be power of 2
Nfft = 2^(nextpow2(Nmin)); %FFT size = the next power of 2 at least as big as Nmin
for count = 1:nsigs
    exampleSig(count,:) = buildWaveformRand(sampPerSym,symRate,Nc,Nsym,guard);
    sigPow(count,:) = abs(ts*fft(exampleSig(count,:),Nfft)).^2;
end
sigPSD = mean(sigPow,1);
%Alternatively, one could also use DFT size equal to the minimum length
%Nfft=Nmin;
%note: fft function in Matlab is just the DFT when Nfft is not a power of 2
%freq domain signal computed using DFT
%fft function of size Nfft automatically zeropads as needed
%fftshift function shifts DC to center of spectrum
signal_freqdomain_centered = fftshift(signal_freqdomain);
fs=1/(Nfft*ts); %actual frequency resolution attained
%set of frequencies for which Fourier transform has been computed using DFT
freqs = ((1:Nfft)-1-Nfft/2)*fs;


figure
hold on;
plot(20*log10(sigPSD))
hold off;
title('spectrum')