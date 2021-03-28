function [] = plotPSD(nsigs,sampPerSym,symRate,Nc,Nsym)

for count = 1:nsigs
    exampleSig(count,:) = buildWaveform(sampPerSym,symRate,Nc,Nsym,0);
    sigPow(count,:) = abs(fftshift(fft(exampleSig(count).IQ))).^2;
end
sigPSD = mean(sigPow,1);
figure
hold on;
plot(sigPSD)
hold off;
title('spectrum')

sampPerSym = 20; %samples per symbol
symRate = 500; %symbols per second
Nsym = 20; %number of symbols in a message
Nc = 5; %number of orthogonal carriers

signal = buildWaveform(sampPerSym,symRate,Nc,Nsym,0);