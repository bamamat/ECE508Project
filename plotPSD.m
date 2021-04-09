function [] = plotPSD(nsigs,sampPerSym,symRate,Nc,Nsym,guard)

for count = 1:nsigs
    exampleSig(count,:) = buildWaveformRand(sampPerSym,symRate,Nc,Nsym,guard);
    sigPow(count,:) = abs(fftshift(fft(exampleSig(count,:)))).^2;
end
sigPSD = mean(sigPow,1);
figure
hold on;
plot(sigPSD)
hold off;
title('spectrum')