function [] = plotPSD(nsigs)

for count = 1:nsigs
    exampleSig(count,:) = buildWaveform(false);
    sigPow(count,:) = abs(fftshift(fft(exampleSig(count,:)))).^2;
end
sigPSD = mean(sigPow,1);
figure
hold on;
plot(sigPSD)
hold off;
title('spectrum')