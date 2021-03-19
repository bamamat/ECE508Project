function [finalSig] = buildWaveform(debug)

sampPerSym = 5; %samples per symbol
symRate = 50; %symbols per second
sampRate = sampPerSym*symRate; %samples per second
Nsym = 10; %number of symbols in a message
Nc = 5; %number of orthogonal carriers
Tofdm = 1/symRate;
a = randi([0,1],Nsym,Nc)*2-1; %content of message
freq_guard = 1;

t = linspace(0,Nsym/symRate,sampRate);
sig = zeros(Nsym,Nc,length(t));

for u = 0:Nsym-1
    for n = 0:Nc-1
        fn = (n+1)/(Tofdm/freq_guard);
        sig(u+1,n+1,:) = a(u+1,n+1)*exp(1i*2*pi*fn*t);
        idx = find(((t-u*Tofdm)/Tofdm)<0 | ((t-u*Tofdm)/Tofdm)>1);
        sig(u+1,n+1,idx) = 0;
    end
end

temp = sum(sig,1);
finalSig = squeeze(sum(temp,2));

if debug
    close all;
    figure
    hold on;
    plot(t,real(finalSig))
    plot(t,imag(finalSig))
    hold off;
    title('IQ signal')

    sigFreq = fftshift(fft(finalSig));
    figure
    hold on;
    plot(t,abs(sigFreq))
    hold off;
    title('spectrum')

    figure
    plot(t,squeeze(sum(sig,1)))
    title('signals')
end

end