function [signal] = buildWaveform(sampPerSym,symRate,Nc,Nsym,debug)
sampRate = sampPerSym*symRate; %samples per second
Tofdm = 1/symRate;
a = randi([0,1],Nsym,Nc)*2-1; %content of message
freq_guard = 1;

t = linspace(0,Nsym/symRate,(Nsym/symRate)*sampRate);
sig = zeros(Nsym,Nc,length(t));

for u = 0:Nsym-1
    for n = 0:Nc-1
        fn = (n+1)/(Tofdm/freq_guard);
        sig(u+1,n+1,:) = a(u+1,n+1)*exp(1i*2*pi*fn*t);
        sig(u+1,n+1,((t-u*Tofdm)/Tofdm)<0 | ((t-u*Tofdm)/Tofdm)>1) = 0;
    end
end

temp = sum(sig,1);
signal.IQ = squeeze(sum(temp,2));
signal.sampRate = sampRate;
signal.symRate = symRate;
signal.numCarrier = Nc;

if debug
    close all;
    figure
    hold on;
    plot(t,real(signal))
    plot(t,imag(signal))
    hold off;
    title('IQ signal')

    sigFreq = fftshift(fft(signal));
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