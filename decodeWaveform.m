function [bits] = decodeWaveform(signal,sampPerSym,symRate,Nc,Nsym,guard)

sampRate = sampPerSym*symRate; %samples per second
Tofdm = 1/symRate;
freq_guard = guard;

t = linspace(0,Nsym/symRate,Nsym*(sampRate/symRate));


for n = 0:Nc-1
    
    fn = (n+1)/(Tofdm/freq_guard);
    %tempBits = bandpass(signal,[fn-0.5/(Tofdm/freq_guard), fn+0.5/(Tofdm/freq_guard)],sampRate);
    tempBits = signal.*exp(-1i*2*pi*fn*t);
    preMean = reshape(tempBits,[],Nsym).';
    bits(:,n+1) = sign(real(mean(preMean,2)));
end

bits = reshape(bits,[],1);

%figure
%plot(bits(1:Nsym))

end