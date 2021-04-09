function [signal] = buildWaveform(bits,sampPerSym,symRate,Nc,Nsym,guard)
debug = 0;
sampRate = sampPerSym*symRate; %samples per second
Tofdm = 1/symRate;
%a = randi([0,1],Nsym,Nc)*2-1; %content of message
if length(bits) ~= Nsym*Nc
    bits = [bits zeros(1,Nsym*Nc-length(bits))];
end
a = reshape(bits,Nsym,Nc);
freq_guard = guard;

t = linspace(0,Nsym/symRate,(Nsym/symRate)*sampRate);

%{
sig = zeros(Nsym,Nc,length(t));

%V1
for u = 0:Nsym-1
    for n = 0:Nc-1
        fn = (n+1)/(Tofdm/freq_guard);
        sig(u+1,n+1,:) = a(u+1,n+1)*exp(1i*2*pi*fn*t);
        sig(u+1,n+1,((t-u*Tofdm)/Tofdm)<0 | ((t-u*Tofdm)/Tofdm)>1) = 0;
    end
end

%V2
for u = 0:Nsym-1
    for n = 0:Nc-1
        fn = (n+1)/(Tofdm/freq_guard);
        tempSig = a(u+1,n+1)*exp(1i*2*pi*fn*t);
        tempSig(((t-u*Tofdm)/Tofdm)<0 | ((t-u*Tofdm)/Tofdm)>1) = 0;
        sig = sig + tempSig;
    end
end

temp = sum(sig,1);
signal = squeeze(sum(temp,2));
%}

signal = zeros(1,length(t));

for n = 0:Nc-1
    fn = (n+1)/(Tofdm/freq_guard);
    tempSig = repelem(a(:,n+1),sampPerSym,1).'.*exp(1i*2*pi*fn*t);
    signal = signal + tempSig;
end


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
end

end