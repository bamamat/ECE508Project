clc;
close all;
clear all;

%plotting
depVar = 'maxRangeDesired';
showPlots = 0;

%constants
c = 299.8E6;%speed of light
boltz = 1.3806485E-23;%boltzman const

%simple Radar parameters
fc=1E9;%Hz
lambda = c/fc;%
txpower=40;%dBm
txgain=0;%dBi
rxgain=0;%dBi
receiverTemp=300;
minRangeDesired = 10e3; %this is an upper limit to the length of the pulse
maxRangeDesired = 50e3;

%antenna postitions
nTxAnt = 256;
nRxAnt = 16;
txSpacing = [0 5 0];
rxSpacing = [0 5 0];
txAntCenter = [10000,0,10];
rxAntCenter = [0,0,0];

%comms parameters
guard = 2;%must be >= 1 (multiplier to the BW separation minimally required)
sampPerSym = 20; %samples per symbol
symRate = 1e6; %symbols per second
sampRate = sampPerSym*symRate;
ts = 1./sampRate;
fs_desired = 1;
Nmin = ceil(1./(fs_desired*ts)); %minimum length DFT for desired frequency granularity
Nfft = 2.^(nextpow2(Nmin)); %FFT size = the next power of 2 at least as big as Nmin
fs=1./(Nfft.*ts);
maxNc = sampPerSym./guard; %consider this an upper bandwidth limit
bandWidth = ((Nfft-1)-1-Nfft./2).*fs;

txAntPos = repmat(txAntCenter,nTxAnt,1).' + txSpacing.'*((-nTxAnt+1)/2:(nTxAnt-1)/2);
rxAntPos = repmat(rxAntCenter,nRxAnt,1).' + rxSpacing.'*((-nRxAnt+1)/2:(nRxAnt-1)/2);

%radar targets
targets(1).rcs = 0;
targets(1).initialPos = [0,0,10];
targets(1).velocity = [-5,0,0];
targets(2).rcs = 5;
targets(2).initialPos = [-20000,0,10];
targets(2).velocity = [10,0,0];
targets(3).rcs = 10;
targets(3).initialPos = [-39000,0,10];
targets(3).velocity = [0,0,0];

%ENCODE DATA
%read an image
RGB = imread('peppers.png');
%add any metadata you want up to 251 unsigned int16s
int16MetaData = [];
binaryIntVect = encodeIntImg(RGB, int16MetaData);

antSelectBits = num2str(binaryIntVect(1:log2(nTxAnt)));
antSelectBits(isspace(antSelectBits)) = '';
antSelect = bin2dec(antSelectBits);
bpskBits = binaryIntVect(log2(nTxAnt)+1:end) * 2 - 1;

%OFDM params
Nsym = floor(minRangeDesired*2./c*symRate); %number of symbols in a message
Nc = min(ceil(length(bpskBits)./Nsym), maxNc);

%more radar params
nPulse = ceil(length(bpskBits)./Nsym./Nc);
pulseWidth=Nsym./symRate;
PRF = floor(1./((2*maxRangeDesired/c)+pulseWidth+1./sampRate));
PRI = 1./PRF;
dataRate = Nc.*Nsym.*PRF;

%basic radar calcs
minUnambigRange=pulseWidth*c./2;
maxUnambigRange=c*(PRI-pulseWidth)./2;

minRange = ceil(minUnambigRange./1000)*1000;
maxRange = floor(maxUnambigRange./1000)*1000;
minTime = minRange*2/c;

disp(['actual min range: ' num2str(minRange(1))])
disp(['actual max range: ' num2str(maxRange(1))])

if showPlots
    if strcmp(depVar,'symRate')
        figure
        plot(symRate/1e3,dataRate/80e6)
        xlabel('Symbol Rate (ksym/s)')
        ylabel('Data Rate (MB/s)')
        title('Data Rate vs Symbol Rate')
        figure
        plot(symRate/1e3,minRange/1e3)
        xlabel('Symbol Rate (ksym/s)')
        ylabel('Min Range (km)')
        title('Min Range vs Symbol Rate')
        figure
        plot(symRate/1e3,maxRange/1e3)
        xlabel('Symbol Rate (ksym/s)')
        ylabel('Max Range (km)')
        title('Max Range vs Symbol Rate')
        figure
        plot(symRate/1e3,nPulse/1e3)
        xlabel('Symbol Rate (ksym/s)')
        ylabel('Number of Pulses')
        title('Number of Pulses vs Symbol Rate')
    end
    if strcmp(depVar,'guard')
        figure
        plot((guard-1)*symRate/1e3,dataRate/8e6)
        xlabel('Freq Guard BW (kHz)')
        ylabel('Data Rate (MB/s)')
        title('Data Rate vs Symbol Rate')
        figure
        plot((guard-1)*symRate/1e3,minRange/1e3)
        xlabel('Freq Guard BW (kHz)')
        ylabel('Min Range (km)')
        title('Min Range vs Symbol Rate')
        figure
        plot((guard-1)*symRate/1e3,maxRange/1e3)
        xlabel('Freq Guard BW (kHz)')
        ylabel('Max Range (km)')
        title('Max Range vs Symbol Rate')
        figure
        plot((guard-1)*symRate/1e3,nPulse/1e3)
        xlabel('Freq Guard BW (kHz)')
        ylabel('Number of Pulses')
        title('Number of Pulses vs Symbol Rate')
    end
    if strcmp(depVar,'minRangeDesired')
        figure
        plot(minRangeDesired/1e3,dataRate/80e6)
        xlabel('Min Range Desired (km)')
        ylabel('Data Rate (MB/s)')
        title('Data Rate vs Min Range Desired')
        figure
        plot(minRangeDesired/1e3,minRange/1e3)
        xlabel('Min Range Desired (km)')
        ylabel('Min Range (km)')
        title('Min Range vs Min Range Desired')
        figure
        plot(minRangeDesired/1e3,maxRange/1e3)
        xlabel('Min Range Desired (km)')
        ylabel('Max Range (km)')
        title('Max Range vs Min Range Desired')
        figure
        plot(minRangeDesired/1e3,nPulse/1e3)
        xlabel('Min Range Desired (km)')
        ylabel('Number of Pulses')
        title('Number of Pulses vs Min Range Desired')
    end
    if strcmp(depVar,'maxRangeDesired')
        figure
        plot(maxRangeDesired/1e3,dataRate/80e6)
        xlabel('Max Range Desired (km)')
        ylabel('Data Rate (MB/s)')
        title('Data Rate vs Max Range Desired')
        figure
        plot(maxRangeDesired/1e3,minRange/1e3)
        xlabel('Max Range Desired (km)')
        ylabel('Min Range (km)')
        title('Min Range vs Max Range Desired')
        figure
        plot(maxRangeDesired/1e3,maxRange/1e3)
        xlabel('Max Range Desired (km)')
        ylabel('Max Range (km)')
        title('Max Range vs Max Range Desired')
        figure
        plot(maxRangeDesired/1e3,nPulse/1e3)
        xlabel('Sample Rate (MS/s)')
        ylabel('Number of Pulses')
        title('Number of Pulses vs Max Range Desired')
    end
else
    
    for pulseCount = 1:128%nPulse

        bits = bpskBits(Nc*Nsym*(pulseCount-1)+1:min(Nc*Nsym*pulseCount,length(bpskBits)));
        %basic signal calcs
        %signal = buildWaveform(bits,sampPerSym,symRate,Nc,Nsym,guard);
        signal = exp(1i*pi*(1e6./pulseWidth)*(1/sampRate:1/sampRate:pulseWidth).^2);
        pulses(:,pulseCount) = signal;
        radarDataLength = floor(((maxRange - minRange)*2/c + pulseWidth)*sampRate);

        %noise calculations
        noisePow = boltz*bandWidth*receiverTemp;
        radarData(:,pulseCount) = noisePow*(randn(1,radarDataLength)+1i*randn(1,radarDataLength));
        radarPos = txAntPos(:,antSelect);

        for tarCount = 1:length(targets)
            tar = targets(tarCount);
            tarPos = tar.initialPos+tar.velocity*PRI*(pulseCount-1);
            pulseTime = PRI*(pulseCount-1);
            tarTime = (norm(tarPos-radarPos.')*2/c) - minTime;
            wholeSampDelay = floor(sampRate*tarTime);
            fracDealy = sampRate*tarTime-floor(sampRate*tarTime);
            tarReturn = ampScale(signal);
            tarReturn = freqShift(tarReturn, sampRate, pulseTime, lambda, [0 0 0], tar.velocity, radarPos, tarPos);
            tarReturn = timeDelay(tarReturn,fracDealy);
            radarData(wholeSampDelay+1:wholeSampDelay+length(tarReturn),pulseCount) = ...
                radarData(wholeSampDelay+1:wholeSampDelay+length(tarReturn),pulseCount) + tarReturn.';
        end

        rangeDoppler(:,pulseCount) = conv(radarData(:,pulseCount),conj(pulses(:,pulseCount)));
        
    end
    
    usedPulses = 128;
    rangeDoppler = rangeDoppler(:,1:usedPulses);
    
    figure
    surf(20*log10(abs(radarData(:,1:usedPulses))));
    
    figure
    surf(20*log10(abs(rangeDoppler(:,1:usedPulses))));
    
    for count = 1:size(rangeDoppler,1)
        matchedFilter(count,1:usedPulses) = fftshift(fft(rangeDoppler(count,1:usedPulses)));
    end

    %RECEIVED DATA
    binaryReceived = uint16((bpskBits + 1) / 2);
    txAntCalc = antSelect;

    RGBrecovered = decodeIntImg(binaryReceived, txAntCalc, nTxAnt);
    figure
    imshow(RGBrecovered)
    
    figure
    surf(20*log10(abs(matchedFilter(:,1:usedPulses))));
    

    %math stuff
    %plotPSD(1000,sampPerSym,symRate,Nc,Nsym,guard)    
    
end
    


function [outSig] = timeDelay(inSig, fracSamp)
    % Design fractional-delay filter
    t = (-(31-1)/2:(31-1)/2)-fracSamp;
    delayFilter = sin(pi*t)./(pi*t);
    outSig = conv([inSig 0],delayFilter,'same');
end

function [outSig] = freqShift(inSig, sampRate, pulseTime, lambda, velRad, velTar, posRad, posTar)
    speed = dot((posTar-posRad.'),(velTar-velRad))/norm((posTar-posRad.'));
    shift = -2*(speed/lambda);
    t = (0:1/sampRate:(length(inSig)-1)/sampRate) + pulseTime;
    outSig = inSig .* exp(1i*2*pi*shift*t);
end

function [outSig] = ampScale(inSig, fc, txPower, rcs, range)
    scale = 1e-13;
    outSig = inSig*scale;
end