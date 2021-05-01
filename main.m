clc;
close all;
clear all;

%plotting
depVar = 'maxRangeDesired';

%flags
showPlots = 0;
showPSD = 0;
useLFM = 0;
useRandAntenna = 0;
runRadar = 1;
runComms = 1;

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
minRangeDesired = 20e3; %this is an upper limit to the length of the pulse
maxRangeDesired = 45e3;

%antenna postitions
nTxAnt = 256;
nRxAnt = 16;
txSpacing = [0 5 0];
rxSpacing = [0 5 0];
txAntCenter = [20000,0,10];
rxAntCenter = [0,0,10];


%define multipath values
directA = 1;
directPhi = 0;
reflectA = 0.98;
reflectPhi = pi;

%comms parameters
guard = 1.25;%must be >= 1 (multiplier to the BW separation minimally required)
sampPerSym = 200; %samples per symbol
symRate = 1000e3; %symbols per second
%calc useful variables
sampRate = sampPerSym*symRate;
ts = 1./sampRate;

Nmin = ceil(1./(ts)); %minimum length DFT for desired frequency granularity
Nfft = 2.^(nextpow2(Nmin)); %FFT size = the next power of 2 at least as big as Nmin
fs=1./(Nfft.*ts);
maxNc = sampPerSym./guard; %consider this an upper bandwidth limit
bandWidth = ((Nfft-1)-1-Nfft./2).*fs;

%place the tx and rx antennas
txAntPos = repmat(txAntCenter,nTxAnt,1).' + txSpacing.'*((-nTxAnt+1)/2:(nTxAnt-1)/2);
rxAntPos = repmat(rxAntCenter,nRxAnt,1).' + rxSpacing.'*((-nRxAnt+1)/2:(nRxAnt-1)/2);

%radar targets
targets(1).rcs = 0;
targets(1).initialPos = [0,0,10];
targets(1).velocity = [-20,0,0];
targets(2).rcs = -3;
targets(2).initialPos = [-20000,0,10];
targets(2).velocity = [40,0,0];
targets(3).rcs = 3;
targets(3).initialPos = [-10000,0,10];
targets(3).velocity = [0,0,0];

%ENCODE DATA
%read an image
RGB = imread('peppers.png');
%add any metadata you want up to 251 unsigned int16s
int16MetaData = [];
binaryIntVect = encodeIntImg(RGB, int16MetaData);
%choose an antenna based on bit data
antSelectBits = num2str(binaryIntVect(1:log2(nTxAnt)));
antSelectBits(isspace(antSelectBits)) = '';
antSelect = bin2dec(antSelectBits);
%remove the bits used in spatial modulation from the message
%also map the binary to symbols [+1/-1]
bpskBits = binaryIntVect(log2(nTxAnt)+1:end) * 2 - 1;

%final OFDM params
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
%round to the nearest Km
minRange = ceil(minUnambigRange./1000)*1000;
maxRange = floor(maxUnambigRange./1000)*1000;

disp(['actual min range: ' num2str(minRange(1))])
disp(['actual max range: ' num2str(maxRange(1))])

if ~showPlots
    
    %calc some useful data
    radarDataLength = ceil(((maxRange - minRange)*2/c + pulseWidth)*sampRate);
    minTime = minRange*2/c;
    noisyBits = [];
    scale = 5e-12;
    pulses = zeros(Nsym*sampPerSym,nPulse);
    %noise calculations
    noisePow = boltz*bandWidth*receiverTemp;
    radarData = noisePow*(randn(radarDataLength,nPulse)+1i*randn(radarDataLength,nPulse));
    
    for pulseCount = 1:nPulse
        
        %grab the bits for the current pulse
        bits = bpskBits(Nc*Nsym*(pulseCount-1)+1:min(Nc*Nsym*pulseCount,length(bpskBits)));
        %basic signal calcs
        if useLFM
            signal = exp(1i*pi*(10e6./pulseWidth)*(1/sampRate:1/sampRate:pulseWidth).^2);
        else
            signal = buildWaveform(bits,sampPerSym,symRate,Nc,Nsym,guard);
        end
        %log the signal for the matched filter
        pulses(:,pulseCount) = signal;
        
        if runRadar
            
            %select the radar position as the antenna position
            if useRandAntenna
                radarPos = txAntPos(:,randi(nTxAnt));
            else
                radarPos = txAntPos(:,antSelect);
            end
            %for each target
            for tarCount = 1:length(targets)
                tar = targets(tarCount);
                %calcualte the position based on the velocity and initial pos
                tarPos = tar.initialPos+tar.velocity*PRI*(pulseCount-1);
                %calculate the pulse time
                pulseTime = PRI*(pulseCount-1);
                %find the time after the radar starts receiving of the return
                tarTime = (norm(tarPos-radarPos.')*2/c) - minTime;
                %take the whole sample delay as an index
                wholeSampDelay = floor(sampRate*tarTime);
                %take the fractional delay for phase delay
                fracDelay = sampRate*tarTime-floor(sampRate*tarTime);
                %send the signal through the channel
                %amplitude scale
                tarReturn = ampScale(signal, scale, tar.rcs);
                %frequency shift
                tarReturn = freqShift(tarReturn, sampRate, pulseTime, ...
                    lambda, [0 0 0], tar.velocity, radarPos, tarPos);
                %phase shift
                tarReturn = timeDelay(tarReturn, fracDelay);
                %multipath
                [response, ~] = tworayResp(radarPos(3),tarPos(3),norm(radarPos(1:2).'-tarPos(1:2)),directPhi,...
                                            reflectPhi,directA,reflectA,fc);
                tarReturn = tarReturn.*response.^2;
                %add the return back to the received data matrix
                radarData(wholeSampDelay+1:wholeSampDelay+length(tarReturn),pulseCount) = ...
                    radarData(wholeSampDelay+1:wholeSampDelay+length(tarReturn),pulseCount) + tarReturn.';
            end
        end
        
        if runComms
            % COMMS SIDE
            noisePowComms = 0.005;
            commsSig = signal + sqrt(noisePowComms).*randn(size(signal));
            noisyBits(end+1:end+Nc*Nsym) = decodeWaveform(commsSig,sampPerSym,symRate,Nc,Nsym,guard);
            %sum(noisyBits(:,pulseCount).' ~= bits)
        end
        
    end
    
    time = linspace(0,pulseWidth, length(pulses));
    figure
    hold on
    plot(time,real(pulses(:,1)))
    plot(time,imag(pulses(:,1)))
    hold off
    xlabel('time (s)')
    ylabel('amplitude')
    title('Time Domain OFDM')
    
    if runRadar
        usedPulses = 64;%2.^floor(log2(nPulse));
        radarData = radarData(:,1:usedPulses);
        pulses = pulses(:,1:usedPulses);

        windowMatch = hamming(size(pulses,1));

        validLen = ceil((maxRange - minRange)*2/c*sampRate);
        rangeDoppler = ifft(fft(radarData,[],1).*conj(fft(pulses.*repmat(windowMatch,1,size(pulses,2)),size(radarData,1),1)));

        windowDoppler = hamming(usedPulses);
        rangeDoppler = rangeDoppler(1:validLen,1:usedPulses);

        matchedFilter = fftshift(fft(windowDoppler.'.*rangeDoppler(:,1:usedPulses),[],2),2);
        freqs = linspace(-PRF/4,PRF/4, usedPulses);
        ranges = linspace(minRange, maxRange, validLen);
        
        figure
        imagesc(freqs*lambda/2,ranges/1e3,300+20*log10(abs(matchedFilter(:,1:usedPulses))));
        xlabel('Speed Towards Radar (m/s)')
        ylabel('Target Range (km)')
        title('Range Doppler Amplitudes')

        figure
        surf(freqs*lambda/2,ranges/1e3,300+20*log10(abs(matchedFilter(:,1:usedPulses))));
        xlabel('Speed Towards Radar (m/s)')
        ylabel('Target Range (km)')
        zlabel('Amplitude (dB)')
        title('Range Doppler')
    end
    
    %RECEIVED DATA
    if runComms
        noisyBits = noisyBits(1:length(bpskBits));
        errors = sum(noisyBits ~= bpskBits);
        binaryReceived = uint16((noisyBits + 1) / 2);
        disp(['error rate: ', num2str(errors/length(bpskBits))])
        disp(['SNR: ', num2str(10*log10(mean(abs(pulses(:,1).^2))/noisePowComms))])
        txAntCalc = antSelect;

        RGBrecovered = decodeIntImg(binaryReceived, txAntCalc, nTxAnt,size(RGB));
        figure
        imshow(RGBrecovered)
    end
    
    %Show OFDM Spectrum
    if showPSD
        plotPSD(1000,sampPerSym,symRate,Nc,Nsym,guard)
    end
    
else
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
end
    
%calculates the distances of the direct path and the reflected path
%inputs: ant1 height, ant2 height, dist between antennas along ground
%outputs: direct path length, reflected path length, angle of reflection
function [Direct,Reflect,theta] = tworayDist(Ant1,Ant2,totalG)
    %calculate distance on ground of two components
    LeftG = totalG*Ant1./(Ant1+Ant2);
    RightG = totalG-LeftG;
    %calculate theta for error checking
    theta = atan(Ant1./LeftG);
    %use pythagorean theorem to calc direct path length
    Direct = sqrt((Ant1-Ant2).^2 + totalG.^2);
    %use pythagorean theorem to calc both parts of the reflection length
    Reflect = sqrt(Ant1.^2 + LeftG.^2) + sqrt(Ant2.^2 + RightG.^2);
end

%calculates the response of the direct path and the reflected path
%inputs: ant1 height, ant2 height, dist between antennas along ground
%direct path phi, reflected path phi, direct A, reflected A, center freq
%outputs: total path response, nominal path response
function [response,responseNom]=tworayResp(h1,h2,distG,dPhi,rPhi,dA,rA,fc)
    %speed of light
    c = 299792458;
    %calculate the distance of the paths
    [directDist,reflectDist,~] = tworayDist(h1, h2, distG);
    %calculate time of the signal paths
    directTau = directDist/c;
    reflectTau = reflectDist/c;
    %calculate the channel responses
    directH = (dA).*exp(-1i*(2*pi*fc*directTau+dPhi));
    reflectH = (rA).*exp(-1i*(2*pi*fc*reflectTau+rPhi));
    %convert to total response and nominal response
    response = directH+reflectH;
    responseNom = directH;
end


function [outSig] = timeDelay(inSig, fracSamp)
    % Design fractional-delay filter
    t = (-(31-1)/2:(31-1)/2)-fracSamp;
    delayFilter = sin(pi*t)./(pi*t);
    outSig = conv([inSig 0],delayFilter,'same');
end

function [outSig] = freqShift(inSig, sampRate, pulseTime, lambda, velRad, velTar, posRad, posTar)
    %find radial speed
    speed = dot((posTar-posRad.'),(velTar-velRad))/norm((posTar-posRad.'));
    %calculate the shift
    shift = -2*(speed/lambda);
    %generate the time vector for the doppler shift
    t = (0:1/sampRate:(length(inSig)-1)/sampRate) + pulseTime;
    %apply the doppler shift
    outSig = inSig .* exp(1i*2*pi*shift*t);
end

function [outSig] = ampScale(inSig, scale, rcs)
    outSig = inSig*scale*10^(rcs/10);
end