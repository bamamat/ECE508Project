clc;
close all;
clear all;

%constants
c = 299.8E6;%speed of light
boltz = 1.3806485E-23;%boltzman const

%simple Radar parameters
PRF=2000;
PRI=1/PRF;%seconds
fc=1E9;%Hz
txpower=40;%dBm
txgain=0;%dBi
rxgain=0;%dBi
receiverTemp=300;
receiverBW=1E6;

%antenna postitions
nTxAnt = 256;
nRxAnt = 16;
txSpacing = [0 5 0];
rxSpacing = [0 5 0];
txAntCenter = [20000,0,10];
rxAntCenter = [0,0,0];

%comms parameters
guard = 1.25;
sampPerSym = 200; %samples per symbol
symRate = 250000; %symbols per second
sampRate = sampPerSym*symRate;

minRangeDesired = 20000; %this is an upper limit to the length of the pulse
maxNc = sampPerSym/guard; %consider this an upper bandwidth limit

txAntPos = repmat(txAntCenter,nTxAnt,1).' + txSpacing.'*((-nTxAnt+1)/2:(nTxAnt-1)/2);
rxAntPos = repmat(rxAntCenter,nRxAnt,1).' + rxSpacing.'*((-nRxAnt+1)/2:(nRxAnt-1)/2);

%radar targets
targets(1).rcs = 0;
targets(1).initialPos = [0,0,0];
targets(1).velocity = [0,0,0];
targets(2).rcs = 5;
targets(2).initialPos = [1000,0,0];
targets(2).velocity = [0,0,0];
targets(3).rcs = 10;
targets(3).initialPos = [-1000,0,0];
targets(3).velocity = [0,0,0];

%ENCODE DATA
RGB = imread('peppers.png');
int16MetaData = [];
binaryIntVect = encodeIntImg(RGB, int16MetaData);

antSelectBits = num2str(binaryIntVect(1:log2(nTxAnt)));
antSelectBits(isspace(antSelectBits)) = '';
antSelect = bin2dec(antSelectBits);
bpskBits = binaryIntVect(log2(nTxAnt)+1:end) * 2 - 1;

%OFDM params

Nsym = floor(minRangeDesired*2/c*symRate); %number of symbols in a message
Nc = min(ceil(length(bpskBits)/Nsym), maxNc);
nPulse = ceil(length(bpskBits)/Nsym/Nc);

pulseWidth=Nsym/symRate;

%basic radar calcs
minUnambigRange=pulseWidth*c/2;
maxUnambigRange=c*(PRI-pulseWidth)/2;

minRange = minRangeDesired;
maxRange = floor(maxUnambigRange/1000)*1000;

for pulseCount = 1:nPulse
    
    bits = bpskBits(Nc*Nsym*(pulseCount-1)+1:min(Nc*Nsym*pulseCount,length(bpskBits)));
    %basic signal calcs
    signal = buildWaveform(bits,sampPerSym,symRate,Nc,Nsym,guard);
    bandWidth = (Nc+1)*Nsym;

    radarDataLength = floor(((maxRange - minRange)/c + pulseWidth)*sampRate);

    %noise calculations
    noisePow = boltz*receiverBW*receiverTemp;
    radarData(:,pulseCount) = noisePow*(randn(1,radarDataLength)+1i*randn(1,radarDataLength));

end

%RECEIVED DATA
binaryReceived = uint16((bpskBits + 1) / 2);
txAntCalc = antSelect;

RGBrecovered = decodeIntImg(binaryReceived, txAntCalc, nTxAnt);
imshow(RGBrecovered)



%math stuff
plotPSD(1000,sampPerSym,symRate,Nc,Nsym,guard)

function [outSig] = timeDelay(inSig, fracSamp)
    % Design fractional-delay filter
    t = (-(31-1)/2:(31-1)/2)-fracSamp;
    delayFilter = sin(pi*t)./(pi*t);
    outSig = conv(inSig,delayFilter,'same');
end

function [outSig] = freqShift(inSig, sampRate, pulseTime, lambda, velRad, velTar, posRad, posTar)
    speed = dot((posTar-posRad),(velTar-velRad))/norm((posTar-posRad));
    shift = -2*(speed/lambda);
    t = (0:1/sampRate:(length(inSig)-1)/sampRate) + pulseTime;
    outSig = inSig .* exp(1i*2*pi*shift*t);
end

function [outSig] = ampScale(inSig, scale)
    outSig = inSig*scale;
end

function [scale] = calcScale(fc, txPower, rcs, range)
    scale = 1;
end