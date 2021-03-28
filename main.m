clc;
close all;
clear all;

%constants
c = 299.8E6;%speed of light
boltz = 1.3806485E-23;%boltzman const

%simple Radar parameters
PRF=2000;
nPulse=16;
PRI=1/PRF;%seconds
fc=1E9;%Hz
txpower=40;%dBm
txgain=0;%dBi
rxgain=0;%dBi
receiverTemp=300;
receiverBW=1000000;

%simple signal parameters
sampPerSym = 30; %samples per symbol
symRate = 250000; %symbols per second
Nsym = 30; %number of symbols in a message
Nc = 8; %number of orthogonal carriers

%basic signal calcs
signal = buildWaveform(sampPerSym,symRate,Nc,Nsym,0);
sampRate = signal.sampRate;
bandWidth = (Nc+1)*Nsym;
pulseWidth=length(signal.IQ)/sampRate;

%basic radar calcs
minUnambigRange=pulseWidth*c/2;
maxUnambigRange=c*(PRI-pulseWidth)/2;

minRange = ceil(minUnambigRange/1000)*1000;
maxRange = floor(maxUnambigRange/1000)*1000;

%comms postitions
txAntPos = [20000,0,10];
rxAntPos = [0,0,0];

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

dataCubeLength = floor(((maxRange - minRange)/c + pulseWidth)*sampRate);

noisePow = boltz*receiverBW*receiverTemp;
dataCube = noisePow*(randn(nPulse,dataCubeLength)+1i*randn(nPulse,dataCubeLength));

%math stuff
%plotPSD(3000,sampPerSym,symRate,Nc,Nsym)