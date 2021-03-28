clc;
close all;
clear all;

c = 299.8E6;%speed of light
boltz = 1.3806485E-23;%boltzman const

sampleRate=10E6;%Hz
PRF = 2000;
PRI=1/PRF;%seconds
fc=1E9;%Hz
txpower=40;%dBm
txgain=0;%dBi
rxgain=0;%dBi

sampPerSym = 30; %samples per symbol
symRate = 250000; %symbols per second
Nsym = 30; %number of symbols in a message
Nc = 8; %number of orthogonal carriers

signal = buildWaveform(sampPerSym,symRate,Nc,Nsym,0);
bandWidth = (Nc+1)*Nsym;
pulseWidth=length(signal.IQ)/signal.sampRate;

minRange=pulseWidth*c/2;
maxRange=c*(PRI-pulseWidth)/2;

receiverTemp=300;
receiverBW=1000;

txAntPos = [20000,0,10];
rxAntPos = [0,0,0];

targets(1).rcs = 0;
targets(1).initialPos = [0,0,0];
targets(1).velocity = [0,0,0];
targets(2).rcs = 5;
targets(2).initialPos = [1000,0,0];
targets(2).velocity = [0,0,0];
targets(3).rcs = 10;
targets(3).initialPos = [-1000,0,0];
targets(3).velocity = [0,0,0];

plotPSD(3000,sampPerSym,symRate,Nc,Nsym)