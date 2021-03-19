clc;
close all;
clear all;

sampleRate=0;
minRange=0;
maxRange=0;
PRI=0;
pulseLength=0;
txpower=0;
txgain=0;
rxgain=0;
fc=0;
receiverTemp=0;
receiverBW=0;
c = 299.8E6;
boltz = 1.3806485E-23;

txAntPos = [];
rxAntPos = [];

targets(1).rcs = 5;
targets(1).initialPos = [0,0,0];
targets(1).velocity = [0,0,0];
targets(2).rcs = 5;
targets(2).initialPos = [0,0,0];
targets(2).velocity = [0,0,0];
targets(3).rcs = 5;
targets(3).initialPos = [0,0,0];
targets(3).velocity = [0,0,0];

plotPSD(10000)