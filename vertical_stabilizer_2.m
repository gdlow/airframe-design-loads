%% General analysis of load case 2 on the vertical stabiliser
clear
clc
close all

%% Enter geometry of tailplane vertical stabilizer
rootChordLen = 5.96; % [m]
tipChordLen = 1.49; % [m]
wingSemiSpan = 22.36; % [m]
N = 100; % use 100 stations, according to lecturer 
n = 1; % load factor = 1
concForce = 50 * 10^3; % [N]
ac_v = 1; % aerodynamic center of vertical stabilizer[m]
% applied torque due to horizontal stabilizer
ac_h = 1; % aerodynamic center of horizontal stabilizer [m]
netForceOnHorizontalPlane = 50 * 10^3;

%% Function call to get shear forces and bending moments
[x, chord, distLoad, shearForce, bendingMoment] = vertical_stabilizer_load(rootChordLen, tipChordLen , wingSemiSpan, N , n, concForce);
bendingMoment = bendingMoment + (netForceOnHorizontalPlane * ac_h);

figure;
plot(x,bendingMoment)

figure;
plot(x,shearForce)
