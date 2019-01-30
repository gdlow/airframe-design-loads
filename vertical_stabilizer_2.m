%% General analysis of load case 2 on the vertical stabiliser
clear
clc
close all

%% Enter geometry of tailplane vertical stabilizer
[vals, names] = xlsread('geometryVariables.xlsx', 'Data', 'B:C');
geoParams = containers.Map(names(2:end,1), vals);
rootChordLen_v = geoParams('rootChordLen_v'); % [m]
tipChordLen_v = geoParams('tipChordLen_v'); % [m]
aspectRatio_v = geoParams('aspectRatio_v');
MAC_v = geoParams('MAC_v'); % [m]
wingSemiSpan_v = aspectRatio_v * MAC_v / 2; % [m]
N = 100; % use 100 stations, according to lecturer 
n = 1; % load factor = 1
concForce = 50 * 10^3; % [N]
ac_v = geoParams('ac_v'); % aerodynamic center of vertical stabilizer[m]
% applied torque due to horizontal stabilizer
ac_h = geoParams('ac_h'); % aerodynamic center of horizontal stabilizer [m]
netForceOnHorizontalPlane = 50 * 10^3; % [N]

%% Function call to get shear forces and bending moments
[x, chord, distLoad, shearForce, bendingMoment] = vertical_stabilizer_load(rootChordLen_v, tipChordLen_v , wingSemiSpan_v, N , n, concForce);
bendingMoment = bendingMoment + (netForceOnHorizontalPlane * ac_h);

figure;
plot(x,bendingMoment)

figure;
plot(x,shearForce)
