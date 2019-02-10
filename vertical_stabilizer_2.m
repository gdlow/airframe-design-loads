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
n = 1.5; % load factor = 1
concForce = 50 * 10^3; % [N]
ac_v = geoParams('ac_v'); % aerodynamic center of vertical stabilizer[m]
% applied torque due to horizontal stabilizer
ac_h = geoParams('ac_h'); % aerodynamic center of horizontal stabilizer [m]
netForceOnHorizontalPlane = 50 * 10^3; % [N]

%% Function call to get shear forces and bending moments
[x, chord, distLift, shearForce, bendingMoment] = vertical_stabilizer_load(rootChordLen_v, tipChordLen_v , wingSemiSpan_v, N , n, concForce);
bendingMoment = bendingMoment + (netForceOnHorizontalPlane * ac_h);

flexAxis = 0.45*chord; % Flex axis is 25 to 65 of chord
ac = 0.25*chord; % aerodynamic center
a = flexAxis - ac; % distence between the flex axis and AC

% distributed lift
figure;
hold on
grid on; grid minor;
title('Distributed Lift along Vertical Stabilizer');
xlabel('Span along Wing [m]');
ylabel('Load per unit meter [N/m]');
plot(x,distLift,'b', 'LineWidth', 2) % lift
legend('Elliptical Lift Distribution');
hold off

% bending moments
figure;
hold on;
grid on; grid minor;
title('Bending Moment Distribution across Vertical Stabilizer');
xlabel('Span along vertical stabilizer [m]');
ylabel('Moment per unit meter [N]');
plot(x, bendingMoment,'k', 'LineWidth', 2) % all
legend('Bending Moment on the Vertical Stabilizer due to lift');
hold off


% shear forces
figure;
hold on;
grid on; grid minor;
title('Shear Force Distribution across Vertical Stabilizer');
xlabel('Span along vertical stabilizer [m]');
ylabel('Shear Force per unit meter [N/m]');
plot(x, shearForce, 'k', 'LineWidth', 2)
legend('Shear Force on the Wing due to lift');
hold off

% torque
T = distLift.*a;
figure;
hold on
grid on; grid minor;
title('Torque Distribution across Vertical Stabilizer');
xlabel('Span along Vertical Stabilizer [m]');
ylabel('Torque per unit meter [N]');
plot(x, T, 'LineWidth', 2) % torque along span
hold off
