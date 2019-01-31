%% General analysis of load case 1 on the horizontal stabiliser
clear
clc
close all

%% Enter geometry of tailplane horizontal stabiliser

% rootChordLen - root chord length [m]
% tipChordLen - tip chord length [m]
% wingSemiSpan - wing semi-span [m]
% takeOffWeight - Take-off weight of entire aircraft [N]
% N - Number of stations [-]
% n - Load factor [-]
% fuelTankLen - Length of fuel tank wrt half span [m]
[vals, names] = xlsread('geometryVariables.xlsx', 'Data', 'B:C');
geoParams = containers.Map(names(2:end,1), vals);
rootChordLen_h = geoParams('rootChordLen_h');
tipChordLen_h = geoParams('tipChordLen_h');
aspectRatio_h = geoParams('aspectRatio_h');
MAC_h = geoParams('MAC_h');
wingSemiSpan = aspectRatio_h*MAC_h/2;
weightStabilizer_h = geoParams('weightStabilizer_h'); % get from AVD report
N = 100; % according to lecturer 
n = 1.5*2.5;

%% Enter operating conditions and pitch aerodyamic moment
CM0_h = geoParams('CM0_h'); % pitch aerodynamic moment from airofoil data
cruiseVelocity = 232.78; % TODO: change this
rho = 1.1; % TODO: change this (at altitude 37000 ft)

%% Obtain lift generated by tailplane
% by considering bending moments around the cg of plane
MAC_wing = geoParams('MAC_wing');
aspectRatio_w = geoParams('aspectRatio_w');
Sref = aspectRatio_w * MAC_wing^2;
x_tail = geoParams('x_tail'); x_cg = geoParams('x_cg'); x_ac = geoParams('x_ac');
momentArmTailPlane = x_tail - x_cg;
momentArmWing = x_cg - x_ac; % assumes x_ac is before x_cg
CM0_w = geoParams('CM0_w');
takeOffWeight = geoParams('takeOffWeight');
% this gets it for half the entire horizontal stabilizer
horizontalStabilizerLift = 0.5 * (1/momentArmTailPlane) * ...
    (-0.5*rho*cruiseVelocity^2*Sref*MAC_wing*CM0_w + n*takeOffWeight*momentArmWing); % TODO: Check direction
%% Function call to get shear forces and bending moments
[x, chord, distLift, distWeightWing, shearForceWing, bendingMomentWing] = ...
    horizontal_stabilizer_load(rootChordLen_h, tipChordLen_h , wingSemiSpan , ...
    horizontalStabilizerLift,weightStabilizer_h, N , n);

flexAxis = 0.45*chord; % Flex axis is 25 to 65 of chord
ac = 0.25*chord; % aerodynamic center
a = flexAxis - ac; % distance between the flex axis and AC
cgWing = 0.55*chord;
b = cgWing - flexAxis; % difference between cg with respect to chord of wing and FA

% load distributions
distLoad = distLift+distWeightWing;

figure;
hold on
plot(x,distLift,'b') % lift
plot(x,distWeightWing,'r') % self weight
plot(x,distLoad,'k') % total load

hold off

% moment distribution
figure;
plot(x,bendingMomentWing,'k') % all


% shear force:
figure;
plot(x, shearForceWing,'k')

M_0 = 0.5*rho*cruiseVelocity^2.*chord.^2*CM0_h;
T = distLift.*a*wingSemiSpan+(distLoad-distLift).*b*wingSemiSpan-M_0;

figure;
hold on
plot(x,T) % torque along span
hold off



