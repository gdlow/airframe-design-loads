%% General analysis of load case 2 on the horizontal stabiliser
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
tipChordLen = geoParams('tipChordLen_h');
aspectRatio_h = geoParams('aspectRatio_h');
MAC_h = geoParams('MAC_h');
wingSemiSpan = aspectRatio_h*MAC_h/2;
weightStabilizer_h = geoParams('weightStabilizer_h'); % get from AVD report
N = 100; % according to lecturer 
n = 1;

%% Enter operating conditions and pitch aerodyamic moment
CM0_h = geoParams('CM0_h'); % pitch aerodynamic moment from airofoil data
cruiseVelocity = 232.78; 
rho = 0.350085776097978;

%% Function call to get shear forces and bending moments (Left Wing)
leftWingLift = 75*10^3; % N
[x_l, chord_l, distLift_l, distWeightWing_l, shearForceWing_l, bendingMomentWing_l] = ...
    horizontal_stabilizer_load(rootChordLen_h, tipChordLen , wingSemiSpan , ...
    leftWingLift,weightStabilizer_h, N , n);

flexAxis = 0.45*chord_l; % Flex axis is 25 to 65 of chord
ac = 0.25*chord_l; % aerodynamic center
a = flexAxis - ac; % distence between the flex axis and AC
cgWing = 0.55*chord_l;
b = cgWing - flexAxis; % difference between cg with respect to chord of wing and FA

% load distributions
distLoad_l = distLift_l+distWeightWing_l;

figure;
hold on
plot(x_l,distLift_l,'b') % lift
plot(x_l,distWeightWing_l,'r') % self weight
plot(x_l,distLoad_l,'k') % total load
hold off

% moment distribution
figure;
plot(x_l,bendingMomentWing_l,'k') % all


% shear force:
figure;
plot(x_l, shearForceWing_l,'k')

M_0_l = 0.5*rho*cruiseVelocity^2.*chord_l.^2*CM0_h;
T_l = distLift_l.*a-(distLoad_l-distLift_l).*b-M_0_l;

figure;
hold on
plot(x_l,T_l) % torque along span
hold off

%% Function call to get shear forces and bending moments (Right Wing)
rightWingLift = 25*10^3; % N
[x_r, chord_r, distLift_r, distWeightWing_r, shearForceWing_r, bendingMomentWing_r] = ...
    horizontal_stabilizer_load(rootChordLen_h, tipChordLen , wingSemiSpan , ...
    rightWingLift,weightStabilizer_h, N , n);

flexAxis = 0.45*chord_r; % Flex axis is 25 to 65 of chord
ac = 0.25*chord_r; % aerodynamic center
a = flexAxis - ac; % distence between the flex axis and AC
cgWing = 0.55*chord_r;
b = cgWing - flexAxis; % difference between cg with respect to chord of wing and FA

% load distributions
distLoad_r = distLift_r+distWeightWing_r;

figure;
hold on
plot(x_r,distLift_r,'b') % lift
plot(x_r,distWeightWing_r,'r') % self weight
plot(x_r,distLoad_r,'k') % total load

hold off

% moment distribution
figure;
plot(x_r,bendingMomentWing_r,'k') % all


% shear force:
figure;
plot(x_r, shearForceWing_r,'k')

M_0_r = 0.5*rho*cruiseVelocity^2.*chord_r.^2*CM0_h;
T_r = distLift_r.*a-(distLoad_r-distLift_r).*b-M_0_r;

figure;
hold on
plot(x_r,T_r) % torque along span
hold off
