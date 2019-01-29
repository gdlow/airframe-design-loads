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

rootChordLen = 3.62;
tipChordLen = 1.45;
aspectRatio = 4;
MACtailPlane = 2.54;
wingSemiSpan = aspectRatio*MACtailPlane/2;
horizontalStabilizerWeight = 0.02 * 78*10^3 * 9.81; % get from AVD report
N = 100; % according to lecturer 
n = 1.5*2.5;

%% Enter operating conditions and pitch aerodyamic moment
Cm0TailPlane = 1.1; % pitch aerodynamic moment from airofoil data
cruiseVelocity = 232.78;
rho = 1.1; % at altitude 37000 ft

%% Obtain lift generated by tailplane
% by considering bending moments around the cg of plane
MACWing = 4.85;
Sref = 6.89*36.41^2;
momentArmTailPlane = 25;
momentArmWing = 3; %xcg-xac
CM0Wing = 1.1;
takeOffWeight = 78*10^3 * 9.81; % get from AVD report
% this gets it for half the entire horizontal stabilizer
horizontalStabilizerLift = 0.5 * (1/momentArmTailPlane) * ...
    (-0.5*rho*cruiseVelocity^2*Sref*MACWing*CM0Wing + n*takeOffWeight*momentArmWing);
%% Function call to get shear forces and bending moments (Left Wing)
leftWingLift = 75*10^3; % N
[x_l, chord_l, distLift_l, distWeightWing_l, shearForceWing_l, bendingMomentWing_l] = ...
    horizontal_stabilizer_load(rootChordLen, tipChordLen , wingSemiSpan , ...
    leftWingLift,horizontalStabilizerWeight, N , n);

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

M_0 = 0.5*rho*cruiseVelocity^2.*chord_l.^2*Cm0TailPlane;
T = distLift_l.*a*wingSemiSpan+(distLoad_l-distLift_l).*b*wingSemiSpan-M_0;

figure;
hold on
plot(x_l,T) % torque along span
hold off

%% Function call to get shear forces and bending moments (Right Wing)
rightWingLift = 25*10^3; % N
[x_r, chord_r, distLift_r, distWeightWing_r, shearForceWing_r, bendingMomentWing_r] = ...
    horizontal_stabilizer_load(rootChordLen, tipChordLen , wingSemiSpan , ...
    rightWingLift,horizontalStabilizerWeight, N , n);

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

M_0 = 0.5*rho*cruiseVelocity^2.*chord_r.^2*Cm0TailPlane;
T = distLift_r.*a*wingSemiSpan+(distLoad_r-distLift_r).*b*wingSemiSpan-M_0;

figure;
hold on
plot(x_r,T) % torque along span
hold off
