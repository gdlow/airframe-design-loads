%% General analysis of load case 1 on the wing group
clear
clc
close all

%% Enter geometry of wing

% rootChordLen - root chord length [m]
% tipChordLen - tip chord length [m]
% wingSemiSpan - wing semi-span [m]
% takeOffWeight - Take-off weight of entire aircraft [N]
% N - Number of stations [-]
% n - Load factor [-]
% fuelTankLen - Length of fuel tank wrt half span [-]
[vals, names] = xlsread('geometryVariables.xlsx', 'Data', 'B:C');
geoParams = containers.Map(names(2:end,1), vals);
rootChordLen_w = geoParams('rootChordLen_w');
tipChordLen_w = geoParams('tipChordLen_w');
wingSemiSpan = geoParams('wingSemiSpan');
takeOffWeight = geoParams('takeOffWeight');
N = 100; % according to lecturer 
n = 1.5*2.5;
fuelTankLen = geoParams('fuelTankLen'); % pct of semi-span from root

%% Enter operating conditions and pitch aerodyamic moment
CM0_w = geoParams('CM0_w'); % pitch aerodynamic moment from airofoil data
cruiseVelocity = 232.78;
rho = 1.1; % at altitude 37000 ft

%% Get values from wing_load and fuel_load

% self weight both lift
[x, chord, distLift, distWeightWing, shearForceWing, bendingMomentWing] = wing_load(rootChordLen_w, tipChordLen_w , wingSemiSpan , takeOffWeight*9.81 , N , n);

% Fuel
[distWeightFuel, shearForceFuel, bendingMomentFuel] = fuel_load(rootChordLen_w, tipChordLen_w , wingSemiSpan , takeOffWeight*9.81 , N , n, fuelTankLen);


flexAxis = 0.45*chord; % Flex axis is 25 to 65 of chord
ac = 0.25*chord; % aerodynamic center
a = flexAxis - ac; % distence between the flex axis and AC
cgWing = 0.55*chord;
b = cgWing - flexAxis; % difference between cg with respect to chord of wing and FA

% Engine: (assumed point load)

% first iteration assumed point load
% 2 entries so it is the same size for plotting
engineLoc_1 = geoParams('engineLoc_1'); % [m]
engineLoc_2 = geoParams('engineLoc_2'); % [m]
ew = geoParams('engineWeight');
engineWeight_1=-ew*9.81; % in [N]
engineWeight_2=-ew*9.81; % in [N]

% load distributions
distLoad = distLift+distWeightWing+distWeightFuel;
panelDist = wingSemiSpan/N;

loc_1 = find(x>engineLoc_1-panelDist,1); % position of engine 1
loc_2 = find(x>engineLoc_2-panelDist,1); % position of engine 2
% add concentrated loads from engines
distLoad(loc_1) = distLoad(loc_1)+engineWeight_1; % total load with engine 1
distLoad(loc_2) = distLoad(loc_2)+engineWeight_2; % total load with engine 2

figure;
hold on
plot(x,distLift,'b') % lift
plot(x,distWeightWing,'r') % self weight
plot(x,distWeightFuel,'g') % fuel weight

plot(x,distLoad,'k') % total load
hold off

totalMoments = bendingMomentWing+bendingMomentFuel;
for i=1:loc_1
    totalMoments(i) = totalMoments(i) - (engineWeight_1*(engineLoc_1-x(i)));
end
for i=1:loc_2
    totalMoments(i) = totalMoments(i) - (engineWeight_2*(engineLoc_2-x(i)));
end
% moment distribution
figure;
hold on
plot(x,bendingMomentWing,'r') % wing self weight and lift 
plot(x,bendingMomentFuel,'b') % fuel
plot(x,totalMoments,'k') % all
hold off

% shear force:
totalShearForce = shearForceWing+shearForceFuel;
% add concentrated loads from engines
totalShearForce(1:loc_1) = totalShearForce(1:loc_1) + engineWeight_1;
totalShearForce(1:loc_2) = totalShearForce(1:loc_2) + engineWeight_2;

figure;
hold on
plot(x,shearForceWing,'r')
plot(x,shearForceFuel,'b')
plot(x,totalShearForce,'k')
hold off

M_0=0.5*rho*cruiseVelocity^2.*chord.^2*CM0_w;
T=distLift.*a*wingSemiSpan+(distLoad-distLift).*b*wingSemiSpan-M_0;

figure;
hold on
plot(x,T) % torque along span
hold off
