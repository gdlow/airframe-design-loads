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
[x, chord, distLift, distWeightWing, shearForceWing, bendingMomentWing] = wing_load(rootChordLen_w, tipChordLen_w , wingSemiSpan , takeOffWeight , N , n);

% Fuel
[distWeightFuel, shearForceFuel, bendingMomentFuel] = fuel_load(rootChordLen_w, tipChordLen_w , wingSemiSpan , takeOffWeight , N , n, fuelTankLen);


flexAxis = 0.45*chord; % Flex axis is 25 to 65 of chord
ac = 0.25*chord; % aerodynamic center
a = flexAxis - ac; % distence between the flex axis and AC
cgWing = 0.55*chord;
b = cgWing - flexAxis; % difference between cg with respect to chord of wing and FA

% Engine: (assumed point load)

% first iteration assumed point load
% 2 entries so it is the same size for plotting
e1 = geoParams('engineLoc_1');
e2 = geoParams('engineLoc_2');
ew = geoParams('engineWeight');
engineLoc_1=[e1,e1]; % in [m]
engineWeight_1=[0,-ew*9.81]; % in [N]

engineLoc_2=[e2,e2]; % in [m]
engineWeight_2=[0,-ew*9.81]; % in [N]

% load distributions
distLoad = distLift+distWeightWing+distWeightFuel;
panelDist = wingSemiSpan/N;

loc_1 = find(x>engineLoc_1(1)-panelDist,1); % position of engine 1
loc_2 = find(x>engineLoc_2(1)-panelDist,1); % position of engine 2

distLoad(loc_1) = distLoad(loc_1)+engineWeight_1(2); % total load with engine 1
distLoad(loc_2) = distLoad(loc_2)+engineWeight_2(2); % total load with engine 2

figure;
hold on
plot(x,distLift,'b') % lift
plot(x,distWeightWing,'r') % self weight
plot(x,distWeightFuel,'g') % fuel weight

plot(x,distLoad,'k') % total load
hold off

totalMoments = bendingMomentWing+bendingMomentFuel-((engineWeight_1(1)*engineLoc_1(2))+(engineWeight_2(1)*engineLoc_2(2)));

% moment distribution
figure;
hold on
plot(x,bendingMomentWing,'r') % wing self weight and lift 
plot(x,bendingMomentFuel,'b') % fuel
plot(x,totalMoments,'k') % all
hold off

% shear force:
totalShearForce = shearForceWing+shearForceFuel;

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
