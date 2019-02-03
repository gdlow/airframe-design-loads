% General analysis for asymmetric loading scenario
% Analysis of cantilever beam constrained at the wing box
% Analysis considered only for aft fuselage
clear;
clc;
close all;
%% Enter geometrical positions
[vals, names] = xlsread('geometryVariables.xlsx', 'Data', 'B:C');
geoParams = containers.Map(names(2:end,1), vals);
radiusFuselage = geoParams('radiusFuselage');
lenFuselage = geoParams('lenFuselage');
x = 0:0.1:lenFuselage;
wingBoxLoc = geoParams('wingBoxLoc');
x_tail = geoParams('x_tail');
% vertical stabiliser
rootChordLen_v = geoParams('rootChordLen_v');
tipChordLen_v = geoParams('tipChordLen_v');
MAC_v = geoParams('MAC_v');
aspectRatio_v = geoParams('aspectRatio_v');
spanVert = aspectRatio_v * MAC_v / 2;
% z location measured from center of fuselage
zLoc = radiusFuselage + spanVert*(rootChordLen_v - MAC_v) / (rootChordLen_v - tipChordLen_v);
% horizontal stabiliser
rootChordLen_h = geoParams('rootChordLen_h');
tipChordLen_h = geoParams('tipChordLen_h');
MAC_h = geoParams('MAC_h');
aspectRatio_h = geoParams('aspectRatio_h');
spanHor = aspectRatio_h * MAC_h / 2;
% on both ends of horizontal stabiliser
yLoc = spanHor*(rootChordLen_h - MAC_h) / (rootChordLen_h - tipChordLen_h);

%% Obtain torque around fuselage and resultant shear flow
T = 50*10^3*zLoc + 50*10^3*yLoc; % CW
areaFuselage = pi*radiusFuselage^2;
q = T/(2*areaFuselage); % extra shear flow due to torque

%% Obtain distributed loads
% Components are either taken to be point loads or evenly distributed loads
% all measurements in S.I units
dataRead = xlsread('loadCaseVariables.xlsx','Data','C:F');
massArray = dataRead(:,1); distStart = dataRead(:,2); distEnd = dataRead(:,3);
conc = dataRead(:,4);
distLoad = zeros(1, length(x));
for i=1:length(massArray)
    if conc(i) == 999
        for j=distStart(i)*10:distEnd(i)*10
            distLoad(j+1) = distLoad(j+1) + massArray(i)/(distEnd(i)-distStart(i));
        end
    else
        distLoad(conc(i)*10+1) = distLoad(conc(i)*10+1) + massArray(i);
    end 
end

%% Multiply by load factor*g
n = 1*9.81;
distLoad = distLoad.*n;

%% Reconfigure distances to start from wingbox
distLoad = distLoad(wingBoxLoc*10:end);
x = 0:0.1:length(distLoad)/10; x = x(1:end-1);
x_tail = x_tail - wingBoxLoc;

%% Obtain reaction forces and moments at wing box
% Force equilibrium
Rf = sum(distLoad) - 100*10^3;
% Moments equilibrium (ACW on wall)
Mf = sum(x.*distLoad) - x_tail*100*10^3;

%% Obtain shear forces and bending moments distribution
% take moments about cut face - moment arm is x(i)-x
shearForce = zeros(1, length(x));
bendingMoments = zeros(1, length(x));
for i=1:length(x)
    if x(i) < x_tail
        shearForce(i) = sum(distLoad(1:i)) - Rf;
        bendingMoments(i) = sum((x(i)-x(1:i)).*distLoad(1:i)) + Mf - Rf*x(i);
    else
        shearForce(i) = sum(distLoad(1:i)) - Rf - 100*10^3;
        bendingMoments(i) = sum((x(i)-x(1:i)).*distLoad(1:i)) + Mf - Rf*x(i) - 100*10^3*(x(i)-x_tail);
    end
end

%% Plot result
figure;
plot(x, distLoad);
figure;
plot(x, shearForce);
figure;
plot(x, bendingMoments);

%% Shear flow
thickness = geoParams('thickness');
radiusFuselage = geoParams('radiusFuselage');
n = 80;
sectorAngle = 360/n;
pitch = pi*2*radiusFuselage/n; % distance between booms
As = geoParams('As'); % stringer x sectional area
y = zeros(1,n);
for i=1:length(y)
    y(i) = radiusFuselage*sind((i-1)*sectorAngle); % height of boom above neutral axis
end
boomArea = As + (thickness*pitch/6)*(2+y(3)/y(2))+(thickness*pitch/6)*(2+y(1)/y(2));
Sy = min(shearForce); % since Sy is defined negative
Ixx = boomArea * sum(y.*y);
cumulative = 0;
qb = zeros(1,n-1);
for i=2:n
    qb(i) = cumulative - Sy/Ixx*sum(boomArea*y(i-1));
    cumulative = qb(i);
end
qs0 = -sum(qb)/n;
qs = qb + qs0 + q;
qsGraph = abs(qs);

%% Plot Shear Flow
theta = 0:0.1:360;
x1 = radiusFuselage*cosd(theta);
y1 = radiusFuselage*sind(theta);
figure;
plot(x1,y1);
hold on;
startPanel = 1;
lenPanel = floor(length(theta)/n);
for i=1:n
    xPanel = (radiusFuselage + qsGraph(i)/(10^7))*cosd(theta(startPanel:startPanel+lenPanel-1));
    yPanel = (radiusFuselage + qsGraph(i)/(10^7))*sind(theta(startPanel:startPanel+lenPanel-1));
    startPanel = startPanel + lenPanel;
    plot(xPanel, yPanel, 'r');
    plot(radiusFuselage*sind((i-1)*sectorAngle), ...
    radiusFuselage*cosd((i-1)*sectorAngle), 'k.', 'MarkerSize', 10);
end
hold off;

%% Bending about x-y plane
pointLoadxy = 50*10^3;
Rfxy = pointLoadxy;
Mfxy = pointLoadxy*x_tail; %ACW
shearForcexy = zeros(1, length(x));
shearForcexy(1:x_tail*10) = - Rfxy;
shearForcexy(x_tail*10+1:end) = pointLoadxy - Rfxy;
bendingMomentsxy = zeros(1,length(x));
bendingMomentsxy(1:x_tail*10) = Mfxy - x(1:x_tail*10).*Rfxy;
bendingMomentsxy(x_tail*10+1:end) = Mfxy - x(x_tail*10+1:end).*Rfxy + (x(x_tail*10+1:end) - x_tail).*pointLoadxy;

figure;
plot(x, shearForcexy);
figure;
plot(x, bendingMomentsxy);