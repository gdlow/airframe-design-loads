% General analysis of a simple beam structure with 3 supported ends
clear;clc; close all;

%% Set RT=0 and initialize positional vars
RT = 0; % Reaction at tail
[vals, names] = xlsread('geometryVariables.xlsx', 'Data', 'B:C');
geoParams = containers.Map(names(2:end,1), vals);
x_tail = geoParams('x_tail');
x_frontSpar = geoParams('x_frontSpar'); % distance to front spar
x_rearSpar = geoParams('x_rearSpar'); % distance to rear spar
%% Get x distribution across fuselage
lenFuselage = geoParams('lenFuselage');
x = 0:0.1:lenFuselage;

%% Get distributed loads
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
n = 1.5*2.5*9.81;
distLoad = distLoad.*n;
%% Obtain reaction forces
RR = (sum(distLoad.*x) - x_tail*RT - x_frontSpar*(sum(distLoad)-RT))/(x_rearSpar-x_frontSpar);
RF = (sum(distLoad.*x) - x_tail*RT - x_rearSpar*(sum(distLoad)-RT))/(x_frontSpar-x_rearSpar);

%% Obtain shear force and bending moment distribution
% take moments about cut face - moment arm is x(i)-x
shearForce = zeros(1, length(x));
bendingMoments = zeros(1, length(x));
for i=1:length(x)
    if x(i) < x_frontSpar
        shearForce(i) = sum(distLoad(1:i));
        bendingMoments(i) = sum((x(i)-x(1:i)).*distLoad(1:i));
    elseif x(i) < x_rearSpar
        shearForce(i) = sum(distLoad(1:i)) - RF;
        bendingMoments(i) = sum((x(i)-x(1:i)).*distLoad(1:i)) - RF*(x(i)-x_frontSpar);
    elseif x(i) < x_tail
        shearForce(i) = sum(distLoad(1:i)) - RF - RR;
        bendingMoments(i) = sum((x(i)-x(1:i)).*distLoad(1:i)) - RF*(x(i)-x_frontSpar) - RR*(x(i)-x_rearSpar);
    else
        shearForce(i) = sum(distLoad(1:i)) - RF - RR - RT;
        bendingMoments(i) = sum((x(i)-x(1:i)).*distLoad(1:i)) - RF*(x(i)-x_frontSpar) - RR*(x(i)-x_rearSpar) - RT*(x(i)-x_tail);
    end
end

%% Plot result
figure;
plot(x, distLoad);
figure;
plot(x, shearForce);
figure;
plot(x, bendingMoments);

%% Shear Flow
thickness = geoParams('thickness');
radiusFuselage = geoParams('radiusFuselage');
n = 80;
sectorAngle = 360/n;
pitch = pi*2*radiusFuselage/n; % distance between booms
As = geoParams('As'); % stringer cross sectional area
y = zeros(1,n);
for i=1:length(y)
    y(i) = radiusFuselage*sind((i-1)*sectorAngle); % height of boom above neutral axis
end
boomArea = As + (thickness*pitch/6)*(2+y(3)/y(2))+(thickness*pitch/6)*(2+y(1)/y(2));
Sy = max(shearForce); 
Ixx = boomArea * sum(y.*y);
cumulative = 0;
qb = zeros(1,n-1);
for i=2:n
    qb(i) = cumulative - Sy/Ixx*sum(boomArea*y(i-1));
    cumulative = qb(i);
end
qs0 = -sum(qb)/n;
qs = qb + qs0;
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

% to show relative shear flow on panels, add qs (/10^7) to current x1 and y1