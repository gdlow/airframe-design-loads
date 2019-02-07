% General analysis of a simple beam structure with 1 reaction on mg
% aircraft is at an AoA to the ground - consider transverse component

clear;clc;

%% Initialize positional vars
[vals, names] = xlsread('geometryVariables.xlsx', 'Data', 'B:C');
geoParams = containers.Map(names(2:end,1), vals);
x_mg = geoParams('x_mg');
aspectRatio_w = geoParams('aspectRatio_w');
MAC_wing = geoParams('MAC_wing');
Sref = aspectRatio_w * MAC_wing^2;
CM0_w = geoParams('CM0_w');
rho = 1; % at 0ft
V_TO = 100; % [ms-1] %TODO: Change this variable
alpha_TO = geoParams('alpha_TO');
%% Get x distribution across fuselage
lenFuselage = geoParams('lenFuselage');
x = 0:0.1:lenFuselage;

%% Re-orientate x-z plane to alpha_TO
x_mg = x_mg * cosd(alpha_TO);
x = x .* cosd(alpha_TO);
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
n = 1*9.81;
distLoad = distLoad.*n;

%% Obtain reaction forces
Rmg = sum(distLoad);

%% Obtain overall Bending moment on fuselage
% +ve ACW

M0_g = 0.5*rho*V_TO^2*Sref*MAC_wing*CM0_w;
overallMoments = sum((x_mg - x).*distLoad) + M0_g;
%% Obtain shear force and bending moment distribution
% take moments about cut face - moment arm is x(i)-x
shearForce = zeros(1, length(x));
bendingMoments = zeros(1, length(x));
for i=1:length(x)
    if x(i) < x_mg
        shearForce(i) = sum(distLoad(1:i));
        bendingMoments(i) = sum((x(i)-x(1:i)).*distLoad(1:i));
    else
        shearForce(i) = sum(distLoad(1:i)) - Rmg;
        bendingMoments(i) = sum((x(i)-x(1:i)).*distLoad(1:i)) - Rmg*(x(i)-x_mg);
    end
end
bendingMoments = bendingMoments - overallMoments; % CW - ACW; bendingMoments defined CW on +ve face
%% Plot result
figure;
plot(x, distLoad);
figure;
hold on;
grid on;
title('Shear Force along fuselage length');
xlabel('x [m]');
ylabel('Shear Force [N]');
plot(x, shearForce, 'LineWidth', 2);
hold off;
figure;
hold on;
grid on;
title('Bending Moment along fuselage length');
xlabel('x [m]');
ylabel('Bending Moment [Nm]');
plot(x, bendingMoments, 'LineWidth', 2);
hold off;