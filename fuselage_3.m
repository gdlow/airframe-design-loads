% General analysis of a simple beam structure with 1 reaction on mg
% aircraft is at an AoA to the ground - consider transverse component

clear;clc;

%% Initialize positional vars
[vals, names] = xlsread('geometryVariables.xlsx', 'Data', 'B:C');
geoParams = containers.Map(names(2:end,1), vals);
x_mg = geoParams('x_mg');
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
            distLoad(j+1) = distLoad(j+1) + massArray(i);
        end
    else
        distLoad(conc(i)*10+1) = distLoad(conc(i)*10+1) + massArray(i);
    end 
end

%% Obtain reaction forces
Rmg = sum(distLoad);

%% Obtain overall Bending moment on fuselage
% +ve ACW
overallMoments = sum((x_mg - x).*distLoad);
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
plot(x, shearForce);
figure;
plot(x, bendingMoments);