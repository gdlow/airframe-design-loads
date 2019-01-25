% General analysis of a simple beam structure with 3 supported ends
clear;clc;

%% Guess RT and initialize positional vars
RT = 0;
xT = 38.0;
xF = 15;
xR = 18;
%% Get x distribution across fuselage
lenFuselage = 42.6;
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
RR = (sum(distLoad.*x) - xT*RT - xF*(sum(distLoad)-RT))/(xR-xF);
RF = (sum(distLoad.*x) - xT*RT - xR*(sum(distLoad)-RT))/(xF-xR);

%% Obtain shear force and bending moment distribution
% take moments about cut face - moment arm is x(i)-x
shearForce = zeros(1, length(x));
bendingMoments = zeros(1, length(x));
for i=1:length(x)
    if x(i) < xF
        shearForce(i) = sum(distLoad(1:i));
        bendingMoments(i) = sum((x(i)-x(1:i)).*distLoad(1:i));
    elseif x(i) < xR
        shearForce(i) = sum(distLoad(1:i)) - RF;
        bendingMoments(i) = sum((x(i)-x(1:i)).*distLoad(1:i)) - RF*(x(i)-xF);
    elseif x(i) < xT
        shearForce(i) = sum(distLoad(1:i)) - RF - RR;
        bendingMoments(i) = sum((x(i)-x(1:i)).*distLoad(1:i)) - RF*(x(i)-xF) - RR*(x(i)-xR);
    else
        shearForce(i) = sum(distLoad(1:i)) - RF - RR - RT;
        bendingMoments(i) = sum((x(i)-x(1:i)).*distLoad(1:i)) - RF*(x(i)-xF) - RR*(x(i)-xR) - RT*(x(i)-xT);
    end
end

%% Plot result
figure;
plot(x, distLoad);
figure;
plot(x, shearForce);
figure;
plot(x, bendingMoments);