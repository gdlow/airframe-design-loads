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
            distLoad(j+1) = distLoad(j+1) + massArray(i);
        end
    else
        distLoad(conc(i)*10+1) = distLoad(conc(i)*10+1) + massArray(i);
    end 
end

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

%% Fuselage structural calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except bendingMoments shearForce x x_frontSpar x_rearSpar x_tail
clc
close all

% Variables
b_f = 0.508; % Frame pitch, from design review [m]
A_s = 500*10^-6; % Stringer area [m^2]
n_stringers = 80; % Number of stringers, from design review
d_fuselage = 4; % Fuselage diameter [m] 
sigma_yield_stringer = 350*10^6; % Yield stress of stringer [Pa]  
t_skin = 0.002; % skin thickness, from design review [m]
E_s = 70*10^9; % Youngs modulus of stringer [Pa]
E_skin = 70*10^9; % Youngs modulus of skin [Pa]
rho_s = 2700; % Density of stringer [kg/m^3]
rho_skin = 2700; % Density of skin [kg/m^3]
R_t = 1;

% Extra thickness due to pressure
delta_P = 

% Stringer positions
b_s = pi*d_fuselage/n_stringers;

stringer_angles = [0:360/n_stringers:360*(1 - 1/n_stringers)]; % angle from horizontal plane, measured around central axis [deg]
x_s = d_fuselage./2.*cosd(stringer_angles); % x-position of stringers
y_s = d_fuselage./2.*sind(stringer_angles); % y-position of stringers

% Stringer sizing due to bending
M_max = max(bendingMoments);
y_max = max(y_s);

A_skin_ideal = t_skin.*b_s./6.*(([y_s(2:end),y_s(1)] + [y_s(end),y_s(1:end - 1)])./y_s + 4);
A_skin_ideal(isnan(A_skin_ideal)) = 0;

I_xx_skin = sum(A_skin_ideal.*y_s.^2);
I_xx_s = A_s*sum(y_s.^2);

I_xx_tot = I_xx_skin + I_xx_s;

sigma_bending = M_max*y_max/I_xx_tot;

% Check that stringers can withstand euler buckling
R_A = A_s/(t_skin*b_f);

h_s = R_A/(1.6*R_t)*b_s;
d_s = 0.3*h_s;
t_s = R_t*t_skin;

I_xx_single_stringer = d_s*t_s^3/12 + (h_s-t_s)^2/4*d_s*t_s + (h_s - 2*t_s)^3*t_s/12;
sigma_euler = pi^2*E_s*I_xx_single_stringer/(A_s*b_f^2);

% Shear Flow
sectorAngle = 360/n_stringers;
y = zeros(1,n_stringers);
for i=1:length(y)
    y(i) = d_fuselage/2*sind((i-1)*sectorAngle); % height of boom above neutral axis
end
A_boom = A_s + A_skin_ideal(2);
Sy = max(shearForce); 
Ixx = A_boom* sum(y.*y);
cumulative = 0;
qb = zeros(1,n_stringers-1);
for i=2:n_stringers
    qb(i) = cumulative - Sy/Ixx*sum(A_boom*y(i));
    cumulative = qb(i);
end
qs0 = -sum(qb)/n_stringers;
qs = qb + qs0;

% Calculate whether or not plate buckling will occur
sigma_ratio = 0.9;
sigma_plate = sigma_ratio*3.62*E_skin*(t_skin/b_s)^2;


% Calculate weight per unit length of skin-stringer configuration
W_skin = d_fuselage*pi*t_skin * rho_skin;
W_s = n_stringers*A_s * rho_s;

W = W_skin + W_s;

% Plot Shear Flow
qsGraph = abs(qs);
theta = 0:0.1:360;
x1 = d_fuselage/2*cosd(theta);
y1 = d_fuselage/2*sind(theta);
figure;
plot(x1,y1);
hold on;
startPanel = 1;
lenPanel = floor(length(theta)/n_stringers);
for i=1:n_stringers
    xPanel = (d_fuselage/2 + qsGraph(i)/(10^7))*cosd(theta(startPanel:startPanel+lenPanel-1));
    yPanel = (d_fuselage/2 + qsGraph(i)/(10^7))*sind(theta(startPanel:startPanel+lenPanel-1));
    startPanel = startPanel + lenPanel;
    plot(xPanel, yPanel, 'r');
    plot(d_fuselage/2*sind((i-1)*sectorAngle), ...
    d_fuselage/2*cosd((i-1)*sectorAngle), 'k.', 'MarkerSize', 10);
end
hold off;

disp(['sigma_bending = ',num2str(sigma_bending/10^6),' MPa'])
disp(['sigma_euler = ',num2str(sigma_euler/10^6),' MPa'])
disp(['sigma_plate = ',num2str(sigma_plate/10^6),' MPa'])