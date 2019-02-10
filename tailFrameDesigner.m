% Script to design an optimal airframe structure for an empennage of T-tail
% configuration. 

% Written by A.Bates for AFD Group 17, Feb 2019:

%                    ** Begin Script **

%% Clear command window and workspace, close all figures:
clear();
clc();

close all;

%% Generic Terms:
N = 100;        %Number of discretised points along surface to consider.

%% Load loading data from excel file:
fileName = 'tailLoading.xlsx';
loadCases = ["H1", "H2_r", "H2_l", "V2"];

loadCaseN = length(loadCases);

for i = 1:loadCaseN
    loadParams(i).case = loadCases(i);
    loadParams(i).Torque = xlsread(fileName, loadParams(i).case, "B5:B104");
    loadParams(i).Moment = xlsread(fileName, loadParams(i).case, "C5:C104");
    loadParams(i).Shear = xlsread(fileName, loadParams(i).case, "D5:D104");
    loadParams(i).Load = xlsread(fileName, loadParams(i).case, "G5:G104");
end

%% Horizontal Stabiliser Geometric Parameters:
hori.b = 10.16;     %[m] - Tip to tip span - not half span!
hori.AR = 4;
hori.taper = 0.4;
hori.lSweep = 35.0;     %[Deg] - Leading edge sweep
hori.cBar = 2.54;       %[m] -  Mean aerodynamic chord
hori.twist = 0;         %[Deg]
hori.dihedral = 0;      %[Deg]
hori.tcr = 0.12;
hori.cR = 3.62;         %[m] - Root chord length
hori.cT = 1.45;         %[m] - Tip chord length
hori.sRef = 25.74;       %[m^2] - Ref area, tip to tip

%% Vertial Stabiliser Geometric Parameters:
vert.b = 5.423;     %[m] - Root to tip span
vert.AR = 1.2;
vert.taper = 0.9;
vert.lSweep = 45.0;     %[Deg] - Leading edge sweep
vert.cBar = 4.52;       %[m] -  Mean aerodynamic chord
vert.twist = 0;         %[Deg]
vert.dihedral = 0;      %[Deg]
vert.tcr = 0.12;
vert.cR = 4.76;         %[m] - Root chord length
vert.cT = 4.28;         %[m] - Tip chord length
vert.sRef = 24.5;       %[m^2] - Ref area, tip to tip

%% Aerofoil Section Data:
hori.lSection = [
  0.000000  -0.000000	;
  0.000150  -0.009560	;
  0.006280  -0.020300	;
  0.018650  -0.031760	;
  0.037300  -0.043240	;
  0.062030  -0.053820	;
  0.092300  -0.062650	;
  0.127320  -0.069150	;
  0.166040  -0.073200	;
  0.207380  -0.075240	;
  0.251310  -0.075970	;
  0.297960  -0.075540	;
  0.346810  -0.074020	;
  0.397330  -0.071500	;
  0.448970  -0.068110	;
  0.501170  -0.063970	;
  0.553350  -0.059240	;
  0.604960  -0.054050	;
  0.655410  -0.048540	;
  0.704170  -0.042850	;
  0.750700  -0.037120	;
  0.794490  -0.031450	;
  0.835060  -0.025970	;
  0.871970  -0.020790	;
  0.904820  -0.016020	;
  0.933240  -0.011760	;
  0.956930  -0.008120	;
  0.975630  -0.005180	;
  0.989140  -0.003020	;
  0.997300  -0.001700	;
  1.000030  -0.001260]   ;

hori.uSection = [
  0.000000 0.000000	;
  0.005330 0.007920	;
  0.015570 0.014010	;
  0.030290 0.018700	;
  0.049150 0.022480	;
  0.071950 0.025860	;
  0.098680 0.029220	;
  0.129540 0.032820	;
  0.164830 0.036600	;
  0.204830 0.040160	;
  0.248690 0.042830	;
  0.295310 0.044460	;
  0.344180 0.045100	;
  0.394760 0.044820	;
  0.446500 0.043710	;
  0.498830 0.041880	;
  0.551170 0.039450	;
  0.602960 0.036550	;
  0.653600 0.033270	;
  0.702570 0.029750	;
  0.749300 0.026070	;
  0.793300 0.022350	;
  0.834070 0.018660	;
  0.871180 0.015120	;
  0.904200 0.011800	;
  0.932790 0.008800	;
  0.956610 0.006210	;
  0.975430 0.004100	;
  0.989010 0.002540	;
  0.997220 0.001580	;
  0.999970 0.001260];

hori.iBox = [0.3 0.75 0.75 0.3; max(hori.uSection(:, 2)) max(hori.uSection(:, 2)) min(hori.lSection(:, 2)) min(hori.lSection(:,2))];
hori.boxT = 0.12;
hori.boxL = 0.75 - 0.3;

hSec = figure;
hold on;
plot(hori.uSection(:,1), hori.uSection(:,2), hori.lSection(:,1), hori.lSection(:,2));
plot(hori.iBox(1, :), hori.iBox(2, :), 'k'); 
axis equal;

vert.uSection = [
  0.000000  0.000000	;
  0.005000  0.009730	;
  0.007500  0.011730	;
  0.012500  0.014920	;
  0.025000  0.020780	;
  0.050000  0.028950	;
  0.075000  0.035040	;
  0.100000  0.039940	;
  0.150000  0.047470	;
  0.200000  0.052870	;
  0.250000  0.056640	;
  0.300000  0.059010	;
  0.350000  0.059950	;
  0.400000  0.059570	;
  0.450000  0.057920	;
  0.500000  0.055170	;
  0.550000  0.051480	;
  0.600000  0.047000	;
  0.650000  0.041860	;
  0.700000  0.036210	;
  0.750000  0.030260	;
  0.800000  0.024260	;
  0.850000  0.018260	;
  0.900000  0.012250	;
  0.950000  0.006250	;
  1.000000  0.000250]   ;	

vert.lSection = [
  0.000000  0.000000	;
  0.005000 -0.009730	;
  0.007500 -0.011730	;
  0.012500 -0.014920	;
  0.025000 -0.020780	;
  0.050000 -0.028950	;
  0.075000 -0.035040	;
  0.100000 -0.039940	;
  0.150000 -0.047470	;
  0.200000 -0.052870	;
  0.250000 -0.056640	;
  0.300000 -0.059010	;
  0.350000 -0.059950	;
  0.400000 -0.059570	;
  0.450000 -0.057920	;
  0.500000 -0.055170	;
  0.550000 -0.051480	;
  0.600000 -0.047000	;
  0.650000 -0.041860	;
  0.700000 -0.036210	;
  0.750000 -0.030260	;
  0.800000 -0.024260	;
  0.850000 -0.018260	;
  0.900000 -0.012250	;
  0.950000 -0.006250	;
  1.000000 -0.000250]   ;

vert.iBox = [0.35 0.68 0.68 0.35; max(vert.uSection(:, 2)) max(vert.uSection(:, 2)) min(vert.lSection(:, 2)) min(vert.lSection(:,2))];
vert.boxT = 0.12;
vert.boxL = 0.68 - 0.35;

vSec = figure;
hold on;
plot(vert.uSection(:,1), vert.uSection(:,2), vert.lSection(:,1), vert.lSection(:,2));
plot(vert.iBox(1, :), vert.iBox(2, :), 'k'); 
axis equal;

%% Model thickness and chord distributions and leading edge location.
%  (Root to tip):
hori.chord = linspace(hori.cR, hori.cT, N);
vert.chord = linspace(vert.cR, vert.cT, N);

hori.maxT = hori.chord .* hori.tcr;
vert.maxT = vert.chord .* vert.tcr;

hori.LE = linspace(0, 3.553, N);
vert.LE = linspace(0, 5.423, N);

hori.span = linspace(0, hori.b / 2 , N);
vert.span = linspace(0, vert.b     , N);

hori.USectionX = hori.uSection(:,1) .* ones([length(hori.uSection), N]);
hori.USectionY = hori.uSection(:,2) .* ones([length(hori.uSection), N]);

hori.LSectionX = hori.lSection(:,1) .* ones([length(hori.lSection), N]);
hori.LSectionY = hori.lSection(:,2) .* ones([length(hori.lSection), N]);

hori.boxH = hori.boxT * hori.chord;
hori.boxC = hori.boxL * hori.chord;

vert.USectionX = vert.uSection(:,1) .* ones([length(vert.uSection), N]);
vert.USectionY = vert.uSection(:,2) .* ones([length(vert.uSection), N]);

vert.LSectionX = vert.lSection(:,1) .* ones([length(vert.lSection), N]);
vert.LSectionY = vert.lSection(:,2) .* ones([length(vert.lSection), N]);

vert.boxH = vert.boxT * vert.chord;
vert.boxC = vert.boxL * vert.chord;

for i = 1:N
    hori.USectionX(:, i) = hori.USectionX(:, i) * hori.chord(i) + hori.LE(i);
    hori.USectionY(:, i) = hori.USectionY(:, i) * hori.chord(i);
    
    hori.LSectionX(:, i) = hori.LSectionX(:, i) * hori.chord(i) + hori.LE(i);
    hori.LSectionY(:, i) = hori.LSectionY(:, i) * hori.chord(i);
    
    vert.USectionX(:, i) = vert.USectionX(:, i) * vert.chord(i) + vert.LE(i);
    vert.USectionY(:, i) = vert.USectionY(:, i) * vert.chord(i);
    
    vert.LSectionX(:, i) = vert.LSectionX(:, i) * vert.chord(i) + vert.LE(i);
    vert.LSectionY(:, i) = vert.LSectionY(:, i) * vert.chord(i);
end

h3D = figure;
hold on;
for i = 1:10:N
    plot3(hori.span(i) * ones(length(hori.USectionY(:, i))), hori.USectionX(:, i), hori.USectionY(:, i), 'b');
    plot3(hori.span(i) * ones(length(hori.LSectionY(:, i))), hori.LSectionX(:, i), hori.LSectionY(:, i), 'b');
end

if i ~= N
    plot3(hori.span(N) * ones(length(hori.USectionY(:, N))), hori.USectionX(:, N), hori.USectionY(:, N), 'b');
    plot3(hori.span(N) * ones(length(hori.LSectionY(:, N))), hori.LSectionX(:, N), hori.LSectionY(:, N), 'b');
end

for i = 1:5:length(hori.lSection)
    plot3(hori.span, hori.LSectionX(i, :), hori.LSectionY(i, :), 'r--');
    plot3(hori.span, hori.USectionX(i, :), hori.USectionY(i, :), 'g');
end

view([22.5, 22.5])
axis equal;

v3D = figure;
hold on;
for i = 1:10:N
    plot3(vert.USectionY(:, i), vert.USectionX(:, i), vert.span(i) * ones(length(vert.USectionY(:, i))), 'b');
    plot3(vert.LSectionY(:, i), vert.LSectionX(:, i), vert.span(i) * ones(length(vert.LSectionY(:, i))), 'b');
end

if i ~= N
    plot3(vert.USectionY(:, N), vert.USectionX(:, N), vert.span(N) * ones(length(vert.USectionY(:, N))), 'b');
    plot3(vert.LSectionY(:, N), vert.LSectionX(:, N), vert.span(N) * ones(length(vert.LSectionY(:, N))), 'b');
end

for i = 1:5:length(vert.lSection)
    plot3(vert.LSectionY(i, :), vert.LSectionX(i, :), vert.span, 'r--');
    plot3(vert.USectionY(i, :), vert.USectionX(i, :), vert.span, 'g');
end

view([22.5, 22.5])
axis equal;

%% Material parameters
E_skin = 71.7 * 10^9; % final E for Al-7075 T6 Aluminum
E_spar = 76.5 * 10^9; % final E for Al-2050 T6 Aluminum
E_rib = 71.7 * 10^9; % final E for Al-2050 T6 Aluminum

%% Spar Sizing:
%Horizontal Stabiliser:

%Spar Webs:
for i = 1:3
    % Shear flow in front and rear spar webs:
    loadParams(i).shearFlowF = loadParams(i).Shear' ./ (2 * hori.boxH) + loadParams(i).Torque' ./(2 * hori.boxH .* hori.boxC);
    loadParams(i).shearFlowR = loadParams(i).Shear' ./ (2 * hori.boxH) - loadParams(i).Torque' ./(2 * hori.boxH .* hori.boxC);
    
    % Required thickness of front and rear spar webs:
    sizing(i).sparTF = (abs(loadParams(i).shearFlowF) .* hori.boxH / (8.1 * E_spar)).^(1/3) * 10e3;
    sizing(i).sparTR = (abs(loadParams(i).shearFlowR) .* hori.boxH / (8.1 * E_spar)).^(1/3) * 10e3;
    
    % Shear Force distribution in front and rear spar webs:
    loadParams(i).shearFS = loadParams(i).Shear' ./ (2 * hori.boxH .* sizing(i).sparTF);
    loadParams(i).shearRS = loadParams(i).Shear' ./ (2 * hori.boxH .* sizing(i).sparTR);
    
    % Max shear stress in front and rear spar webs:
    loadParams(i).maxShearFS = max(loadParams(i).shearFS);
    loadParams(i).maxShearRS = max(loadParams(i).shearRS);
    
    % Critical shear stress in front and rear spar webs:
    loadParams(i).critShearFS = 8.98 .* E_spar .* (sizing(i).sparTF ./ hori.boxH).^2;
    loadParams(i).critShearRS = 8.98 .* E_spar .* (sizing(i).sparTR ./ hori.boxH).^2;
end
