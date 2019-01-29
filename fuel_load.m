function [distWeightFuel, shearForceFuel, bendingMomentFuel] = ...
    fuel_load(rootChordLen, tipChordLen , wingSemiSpan , takeOffWeight , N , n, fuelTankLen)

% INPUT:
% rootChordLen - root chord length [m]
% tipChordLen - tip chord length [m]
% wingSemiSpan - half-spanof wing [m]
% takeOffWeight - Take-off weight of entire aircraft [N]
% N - Number of stations [-]
% n - Load factor [-]
% fuelTankLen - Length of fuel tank wrt half span [-]
% OUTPUT:
% distWeightFuel - Distribution of weight of fuel [N/m]
% shearForceFuel - Shear force due to fuel [N]
% bendingMomentFuel - Moment due to fuel [Nm]

x = linspace(0, wingSemiSpan, N);
weightFuel= 0.83*0.35*takeOffWeight; % weight of fuel (35%)*0.83 stored in wing
% Length of fuel tank from root (currently in terms of % of span)
Lf = fuelTankLen*wingSemiSpan;
C0f = 0.45*rootChordLen;
Ctf = 0.45*(rootChordLen- (rootChordLen-tipChordLen)*Lf/wingSemiSpan);
distWeightFuel = zeros(1,length(x));
distWeightFuel(x<Lf) = -(weightFuel*n)/(Lf*(C0f+Ctf))*(C0f- (C0f-Ctf)/Lf.*x(x<Lf) );
qf = @(x)(-(weightFuel*n)/(Lf*(C0f+Ctf))*(C0f- (C0f-Ctf)/Lf.*x));
% Shear Force
shearForceFuel= zeros(1,length(x));
for i=1:length(x)-1
    if x(i)< Lf
        shearForceFuel(i) = integral(qf,x(i), max(x(x<=Lf)));
    end
end
% Moment
bendingMomentFuel = zeros(1,length(x));
for i=1:length(x)-1
    bendingMomentFuel(i) = -trapz(x(i:length(x)), shearForceFuel(i:length(x)));
end