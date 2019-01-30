function [x, chord, distLoad, shearForce, bendingMoment] = ...
    vertical_stabilizer_load(rootChordLen, tipChordLen , wingSemiSpan, N , n, concForce)

% INPUT:
% rootChordLen - root chord length [m]
% tipChordLen - tip chord length [m]
% wingSemiSpan - half-spanof wing [m]
% N - Number of stations [-]
% n - Load factor [-]
% concForce - concentrated Force [N]


% OUTPUT:
% x - Span along wing [m]
% chord - Chord along span (array) [m]
% distLoad - distribution of lift [N/m]
% shearForce - shear force [N]
% bendingMoment - Moment [Nm]

x = linspace(0,wingSemiSpan, N); % along wing span [m]
chord = rootChordLen-((rootChordLen-tipChordLen)/wingSemiSpan)*x; %chord length [m]

% Wing load due to lift
ql = @(x)(((4*n*concForce)*sqrt(wingSemiSpan^2-x.^2))/(pi*wingSemiSpan^2)); % lift function
distLoad = ql(x);


% Calculate shear force
shearForce = zeros(1,length(x));
for i=1:length(x)-1
    shearForce(i) = integral(ql, x(i), wingSemiSpan);
end


% Calculate moment
bendingMoment = zeros(1,length(x));
for i=1:length(x)-1
    bendingMoment(i) = trapz(x(i:length(x)), shearForce(i:length(x)));
end
