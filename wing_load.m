function [x, chord, distLift, distWeight ,shearForce, bendingMoment] = ...
    wing_load(rootChordLen, tipChordLen, wingSemiSpan, takeOffWeight, N ,n)

% INPUT:
% rootChordLen - root chord length [m]
% tipChordLen - tip chord length [m]
% wingSemiSpan - half-spanof wing [m]
% takeOffWeight - Take-off weight of entire aircraft [N]
% N - Number of stations [-]
% n - Load factor [-]


% OUTPUT:
% x - Span along wing [m]
% chord - Chord along span (array) [m]
% distLift - distribution of lift [N/m]
% distWeight - distribution of structural weight ` [N/m]
% shearForce - shear force [N]
% bendingMoment - Moment [Nm]



x = linspace(0, wingSemiSpan, N); % along wing span [m]
chord = rootChordLen - ((rootChordLen - tipChordLen)/wingSemiSpan) * x; %chord length [m]

wingStructureWeight = 0.11 * takeOffWeight; % weight of wing structure (11%) [N]


% Wing load due to lift
%q_l= (2*W_TO*n)/(pi*L^2)*sqrt(L^2-x.^2);
ql = @(x)(2*takeOffWeight*n)/(pi*wingSemiSpan^2)*sqrt(wingSemiSpan^2-x.^2); % lift function
distLift= ql(x);


% Wing load due to wing structure
qw= @(x)(-(wingStructureWeight*n)/(wingSemiSpan*(rootChordLen+tipChordLen))*(rootChordLen- (rootChordLen-tipChordLen)/wingSemiSpan.*x )); % weight structure function
distWeight= qw(x);


% Calculate shear force
shearForce = zeros(1,length(x));
for i=1:length(x)-1
shearForce(i) = integral(ql, x(i), wingSemiSpan )+ integral(qw, x(i), wingSemiSpan );
end


%Calculate moment
bendingMoment = zeros(1,length(x));
for i=1:length(x)-1
bendingMoment(i) = trapz(x(i:length(x)), shearForce(i:length(x)));
end
