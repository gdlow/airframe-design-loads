function [x, chord, distLift, distWeight ,shearForce, bendingMoment] = ...
    horizontal_stabilizer_load(rootChordLen, tipChordLen, wingSemiSpan,...
    horizontalStabilizerLift,horizontalStabilizerWeight , N ,n)

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

% Wing load due to lift
%q_l= (2*W_TO*n)/(pi*L^2)*sqrt(L^2-x.^2);
ql = @(x)(-4*horizontalStabilizerLift*n)*sqrt(wingSemiSpan^2-x.^2)/(pi*wingSemiSpan^2); % lift function
distLift= ql(x);


% Wing load due to wing structure
qw= @(x)(-(horizontalStabilizerWeight*n)/(wingSemiSpan*(rootChordLen+tipChordLen))*(rootChordLen- (rootChordLen-tipChordLen)/wingSemiSpan.*x ));
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
