%% General analysis of load case 1 on the wing group
clear
clc
close all

%% Enter geometry of wing

% rootChordLen - root chord length [m]
% tipChordLen - tip chord length [m]
% wingSemiSpan - wing semi-span [m]
% takeOffWeight - Take-off weight of entire aircraft [N]
% N - Number of stations [-]
% n - Load factor [-]
% fuelTankLen - Length of fuel tank wrt half span [-]
[vals, names] = xlsread('geometryVariables.xlsx', 'Data', 'B:C');
geoParams = containers.Map(names(2:end,1), vals);
rootChordLen_w = geoParams('rootChordLen_w');
tipChordLen_w = geoParams('tipChordLen_w');
wingSemiSpan = geoParams('wingSemiSpan');
takeOffWeight = geoParams('takeOffWeight');
N = 100; % according to lecturer 
n = 1.5*2.5;
fuelTankLen = geoParams('fuelTankLen'); % pct of semi-span from root

%% Enter operating conditions and pitch aerodyamic moment
CM0_w = geoParams('CM0_w'); % pitch aerodynamic moment from airofoil data
cruiseVelocity = 232.78;
rho = 1.1; % at altitude 37000 ft

%% Get values from wing_load and fuel_load

% self weight both lift
[x, chord, distLift, distWeightWing, shearForceWing, bendingMomentWing] = wing_load(rootChordLen_w, tipChordLen_w , wingSemiSpan , takeOffWeight*9.81 , N , n);

% Fuel
[distWeightFuel, shearForceFuel, bendingMomentFuel] = fuel_load(rootChordLen_w, tipChordLen_w , wingSemiSpan , takeOffWeight*9.81 , N , n, fuelTankLen);


flexAxis = 0.45*chord; % Flex axis is 25 to 65 of chord
ac = 0.25*chord; % aerodynamic center
a = flexAxis - ac; % distence between the flex axis and AC
cgWing = 0.55*chord;
b = cgWing - flexAxis; % difference between cg with respect to chord of wing and FA

% Engine: (assumed point load)

% first iteration assumed point load
% 2 entries so it is the same size for plotting
engineLoc_1 = geoParams('engineLoc_1'); % [m]
engineLoc_2 = geoParams('engineLoc_2'); % [m]
ew = geoParams('engineWeight');
engineWeight_1=-ew*9.81; % in [N]
engineWeight_2=-ew*9.81; % in [N]

% load distributions
distLoad = distLift+distWeightWing+distWeightFuel;
panelDist = wingSemiSpan/N;

loc_1 = find(x>engineLoc_1-panelDist,1); % position of engine 1
loc_2 = find(x>engineLoc_2-panelDist,1); % position of engine 2
% add concentrated loads from engines
distLoad(loc_1) = distLoad(loc_1)+engineWeight_1; % total load with engine 1
distLoad(loc_2) = distLoad(loc_2)+engineWeight_2; % total load with engine 2

figure;
hold on
plot(x,distLift,'b') % lift
plot(x,distWeightWing,'r') % self weight
plot(x,distWeightFuel,'g') % fuel weight

plot(x,distLoad,'k') % total load
hold off

totalMoments = bendingMomentWing+bendingMomentFuel;
for i=1:loc_1
    totalMoments(i) = totalMoments(i) - (engineWeight_1*(engineLoc_1-x(i)));
end
for i=1:loc_2
    totalMoments(i) = totalMoments(i) - (engineWeight_2*(engineLoc_2-x(i)));
end
% moment distribution
figure;
hold on
plot(x,bendingMomentWing,'r') % wing self weight and lift 
plot(x,bendingMomentFuel,'b') % fuel
plot(x,totalMoments,'k') % all
hold off

% shear force:
totalShearForce = shearForceWing+shearForceFuel;
% add concentrated loads from engines
totalShearForce(1:loc_1) = totalShearForce(1:loc_1) + engineWeight_1;
totalShearForce(1:loc_2) = totalShearForce(1:loc_2) + engineWeight_2;

figure;
hold on
plot(x,shearForceWing,'r')
plot(x,shearForceFuel,'b')
plot(x,totalShearForce,'k')
hold off

M_0=0.5*rho*cruiseVelocity^2.*chord.^2*CM0_w;
T=distLift.*a*wingSemiSpan+(distLoad-distLift).*b*wingSemiSpan-M_0;

figure;
hold on
plot(x,T) % torque along span
hold off
%% nic's wing structure sizing

% fixed wing parameters determined during AVD projet
b = 37.40; % span
AR = 7.5;
taperRatio = 0.20; 
lead_sweep = 27; % degrees
c_bar = 5.57; % mean aerodynamic chord
twist = -5; % degrees
dihedral = -2; % degrees
t_over_c = 0.12; % thickness to chord ratio
c_Root = 8.10;
c_tip = 1.62;
S_ref = 176.71; % m^2

%% material parameters
E_skin = 76 * 10^9; % look up proper value Al-2090
E_spar = 76.5 * 10^9; % look up proper value Al-7249
E_rib = 76 * 10^9;

%% skin-stringer panel geometry - standard z stringer
% close all
% % ratios to decide
% d_over_h = 0.3; % fixed due to data avaiable
% ts_over_t = [0.40; 0.50; 0.60; 0.70; 0.80; 0.90; 1.0; 1.2; 1.4];
% 
% % sigma_crit has to be greater than zero or the skin thickness without the
% % stringer will be smaller than with stringer i.e. inefficient
% 
% % trial 1
% A_s_over_bt = 0.5;
% sigma_ratio = [0.3; 0.7; 1.05; 1.15; 1.25; 1.3; 1.38; 1.48; 1.55];
% % 
% % % trial 2
% % A_s_over_bt = 0.8;
% % sigma_ratio = [0.175; 0.3; 0.6; 0.975; 1.15; 1.25; 1.35; 1.475; 1.575];
% 
% 
% L = 0.5; % rib spacing guess
% 
% % cover load per unit length N
% h_wingBox = 0.10 * chord; % height of wing box assumed 10% chord
% width_wingBox = (0.65-0.15) * chord; % width of wing box
% N_cover = (totalMoments ./ (h_wingBox .* width_wingBox))';
% 
% % optimising ts_over_t ratio
% 
% for b1 = 0.01:0.01:0.1
%     figure
%     hold on
%     for i=1:9
%         % WITH STRINGERS
%         t1(:,i) = ((N_cover.*b1^2) ./ (sigma_ratio(i) .* 3.62 * E_skin)).^(1/3);
%         A_s_stringer(:,i) = b1 .* ts_over_t(i) * t1(:,i);
%         %h_stringer(:,i) = A_s_stringer(:,i)./(1.6 .* ts_over_t(i) .* t1(:,i));
%         sigma_stringer(:,i) = N_cover./t1(:,i);
%         plot(x,t1(:,i),'b');
%         
%         % WITHOUT STRINGERS
%         t1_nada(:,i) = ((N_cover.*b1^2) ./ (3.62 * E_skin)).^(1/3);
%         %A_s_nada(:,i) = b1 .* t1_nada(:,i);
%         sigma_nada(:,i) = N_cover./t1_nada(:,i);        
%         plot(x,t1_nada(:,i),'r');
%         
%         % effective calculations with  stringers used
%         effective_T(i,:) = (t1(:,i) + A_s_stringer(:,i)./b1);
%         sigma_z(:,i) = N_cover./effective_T(i,:)';
%         F(:,i)  = sigma_stringer(:,i) .* sqrt(L ./ (N_cover.*E_skin)); % farrar efficiency, maximise at root
%         %plot(x,effective_T,'m')
%         
%     end
%         % failure thickness
%         t_failure  = N_cover./(520*10^6);
%         plot(x,t_failure,'g')
% %     % plotting yield stress
% %     plot([x(1),x(100)],[520*10^6,520*10^6],'k')
% %     plot(x,t1(:,i),'o')
%     hold off
% end

% % optimising for stringer pitch
% counter=0;
% figure
% hold on
% for b1=0.01:0.01:0.1
%     counter=counter+1;
%     
%     % skin thickness
%     t1(counter,:) = ((N_cover.*b1^2)/(sigma_ratio*3.62*E_skin)).^(1/3);
%     t1_max(counter,1) = max(t1(counter,:));
%     
%     % total cross-sectional area of stringers
%     A_s(counter,:) = b1.*t1(counter,:);
%     total_A(counter,1) = sum(A_s(counter,:));
%     
%     % stringer height
%     h_stringer(counter,:) = A_s(counter,:)./(1.6*t1(counter,:)); % assumed t/t_s = 1
%     h_max(counter,1) = max(h_stringer(counter,:));
%     
%     % cover stress
%     sigma(counter,:) = N_cover./t1(counter,:);
%     
%     % effective thickness
%     T(counter,:) = t1(counter,:) + A_s(counter,:)./b1;
%     
%     % longitudinal stress
%     sigma_z(counter,:) = N_cover./T(counter,:);
%     plot(x,sigma_z(counter,:))
%     
%     % checking ratios
%     ratio1(counter,:) = A_s(counter,:)/(b1*t1(counter,:));
% end
% plot([x(1),x(100)],[520*10^6,520*10^6],'k')
% hold off


%% spar sizing
close all
h_wingBox = 0.10 * chord; % height of wing box assumed 10% chord
width_wingBox = (0.65-0.15) * chord; % width of wing box

shear_flow_front = totalShearForce./(2 .* h_wingBox) + T./(2 .* h_wingBox .* width_wingBox);
shear_flow_rear = totalShearForce./(2 .* h_wingBox) - T./(2 .* h_wingBox .* width_wingBox);

t_front = (abs(shear_flow_front) .* h_wingBox / (8.1 * E_spar)).^(1/3);
t_rear = (abs(shear_flow_rear) .* h_wingBox / (8.1 * E_spar)).^(1/3);

figure
hold on
plot(x,abs(shear_flow_front),'b')
plot(x,abs(shear_flow_rear),'r')
title('shear flow')
hold off

figure
hold on
plot(x,t_front,'b')
plot(x,t_rear,'r')
title('spar web thickness')
legend('front','rear')
hold off

% incorporating manufacturing constraints into design
for i = 1:100
    
    % manufacturing constraint of 1mm minimum thickness
    % convert to mm
    t_front(i) = t_front(i) * 10^3;
    t_rear(i) = t_rear(i) * 10^3;
    
    % front spar
    if t_front(i)<1
        t_front(i) = 1;
    else
        t_front(i) = t_front(i);
    end
    
    % rear spar
    if t_rear(i)<1
        t_rear(i) = 1;
    else
        t_rear(i) = t_rear(i);
    end
    
    % rounding up to nearest half milimeter
    t_front_min(i) = ceil(t_front(i));
    t_rear_min(i) = ceil(t_rear(i));
        
    % front spar
    if t_front_min(i)-t_front(i)>0.5
        
        t_front_min_1(i) = (t_front_min(i) - 0.5) * 10^-3;
        
    elseif t_front_min(i)-t_front(i)<0.5
        
        t_front_min_1(i) = t_front_min(i) * 10^-3;
    end
    
    % rear spar
    if t_rear_min(i)-t_rear(i)>0.5
        
        t_rear_min_1(i) = (t_rear_min(i) - 0.5) * 10^-3; % back to mm
        
    elseif t_rear_min(i)-t_rear(i)<0.5
        
        t_rear_min_1(i) = t_rear_min(i) * 10^-3; % back to mm
    end
end

figure
hold on
plot(x,t_front_min_1)
plot(x,t_rear_min_1)
hold off

% optimising spar thickness to give smooth distribution to have smooth
% decreasing variation in spar web thickness along half span

t_front_min_2 = t_front_min_1;
t_rear_min_2 =t_rear_min_1;

% front spar
for j=1:99
    if t_front_min_2(j+1)>t_front_min_2(j)
       t_front_min_2(j) = t_front_min_2(j+1);
    end
end

for j = 1:80
   % rear spar
   if t_rear_min_2(j) < 6.5*10^-3
      t_rear_min_2(j) = 6.5*10^-3;
   end    
end

figure
hold on
plot(x,t_front_min_1,'-.','LineWidth',1.5)
plot(x,t_rear_min_1,'-.','LineWidth',1.5)
plot(x,t_front_min_2,'--','LineWidth',1.5)
plot(x,t_rear_min_2,'--','LineWidth',1.5)
legend('Manufacturing constrained front spar','Manufacturing constrained rear spar','Optimised front spar','Optimised rear spar')
hold off
%%
% calculating shear stress in spar webs
shear_front_spar = totalShearForce ./ (2 .* h_wingBox .* t_front);
shear_rear_spar = totalShearForce ./ (2 .* h_wingBox .* t_rear);

% max shear stress must be below yield stress
shear_front_max = max(shear_front_spar);
shear_rear_max = max(shear_rear_spar);

% critical stress
shear_crit_front = 8.98 .* E_spar .* (t_front./h_wingBox).^2;
shear_crit_rear = 8.98 .* E_spar .* (t_rear./h_wingBox).^2;

max(shear_crit_front)
max(shear_crit_rear)

% yield stress at any point

% spar cap area THIS IS WRONG
figure
yield_spar = 550 * 10^6; % UPDATE THIS VALUE
A = totalMoments ./ ((yield_spar)*h_wingBox);
plot(x,A)

% plotting aerofoil profile

% upper surface
fid = fopen('NACA63(1)-412UPPER.txt');
vals1 = fscanf(fid, '%f %f ',[2 Inf]);
coordinates1 = transpose(vals1);
fclose(fid);

% lower surface
fid = fopen('NACA63(1)-412LOWER.txt');
vals2 = fscanf(fid, '%f %f ',[2 Inf]);
coordinates2 = transpose(vals2);
fclose(fid);

figure
hold on
leg1 = plot(coordinates1(:,1),coordinates1(:,2),'k');
plot(coordinates2(:,1),coordinates2(:,2),'k');
ylim([-0.35,0.35]);
xlabel('x/c','FontSize',18)
ylabel('y/c','FontSize',18)

% plotting wing box, actual shape of a rectangle
% front spar is at 15% chord and rear spar is at 65% chord
leg2 = plot([0.15,0.15],[0.061795,-0.034220],'r-','LineWidth',1.5);
plot([0.65,0.65],[-0.0179087,0.0590905],'r-','LineWidth',1.5);
plot(coordinates1([9:19],1),coordinates1([9:19],2),'r-','LineWidth',1.5);
plot(coordinates2([9:19],1),coordinates2([9:19],2),'r-','LineWidth',1.5);

% plotting idealised rectangular wing box
leg3 = plot([0.15,0.65,0.65,0.15,0.15],[-0.035,-0.035,0.065,0.065,-0.035],'bo-.','LineWidth',1.5)
legend([leg1,leg2,leg3],'NACA profile','Real wing box section','Idealised rectangular wing box section')
hold off


%% integrally stiffened panels ideal
close all
% cover load per unit length N
h_wingBox = 0.10 * chord; % height of wing box assumed 10% chord
width_wingBox = (0.65-0.15) * chord; % width of wing box
N_cover = (totalMoments ./ (h_wingBox .* width_wingBox))';

% integrally stiffened panels - realistic approach
close all
% trial one, fixing values of b, R_t and R_b
R_b_1 = 0.65; R_t_1 = 2.25; sigma_ratio_1 = 1.6;

F_new_1 = 1.314 * ((R_b_1^3 * R_t_1 * (4 + R_b_1 * R_t_1))^0.25)/(1 + R_b_1 * R_t_1) * sigma_ratio_1^0.25;

%rib sizing
% rib thickness
t_rib = ((abs(distLoad) .* (h_wingBox).^2)./(3.62 .* chord .* E_rib)).^(1/3);
figure
plot(x,t_rib)

% rib spacing % farrar efficieny not constant!
L = (4 .* F_new_1^2 .* width_wingBox.^2 .* t_rib.^2 .* E_rib ./ N_cover').^(1/3);
figure
plot(x,L)
L_min = min(L)

% stringer pitch optimised for the root
b_new_1 = 1.103 * ((sqrt(F_new_1)*(1 + R_b_1 * R_t_1)/((R_b_1^3 * R_t_1 * (4 + R_b_1 * R_t_1))^0.5))) .* (N_cover(1)' .* L_min^3 ./ E_skin).^0.25;

% WITHOUT STRINGERS
t1_nada = ((N_cover.*L_min^2) ./ (3.62 * E_skin)).^(1/3);
A_s_nada = L_min .* t1_nada;
sigma_nada = N_cover./t1_nada;

% WITH STRINGERS
t_new_1 = ((N_cover .* b_new_1^2) ./ (sigma_ratio_1 * 3.62 * E_skin * (1 + R_b_1 * R_t_1))).^(1/3);
A_s = t_new_1 * b_new_1;
T_effective = (t_new_1 + A_s./b_new_1);
sigma_new_1 = N_cover./(t_new_1 * (1+R_t_1 * R_b_1));
sigma_z = N_cover./T_effective;

% yield skin thickness
t_failure  = N_cover ./ (520*10^6);

% plotting thickness distriutions
figure
hold on
stairs(x,t1_nada,'r');
stairs(x,t_new_1,'b');
stairs(x,t_failure,'g');
hold off

figure
hold on
plot(x,sigma_nada);
plot(x,sigma_new_1);
plot(x,sigma_z);
plot([x(1), x(100)], [434.5*10^6 434.5*10^6])
legend('no stringers','stringers included','sigma_z','sigma yield')
title('stress variation along half span')
hold off
% SANITY CHECK
% t_new_1 = ((900000 .* b_new_1^2) ./ (1.75 * 3.62 * (70 * 10^9) * (1 + R_b_1 * R_t_1))).^(1/3);

% buckling stress
sigma_crit = N_cover ./ (t_new_1 * (1+R_t_1 * R_b_1));

% SANITY CHECK
% sigma_crit = 900000 ./ (t_new_1 * (1+R_t_1 * R_b_1));

% 
% figure
% hold on
% plot(x,sigma_crit)
% % plotting yield stress
% plot([x(1),x(100)],[520*10^6,520*10^6],'k')
% hold off

% figure 
% hold on
% plot(x,t_new_1)
% % failure thickness
% t_failure  = N_cover./(520*10^6);
% plot(x,t_failure,'g')
% hold off
%% rib sizing
% close all

h_wingBox = 0.10 * chord; % height of wing box assumed 10% chord
width_wingBox = (0.65-0.15) * chord; % width of wing box

% rib thickness
t_rib = ((abs(distLoad) .* (h_wingBox).^2)./(3.62 .* chord .* E_rib)).^(1/3);
figure
plot(x,t_rib)

% rib spacing
L = (4 .* 0.81^2 .* width_wingBox.^2 .* t_rib.^2 .* E_rib ./ N_cover').^(1/3);
figure
plot(x,L)
L_min = min(L)


% stress in ribs
sigma_ribs = distLoad ./ (t_rib .* chord);

% crushing force
