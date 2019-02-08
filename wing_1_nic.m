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
rho = 0.350085776097978; % at altitude 37000 ft

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
T=distLift.*a-(distLoad-distLift).*b-M_0;

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
stairs(x,t_front,'b')
stairs(x,t_rear,'r')
legend('Front spar','Rear spar')
xlabel('Span-wise distance from root (m)')
ylabel('Spar web (m)')
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

% rear spar
for j=1:99
    if t_rear_min_2(j+1)>t_rear_min_2(j)
        t_rear_min_2(j) = t_rear_min_2(j+1);
    end
end
for j=1:98
    if t_rear_min_2(j+2)>t_rear_min_2(j)
        t_rear_min_2(j) = t_rear_min_2(j+2);
    end
end

% for j = 1:80
%    % rear spar
%    if t_rear_min_2(j) < 6.5*10^-3
%       t_rear_min_2(j) = 6.5*10^-3;
%    end    
% end

figure
hold on
plot(x,t_front_min_1,':','LineWidth',2.5)
plot(x,t_rear_min_1,':','LineWidth',2.5)
plot(x,t_front_min_2,'--','LineWidth',2)
plot(x,t_rear_min_2,'--','LineWidth',2)
xlabel('Span-wise distance from root (m)')
ylabel('Spar web (m)')
legend('Manufacturing constrained front spar','Manufacturing constrained rear spar','Optimised front spar','Optimised rear spar')
hold off

% calculating shear stress in spar webs
shear_front_spar = totalShearForce ./ (2 .* h_wingBox .* t_front_min_2);
shear_rear_spar = totalShearForce ./ (2 .* h_wingBox .* t_rear_min_2);

% max shear stress must be below yield stress
shear_front_max = max(shear_front_spar);
shear_rear_max = max(shear_rear_spar);

% critical stress
shear_crit_front = 8.98 .* E_spar .* (t_front_min_2./h_wingBox).^2;
shear_crit_rear = 8.98 .* E_spar .* (t_rear_min_2./h_wingBox).^2;

max(shear_crit_front)
max(shear_crit_rear)


% spar cap area 
figure
yield_spar = 529.5 * 10^6; % UPDATE THIS VALUE
A = totalMoments ./ ((yield_spar).*h_wingBox);
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

%% Digitalising Catchpole Diagram
close all

% R_t = 1.0
fid = fopen('R_t_1_0.txt');
vals_1 = fscanf(fid, '%f %f ',[2 Inf]);
catchpole_1 = transpose(vals_1);
fclose(fid);

% R_t = 1.5
fid = fopen('R_t_1_5.txt');
vals_2 = fscanf(fid, '%f %f ',[2 Inf]);
catchpole_2 = transpose(vals_2);
fclose(fid);

% R_t = 2.0
fid = fopen('R_t_2_0.txt');
vals_3 = fscanf(fid, '%f %f ',[2 Inf]);
catchpole_3 = transpose(vals_3);
fclose(fid);

% R_t = 2.5
fid = fopen('R_t_2_5.txt');
vals_4 = fscanf(fid, '%f %f ',[2 Inf]);
catchpole_4 = transpose(vals_4);
fclose(fid);

% R_t = 5.0
fid = fopen('R_t_5_0.txt');
vals_5 = fscanf(fid, '%f %f ',[2 Inf]);
catchpole_5 = transpose(vals_5);
fclose(fid);

figure
hold on
plot(catchpole_1(:,1),catchpole_1(:,2),'r-')
plot(catchpole_2(:,1),catchpole_2(:,2),'b-')
plot(catchpole_3(:,1),catchpole_3(:,2),'r-.')
plot(catchpole_4(:,1),catchpole_4(:,2),'b-.')
plot(catchpole_5(:,1),catchpole_5(:,2),'r--')
xlabel('$\frac{h}{b}$','Interpreter','latex','FontSize',20)
ylabel('$\frac{\sigma_{cr}}{\sigma_0}$','Interpreter','latex','FontSize',20)
legend('$\frac{t_s}{t} = 1.0$','$\frac{t_s}{t} = 1.5$','$\frac{t_s}{t} = 2.0$','$\frac{t_s}{t} = 2.5$','$\frac{t_s}{t} = 5.0$','Interpreter','latex','FontSize',12)
hold off

% possible parameter values
R_t = [1.0; 1.5; 2.0; 2.5; 5.0];
R_b = [0.2; 0.4; 0.6; 0.8; 1.0; 1.2];
% sigma ratio at R_b values in previous line for different R_t curves
sigma_ratio(1,[1:6]) = [catchpole_1(9,2), catchpole_1(19,2), catchpole_1(31,2), catchpole_1(42,2), catchpole_1(53,2), catchpole_1(63,2)]; % R_t = 1.0
sigma_ratio(2,[1:6]) = [catchpole_2(7,2), catchpole_2(10,2), catchpole_2(20,2), catchpole_2(33,2), catchpole_2(44,2), catchpole_2(52,2)]; % R_t = 1.5
sigma_ratio(3,[1:6]) = [catchpole_3(16,2), catchpole_3(26,2), catchpole_3(35,2), catchpole_3(51,2), catchpole_3(63,2), catchpole_3(73,2)]; % R_t = 2.0
sigma_ratio(4,[1:6]) = [catchpole_4(16,2), catchpole_4(27,2), catchpole_4(37,2), catchpole_4(49,2), catchpole_4(65,2), catchpole_4(76,2)]; % R_t = 2.5
sigma_ratio(5,[1:6]) = [catchpole_5(16,2), catchpole_5(27,2), catchpole_5(37,2), catchpole_5(48,2), catchpole_5(63,2), catchpole_5(73,2)]; % R_t = 5.0


%% integrally stiffened panels ideal - LOOK HERE
close all
clc
% cover load per unit length N
h_wingBox = 0.10 * chord; % height of wing box assumed 10% chord
width_wingBox = (0.65-0.15) * chord; % width of wing box
N_cover = (totalMoments ./ (h_wingBox .* width_wingBox))';

% for l = 1: 100 % loop to vary rib spacing L
%     for k = 1:100 % loop to vary stringer pitch b
%         for i = 1:6 % loop iterating R_b
%             for j = 1:5 % loop iterating R_t
%                 L_range = linspace(0.1,2,100);
%                 L = L_range(l);
%                 b_range = linspace(0.01,0.4,100);
%                 b(i,j,k) = b_range(k);
%                 R_b1 = R_b(i);
%                 R_t1 = R_t(j);
%                 sigma_ratio1 = sigma_ratio(j,i);
%                 t_with(i,j,k) = ((N_cover(1) * b(i,j,k)^2) / (sigma_ratio1 * 3.62 * E_skin * (1 + R_b1 * R_t1))).^(1/3);
%                 F(i,j,k) = 1.314 * (((R_b1^3 * R_t1 * (4 + R_b1 * R_t1))^0.25) / (1 + R_b1 * R_t1)) * sigma_ratio1^0.25;
% %                 t_s(i,j,k) = R_t * t_with(i,j,k)
% %                 A_s(i,j,k) = b(i,j,k) * t(i,j,k) 
%             end
%         end
%     end
% end

% Other code attempts found here aka code graveyard 
% L = 0.1; % intial guess, needs to be optimised as well

% % optimisation time
% for i = 1:5
%     R_t1 = R_t(i);
%     for j = 1:6
%         R_b1 = R_b(j);
%         F(i,j) = 1.314 * ((R_b1^3 * R_t1 * (4 + R_b1 * R_t1))^0.25) / (1 + R_b1 * R_t1) * sigma_ratio(i,j)^0.25;
%         % calculating ideal pitch at the root and setting this as constant
%         %for the rest of the wing. 
%         b(i,j) = 1.103 * ((sqrt(F(i,j))*(1 + R_b1 * R_t1)/((R_b1^3 * R_t1 * (4 + R_b1 * R_t1))^0.5))) * (N_cover(1)' * L^3 / E_skin)^0.25;
%         
%         % skin thickness and stress at root WITHOUT stringers
%         t1_nada(i,j) = ((N_cover(1) * L^2) / (3.62 * E_skin))^(1/3);
%         sigma_nada(i,j) = N_cover(1)/t1_nada(i,j);
%         
%         % skin thickness and stress at root WITH stringers
%         t1_with(i,j) = ((N_cover(1) * b(i,j)^2) / (sigma_ratio(i,j) * 3.62 * E_skin * (1 + R_b1 * R_t1))).^(1/3);
%         % t1_with(i,j) = (1 / (F(i,j) * (1 + R_t1 * R_b1))) * (N_cover(1) * L / E_skin)^(1/2);
%         sigma_with(i,j) = N_cover(1)/(t1_with(i,j) * (1+R_t1 * R_b1));
%         h_stringer(i,j) = b(i,j) * R_b1;
%         t_stringer(i,j) = R_t1 * t1_with(i,j);
%         A_s(i,j) = t_stringer(i,j) * b(i,j);
%         t_effective(i,j) = (t1_with(i,j) + A_s(i,j) / b(i,j));
%         sigma_z(i,j) = N_cover(1) / t_effective(i,j);
%         
%         % euler buckling
% %         y(i,j) = b(i,j) - 
% %         I = (1/3) * (t_with(i,j) * y(i,j)^3 + b(i,j) * (b(i,j) - y(i,j))^3 - (b(i,j) - y(i,j)) * (b(i,j) - y(i,j) - t_with(i,j))^3)
% %         A_column(i,j) = t_stringer(i,j) * h_stringer(i,j); % area of column
% %         I(i,j) = t_stringer(i,j) * h_stringer(i,j)^3 / 12;
% %         sigma_euler(i,j) = (pi^2 * E_skin * I(i,j)) / (4 * h_stringer(i,j)^2 * A_column(i,j));
%         rho_squared(i,j) = (b(i,j)^2 * R_t1 * R_b1^3 * (4 + R_b1 * R_t1)) / (12 * (1 + R_t1 * R_b1)^2);
%         sigma_euler(i,j) = pi^2 * E_skin  * rho_squared(i,j) / L^2;
%     end
% end

% % decide on ratio selection
% R_b_trial = R_b(3); R_t_trial = R_t(3); sigma_ratio_trial = sigma_ratio(3,3);
% 
% F_trial = 1.314 * ((R_b_trial^3 * R_t_trial * (4 + R_b_trial * R_t_trial))^0.25) / (1 + R_b_trial * R_t_trial) * sigma_ratio_trial^0.25;
% b_trial = 1.103 * ((sqrt(F_trial)*(1 + R_b_trial * R_t_trial)/((R_b_trial^3 * R_t_trial * (4 + R_b_trial * R_t_trial))^0.5))) * (N_cover(1) * L^3 / E_skin)^0.25;
% t1_with_trial = ((N_cover .* b_trial.^2) ./ (sigma_ratio_trial .* 3.62 * E_skin .* (1 + R_b_trial .* R_t_trial))).^(1/3);
% t_stringer_trial = R_t_trial * t1_with_trial;
% h_stringer_trial = R_b_trial * b_trial;
% sigma_with_trial = N_cover ./ (t1_with_trial * (1 + R_t_trial * R_b_trial));
% 
% 
% figure
% hold on
% stairs(x,t1_with_trial)
% stairs(x,t_stringer_trial)
% stairs(x,ones(100,1) .* h_stringer_trial)
% hold off
% 
% figure
% hold on
% plot(x,sigma_with_trial)
% plot([x(1),x(100)],[417*10^6, 417*10^6])
% hold off

L = 0.5; % intial guess, needs to be optimised as well

% optimisation time
for k = 1:100
    
    for i = 1:5
        R_t1 = R_t(i);
        
        for j = 1:6
            R_b1 = R_b(j);
            
            for L = 0.1:0.001:5
                % Farrar Efficiency
                F(i,j) = 1.314 * ((R_b1^3 * R_t1 * (4 + R_b1 * R_t1))^0.25) / (1 + R_b1 * R_t1) * sigma_ratio(i,j)^0.25;
                
                % calculating ideal pitch at the root and setting this as constant
                %for the rest of the wing.
                b(i,j) = 1.103 * ((sqrt(F(i,j))*(1 + R_b1 * R_t1)/((R_b1^3 * R_t1 * (4 + R_b1 * R_t1))^0.5))) * (N_cover(1) * L^3 / E_skin)^0.25;
                
                % calculating how Farrar efficiency varies along Span with
                % pitch set constant to value at root
                F_variable(i,j,k) = ((b(i,j)^2 * R_b1^3 * R_t1 * (4 + R_b1 * R_t1)) / (1.103^2 * (1 + R_b1 * R_t1)^2)) * (E_skin / (N_cover(k) * L^3))^(1/2);
                
                % skin thickness and stress at root WITHOUT stringers
                t1_nada(i,j,k) = ((N_cover(k) * L^2) / (3.62 * E_skin))^(1/3);
                
                % skin thickness and stress at root WITH stringers
                t1_with(i,j,k) = ((N_cover(k) * b(i,j)^2) / (sigma_ratio(i,j) * 3.62 * E_skin * (1 + R_b1 * R_t1)))^(1/3);
                
                % initial buckling stress
                sigma_with(i,j,k) = N_cover(k)/(t1_with(i,j,k) * (1+R_t1 * R_b1));
                
                
                % Euler (flexural) buckling stress
                rho_squared(i,j,k) = (b(i,j)^2 * R_t1 * R_b1^3 * (4 + R_b1 * R_t1)) / (12 * (1 + R_t1 * R_b1)^2);
                sigma_euler(i,j,k) = pi^2 * E_skin  * rho_squared(i,j) / L^2;
                
                t_stringer(i,j,k) = R_t1 * t1_with(i,j,k);
                h_stringer(i,j,k) = R_b1 * b(i,j);
                
                t_e(i,j,k) = (t1_with(i,j,k) * (1 + R_t1 * R_b1));
                
                I(i,j,k) = 2 * t_e(i,j,k) * chord(k) * (h_wingBox(k) / 2)^2;
                F_crush(i,j,k) = (totalMoments(k)^2 * h_wingBox(k) * L * t_e(i,j,k)* chord(k)) / (2 * E_skin * I(i,j,k)^2);
                t_rib(i,j,k) = ((F_crush(i,j,k) * h_wingBox(k)^2) / (3.62 * E_skin))^(1/3);
                
                L_iteration = (4 .* F_variable(i,j,k)^2 .* h_wingBox(k).^2 .* t_rib(i,j,k).^2 .* E_rib ./ N_cover(k)).^(1/3);
                
                eps = L_iteration - L;
                
                if eps<0.001
                    L_final(i,j,k) = L_iteration;
                    
                    % calculating Farrar efficiency with stress
                    F_with_stress(i,j,k)  = sigma_with(i,j,k) * sqrt(L_final/(N_cover(1) * E_skin));
                    break
                end
            end
            
        end
    end
end


% spar_Cap = 3 * A_s % rule of thumb :)

% R_b = 0.6; R_t = 2.0; F = 0.8059; sigma_ratio = 1.47566;
% a1 = ((4 * F^2 * h_wingBox(1)^2 * h_wingBox(1)^(2/3) * E_skin)/(N_cover(1) * (3.62 * E_skin)^(2/3)))^(2/3);
% a2 = ((totalMoments(1)^2 * h_wingBox(1) * chord(1))/(8 * width_wingBox(1)^2 * E_skin * (h_wingBox(1)/2)^4)) ^ (2/9);
% a3 = (1+R_b * R_t)^(-2/9);
% a4 = ((N_cover(1))/(sigma_ratio * 3.62 * E_skin * (1 + R_b * R_t)))^(-2/27);
% a5  = (1.103 * ((sqrt(F)*(1 + R_b * R_t)/((R_b^3 * R_t * (4 + R_b * R_t))^0.5))) * (N_cover(1) / E_skin)^0.25)^(-2/27);
% 
% L_niklas = (a1*a2*a3*a4*a5)^(1 / (1 - 2/9 + 6/108))

% % corresponding rib thickness for trial
% % crushing force
% t_e_trial = (t1_with_trial * (1 + R_t_trial * R_b_trial));
% I_trial = 2 .* t_e_trial' .* chord .* h_wingBox ./ 2;
% F_crush_trial = (totalMoments.^2 .* h_wingBox .* L .* t_e_trial' .* chord) ./ (2 .* E_skin * I_trial.^2);
% t_rib = ((F_crush_trial .* h_wingBox.^2) ./ (3.62 * E_skin)).^(1/3);
% 
% figure
% plot(x,t_rib)
% 
% % optimal rib spacing
% L_trial = (4 .* F_trial^2 .* width_wingBox.^2 .* t_rib.^2 .* E_rib ./ N_cover').^(1/3);
% figure
% plot(x,L_trial)

% figure
% hold on
% plot(sigma_with,'rx')
% plot(sigma_nada,'bo')


%%
% % integrally stiffened panels - realistic approach
% close all
% % trial one, fixing values of b, R_t and R_b
% R_b_1 = 0.65; R_t_1 = 2.25; sigma_ratio_1 = 1.6;
% 
% F_new_1 = 1.314 * ((R_b_1^3 * R_t_1 * (4 + R_b_1 * R_t_1))^0.25)/(1 + R_b_1 * R_t_1) * sigma_ratio_1^0.25;
% 
% L = 0.5;
% 
% % stringer pitch optimised for the root
% b_new_1 = 1.103 * ((sqrt(F_new_1)*(1 + R_b_1 * R_t_1)/((R_b_1^3 * R_t_1 * (4 + R_b_1 * R_t_1))^0.5))) .* (N_cover(1)' .* L^3 ./ E_skin).^0.25;
% 
% % WITHOUT STRINGERS
% t1_nada = ((N_cover.*L^2) ./ (3.62 * E_skin)).^(1/3);
% A_s_nada = L .* t1_nada;
% sigma_nada = N_cover./t1_nada;
% 
% % WITH STRINGERS
% t_new_1 = ((N_cover .* b_new_1^2) ./ (sigma_ratio_1 * 3.62 * E_skin * (1 + R_b_1 * R_t_1))).^(1/3);
% A_s = t_new_1 * b_new_1;
% T_effective = (t_new_1 + A_s./b_new_1);
% sigma_new_1 = N_cover./(t_new_1 * (1+R_t_1 * R_b_1));
% sigma_z = N_cover./T_effective;
% 
% % yield skin thickness
% t_failure  = N_cover ./ (520*10^6);
% 
% % plotting thickness distriutions
% figure
% hold on
% stairs(x,t1_nada,'r');
% stairs(x,t_new_1,'b');
% stairs(x,t_failure,'g');
% hold off
% 
% figure
% hold on
% plot(x,sigma_nada);
% plot(x,sigma_new_1);
% plot(x,sigma_z);
% plot([x(1), x(100)], [434.5*10^6 434.5*10^6])
% legend('no stringers','stringers included','sigma_z','sigma yield')
% title('stress variation along half span')
% hold off
% SANITY CHECK
% t_new_1 = ((900000 .* b_new_1^2) ./ (1.75 * 3.62 * (70 * 10^9) * (1 + R_b_1 * R_t_1))).^(1/3);

% % buckling stress
% sigma_crit = N_cover ./ (t_new_1 * (1+R_t_1 * R_b_1));

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
% %% rib sizing
% % close all
% 
% h_wingBox = 0.10 * chord; % height of wing box assumed 10% chord
% width_wingBox = (0.65-0.15) * chord; % width of wing box
% 
% % rib thickness
% t_rib = ((abs(distLoad) .* (h_wingBox).^2)./(3.62 .* chord .* E_rib)).^(1/3);
% figure
% plot(x,t_rib)
% 
% % rib spacing
% L = (4 .* 0.81^2 .* width_wingBox.^2 .* t_rib.^2 .* E_rib ./ N_cover').^(1/3);
% figure
% plot(x,L)
% L_min = min(L)
% 
% 
% % stress in ribs
% sigma_ribs = distLoad ./ (t_rib .* chord);
% 
% % crushing force
% 
% % %%
% % %rib sizing
% % % determining crushing force F_crush
% % L_trial = 0.5
% % 
% % % rib thickness
% % t_rib = ((abs(distLoad) .* (h_wingBox).^2)./(3.62 .* chord .* E_rib)).^(1/3);
% % figure
% % plot(x,t_rib)
% % 
% % % rib spacing % farrar efficieny not constant!
% % L = (4 .* F_new_1^2 .* width_wingBox.^2 .* t_rib.^2 .* E_rib ./ N_cover').^(1/3);
% % figure
% % plot(x,L)
% % L_min = min(L)
