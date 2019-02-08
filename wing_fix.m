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
grid on; grid minor;
title('Distributed Loads on Main Wing');
xlabel('Span along Wing [m]');
ylabel('Load per unit meter [N/m]');
plot(x,distLift,'b', 'LineWidth', 2) % lift
plot(x,distWeightWing,'r', 'LineWidth', 2) % self weight
plot(x,distWeightFuel,'g', 'LineWidth', 2) % fuel weight
plot(x,distLoad,'k','LineWidth', 2) % total load
legend('Elliptical Lift Distribution', 'Distributed Weight due to Self-Weight', ...
    'Fuel Weight', 'Total Load (including Engine Loads)');
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
grid on; grid minor;
title('Bending Moment Distribution across Wing');
xlabel('Span along wing [m]');
ylabel('Moment per unit meter [N]');
plot(x,bendingMomentWing,'r', 'LineWidth', 2) % wing self weight and lift 
plot(x,bendingMomentFuel,'b', 'LineWidth', 2) % fuel
plot(x,totalMoments,'k', 'LineWidth', 2) % all
legend('Bending Moment on the Wing due to self-weight and lift', ...
    'Bending Moment due to fuel weight distribution', ...
    'Overall Bending Moments (including Engine loads)');
hold off

% shear force:
totalShearForce = shearForceWing+shearForceFuel;
% add concentrated loads from engines
totalShearForce(1:loc_1) = totalShearForce(1:loc_1) + engineWeight_1;
totalShearForce(1:loc_2) = totalShearForce(1:loc_2) + engineWeight_2;

figure;
hold on
grid on; grid minor;
title('Shear Force Distribution across Wing');
xlabel('Span along wing [m]');
ylabel('Shear Force per unit meter [N/m]');
plot(x,shearForceWing,'r', 'LineWidth', 2)
plot(x,shearForceFuel,'b', 'LineWidth', 2)
plot(x,totalShearForce,'k', 'LineWidth', 2)
legend('Shear Force on the Wing due to self-weight and lift', ...
    'Shear Force due to fuel weight distribution', ...
    'Overall Shear Force (including Engine loads)');
hold off

M_0=0.5*rho*cruiseVelocity^2.*chord.^2*CM0_w;
distLoad(1:loc_1-1) = distLoad(1:loc_1-1) + engineWeight_1;
distLoad(1:loc_2-1) = distLoad(1:loc_2-1) + engineWeight_2;
T=distLift.*a-(distLoad-distLift).*b-M_0;

figure;
hold on
grid on; grid minor;
title('Torque Distribution across Wing');
xlabel('Span along wing [m]');
ylabel('Torque per unit meter [N]');
plot(x,T, 'LineWidth', 2) % torque along span
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
E_skin = 71.7 * 10^9; % final E for Al-7075 T6 Aluminum
E_spar = 76.5 * 10^9; % final E for Al-2050 T6 Aluminum
E_rib = 71.7 * 10^9; % final E for Al-2050 T6 Aluminum

%% spar sizing

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
leg3 = plot([0.15,0.65,0.65,0.15,0.15],[-0.035,-0.035,0.065,0.065,-0.035],'bo-.','LineWidth',1.5);
legend([leg1,leg2,leg3],'NACA profile','Real wing box section','Idealised rectangular wing box section')
hold off

%% Digitalising Catchpole Diagram

% R_t = 0.75
fid = fopen('R_t_0.75.txt');
vals_1 = fscanf(fid, '%f %f ',[2 Inf]);
catchpole_1 = transpose(vals_1);
fclose(fid);

% R_t = 1.00
fid = fopen('R_t_1.00.txt');
vals_2 = fscanf(fid, '%f %f ',[2 Inf]);
catchpole_2 = transpose(vals_2);
fclose(fid);

% R_t = 1.25
fid = fopen('R_t_1.25.txt');
vals_3 = fscanf(fid, '%f %f ',[2 Inf]);
catchpole_3 = transpose(vals_3);
fclose(fid);

% R_t = 1.50
fid = fopen('R_t_1.50.txt');
vals_4 = fscanf(fid, '%f %f ',[2 Inf]);
catchpole_4 = transpose(vals_4);
fclose(fid);

% R_t = 1.75
fid = fopen('R_t_1.75.txt');
vals_5 = fscanf(fid, '%f %f ',[2 Inf]);
catchpole_5 = transpose(vals_5);
fclose(fid);

% R_t = 1.875
fid = fopen('R_t_1.875.txt');
vals_6 = fscanf(fid, '%f %f ',[2 Inf]);
catchpole_6 = transpose(vals_6);
fclose(fid);

% R_t = 2.00
fid = fopen('R_t_2.00.txt');
vals_7 = fscanf(fid, '%f %f ',[2 Inf]);
catchpole_7 = transpose(vals_7);
fclose(fid);

% R_t = 2.25
fid = fopen('R_t_2.25.txt');
vals_8 = fscanf(fid, '%f %f ',[2 Inf]);
catchpole_8 = transpose(vals_8);
fclose(fid);

% R_t = 2.50
fid = fopen('R_t_2.50.txt');
vals_9 = fscanf(fid, '%f %f ',[2 Inf]);
catchpole_9 = transpose(vals_9);
fclose(fid);

% R_t = 3.00
fid = fopen('R_t_3.00.txt');
vals_10 = fscanf(fid, '%f %f ',[2 Inf]);
catchpole_10 = transpose(vals_10);
fclose(fid);


%% integrally stiffened panels ideal - LOOK HERE
close all
clc

% cover load per unit length N

W_optimal = inf;
rho_ribs = 2700; % Al 7075
rho_panels = 2810; % Al 2050

for n_stringers = 25:1:55

    for h = 0.02:0.001:0.05

        for R_t = [0.75]

            for L = 0.4:0.01:1     
                
                chord1 = (c_tip - c_Root)./(wingSemiSpan).*[0:L:wingSemiSpan] + c_Root;
                h_wingBox1 = 0.10 .* chord1; % height of wing box assumed 10% chord
                width_wingBox1 = (0.65-0.15) .* chord1; % width of wing box
                TotalMoments = interp1(x,totalMoments,[0:L:wingSemiSpan]);
                N_cover1 = (TotalMoments ./ (h_wingBox1 .* width_wingBox1));

                b = width_wingBox1(1)/(n_stringers + 1);
                R_b = h/b;

                % Find sigma_ratio

                if (R_t == 0.75)
                    sigma_ratio = interp1q(catchpole_1(:,1),catchpole_1(:,2),R_b)/3.62;
                elseif(R_t == 1)
                    sigma_ratio = interp1q(catchpole_2(:,1),catchpole_2(:,2),R_b)/3.62;
                elseif(R_t == 1.25)
                    sigma_ratio = interp1q(catchpole_3(:,1),catchpole_3(:,2),R_b)/3.62;
                elseif(R_t == 1.5)
                    sigma_ratio = interp1q(catchpole_4(:,1),catchpole_4(:,2),R_b)/3.62;
                elseif(R_t == 1.75)
                    sigma_ratio = interp1q(catchpole_5(:,1),catchpole_5(:,2),R_b)/3.62;
                elseif(R_t == 1.85)
                    sigma_ratio = interp1q(catchpole_6(:,1),catchpole_6(:,2),R_b)/3.62;
                elseif(R_t == 2)
                    sigma_ratio = interp1q(catchpole_7(:,1),catchpole_7(:,2),R_b)/3.62;
                elseif(R_t == 2.25)
                    sigma_ratio = interp1q(catchpole_8(:,1),catchpole_8(:,2),R_b)/3.62;
                elseif(R_t == 2.5)
                    sigma_ratio = interp1q(catchpole_9(:,1),catchpole_9(:,2),R_b/3.62);
                elseif(R_t == 3)
                    sigma_ratio = interp1q(catchpole_10(:,1),catchpole_10(:,2),R_b)/3.62;
                end

                t_plate = (N_cover1.*b.^2./(sigma_ratio.*3.62.*E_skin.*(1 + R_b.*R_t))).^(1/3);
                rho_squared = b^2*R_t*R_b^3*(4 + R_b*R_t)/(12*(1 + R_b*R_t));
                t_euler = N_cover1.*L.^2./((1 + R_b.*R_t)*pi^2*E_skin*rho_squared);

                if(t_plate(1) < t_euler(1))
                    t = t_euler;
                else
                    t = t_plate;
                end

                farrar = 1.314 * ((R_b^3 * R_b * (4 + R_b * R_t))^0.25) / (1 + R_b * R_t) * sigma_ratio^0.25;
                sigma_crit = N_cover1(1)/(t(1) * (1+R_t * R_b));

                if(farrar > 0.81)
                    break
                end

                t_s = R_t .* t;
                A_s = h .* t;
                
                if(h/t_s(1) > 8)
                    continue
                end
                
                t_e = (t .* (1 + R_t * R_b));

                I = 2 .* t_e .* chord1 .* (h_wingBox1 ./ 2).^2;
                F_crush = (TotalMoments.^2 .* h_wingBox1 .* L .* t_e.* chord1) ./ (2 .* E_skin * I.^2);
                t_rib = ((F_crush .* h_wingBox1.^2) ./ (3.62 * E_skin)).^(1/3);

                W_ribs = sum(t_rib.*h_wingBox1.*width_wingBox1)*rho_ribs;
                W_panels = sum((t(2:end) + t(1:end-1))/2.*(1 + R_b*R_t).*(width_wingBox1(2:end) + width_wingBox1(1,end-1))*L/2)*rho_panels;
                W = W_ribs + W_panels;
                
                if(W < W_optimal)
                    R_b_optimal = R_b;
                    R_t_optimal = R_t;
                    b_optimal = b;
                    n_stringers_optimal = n_stringers; 
                    L_optimal = L;
                    h_optimal = h;
                    t_optimal = t;
                    t_s_optimal = t_s;
                    sigma_crit_optimal = sigma_crit;
                    farrar_optimal = farrar;
                    W_optimal = W;

                end

            end

        end

    end

end

disp(['R_b = ',num2str(R_b_optimal)])
disp(['R_t = ',num2str(R_t_optimal)])
disp(['b = ',num2str(b_optimal*100),' cm'])
disp(['Stringers = ',num2str(n_stringers_optimal)])
disp(['L = ',num2str(L_optimal),' m'])
disp(['h = ',num2str(h_optimal*100),' cm'])
disp(['t_root = ',num2str(t_optimal(1)*10^3),' mm'])
disp(['t_s_root = ',num2str(t_s_optimal(1)*10^3),' mm'])
disp(['sigma_crit = ',num2str(sigma_crit_optimal/10^6),' MPa'])
disp(['Farrar efficiency = ',num2str(farrar_optimal)])
disp(['Weight = ',num2str(W_optimal),' kg'])
disp(['Stringer area = ',num2str(A_s(1)*10^6),'mm^2'])


% % optimisation time
% sigma_yield = 300.*10.^6;
% W_optimal = inf;
% tic
% % for b = 0.02:0.01:0.3
% % disp(['Running for b = ',num2str(b)])
%     for R_b = 0.1:0.01:1.1
% 
%         for R_t = [1,1.5,2,2.5,5]
% 
%             for L = 0.1:0.01:1.5
% 
%                 span_loc = [0:L:wingSemiSpan];
%                 panel_loc = (span_loc(2:end) + span_loc(1:end - 1))/2;
%                 rib_loc = span_loc(2:end - 1);
%                 chord_panel = (c_tip - c_Root)./(wingSemiSpan).*panel_loc + c_Root;
%                 chord_rib = (c_tip - c_Root)./(wingSemiSpan).*rib_loc + c_Root;
%                 h_wingBox_rib = 0.10 .* chord_rib; % height of wing box assumed 10% chord
%                 h_wingBox_panel = 0.10 .* chord_panel; % height of wing box assumed 10% chord
%                 width_wingBox_rib = (0.65-0.15) .* chord_rib; % width of wing box
%                 width_wingBox_panel = (0.65-0.15) .* chord_panel; % width of wing box
%                 TotalMoments_rib = totalMoments([ceil(rib_loc./(wingSemiSpan./100))]); 
%                 TotalMoments_panel = totalMoments([ceil(panel_loc./(wingSemiSpan./100))]); 
%                 N_cover = (TotalMoments_panel ./ (h_wingBox_panel .* width_wingBox_panel))';
% 
%                 if (R_t == 1)
%                     sigma_ratio = interp1q(catchpole_1(:,1),catchpole_1(:,2),R_b);
%                 elseif(R_t == 1.5)
%                     sigma_ratio = interp1q(catchpole_2(:,1),catchpole_2(:,2),R_b);
%                 elseif(R_t == 2)
%                     sigma_ratio = interp1q(catchpole_3(:,1),catchpole_3(:,2),R_b);
%                 elseif(R_t == 2.5)
%                     sigma_ratio = interp1q(catchpole_4(:,1),catchpole_4(:,2),R_b);
%                 elseif(R_t == 5)
%                     sigma_ratio = interp1q(catchpole_5(:,1),catchpole_5(:,2),R_b);
%                 end
% 
%                 
%                 rho_squared = sigma_yield.*L.^2./(E_skin.*pi.^2);
%                 b = (rho_squared.*12.*(1 + R_t.*R_b)./(R_t.*R_b.^3.*(4 + R_t.*R_b))).^(1/2);
%                 
%                 b = 0.05;
%                 
%                 t = (N_cover.*b.^2./(sigma_ratio.*3.62.*E_skin.*(1 + R_b.*R_t))).^(1/3);
% 
% 
%                 t = t';
%                 
%                 t_s = R_t.*t;
%                 h = b.*R_b;
% 
%                 I = 2 .* (t(2:end) + t(1:end - 1))./2.*(1 + R_b.*R_t) .* chord_rib .* (h_wingBox_rib ./ 2).^2;
%                 F_crush = (TotalMoments_rib.^2 .* h_wingBox_rib .* L .* (t(2:end) + t(1:end - 1))./2.*(1 + R_b.*R_t).* chord_rib) ./ (2 .* E_skin * I.^2);
%                 t_rib = ((F_crush .* h_wingBox_rib.^2) ./ (3.62 * E_skin)).^(1/3);
% 
%                 W_panels = sum(t.*(1 + R_b.*R_t).*width_wingBox_panel.*L);
%                 W_ribs = sum(t_rib.*h_wingBox_rib.*width_wingBox_rib);
% 
%                 W = W_ribs + W_panels;
% 
%                 if(W < W_optimal)
%                     R_b_optimal = R_b;
%                     R_t_optimal = R_t;
%                     b_optimal = b;
%                     L_optimal = L;
%                     t_optimal = t;
%                     t_s_optimal = t_s;
%                     h_optimal = h;
%                     t_rib_optimal = t_rib;
%                     W_optimal = W;
%                 end
% 
%             end
% 
%         end
% 
%     end
% 
% % end
% toc
% disp(['W_optimal = ',num2str(W_optimal)])
% disp(['R_b_optimal = ',num2str(R_b_optimal)])
% disp(['R_t_optimal = ',num2str(R_t_optimal)])
% disp(['b_optimal = ',num2str(b_optimal)])
% disp(['L_optimal = ',num2str(L_optimal)])
% disp(['h_optimal = ',num2str(h_optimal)])
% disp(['t_rib_optimal = ',num2str(t_rib_optimal)])
% disp(['t_optimal = ',num2str(t_optimal)])
