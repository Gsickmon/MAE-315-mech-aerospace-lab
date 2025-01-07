%lab 0 mech aero
clear
close
clc
%% DATA
syms d F w0 t0 u_c u_tm u_d u_f L 
% d= displacment
% F= force
% t0 = thickness
% w0 = width
% u_c = caliper uncertianty
% u_tm = tapemeasure uncertianty
% u_d = displacment uncertianty
% u_f = force uncertianty
% L = length

A = t0.*w0; % area
stress = F./A;
strain = d./L;
E = stress./strain;

u_stress = sqrt((diff(stress,F).*u_f).^2 + (diff(stress,t0).*u_c).^2 + (diff(stress,w0).*u_c).^2);
u_strain = sqrt((diff(strain,d).*u_d).^2 + (diff(strain,L).*u_tm).^2);
u_E = sqrt((diff(E,F).*u_f).^2 + (diff(E,t0).*u_c).^2 + (diff(E,w0).*u_c).^2 + ...
    (diff(E,d).*u_d).^2 + (diff(E,L).*u_tm).^2);

% manual differentiation
%u_stress = sqrt(((1/(t0.*w0)).*u_f).^2 + ((avgF/((t0.^2).*w0).*u_t).^2) + (-avgF/(t0.*w0.^2).*u_w).^2);
%u_strain = sqrt((((1./L).*u_d).^2)+(((avgd./(L.^2)).*u_L).^2));

%Data input
rawdata = importdata('Intro_specimen.txt','\t',2);
info = rawdata.data;
d = info(:,1); % displacment data
F = info(:,2); % force data
s_E = info(:,3); % strain data for E(not used in this but avaliable)
w0 = 0.412; %in
t0 = 0.057; %in
u_c = 0.0005; %in
u_tm = 1/64; %in
u_d = 0.012;
u_f = 0.012;
L=6.125; %in
A = eval(A);
u_stress = eval(u_stress);
u_strain = eval(u_strain);
u_E = eval(u_E);
u_E = mean(u_E(1:50)); % to find definite number for youngs modulus
stress = eval(stress);
strain = eval(strain);


%% Finding Values

% Finding ultimite stress and strain
ultimite_stress = max(stress); % could also be stress(find(max(stress) == stress))
ultimite_strain = strain(find(max(stress) == stress));

%Rupture stress and strain
Rupture_stress = stress(end);
Rupture_strain = strain(end);

%yield point 0.2 percent offset and Youngs Modulus(E)

p = polyfit(strain(1:50),stress(1:50),1);
E = p(1);
b = stress(1)-p(1).*(strain(1) + 0.002);
p(2) = b;
offset = polyval(p,strain);
yield_intercept = find(offset >= stress,1);
yield_x = strain(yield_intercept);
yield_y = stress(yield_intercept);
start_index = find(offset >= stress(1),1);




%% plotting 

figure(1)
plot(strain,stress,'blue')
hold on
plot(yield_x,yield_y,'diamond','MarkerSize',7,'MarkerFaceColor',[0,0,1])
hold on
plot(strain(start_index:yield_intercept),offset(start_index:yield_intercept),'green')
hold on
plot(ultimite_strain,ultimite_stress,'diamond','MarkerSize',7,'MarkerFaceColor',[0,1,0])
hold on
plot(Rupture_strain,Rupture_stress,'diamond','MarkerSize',7,'MarkerFaceColor',[1,0,0])
title('strain vs stress of carbon steel')
xlabel('strain')
ylabel('stress(lb/in^2)')
errorbar(strain(1:100:end),stress(1:100:end),u_stress(1:100:end),'vertical','.') % vertical errorbars for stress
errorbar(strain(1:100:end),stress(1:100:end),u_strain(1:100:end),'horizontal','.') % horizontal errorbars for strain
legend('strain vs stress curve','yield point','0.2% offset line','ultimite stress and strain', ...
    'rupture point','uncertianty of stress','uncertianty of strain','Location','southeast')

fprintf('Ultimate Stress: %g +- %g\n',ultimite_stress, u_stress(find(max(stress) == stress)))
fprintf('Ultimate Strain: %g +- %g\n',ultimite_strain, u_strain(find(max(stress) == stress)))
fprintf('Rupture Stress: %g +- %g\n',Rupture_stress, u_stress(end))
fprintf('Rupture Strain: %g +- %g\n',Rupture_strain, u_strain(end))
fprintf('Youngs Modulus: %g +- %g\n',E, u_E)
fprintf('Yield Point: (%g,%g)\n',yield_x,yield_y)