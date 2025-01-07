%Lab 1
%MAE 315
%Group 1
clear
close
clc

%% Q1
%% Steel with Extensometer
%DATA
syms d F w0 t0 u_c u_tm u_d u_f L tf wf strain_yg_exten u_extens F2 u_f_2
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
Q = stress .* (1 + strain);
Af = tf.*wf ; % final area
Area_reduction = (1 - (Af./A)).*100 ; % area reduction percent
stress2 = F2 ./A;
yung = stress2./ strain_yg_exten;

u_stress = sqrt((diff(stress,F).*u_f).^2 + (diff(stress,t0).*u_c).^2 + (diff(stress,w0).*u_c).^2);
u_strain = sqrt((diff(strain,d).*u_d).^2 + (diff(strain,L).*u_tm).^2);
u_E = sqrt((diff(E,F).*u_f).^2 + (diff(E,t0).*u_c).^2 + (diff(E,w0).*u_c).^2 + ...
    (diff(E,d).*u_d).^2 + (diff(E,L).*u_tm).^2);
U_Stress_true = sqrt((diff(Q,F).*u_f).^2 + (diff(Q,t0).*u_c).^2 + (diff(Q,w0).*u_c).^2 + ...
    (diff(Q,d).*u_d).^2 + (diff(Q,L).*u_tm).^2);
u_young_MMTS = sqrt((diff(E,F).*u_f).^2 + (diff(E,t0).*u_c).^2 + (diff(E,w0).*u_c).^2 + ...
    (diff(E,d).*u_d).^2 + (diff(E,L).*u_tm).^2);

u_young_exten = sqrt((diff(yung,F2).*u_f_2).^2 + (diff(yung,t0).* u_c).^2 + (diff(yung,w0).*u_c).^2 + (diff(yung,strain_yg_exten) .* u_extens).^2);
U_Area_reduction = sqrt((diff(Af,tf).*u_c).^2 + (diff(Af,wf).*u_c).^2 + (diff(A,t0).*u_c).^2 + (diff(A,w0).*u_c).^2);
poissons = -((tf- t0)./t0)./(d./L);
u_poisson = sqrt((diff(poissons,t0) .* u_c).^2 + (diff(poissons,tf) .* u_c).^2 + (diff(poissons,d) .* u_d).^2 + (diff(poissons,L) .* u_tm).^2);

% manual differentiation
%u_stress = sqrt(((1/(t0.*w0)).*u_f).^2 + ((avgF/((t0.^2).*w0).*u_t).^2) + (-avgF/(t0.*w0.^2).*u_w).^2);
%u_strain = sqrt((((1./L).*u_d).^2)+(((avgd./(L.^2)).*u_L).^2));
disp('steel')
%Data input
rawdata = importdata('M008_Steel.dat','\t',5);
info = rawdata.data;
d = info(:,2); % displacment data
F = info(:,1); % force data
F2 = F(1:50);
s_E = info(:,3); % strain data for E(not used in this but avaliable)

w0 = 0.523; %in change these values for each sample
t0 = 0.062; %in change these values for each sample
L = 6.5; %in change these values for each sample
tf = 0.056;  % final thickness (change)
Lf = 6.9; % final length 
wf = 0.426; %final width
u_extens = 5 .* 10 .^-8; % Reading on axial strain is x e^-7, half of that

u_c = 0.0005; %in 
u_tm = 1/64; %in
u_d = 0.012 .* d;
u_f = 0.012 .* F;
A = eval(A);
Af = eval(Af);
u_stress = eval(u_stress);
u_strain = eval(u_strain);
u_E = eval(u_E);
u_E = mean(u_E(1:50)); % to find definite number for youngs modulus
stress = eval(stress);
strain = eval(strain);
U_strain_true = u_strain;
U_Stress_true = eval(U_Stress_true);
Area_reduction = eval(Area_reduction);
U_Area_reduction = eval(U_Area_reduction);

% Area reduction

% t0 = intial thickness
% w0 = intial thickness
% tf = final thickness 
% L = intial length
% wf =final thickness 

%Finding Values

% Finding ultimite stress and strain
ultimite_stress = max(stress); % could also be stress(find(max(stress) == stress))
ultimite_strain = strain(find(max(stress) == stress));

%Rupture stress and strain
Rupture_stress = stress(end);
Rupture_strain = strain(end);


% true values 
true_stress = stress .* (1+strain);
true_strain = log(1+strain);

True_Rupture_stress = true_stress(end);
True_Rupture_strain = true_strain(end);



fprintf('Ultimate Stress: %g +- %g\n',ultimite_stress, u_stress(find(max(stress) == stress)))
fprintf('Ultimate Strain: %g +- %g\n',ultimite_strain, u_strain(find(max(stress) == stress)))
fprintf('Rupture Stress: %g +- %g\n',Rupture_stress, u_stress(end))
fprintf('Rupture Strain: %g +- %g\n',Rupture_strain, u_strain(end))
fprintf('Rupture Strain: %g +- %g\n',Area_reduction, U_Area_reduction)



% Resiliance and toughness

i = 1;
area_resilience = 0;
U_resilience = 0;
while i <=50
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain(i+1) - strain(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 
area_resilience = area_resilience + area_trap;

d_u_stress_r = u_stress(i) .* (strain(i+1) - strain(i));
U_resilience = d_u_stress_r + U_resilience;

i = i + 1; 
end

i = 1;
area_toughness = 0;
U_toughness = 0;

while i <= length(stress) - 1
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain(i+1) - strain(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 

area_toughness = area_toughness + area_trap;
d_u_stress_t = u_stress(i) .* (strain(i+1) - strain(i));
U_toughness = d_u_stress_t + U_toughness;

i = i + 1; 
end




% True Uncertainty

rupture_index = length(true_stress);
U_Stress_True_rupture = U_Stress_true(rupture_index);
U_Strain_True_rupture = U_strain_true(rupture_index);

fprintf('Area Resilience: %g +- %g\n',area_resilience, U_resilience)
fprintf('Area Toughness: %g +- %g\n',area_toughness, U_toughness)



fprintf('True Ruptur Stress: %g +- %g\n',True_Rupture_stress, U_Stress_True_rupture)
fprintf('True Ruptur Strain: %g +- %g\n',True_Rupture_strain, U_Strain_True_rupture)

% MMTS

stress_yg_MMTS = stress(1:1:50);
strain_yg_MMTS = strain(1:1:50);


estima = fitlm(strain_yg_MMTS,stress_yg_MMTS);
x1=estima.Coefficients;
youngs_modulus_MMTS=x1{'x1', 'Estimate'};


strain_m_MMTS = strain * 1.002;
strain_m_MMTS = strain_m_MMTS(1:100);
vertical_offset_MMTS = stress(1) - youngs_modulus_MMTS .* strain_m_MMTS(1);

stress_m_MMTS = zeros(100,1);
i = 1;
while i <= 100
stress_m_MMTS(i) = youngs_modulus_MMTS .* strain_m_MMTS(i) + vertical_offset_MMTS;
i = i+1;
end

stress_test = stress (1:100);

index_MMTS = find(stress_m_MMTS > stress_test, 1);
yield_stress_MMTS = stress(index_MMTS);
yield_strain_MMTS = strain(index_MMTS);


oset_plot_MMTS = stress_m_MMTS(1:index_MMTS);
strain_oset_MMTS = strain_m_MMTS(1:index_MMTS);

temp = eval(u_young_MMTS);
temp2 = temp(1:index_MMTS);
u_young_MMTS = mean(temp2);

% Exten
stress_yg_exten = stress(1:1:50);
s_E = info(:,3);
strain_yg_exten = s_E(1:1:50);


estima = fitlm(strain_yg_exten,stress_yg_exten);
x1=estima.Coefficients;
youngs_modulus_extensometer=x1{'x1', 'Estimate'};


strain_m_exten = strain * 1.002;
strain_m_exten = strain_m_exten(1:100);
vertical_offset_exten = stress(1) - youngs_modulus_extensometer .* strain_m_exten(1);

stress_m_exten = zeros(100,1);
i = 1;
while i <= 100
stress_m_exten(i) = youngs_modulus_extensometer .* strain_m_exten(i) + vertical_offset_exten;
i = i+1;
end

stress_test = stress (1:100);

index_exten = find(stress_m_exten > stress_test, 1);
yield_stress_Extensometer = stress(index_exten);
yield_strain_Extensometer = strain(index_exten);


oset_plot_exten = stress_m_exten(1:index_exten);
strain_oset_exten = strain_m_exten(1:index_exten);




fprintf('Yield Stress: %g\n',yield_stress_MMTS)
fprintf('Yield Strain: %g\n',yield_strain_MMTS)
fprintf('Youngs Modulus: %g +- %g\n', youngs_modulus_MMTS, u_young_MMTS)
% Poisson Ratio


lateral_strain = ((tf-t0)/t0);
poissons = -(lateral_strain)/(strain(end));
u_poisson = eval(u_poisson);
uid = length(u_poisson);
U_poissons = u_poisson(uid);
fprintf('Poissons Ratio: %g +- %g\n',poissons, U_poissons)



% Plotting 

figure(1)
plot(strain,stress,'blue')
hold on

plot(ultimite_strain,ultimite_stress,'diamond','MarkerSize',7,'MarkerFaceColor',[0,1,0])
hold on
plot(Rupture_strain,Rupture_stress,'diamond','MarkerSize',7,'MarkerFaceColor',[1,0,0])
title('Strain vs Stress of Steel')
xlabel('Strain')
ylabel('Stress(lb/in^2)')
errorbar(strain(1:40:end),stress(1:40:end),u_stress(1:40:end),'vertical','.') % vertical errorbars for stress
errorbar(strain(1:40:end),stress(1:40:end),u_strain(1:40:end),'horizontal','.') % horizontal errorbars for strain
plot(true_strain, true_stress,'red')
hold on
plot(True_Rupture_strain,True_Rupture_stress,'diamond','MarkerSize',7,'MarkerFaceColor',[0.5,0.5,0])
hold on
plot(strain_oset_MMTS,oset_plot_MMTS,'yellow')
hold on
plot(yield_strain_Extensometer,yield_stress_Extensometer,'diamond','MarkerSize',7,'MarkerFaceColor',[1,1,1])
hold on
plot(strain_oset_exten,oset_plot_exten,'green')
hold on
plot(yield_strain_MMTS,yield_stress_MMTS,'diamond','MarkerSize',7,'MarkerFaceColor',[0,1,1])
legend('strain vs stress curve','ultimite stress and strain', 'rupture point', ...
   'uncertianty of stress','uncertianty of strain','true strain vs stress curve','true rupture','0.2% offset line','yield point','Location','south')

fprintf('Youngs Modulus: %g +- %g\n',youngs_modulus_MMTS, u_E)

% True Uncertainty

rupture_index = length(true_stress);
U_Stress_True_rupture = U_Stress_true(rupture_index);
U_Strain_True_rupture = U_strain_true(rupture_index);
%% Steel
% DATA
syms d F w0 t0 u_c u_tm u_d u_f L tf wf strain_yg_exten u_extens F2 u_f_2
% d= displacment
% F= force
% t0 = thickness
% w0 = width
% wf = final width
% tf = final thickness
% u_c = caliper uncertianty
% u_tm = tapemeasure uncertianty
% u_d = displacment uncertianty
% u_f = force uncertianty
% L = length

A = t0.*w0; % area
stress = F./A;
strain = d./L;
Af = tf.*wf ; % final area
Area_reduction = (1 - (Af./A)).*100 ; % area reduction percent 
E = stress./strain;
stress2 = F2 ./A;
yung = stress2./ strain_yg_exten;
Q = stress .* (1 + strain); % true stress
u_stress = sqrt((diff(stress,F).*u_f).^2 + (diff(stress,t0).*u_c).^2 + (diff(stress,w0).*u_c).^2);
u_strain = sqrt((diff(strain,d).*u_d).^2 + (diff(strain,L).*u_tm).^2);
u_E = sqrt((diff(E,F).*u_f).^2 + (diff(E,t0).*u_c).^2 + (diff(E,w0).*u_c).^2 + ...
    (diff(E,d).*u_d).^2 + (diff(E,L).*u_tm).^2);
u_young_MMTS = sqrt((diff(E,F).*u_f).^2 + (diff(E,t0).*u_c).^2 + (diff(E,w0).*u_c).^2 + ...
    (diff(E,d).*u_d).^2 + (diff(E,L).*u_tm).^2);

u_young_exten = sqrt((diff(yung,F2).*u_f_2).^2 + (diff(yung,t0).* u_c).^2 + (diff(yung,w0).*u_c).^2 + (diff(yung,strain_yg_exten) .* u_extens).^2);


U_Stress_true = sqrt((diff(Q,F).*u_f).^2 + (diff(Q,t0).*u_c).^2 + (diff(Q,w0).*u_c).^2 + ...
    (diff(Q,d).*u_d).^2 + (diff(Q,L).*u_tm).^2);
U_Area_reduction = sqrt((diff(Af,tf).*u_c).^2 + (diff(Af,wf).*u_c).^2 + (diff(A,t0).*u_c).^2 + (diff(A,w0).*u_c).^2);
poissons = -((tf- t0)./t0)./(d./L);
u_poisson = sqrt((diff(poissons,t0) .* u_c).^2 + (diff(poissons,tf) .* u_c).^2 + (diff(poissons,d) .* u_d).^2 + (diff(poissons,L) .* u_tm).^2);
% manual differentiation
%u_stress = sqrt(((1/(t0.*w0)).*u_f).^2 + ((avgF/((t0.^2).*w0).*u_t).^2) + (-avgF/(t0.*w0.^2).*u_w).^2);
%u_strain = sqrt((((1./L).*u_d).^2)+(((avgd./(L.^2)).*u_L).^2));
disp('steel')
%Data input
rawdata = importdata('M008_Steel.dat','\t',5);
info = rawdata.data;
d = info(:,2); % displacment data
F = info(:,1); % force data
F2 = F(1:50);
s_E = info(:,3); % strain data for E(not used in this but avaliable)

L=6.5; %in change these values for each sample
w0 = 0.523; %in change these values for each sample
t0 = 0.062; %in change these values for each sample
Lf = 6.9; % final length 
wf = 0.426; %final width
tf = 0.056;  % final thickness (change)

u_extens = 5 .* 10 .^-8; % Reading on axial strain is x e^-7, half of that
u_c = 0.0005; %in 
u_tm = 1/64; %in
u_d = 0.012 .* d;
u_f = 0.012 .* F;
u_f_2 = .012 .*F2;
A = eval(A);
Af = eval(Af);
u_stress = eval(u_stress);
u_strain = eval(u_strain);

stress = eval(stress);
strain = eval(strain);
U_strain_true = u_strain;
U_Stress_true = eval(U_Stress_true);
Area_reduction = eval(Area_reduction);
U_Area_reduction = eval(U_Area_reduction);

% Area reduction

% t0 = intial thickness
% w0 = intial thickness
% tf = final thickness 
% L = intial length
% wf =final thickness 




% Finding Values

% Finding ultimite stress and strain
ultimite_stress = max(stress); % could also be stress(find(max(stress) == stress))
ultimite_strain = strain(find(max(stress) == stress));

%Rupture stress and strain
Rupture_stress = stress(end);
Rupture_strain = strain(end);

%yield point 0.2 percent offset and Youngs Modulus(E)

p = polyfit(strain(1:50),stress(1:50),1);
E = p(1);
b = stress(1)-E.*(strain(1) .* 1.002);
p(2) = b;
offset = polyval(p,strain(1));
yield_intercept = find(offset >= stress,1);
yield_x = strain(yield_intercept);
yield_y = stress(yield_intercept);
start_index = find(offset >= stress(1),1);

% true values 
true_stress = stress .* (1+strain);
true_strain = log(1+strain);

True_Rupture_stress = true_stress(end);
True_Rupture_strain = true_strain(end);



fprintf('Ultimate Stress: %g +- %g\n',ultimite_stress, u_stress(find(max(stress) == stress)))
fprintf('Ultimate Strain: %g +- %g\n',ultimite_strain, u_strain(find(max(stress) == stress)))
fprintf('Rupture Stress: %g +- %g\n',Rupture_stress, u_stress(end))
fprintf('Rupture Strain: %g +- %g\n',Rupture_strain, u_strain(end))
fprintf('Yield Point: (%g,%g)\n',yield_x,yield_y)



% Resiliance and toughness

i = 1;
area_resilience = 0;
U_resilience = 0;
while i <=50
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain(i+1) - strain(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 
area_resilience = area_resilience + area_trap;

d_u_stress_r = u_stress(i) .* (strain(i+1) - strain(i));
U_resilience = d_u_stress_r + U_resilience;

i = i + 1; 
end

i = 1;
area_toughness = 0;
U_toughness = 0;

while i <= length(stress) - 1
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain(i+1) - strain(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 

area_toughness = area_toughness + area_trap;
d_u_stress_t = u_stress(i) .* (strain(i+1) - strain(i));
U_toughness = d_u_stress_t + U_toughness;

i = i + 1; 
end

% True Uncertainty

rupture_index = length(true_stress);
U_Stress_True_rupture = U_Stress_true(rupture_index);
U_Strain_True_rupture = U_strain_true(rupture_index);

fprintf('Area Resilience: %g +- %g\n',area_resilience, U_resilience)
fprintf('Area Toughness: %g +- %g\n',area_toughness, U_toughness)

fprintf('True Rupture Stress: %g\n',True_Rupture_stress)
fprintf('True Rupture Strain: %g\n',True_Rupture_strain)

fprintf('True Ruptur Stress: %g +- %g\n',True_Rupture_stress, U_Stress_True_rupture)
fprintf('True Ruptur Strain: %g +- %g\n',True_Rupture_strain, U_Strain_True_rupture)





% MMTS

stress_yg_MMTS = stress(1:1:50);
strain_yg_MMTS = strain(1:1:50);


estima = fitlm(strain_yg_MMTS,stress_yg_MMTS);
x1=estima.Coefficients;
youngs_modulus_MMTS=x1{'x1', 'Estimate'};


strain_m_MMTS = strain * 1.002;
strain_m_MMTS = strain_m_MMTS(1:100);
vertical_offset_MMTS = stress(1) - youngs_modulus_MMTS .* strain_m_MMTS(1);

stress_m_MMTS = zeros(100,1);
i = 1;
while i <= 100
stress_m_MMTS(i) = youngs_modulus_MMTS .* strain_m_MMTS(i) + vertical_offset_MMTS;
i = i+1;
end

stress_test = stress (1:100);

index_MMTS = find(stress_m_MMTS > stress_test, 1);
yield_stress_MMTS = stress(index_MMTS);
yield_strain_MMTS = strain(index_MMTS);


oset_plot_MMTS = stress_m_MMTS(1:index_MMTS);
strain_oset_MMTS = strain_m_MMTS(1:index_MMTS);


i = 1;
area_resilience_MMTS = 0;
U_resilience_MMTS = 0;
h1 = 0;
h2 = 0;
b = 0;
area_rect = 0;
area_trap = 0;

while i <= index_MMTS
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain_m_MMTS(i+1) - strain_m_MMTS(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 
area_resilience_MMTS = area_resilience_MMTS + area_trap;

d_u_stress_r_MMTS = u_stress(i) .* (strain_m_MMTS(i+1) - strain_m_MMTS(i));
U_resilience_MMTS = d_u_stress_r_MMTS + U_resilience_MMTS;

i = i + 1; 
end
temp = eval(u_young_MMTS);
temp2 = temp(1:index_MMTS);
u_young_MMTS = mean(temp2)
% Exten
stress_yg_exten = stress(1:1:50);
s_E = info(:,3);
strain_yg_exten = s_E(1:1:50);


estima = fitlm(strain_yg_exten,stress_yg_exten);
x1=estima.Coefficients;
youngs_modulus_extensometer=x1{'x1', 'Estimate'};


strain_m_exten = strain * 1.002;
strain_m_exten = strain_m_exten(1:100);
vertical_offset_exten = stress(1) - youngs_modulus_extensometer .* strain_m_exten(1);

stress_m_exten = zeros(100,1);
i = 1;
while i <= 100
stress_m_exten(i) = youngs_modulus_extensometer .* strain_m_exten(i) + vertical_offset_exten;
i = i+1;
end

stress_test = stress (1:100);

index_exten = find(stress_m_exten > stress_test, 1);
yield_stress_Extensometer = stress(index_exten);
yield_strain_Extensometer = strain(index_exten);


oset_plot_exten = stress_m_exten(1:index_exten);
strain_oset_exten = strain_m_exten(1:index_exten);

i = 1;
area_resilience_EXTEN = 0;
U_resilience_EXTEN = 0;
h1 = 0;
h2 = 0;
b = 0;
area_rect = 0;
area_trap = 0;

while i <= index_exten
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain_m_exten(i+1) - strain_m_exten(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 
area_resilience_EXTEN = area_resilience_EXTEN + area_trap;

d_u_stress_r_EXTEN = u_stress(i) .* (strain_m_exten(i+1) - strain_m_exten(i));
U_resilience_EXTEN = d_u_stress_r_EXTEN + U_resilience_EXTEN;

i = i + 1; 
end
temp3 = eval(u_young_exten);
temp4 = temp3(1:index_exten);
u_young_Exten = mean(temp4)

fprintf('Yield Stress: %g\n',yield_stress_MMTS)
fprintf('Yield Strain: %g\n',yield_strain_MMTS)
fprintf('Area Reduction Percent: %g +- %g\n',Area_reduction, U_Area_reduction)

% poisson Ratio


lateral_strain = ((tf-t0)/t0);
poissons = -(lateral_strain)/(strain(end));
u_poisson = eval(u_poisson);
uid = length(u_poisson);
U_poissons = u_poisson(uid);
fprintf('Poissons Ratio: %g +- %g\n',poissons, U_poissons)

% plotting 

figure(2)
plot(strain,stress,'blue')
hold on
plot(yield_x,yield_y,'diamond','MarkerSize',7,'MarkerFaceColor',[0,0,1])
hold on
plot(strain(start_index:yield_intercept),offset(start_index:yield_intercept),'green')
hold on
plot(ultimite_strain,ultimite_stress,'diamond','MarkerSize',7,'MarkerFaceColor',[0,1,0])
hold on
plot(Rupture_strain,Rupture_stress,'diamond','MarkerSize',7,'MarkerFaceColor',[1,0,0])
title('strain vs stress of Steel')
xlabel('strain')
ylabel('stress(lb/in^2)')
errorbar(strain(1:40:end),stress(1:40:end),u_stress(1:40:end),'vertical','.') % vertical errorbars for stress
errorbar(strain(1:40:end),stress(1:40:end),u_strain(1:40:end),'horizontal','.') % horizontal errorbars for strain
plot(true_strain, true_stress,'red')
hold on
plot(True_Rupture_strain,True_Rupture_stress,'diamond','MarkerSize',7,'MarkerFaceColor',[0.5,0.5,0])
hold on
plot(strain_oset_MMTS,oset_plot_MMTS,'green')
hold on
%plot(yield_strain_Extensometer,yield_stress_Extensometer,'diamond','MarkerSize',7,'MarkerFaceColor',[1,1,1])
hold on
%plot(strain_oset_exten,oset_plot_exten)
hold on
plot(yield_strain_MMTS,yield_stress_MMTS,'diamond','MarkerSize',7,'MarkerFaceColor',[0,1,1])
legend('strain vs stress curve','ultimite stress and strain', 'rupture point', ...
   'uncertianty of stress','uncertianty of strain','true strain vs stress curve','true rupture','0.2% offset line','yield point','Location','southeast')

%fprintf('Youngs Modulus: %g +- %g\n',youngs_modulus_MMTS, u_E)

%% Aluminum
clear
% DATA
syms d F w0 t0 u_c u_tm u_d u_f L tf wf strain_yg_exten u_extens F2 u_f_2
% d= displacment
% F= force
% t0 = thickness
% w0 = width
% wf = final width
% tf = final thickness
% u_c = caliper uncertianty
% u_tm = tapemeasure uncertianty
% u_d = displacment uncertianty
% u_f = force uncertianty
% L = length

A = t0.*w0; % area
stress = F./A;
strain = d./L;
Af = tf.*wf ; % final area
Area_reduction = (1 - (Af./A)).*100 ; % area reduction percent 
E = stress./strain;
stress2 = F2 ./A;
yung = stress2./ strain_yg_exten;
Q = stress .* (1 + strain); % true stress
u_stress = sqrt((diff(stress,F).*u_f).^2 + (diff(stress,t0).*u_c).^2 + (diff(stress,w0).*u_c).^2);
u_strain = sqrt((diff(strain,d).*u_d).^2 + (diff(strain,L).*u_tm).^2);
u_E = sqrt((diff(E,F).*u_f).^2 + (diff(E,t0).*u_c).^2 + (diff(E,w0).*u_c).^2 + ...
    (diff(E,d).*u_d).^2 + (diff(E,L).*u_tm).^2);
u_young_MMTS = sqrt((diff(E,F).*u_f).^2 + (diff(E,t0).*u_c).^2 + (diff(E,w0).*u_c).^2 + ...
    (diff(E,d).*u_d).^2 + (diff(E,L).*u_tm).^2);

u_young_exten = sqrt((diff(yung,F2).*u_f_2).^2 + (diff(yung,t0).* u_c).^2 + (diff(yung,w0).*u_c).^2 + (diff(yung,strain_yg_exten) .* u_extens).^2);


U_Stress_true = sqrt((diff(Q,F).*u_f).^2 + (diff(Q,t0).*u_c).^2 + (diff(Q,w0).*u_c).^2 + ...
    (diff(Q,d).*u_d).^2 + (diff(Q,L).*u_tm).^2);
U_Area_reduction = sqrt((diff(Af,tf).*u_c).^2 + (diff(Af,wf).*u_c).^2 + (diff(A,t0).*u_c).^2 + (diff(A,w0).*u_c).^2);
poissons = -((tf- t0)./t0)./(d./L);
u_poisson = sqrt((diff(poissons,t0) .* u_c).^2 + (diff(poissons,tf) .* u_c).^2 + (diff(poissons,d) .* u_d).^2 + (diff(poissons,L) .* u_tm).^2);
% manual differentiation
%u_stress = sqrt(((1/(t0.*w0)).*u_f).^2 + ((avgF/((t0.^2).*w0).*u_t).^2) + (-avgF/(t0.*w0.^2).*u_w).^2);
%u_strain = sqrt((((1./L).*u_d).^2)+(((avgd./(L.^2)).*u_L).^2));
disp('Aluminum')
%Data input
rawdata = importdata('M008_Aluminum.dat','\t',5);
info = rawdata.data;
d = info(:,2); % displacment data
F = info(:,1); % force data
F2 = F(1:50);
s_E = info(:,3); % strain data for E(not used in this but avaliable)

L=6.5; %in change these values for each sample
w0 = 0.518; %in change these values for each sample
t0 = 0.062; %in change these values for each sample
Lf = 6.625; % final length 
wf = 0.462; %final width
tf = 0.0615;  % final thickness (change)

u_extens = 5 .* 10 .^-8; % Reading on axial strain is x e^-7, half of that
u_c = 0.0005; %in 
u_tm = 1/64; %in
u_d = 0.012 .* d;
u_f = 0.012 .* F;
u_f_2 = .012 .*F2;
A = eval(A);
Af = eval(Af);
u_stress = eval(u_stress);
u_strain = eval(u_strain);

stress = eval(stress);
strain = eval(strain);
U_strain_true = u_strain;
U_Stress_true = eval(U_Stress_true);
Area_reduction = eval(Area_reduction);
U_Area_reduction = eval(U_Area_reduction);

% Area reduction

% t0 = intial thickness
% w0 = intial thickness
% tf = final thickness 
% L = intial length
% wf =final thickness 




% Finding Values

% Finding ultimite stress and strain
ultimite_stress = max(stress); % could also be stress(find(max(stress) == stress))
ultimite_strain = strain(find(max(stress) == stress));

%Rupture stress and strain
Rupture_stress = stress(end);
Rupture_strain = strain(end);

%yield point 0.2 percent offset and Youngs Modulus(E)

p = polyfit(strain(1:50),stress(1:50),1);
E = p(1);
b = stress(1)-E.*(strain(1) .* 1.002);
p(2) = b;
offset = polyval(p,strain(1));
yield_intercept = find(offset >= stress,1);
yield_x = strain(yield_intercept);
yield_y = stress(yield_intercept);
start_index = find(offset >= stress(1),1);

% true values 
true_stress = stress .* (1+strain);
true_strain = log(1+strain);

True_Rupture_stress = true_stress(end);
True_Rupture_strain = true_strain(end);



fprintf('Ultimate Stress: %g +- %g\n',ultimite_stress, u_stress(find(max(stress) == stress)))
fprintf('Ultimate Strain: %g +- %g\n',ultimite_strain, u_strain(find(max(stress) == stress)))
fprintf('Rupture Stress: %g +- %g\n',Rupture_stress, u_stress(end))
fprintf('Rupture Strain: %g +- %g\n',Rupture_strain, u_strain(end))
fprintf('Yield Point: (%g,%g)\n',yield_x,yield_y)



% Resiliance and toughness

i = 1;
area_resilience = 0;
U_resilience = 0;
while i <=50
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain(i+1) - strain(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 
area_resilience = area_resilience + area_trap;

d_u_stress_r = u_stress(i) .* (strain(i+1) - strain(i));
U_resilience = d_u_stress_r + U_resilience;

i = i + 1; 
end

i = 1;
area_toughness = 0;
U_toughness = 0;

while i <= length(stress) - 1
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain(i+1) - strain(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 

area_toughness = area_toughness + area_trap;
d_u_stress_t = u_stress(i) .* (strain(i+1) - strain(i));
U_toughness = d_u_stress_t + U_toughness;

i = i + 1; 
end

% True Uncertainty

rupture_index = length(true_stress);
U_Stress_True_rupture = U_Stress_true(rupture_index);
U_Strain_True_rupture = U_strain_true(rupture_index);

fprintf('Area Resilience: %g +- %g\n',area_resilience, U_resilience)
fprintf('Area Toughness: %g +- %g\n',area_toughness, U_toughness)

fprintf('True Rupture Stress: %g\n',True_Rupture_stress)
fprintf('True Rupture Strain: %g\n',True_Rupture_strain)

fprintf('True Ruptur Stress: %g +- %g\n',True_Rupture_stress, U_Stress_True_rupture)
fprintf('True Ruptur Strain: %g +- %g\n',True_Rupture_strain, U_Strain_True_rupture)





% MMTS

stress_yg_MMTS = stress(1:1:50);
strain_yg_MMTS = strain(1:1:50);


estima = fitlm(strain_yg_MMTS,stress_yg_MMTS);
x1=estima.Coefficients;
youngs_modulus_MMTS=x1{'x1', 'Estimate'};


strain_m_MMTS = strain * 1.002;
strain_m_MMTS = strain_m_MMTS(1:100);
vertical_offset_MMTS = stress(1) - youngs_modulus_MMTS .* strain_m_MMTS(1);

stress_m_MMTS = zeros(100,1);
i = 1;
while i <= 100
stress_m_MMTS(i) = youngs_modulus_MMTS .* strain_m_MMTS(i) + vertical_offset_MMTS;
i = i+1;
end

stress_test = stress (1:100);

index_MMTS = find(stress_m_MMTS > stress_test, 1);
yield_stress_MMTS = stress(index_MMTS);
yield_strain_MMTS = strain(index_MMTS);


oset_plot_MMTS = stress_m_MMTS(1:index_MMTS);
strain_oset_MMTS = strain_m_MMTS(1:index_MMTS);


i = 1;
area_resilience_MMTS = 0;
U_resilience_MMTS = 0;
h1 = 0;
h2 = 0;
b = 0;
area_rect = 0;
area_trap = 0;

while i <= index_MMTS
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain_m_MMTS(i+1) - strain_m_MMTS(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 
area_resilience_MMTS = area_resilience_MMTS + area_trap;

d_u_stress_r_MMTS = u_stress(i) .* (strain_m_MMTS(i+1) - strain_m_MMTS(i));
U_resilience_MMTS = d_u_stress_r_MMTS + U_resilience_MMTS;

i = i + 1; 
end
temp = eval(u_young_MMTS);
temp2 = temp(1:index_MMTS);
u_young_MMTS = mean(temp2)
% Exten
stress_yg_exten = stress(1:1:50);
s_E = info(:,3);
strain_yg_exten = s_E(1:1:50);


estima = fitlm(strain_yg_exten,stress_yg_exten);
x1=estima.Coefficients;
youngs_modulus_extensometer=x1{'x1', 'Estimate'};


strain_m_exten = strain * 1.002;
strain_m_exten = strain_m_exten(1:100);
vertical_offset_exten = stress(1) - youngs_modulus_extensometer .* strain_m_exten(1);

stress_m_exten = zeros(100,1);
i = 1;
while i <= 100
stress_m_exten(i) = youngs_modulus_extensometer .* strain_m_exten(i) + vertical_offset_exten;
i = i+1;
end

stress_test = stress (1:100);

index_exten = find(stress_m_exten > stress_test, 1);
yield_stress_Extensometer = stress(index_exten);
yield_strain_Extensometer = strain(index_exten);


oset_plot_exten = stress_m_exten(1:index_exten);
strain_oset_exten = strain_m_exten(1:index_exten);

i = 1;
area_resilience_EXTEN = 0;
U_resilience_EXTEN = 0;
h1 = 0;
h2 = 0;
b = 0;
area_rect = 0;
area_trap = 0;

while i <= index_exten
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain_m_exten(i+1) - strain_m_exten(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 
area_resilience_EXTEN = area_resilience_EXTEN + area_trap;

d_u_stress_r_EXTEN = u_stress(i) .* (strain_m_exten(i+1) - strain_m_exten(i));
U_resilience_EXTEN = d_u_stress_r_EXTEN + U_resilience_EXTEN;

i = i + 1; 
end
temp3 = eval(u_young_exten);
temp4 = temp3(1:index_exten);
u_young_Exten = mean(temp4)

fprintf('Yield Stress: %g\n',yield_stress_MMTS)
fprintf('Yield Strain: %g\n',yield_strain_MMTS)
fprintf('Area Reduction Percent: %g +- %g\n',Area_reduction, U_Area_reduction)

% poisson Ratio


lateral_strain = ((tf-t0)/t0);
poissons = -(lateral_strain)/(strain(end));
u_poisson = eval(u_poisson);
uid = length(u_poisson);
U_poissons = u_poisson(uid);
fprintf('Poissons Ratio: %g +- %g\n',poissons, U_poissons)

% plotting 

figure(3)
plot(strain,stress,'blue')
hold on
plot(yield_x,yield_y,'diamond','MarkerSize',7,'MarkerFaceColor',[0,0,1])
hold on
plot(strain(start_index:yield_intercept),offset(start_index:yield_intercept),'green')
hold on
plot(ultimite_strain,ultimite_stress,'diamond','MarkerSize',7,'MarkerFaceColor',[0,1,0])
hold on
plot(Rupture_strain,Rupture_stress,'diamond','MarkerSize',7,'MarkerFaceColor',[1,0,0])
title('strain vs stress of Aluminum')
xlabel('strain')
ylabel('stress(lb/in^2)')
errorbar(strain(1:40:end),stress(1:40:end),u_stress(1:40:end),'vertical','.') % vertical errorbars for stress
errorbar(strain(1:40:end),stress(1:40:end),u_strain(1:40:end),'horizontal','.') % horizontal errorbars for strain
plot(true_strain, true_stress,'red')
hold on
plot(True_Rupture_strain,True_Rupture_stress,'diamond','MarkerSize',7,'MarkerFaceColor',[0.5,0.5,0])
hold on
plot(strain_oset_MMTS,oset_plot_MMTS,'green')
hold on
%plot(yield_strain_Extensometer,yield_stress_Extensometer,'diamond','MarkerSize',7,'MarkerFaceColor',[1,1,1])
hold on
%plot(strain_oset_exten,oset_plot_exten)
hold on
plot(yield_strain_MMTS,yield_stress_MMTS,'diamond','MarkerSize',7,'MarkerFaceColor',[0,1,1])
legend('strain vs stress curve','ultimite stress and strain', 'rupture point', ...
   'uncertianty of stress','uncertianty of strain','true strain vs stress curve','true rupture','0.2% offset line','yield point','Location','southeast')

%fprintf('Youngs Modulus: %g +- %g\n',youngs_modulus_MMTS, u_E)


%% Carbon Fiber 45
clear
% DATA
syms d F w0 t0 u_c u_tm u_d u_f L tf wf strain_yg_exten u_extens F2 u_f_2
% d= displacment
% F= force
% t0 = thickness
% w0 = width
% wf = final width
% tf = final thickness
% u_c = caliper uncertianty
% u_tm = tapemeasure uncertianty
% u_d = displacment uncertianty
% u_f = force uncertianty
% L = length

A = t0.*w0; % area
stress = F./A;
strain = d./L;
Af = tf.*wf ; % final area
Area_reduction = (1 - (Af./A)).*100 ; % area reduction percent 
E = stress./strain;
stress2 = F2 ./A;
yung = stress2./ strain_yg_exten;
Q = stress .* (1 + strain); % true stress
u_stress = sqrt((diff(stress,F).*u_f).^2 + (diff(stress,t0).*u_c).^2 + (diff(stress,w0).*u_c).^2);
u_strain = sqrt((diff(strain,d).*u_d).^2 + (diff(strain,L).*u_tm).^2);
u_E = sqrt((diff(E,F).*u_f).^2 + (diff(E,t0).*u_c).^2 + (diff(E,w0).*u_c).^2 + ...
    (diff(E,d).*u_d).^2 + (diff(E,L).*u_tm).^2);
u_young_MMTS = sqrt((diff(E,F).*u_f).^2 + (diff(E,t0).*u_c).^2 + (diff(E,w0).*u_c).^2 + ...
    (diff(E,d).*u_d).^2 + (diff(E,L).*u_tm).^2);

u_young_exten = sqrt((diff(yung,F2).*u_f_2).^2 + (diff(yung,t0).* u_c).^2 + (diff(yung,w0).*u_c).^2 + (diff(yung,strain_yg_exten) .* u_extens).^2);


U_Stress_true = sqrt((diff(Q,F).*u_f).^2 + (diff(Q,t0).*u_c).^2 + (diff(Q,w0).*u_c).^2 + ...
    (diff(Q,d).*u_d).^2 + (diff(Q,L).*u_tm).^2);
U_Area_reduction = sqrt((diff(Af,tf).*u_c).^2 + (diff(Af,wf).*u_c).^2 + (diff(A,t0).*u_c).^2 + (diff(A,w0).*u_c).^2);
poissons = -((tf- t0)./t0)./(d./L);
u_poisson = sqrt((diff(poissons,t0) .* u_c).^2 + (diff(poissons,tf) .* u_c).^2 + (diff(poissons,d) .* u_d).^2 + (diff(poissons,L) .* u_tm).^2);
% manual differentiation
%u_stress = sqrt(((1/(t0.*w0)).*u_f).^2 + ((avgF/((t0.^2).*w0).*u_t).^2) + (-avgF/(t0.*w0.^2).*u_w).^2);
%u_strain = sqrt((((1./L).*u_d).^2)+(((avgd./(L.^2)).*u_L).^2));
disp('Carbon Fiber 45')
%Data input
rawdata = importdata('M008_Carbon_Fiber_45.dat','\t',5);
info = rawdata.data;
d = info(:,2); % displacment data
F = info(:,1); % force data
F2 = F(1:50);
s_E = info(:,3); % strain data for E(not used in this but avaliable)

L=6.5; %in change these values for each sample
w0 = 0.518; %in change these values for each sample
t0 = 0.060; %in change these values for each sample
Lf = 7.375; % final length 
wf = 0.481; %final width
tf = 0.060;  % final thickness (change)

u_extens = 5 .* 10 .^-8; % Reading on axial strain is x e^-7, half of that
u_c = 0.0005; %in 
u_tm = 1/64; %in
u_d = 0.012 .* d;
u_f = 0.012 .* F;
u_f_2 = .012 .*F2;
A = eval(A);
Af = eval(Af);
u_stress = eval(u_stress);
u_strain = eval(u_strain);

stress = eval(stress);
strain = eval(strain);
U_strain_true = u_strain;
U_Stress_true = eval(U_Stress_true);
Area_reduction = eval(Area_reduction);
U_Area_reduction = eval(U_Area_reduction);

% Area reduction

% t0 = intial thickness
% w0 = intial thickness
% tf = final thickness 
% L = intial length
% wf =final thickness 




% Finding Values

% Finding ultimite stress and strain
ultimite_stress = max(stress); % could also be stress(find(max(stress) == stress))
ultimite_strain = strain(find(max(stress) == stress));

%Rupture stress and strain
Rupture_stress = stress(end);
Rupture_strain = strain(end);

%yield point 0.2 percent offset and Youngs Modulus(E)

p = polyfit(strain(1:50),stress(1:50),1);
E = p(1);
b = stress(1)-E.*(strain(1) .* 1.002);
p(2) = b;
offset = polyval(p,strain(1));
yield_intercept = find(offset >= stress,1);
yield_x = strain(yield_intercept);
yield_y = stress(yield_intercept);
start_index = find(offset >= stress(1),1);

% true values 
true_stress = stress .* (1+strain);
true_strain = log(1+strain);

True_Rupture_stress = true_stress(end);
True_Rupture_strain = true_strain(end);



fprintf('Ultimate Stress: %g +- %g\n',ultimite_stress, u_stress(find(max(stress) == stress)))
fprintf('Ultimate Strain: %g +- %g\n',ultimite_strain, u_strain(find(max(stress) == stress)))
fprintf('Rupture Stress: %g +- %g\n',Rupture_stress, u_stress(end))
fprintf('Rupture Strain: %g +- %g\n',Rupture_strain, u_strain(end))
fprintf('Yield Point: (%g,%g)\n',yield_x,yield_y)



% Resiliance and toughness

i = 1;
area_resilience = 0;
U_resilience = 0;
while i <=50
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain(i+1) - strain(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 
area_resilience = area_resilience + area_trap;

d_u_stress_r = u_stress(i) .* (strain(i+1) - strain(i));
U_resilience = d_u_stress_r + U_resilience;

i = i + 1; 
end

i = 1;
area_toughness = 0;
U_toughness = 0;

while i <= length(stress) - 1
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain(i+1) - strain(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 

area_toughness = area_toughness + area_trap;
d_u_stress_t = u_stress(i) .* (strain(i+1) - strain(i));
U_toughness = d_u_stress_t + U_toughness;

i = i + 1; 
end

% True Uncertainty

rupture_index = length(true_stress);
U_Stress_True_rupture = U_Stress_true(rupture_index);
U_Strain_True_rupture = U_strain_true(rupture_index);

fprintf('Area Resilience: %g +- %g\n',area_resilience, U_resilience)
fprintf('Area Toughness: %g +- %g\n',area_toughness, U_toughness)

fprintf('True Rupture Stress: %g\n',True_Rupture_stress)
fprintf('True Rupture Strain: %g\n',True_Rupture_strain)

fprintf('True Ruptur Stress: %g +- %g\n',True_Rupture_stress, U_Stress_True_rupture)
fprintf('True Ruptur Strain: %g +- %g\n',True_Rupture_strain, U_Strain_True_rupture)





% MMTS

stress_yg_MMTS = stress(1:1:50);
strain_yg_MMTS = strain(1:1:50);


estima = fitlm(strain_yg_MMTS,stress_yg_MMTS);
x1=estima.Coefficients;
youngs_modulus_MMTS=x1{'x1', 'Estimate'};


strain_m_MMTS = strain * 1.002;
strain_m_MMTS = strain_m_MMTS(1:100);
vertical_offset_MMTS = stress(1) - youngs_modulus_MMTS .* strain_m_MMTS(1);

stress_m_MMTS = zeros(100,1);
i = 1;
while i <= 100
stress_m_MMTS(i) = youngs_modulus_MMTS .* strain_m_MMTS(i) + vertical_offset_MMTS;
i = i+1;
end

stress_test = stress (1:100);

index_MMTS = find(stress_m_MMTS > stress_test, 1);
yield_stress_MMTS = stress(index_MMTS);
yield_strain_MMTS = strain(index_MMTS);


oset_plot_MMTS = stress_m_MMTS(1:index_MMTS);
strain_oset_MMTS = strain_m_MMTS(1:index_MMTS);


i = 1;
area_resilience_MMTS = 0;
U_resilience_MMTS = 0;
h1 = 0;
h2 = 0;
b = 0;
area_rect = 0;
area_trap = 0;

while i <= index_MMTS
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain_m_MMTS(i+1) - strain_m_MMTS(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 
area_resilience_MMTS = area_resilience_MMTS + area_trap;

d_u_stress_r_MMTS = u_stress(i) .* (strain_m_MMTS(i+1) - strain_m_MMTS(i));
U_resilience_MMTS = d_u_stress_r_MMTS + U_resilience_MMTS;

i = i + 1; 
end
temp = eval(u_young_MMTS);
temp2 = temp(1:index_MMTS);
u_young_MMTS = mean(temp2)
% Exten
stress_yg_exten = stress(1:1:50);
s_E = info(:,3);
strain_yg_exten = s_E(1:1:50);


estima = fitlm(strain_yg_exten,stress_yg_exten);
x1=estima.Coefficients;
youngs_modulus_extensometer=x1{'x1', 'Estimate'};


strain_m_exten = strain * 1.002;
strain_m_exten = strain_m_exten(1:100);
vertical_offset_exten = stress(1) - youngs_modulus_extensometer .* strain_m_exten(1);

stress_m_exten = zeros(100,1);
i = 1;
while i <= 100
stress_m_exten(i) = youngs_modulus_extensometer .* strain_m_exten(i) + vertical_offset_exten;
i = i+1;
end

stress_test = stress (1:100);

index_exten = find(stress_m_exten > stress_test, 1);
yield_stress_Extensometer = stress(index_exten);
yield_strain_Extensometer = strain(index_exten);


oset_plot_exten = stress_m_exten(1:index_exten);
strain_oset_exten = strain_m_exten(1:index_exten);

i = 1;
area_resilience_EXTEN = 0;
U_resilience_EXTEN = 0;
h1 = 0;
h2 = 0;
b = 0;
area_rect = 0;
area_trap = 0;

while i <= index_exten
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain_m_exten(i+1) - strain_m_exten(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 
area_resilience_EXTEN = area_resilience_EXTEN + area_trap;

d_u_stress_r_EXTEN = u_stress(i) .* (strain_m_exten(i+1) - strain_m_exten(i));
U_resilience_EXTEN = d_u_stress_r_EXTEN + U_resilience_EXTEN;

i = i + 1; 
end
temp3 = eval(u_young_exten);
temp4 = temp3(1:index_exten);
u_young_Exten = mean(temp4)

fprintf('Yield Stress: %g\n',yield_stress_MMTS)
fprintf('Yield Strain: %g\n',yield_strain_MMTS)
fprintf('Area Reduction Percent: %g +- %g\n',Area_reduction, U_Area_reduction)

% poisson Ratio


lateral_strain = ((tf-t0)/t0);
poissons = -(lateral_strain)/(strain(end));
u_poisson = eval(u_poisson);
uid = length(u_poisson);
U_poissons = u_poisson(uid);
fprintf('Poissons Ratio: %g +- %g\n',poissons, U_poissons)

% plotting 

figure(4)
plot(strain,stress,'blue')
hold on
plot(yield_x,yield_y,'diamond','MarkerSize',7,'MarkerFaceColor',[0,0,1])
hold on
plot(strain(start_index:yield_intercept),offset(start_index:yield_intercept),'green')
hold on
plot(ultimite_strain,ultimite_stress,'diamond','MarkerSize',7,'MarkerFaceColor',[0,1,0])
hold on
plot(Rupture_strain,Rupture_stress,'diamond','MarkerSize',7,'MarkerFaceColor',[1,0,0])
title('strain vs stress of Carbon Fiber 45')
xlabel('strain')
ylabel('stress(lb/in^2)')
errorbar(strain(1:40:end),stress(1:40:end),u_stress(1:40:end),'vertical','.') % vertical errorbars for stress
errorbar(strain(1:40:end),stress(1:40:end),u_strain(1:40:end),'horizontal','.') % horizontal errorbars for strain
plot(true_strain, true_stress,'red')
hold on
plot(True_Rupture_strain,True_Rupture_stress,'diamond','MarkerSize',7,'MarkerFaceColor',[0.5,0.5,0])
hold on
plot(strain_oset_MMTS,oset_plot_MMTS,'green')
hold on
%plot(yield_strain_Extensometer,yield_stress_Extensometer,'diamond','MarkerSize',7,'MarkerFaceColor',[1,1,1])
hold on
%plot(strain_oset_exten,oset_plot_exten)
hold on
plot(yield_strain_MMTS,yield_stress_MMTS,'diamond','MarkerSize',7,'MarkerFaceColor',[0,1,1])
legend('strain vs stress curve','ultimite stress and strain', 'rupture point', ...
   'uncertianty of stress','uncertianty of strain','true strain vs stress curve','true rupture','0.2% offset line','yield point','Location','southeast')

%fprintf('Youngs Modulus: %g +- %g\n',youngs_modulus_MMTS, u_E)

%% Carbon Fiber 90
clear
% DATA
syms d F w0 t0 u_c u_tm u_d u_f L tf wf strain_yg_exten u_extens F2 u_f_2
% d= displacment
% F= force
% t0 = thickness
% w0 = width
% wf = final width
% tf = final thickness
% u_c = caliper uncertianty
% u_tm = tapemeasure uncertianty
% u_d = displacment uncertianty
% u_f = force uncertianty
% L = length

A = t0.*w0; % area
stress = F./A;
strain = d./L;
Af = tf.*wf ; % final area
Area_reduction = (1 - (Af./A)).*100 ; % area reduction percent 
E = stress./strain;
stress2 = F2 ./A;
yung = stress2./ strain_yg_exten;
Q = stress .* (1 + strain); % true stress
u_stress = sqrt((diff(stress,F).*u_f).^2 + (diff(stress,t0).*u_c).^2 + (diff(stress,w0).*u_c).^2);
u_strain = sqrt((diff(strain,d).*u_d).^2 + (diff(strain,L).*u_tm).^2);
u_E = sqrt((diff(E,F).*u_f).^2 + (diff(E,t0).*u_c).^2 + (diff(E,w0).*u_c).^2 + ...
    (diff(E,d).*u_d).^2 + (diff(E,L).*u_tm).^2);
u_young_MMTS = sqrt((diff(E,F).*u_f).^2 + (diff(E,t0).*u_c).^2 + (diff(E,w0).*u_c).^2 + ...
    (diff(E,d).*u_d).^2 + (diff(E,L).*u_tm).^2);

u_young_exten = sqrt((diff(yung,F2).*u_f_2).^2 + (diff(yung,t0).* u_c).^2 + (diff(yung,w0).*u_c).^2 + (diff(yung,strain_yg_exten) .* u_extens).^2);


U_Stress_true = sqrt((diff(Q,F).*u_f).^2 + (diff(Q,t0).*u_c).^2 + (diff(Q,w0).*u_c).^2 + ...
    (diff(Q,d).*u_d).^2 + (diff(Q,L).*u_tm).^2);
U_Area_reduction = sqrt((diff(Af,tf).*u_c).^2 + (diff(Af,wf).*u_c).^2 + (diff(A,t0).*u_c).^2 + (diff(A,w0).*u_c).^2);
poissons = -((tf- t0)./t0)./(d./L);
u_poisson = sqrt((diff(poissons,t0) .* u_c).^2 + (diff(poissons,tf) .* u_c).^2 + (diff(poissons,d) .* u_d).^2 + (diff(poissons,L) .* u_tm).^2);
% manual differentiation
%u_stress = sqrt(((1/(t0.*w0)).*u_f).^2 + ((avgF/((t0.^2).*w0).*u_t).^2) + (-avgF/(t0.*w0.^2).*u_w).^2);
%u_strain = sqrt((((1./L).*u_d).^2)+(((avgd./(L.^2)).*u_L).^2));
disp('Carbon Fiber 90')
%Data input
rawdata = importdata('M008_CF90.dat','\t',5);
info = rawdata.data;
d = info(:,2); % displacment data
F = info(:,1); % force data
F2 = F(1:50);
s_E = info(:,3); % strain data for E(not used in this but avaliable)

L=6.5; %in change these values for each sample
w0 = 0.510; %in change these values for each sample
t0 = 0.072; %in change these values for each sample
Lf = 6.5; % final length 
wf = 0.510; %final width
tf = 0.064;  % final thickness (change)

u_extens = 5 .* 10 .^-8; % Reading on axial strain is x e^-7, half of that
u_c = 0.0005; %in 
u_tm = 1/64; %in
u_d = 0.012 .* d;
u_f = 0.012 .* F;
u_f_2 = .012 .*F2;
A = eval(A);
Af = eval(Af);
u_stress = eval(u_stress);
u_strain = eval(u_strain);

stress = eval(stress);
strain = eval(strain);
U_strain_true = u_strain;
U_Stress_true = eval(U_Stress_true);
Area_reduction = eval(Area_reduction);
U_Area_reduction = eval(U_Area_reduction);

% Area reduction

% t0 = intial thickness
% w0 = intial thickness
% tf = final thickness 
% L = intial length
% wf =final thickness 




% Finding Values

% Finding ultimite stress and strain
ultimite_stress = max(stress); % could also be stress(find(max(stress) == stress))
ultimite_strain = strain(find(max(stress) == stress));

%Rupture stress and strain
Rupture_stress = stress(end);
Rupture_strain = strain(end);

%yield point 0.2 percent offset and Youngs Modulus(E)

p = polyfit(strain(1:50),stress(1:50),1);
E = p(1);
b = stress(1)-E.*(strain(1) .* 1.002);
p(2) = b;
offset = polyval(p,strain(1));
yield_intercept = find(offset >= stress,1);
yield_x = strain(yield_intercept);
yield_y = stress(yield_intercept);
start_index = find(offset >= stress(1),1);

% true values 
true_stress = stress .* (1+strain);
true_strain = log(1+strain);

True_Rupture_stress = true_stress(end);
True_Rupture_strain = true_strain(end);



fprintf('Ultimate Stress: %g +- %g\n',ultimite_stress, u_stress(find(max(stress) == stress)))
fprintf('Ultimate Strain: %g +- %g\n',ultimite_strain, u_strain(find(max(stress) == stress)))
fprintf('Rupture Stress: %g +- %g\n',Rupture_stress, u_stress(end))
fprintf('Rupture Strain: %g +- %g\n',Rupture_strain, u_strain(end))
fprintf('Yield Point: (%g,%g)\n',yield_x,yield_y)



% Resiliance and toughness

i = 1;
area_resilience = 0;
U_resilience = 0;
while i <=50
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain(i+1) - strain(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 
area_resilience = area_resilience + area_trap;

d_u_stress_r = u_stress(i) .* (strain(i+1) - strain(i));
U_resilience = d_u_stress_r + U_resilience;

i = i + 1; 
end

i = 1;
area_toughness = 0;
U_toughness = 0;

while i <= length(stress) - 1
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain(i+1) - strain(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 

area_toughness = area_toughness + area_trap;
d_u_stress_t = u_stress(i) .* (strain(i+1) - strain(i));
U_toughness = d_u_stress_t + U_toughness;

i = i + 1; 
end

% True Uncertainty

rupture_index = length(true_stress);
U_Stress_True_rupture = U_Stress_true(rupture_index);
U_Strain_True_rupture = U_strain_true(rupture_index);

fprintf('Area Resilience: %g +- %g\n',area_resilience, U_resilience)
fprintf('Area Toughness: %g +- %g\n',area_toughness, U_toughness)

fprintf('True Rupture Stress: %g\n',True_Rupture_stress)
fprintf('True Rupture Strain: %g\n',True_Rupture_strain)

fprintf('True Ruptur Stress: %g +- %g\n',True_Rupture_stress, U_Stress_True_rupture)
fprintf('True Ruptur Strain: %g +- %g\n',True_Rupture_strain, U_Strain_True_rupture)





% MMTS

stress_yg_MMTS = stress(1:1:50);
strain_yg_MMTS = strain(1:1:50);


estima = fitlm(strain_yg_MMTS,stress_yg_MMTS);
x1=estima.Coefficients;
youngs_modulus_MMTS=x1{'x1', 'Estimate'};


strain_m_MMTS = strain * 1.002;
strain_m_MMTS = strain_m_MMTS(1:100);
vertical_offset_MMTS = stress(1) - youngs_modulus_MMTS .* strain_m_MMTS(1);

stress_m_MMTS = zeros(100,1);
i = 1;
while i <= 100
stress_m_MMTS(i) = youngs_modulus_MMTS .* strain_m_MMTS(i) + vertical_offset_MMTS;
i = i+1;
end

stress_test = stress (1:100);

index_MMTS = find(stress_m_MMTS > stress_test, 1);
yield_stress_MMTS = stress(index_MMTS);
yield_strain_MMTS = strain(index_MMTS);


oset_plot_MMTS = stress_m_MMTS(1:index_MMTS);
strain_oset_MMTS = strain_m_MMTS(1:index_MMTS);


i = 1;
area_resilience_MMTS = 0;
U_resilience_MMTS = 0;
h1 = 0;
h2 = 0;
b = 0;
area_rect = 0;
area_trap = 0;

while i <= index_MMTS
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain_m_MMTS(i+1) - strain_m_MMTS(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 
area_resilience_MMTS = area_resilience_MMTS + area_trap;

d_u_stress_r_MMTS = u_stress(i) .* (strain_m_MMTS(i+1) - strain_m_MMTS(i));
U_resilience_MMTS = d_u_stress_r_MMTS + U_resilience_MMTS;

i = i + 1; 
end
temp = eval(u_young_MMTS);
temp2 = temp(1:index_MMTS);
u_young_MMTS = mean(temp2)
% Exten
stress_yg_exten = stress(1:1:50);
s_E = info(:,3);
strain_yg_exten = s_E(1:1:50);


estima = fitlm(strain_yg_exten,stress_yg_exten);
x1=estima.Coefficients;
youngs_modulus_extensometer=x1{'x1', 'Estimate'};


strain_m_exten = strain * 1.002;
strain_m_exten = strain_m_exten(1:100);
vertical_offset_exten = stress(1) - youngs_modulus_extensometer .* strain_m_exten(1);

stress_m_exten = zeros(100,1);
i = 1;
while i <= 100
stress_m_exten(i) = youngs_modulus_extensometer .* strain_m_exten(i) + vertical_offset_exten;
i = i+1;
end

stress_test = stress (1:100);

index_exten = find(stress_m_exten > stress_test, 1);
yield_stress_Extensometer = stress(index_exten);
yield_strain_Extensometer = strain(index_exten);


oset_plot_exten = stress_m_exten(1:index_exten);
strain_oset_exten = strain_m_exten(1:index_exten);

i = 1;
area_resilience_EXTEN = 0;
U_resilience_EXTEN = 0;
h1 = 0;
h2 = 0;
b = 0;
area_rect = 0;
area_trap = 0;

while i <= index_exten
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain_m_exten(i+1) - strain_m_exten(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 
area_resilience_EXTEN = area_resilience_EXTEN + area_trap;

d_u_stress_r_EXTEN = u_stress(i) .* (strain_m_exten(i+1) - strain_m_exten(i));
U_resilience_EXTEN = d_u_stress_r_EXTEN + U_resilience_EXTEN;

i = i + 1; 
end
temp3 = eval(u_young_exten);
temp4 = temp3(1:index_exten);
u_young_Exten = mean(temp4)

fprintf('Yield Stress: %g\n',yield_stress_MMTS)
fprintf('Yield Strain: %g\n',yield_strain_MMTS)
fprintf('Area Reduction Percent: %g +- %g\n',Area_reduction, U_Area_reduction)

% poisson Ratio


lateral_strain = ((tf-t0)/t0);
poissons = -(lateral_strain)/(strain(end));
u_poisson = eval(u_poisson);
uid = length(u_poisson);
U_poissons = u_poisson(uid);
fprintf('Poissons Ratio: %g +- %g\n',poissons, U_poissons)

% plotting 

figure(5)
plot(strain,stress,'blue')
hold on
plot(yield_x,yield_y,'diamond','MarkerSize',7,'MarkerFaceColor',[0,0,1])
hold on
plot(strain(start_index:yield_intercept),offset(start_index:yield_intercept),'green')
hold on
plot(ultimite_strain,ultimite_stress,'diamond','MarkerSize',7,'MarkerFaceColor',[0,1,0])
hold on
plot(Rupture_strain,Rupture_stress,'diamond','MarkerSize',7,'MarkerFaceColor',[1,0,0])
title('strain vs stress of Carbon Fiber 90')
xlabel('strain')
ylabel('stress(lb/in^2)')
errorbar(strain(1:40:end),stress(1:40:end),u_stress(1:40:end),'vertical','.') % vertical errorbars for stress
errorbar(strain(1:40:end),stress(1:40:end),u_strain(1:40:end),'horizontal','.') % horizontal errorbars for strain
plot(true_strain, true_stress,'red')
hold on
plot(True_Rupture_strain,True_Rupture_stress,'diamond','MarkerSize',7,'MarkerFaceColor',[0.5,0.5,0])
hold on
plot(strain_oset_MMTS,oset_plot_MMTS,'green')
hold on
%plot(yield_strain_Extensometer,yield_stress_Extensometer,'diamond','MarkerSize',7,'MarkerFaceColor',[1,1,1])
hold on
%plot(strain_oset_exten,oset_plot_exten)
hold on
plot(yield_strain_MMTS,yield_stress_MMTS,'diamond','MarkerSize',7,'MarkerFaceColor',[0,1,1])
legend('strain vs stress curve','ultimite stress and strain', 'rupture point', ...
   'uncertianty of stress','uncertianty of strain','true strain vs stress curve','true rupture','0.2% offset line','yield point','Location','southeast')

%fprintf('Youngs Modulus: %g +- %g\n',youngs_modulus_MMTS, u_E)


%% Plastic 1
clear
% DATA
syms d F w0 t0 u_c u_tm u_d u_f L tf wf strain_yg_exten u_extens F2 u_f_2
% d= displacment
% F= force
% t0 = thickness
% w0 = width
% wf = final width
% tf = final thickness
% u_c = caliper uncertianty
% u_tm = tapemeasure uncertianty
% u_d = displacment uncertianty
% u_f = force uncertianty
% L = length

A = t0.*w0; % area
stress = F./A;
strain = d./L;
Af = tf.*wf ; % final area
Area_reduction = (1 - (Af./A)).*100 ; % area reduction percent 
E = stress./strain;
stress2 = F2 ./A;
yung = stress2./ strain_yg_exten;
Q = stress .* (1 + strain); % true stress
u_stress = sqrt((diff(stress,F).*u_f).^2 + (diff(stress,t0).*u_c).^2 + (diff(stress,w0).*u_c).^2);
u_strain = sqrt((diff(strain,d).*u_d).^2 + (diff(strain,L).*u_tm).^2);
u_E = sqrt((diff(E,F).*u_f).^2 + (diff(E,t0).*u_c).^2 + (diff(E,w0).*u_c).^2 + ...
    (diff(E,d).*u_d).^2 + (diff(E,L).*u_tm).^2);
u_young_MMTS = sqrt((diff(E,F).*u_f).^2 + (diff(E,t0).*u_c).^2 + (diff(E,w0).*u_c).^2 + ...
    (diff(E,d).*u_d).^2 + (diff(E,L).*u_tm).^2);

u_young_exten = sqrt((diff(yung,F2).*u_f_2).^2 + (diff(yung,t0).* u_c).^2 + (diff(yung,w0).*u_c).^2 + (diff(yung,strain_yg_exten) .* u_extens).^2);


U_Stress_true = sqrt((diff(Q,F).*u_f).^2 + (diff(Q,t0).*u_c).^2 + (diff(Q,w0).*u_c).^2 + ...
    (diff(Q,d).*u_d).^2 + (diff(Q,L).*u_tm).^2);
U_Area_reduction = sqrt((diff(Af,tf).*u_c).^2 + (diff(Af,wf).*u_c).^2 + (diff(A,t0).*u_c).^2 + (diff(A,w0).*u_c).^2);
poissons = -((tf- t0)./t0)./(d./L);
u_poisson = sqrt((diff(poissons,t0) .* u_c).^2 + (diff(poissons,tf) .* u_c).^2 + (diff(poissons,d) .* u_d).^2 + (diff(poissons,L) .* u_tm).^2);
% manual differentiation
%u_stress = sqrt(((1/(t0.*w0)).*u_f).^2 + ((avgF/((t0.^2).*w0).*u_t).^2) + (-avgF/(t0.*w0.^2).*u_w).^2);
%u_strain = sqrt((((1./L).*u_d).^2)+(((avgd./(L.^2)).*u_L).^2));
disp('Plastic 1')
%Data input
rawdata = importdata('M008_Plastic_1.dat','\t',5);
info = rawdata.data;
d = info(:,2); % displacment data
F = info(:,1); % force data
F2 = F(1:50);
s_E = info(:,3); % strain data for E(not used in this but avaliable)

L=6.5; %in change these values for each sample
w0 = 0.506; %in change these values for each sample
t0 = 0.137; %in change these values for each sample
Lf = 6.625; % final length 
wf = 0.505; %final width
tf = 0.135;  % final thickness (change)

u_extens = 5 .* 10 .^-8; % Reading on axial strain is x e^-7, half of that
u_c = 0.0005; %in 
u_tm = 1/64; %in
u_d = 0.012 .* d;
u_f = 0.012 .* F;
u_f_2 = .012 .*F2;
A = eval(A);
Af = eval(Af);
u_stress = eval(u_stress);
u_strain = eval(u_strain);

stress = eval(stress);
strain = eval(strain);
U_strain_true = u_strain;
U_Stress_true = eval(U_Stress_true);
Area_reduction = eval(Area_reduction);
U_Area_reduction = eval(U_Area_reduction);

% Area reduction

% t0 = intial thickness
% w0 = intial thickness
% tf = final thickness 
% L = intial length
% wf =final thickness 




% Finding Values

% Finding ultimite stress and strain
ultimite_stress = max(stress); % could also be stress(find(max(stress) == stress))
ultimite_strain = strain(find(max(stress) == stress));

%Rupture stress and strain
Rupture_stress = stress(end);
Rupture_strain = strain(end);

%yield point 0.2 percent offset and Youngs Modulus(E)

p = polyfit(strain(1:50),stress(1:50),1);
E = p(1);
b = stress(1)-E.*(strain(1) .* 1.002);
p(2) = b;
offset = polyval(p,strain(1));
yield_intercept = find(offset >= stress,1);
yield_x = strain(yield_intercept);
yield_y = stress(yield_intercept);
start_index = find(offset >= stress(1),1);

% true values 
true_stress = stress .* (1+strain);
true_strain = log(1+strain);

True_Rupture_stress = true_stress(end);
True_Rupture_strain = true_strain(end);



fprintf('Ultimate Stress: %g +- %g\n',ultimite_stress, u_stress(find(max(stress) == stress)))
fprintf('Ultimate Strain: %g +- %g\n',ultimite_strain, u_strain(find(max(stress) == stress)))
fprintf('Rupture Stress: %g +- %g\n',Rupture_stress, u_stress(end))
fprintf('Rupture Strain: %g +- %g\n',Rupture_strain, u_strain(end))
fprintf('Yield Point: (%g,%g)\n',yield_x,yield_y)



% Resiliance and toughness

i = 1;
area_resilience = 0;
U_resilience = 0;
while i <=50
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain(i+1) - strain(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 
area_resilience = area_resilience + area_trap;

d_u_stress_r = u_stress(i) .* (strain(i+1) - strain(i));
U_resilience = d_u_stress_r + U_resilience;

i = i + 1; 
end

i = 1;
area_toughness = 0;
U_toughness = 0;

while i <= length(stress) - 1
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain(i+1) - strain(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 

area_toughness = area_toughness + area_trap;
d_u_stress_t = u_stress(i) .* (strain(i+1) - strain(i));
U_toughness = d_u_stress_t + U_toughness;

i = i + 1; 
end

% True Uncertainty

rupture_index = length(true_stress);
U_Stress_True_rupture = U_Stress_true(rupture_index);
U_Strain_True_rupture = U_strain_true(rupture_index);

fprintf('Area Resilience: %g +- %g\n',area_resilience, U_resilience)
fprintf('Area Toughness: %g +- %g\n',area_toughness, U_toughness)

fprintf('True Rupture Stress: %g\n',True_Rupture_stress)
fprintf('True Rupture Strain: %g\n',True_Rupture_strain)

fprintf('True Ruptur Stress: %g +- %g\n',True_Rupture_stress, U_Stress_True_rupture)
fprintf('True Ruptur Strain: %g +- %g\n',True_Rupture_strain, U_Strain_True_rupture)





% MMTS

stress_yg_MMTS = stress(1:1:50);
strain_yg_MMTS = strain(1:1:50);


estima = fitlm(strain_yg_MMTS,stress_yg_MMTS);
x1=estima.Coefficients;
youngs_modulus_MMTS=x1{'x1', 'Estimate'};


strain_m_MMTS = strain * 1.002;
strain_m_MMTS = strain_m_MMTS(1:100);
vertical_offset_MMTS = stress(1) - youngs_modulus_MMTS .* strain_m_MMTS(1);

stress_m_MMTS = zeros(100,1);
i = 1;
while i <= 100
stress_m_MMTS(i) = youngs_modulus_MMTS .* strain_m_MMTS(i) + vertical_offset_MMTS;
i = i+1;
end

stress_test = stress (1:100);

index_MMTS = find(stress_m_MMTS > stress_test, 1);
yield_stress_MMTS = stress(index_MMTS);
yield_strain_MMTS = strain(index_MMTS);


oset_plot_MMTS = stress_m_MMTS(1:index_MMTS);
strain_oset_MMTS = strain_m_MMTS(1:index_MMTS);


i = 1;
area_resilience_MMTS = 0;
U_resilience_MMTS = 0;
h1 = 0;
h2 = 0;
b = 0;
area_rect = 0;
area_trap = 0;

while i <= index_MMTS
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain_m_MMTS(i+1) - strain_m_MMTS(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 
area_resilience_MMTS = area_resilience_MMTS + area_trap;

d_u_stress_r_MMTS = u_stress(i) .* (strain_m_MMTS(i+1) - strain_m_MMTS(i));
U_resilience_MMTS = d_u_stress_r_MMTS + U_resilience_MMTS;

i = i + 1; 
end
temp = eval(u_young_MMTS);
temp2 = temp(1:index_MMTS);
u_young_MMTS = mean(temp2)
% Exten
stress_yg_exten = stress(1:1:50);
s_E = info(:,3);
strain_yg_exten = s_E(1:1:50);


estima = fitlm(strain_yg_exten,stress_yg_exten);
x1=estima.Coefficients;
youngs_modulus_extensometer=x1{'x1', 'Estimate'};


strain_m_exten = strain * 1.002;
strain_m_exten = strain_m_exten(1:100);
vertical_offset_exten = stress(1) - youngs_modulus_extensometer .* strain_m_exten(1);

stress_m_exten = zeros(100,1);
i = 1;
while i <= 100
stress_m_exten(i) = youngs_modulus_extensometer .* strain_m_exten(i) + vertical_offset_exten;
i = i+1;
end

stress_test = stress (1:100);

index_exten = find(stress_m_exten > stress_test, 1);
yield_stress_Extensometer = stress(index_exten);
yield_strain_Extensometer = strain(index_exten);


oset_plot_exten = stress_m_exten(1:index_exten);
strain_oset_exten = strain_m_exten(1:index_exten);

i = 1;
area_resilience_EXTEN = 0;
U_resilience_EXTEN = 0;
h1 = 0;
h2 = 0;
b = 0;
area_rect = 0;
area_trap = 0;

while i <= index_exten
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain_m_exten(i+1) - strain_m_exten(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 
area_resilience_EXTEN = area_resilience_EXTEN + area_trap;

d_u_stress_r_EXTEN = u_stress(i) .* (strain_m_exten(i+1) - strain_m_exten(i));
U_resilience_EXTEN = d_u_stress_r_EXTEN + U_resilience_EXTEN;

i = i + 1; 
end
temp3 = eval(u_young_exten);
temp4 = temp3(1:index_exten);
u_young_Exten = mean(temp4)

fprintf('Yield Stress: %g\n',yield_stress_MMTS)
fprintf('Yield Strain: %g\n',yield_strain_MMTS)
fprintf('Area Reduction Percent: %g +- %g\n',Area_reduction, U_Area_reduction)

% poisson Ratio


lateral_strain = ((tf-t0)/t0);
poissons = -(lateral_strain)/(strain(end));
u_poisson = eval(u_poisson);
uid = length(u_poisson);
U_poissons = u_poisson(uid);
fprintf('Poissons Ratio: %g +- %g\n',poissons, U_poissons)

% plotting 

figure(6)
plot(strain,stress,'blue')
hold on
plot(yield_x,yield_y,'diamond','MarkerSize',7,'MarkerFaceColor',[0,0,1])
hold on
plot(strain(start_index:yield_intercept),offset(start_index:yield_intercept),'green')
hold on
plot(ultimite_strain,ultimite_stress,'diamond','MarkerSize',7,'MarkerFaceColor',[0,1,0])
hold on
plot(Rupture_strain,Rupture_stress,'diamond','MarkerSize',7,'MarkerFaceColor',[1,0,0])
title('strain vs Plastic 1')
xlabel('strain')
ylabel('stress(lb/in^2)')
errorbar(strain(1:40:end),stress(1:40:end),u_stress(1:40:end),'vertical','.') % vertical errorbars for stress
errorbar(strain(1:40:end),stress(1:40:end),u_strain(1:40:end),'horizontal','.') % horizontal errorbars for strain
plot(true_strain, true_stress,'red')
hold on
plot(True_Rupture_strain,True_Rupture_stress,'diamond','MarkerSize',7,'MarkerFaceColor',[0.5,0.5,0])
hold on
plot(strain_oset_MMTS,oset_plot_MMTS,'green')
hold on
%plot(yield_strain_Extensometer,yield_stress_Extensometer,'diamond','MarkerSize',7,'MarkerFaceColor',[1,1,1])
hold on
%plot(strain_oset_exten,oset_plot_exten)
hold on
plot(yield_strain_MMTS,yield_stress_MMTS,'diamond','MarkerSize',7,'MarkerFaceColor',[0,1,1])
legend('strain vs stress curve','ultimite stress and strain', 'rupture point', ...
   'uncertianty of stress','uncertianty of strain','true strain vs stress curve','true rupture','0.2% offset line','yield point','Location','southeast')

%fprintf('Youngs Modulus: %g +- %g\n',youngs_modulus_MMTS, u_E)

%% Plastic 2
clear
% DATA
syms d F w0 t0 u_c u_tm u_d u_f L tf wf strain_yg_exten u_extens F2 u_f_2
% d= displacment
% F= force
% t0 = thickness
% w0 = width
% wf = final width
% tf = final thickness
% u_c = caliper uncertianty
% u_tm = tapemeasure uncertianty
% u_d = displacment uncertianty
% u_f = force uncertianty
% L = length

A = t0.*w0; % area
stress = F./A;
strain = d./L;
Af = tf.*wf ; % final area
Area_reduction = (1 - (Af./A)).*100 ; % area reduction percent 
E = stress./strain;
stress2 = F2 ./A;
yung = stress2./ strain_yg_exten;
Q = stress .* (1 + strain); % true stress
u_stress = sqrt((diff(stress,F).*u_f).^2 + (diff(stress,t0).*u_c).^2 + (diff(stress,w0).*u_c).^2);
u_strain = sqrt((diff(strain,d).*u_d).^2 + (diff(strain,L).*u_tm).^2);
u_E = sqrt((diff(E,F).*u_f).^2 + (diff(E,t0).*u_c).^2 + (diff(E,w0).*u_c).^2 + ...
    (diff(E,d).*u_d).^2 + (diff(E,L).*u_tm).^2);
u_young_MMTS = sqrt((diff(E,F).*u_f).^2 + (diff(E,t0).*u_c).^2 + (diff(E,w0).*u_c).^2 + ...
    (diff(E,d).*u_d).^2 + (diff(E,L).*u_tm).^2);

u_young_exten = sqrt((diff(yung,F2).*u_f_2).^2 + (diff(yung,t0).* u_c).^2 + (diff(yung,w0).*u_c).^2 + (diff(yung,strain_yg_exten) .* u_extens).^2);


U_Stress_true = sqrt((diff(Q,F).*u_f).^2 + (diff(Q,t0).*u_c).^2 + (diff(Q,w0).*u_c).^2 + ...
    (diff(Q,d).*u_d).^2 + (diff(Q,L).*u_tm).^2);
U_Area_reduction = sqrt((diff(Af,tf).*u_c).^2 + (diff(Af,wf).*u_c).^2 + (diff(A,t0).*u_c).^2 + (diff(A,w0).*u_c).^2);
poissons = -((tf- t0)./t0)./(d./L);
u_poisson = sqrt((diff(poissons,t0) .* u_c).^2 + (diff(poissons,tf) .* u_c).^2 + (diff(poissons,d) .* u_d).^2 + (diff(poissons,L) .* u_tm).^2);
% manual differentiation
%u_stress = sqrt(((1/(t0.*w0)).*u_f).^2 + ((avgF/((t0.^2).*w0).*u_t).^2) + (-avgF/(t0.*w0.^2).*u_w).^2);
%u_strain = sqrt((((1./L).*u_d).^2)+(((avgd./(L.^2)).*u_L).^2));
disp('Plastic 2')
%Data input
rawdata = importdata('M008_22_plastic2.dat','\t',5);
info = rawdata.data;
d = info(:,2); % displacment data
F = info(:,1); % force data
F2 = F(1:50);
s_E = info(:,3); % strain data for E(not used in this but avaliable)

L=6.5; %in change these values for each sample
w0 = 0.506; %in change these values for each sample
t0 = 0.137; %in change these values for each sample
Lf = 6.625; % final length 
wf = 0.505; %final width
tf = 0.135;  % final thickness (change)

u_extens = 5 .* 10 .^-8; % Reading on axial strain is x e^-7, half of that
u_c = 0.0005; %in 
u_tm = 1/64; %in
u_d = 0.012 .* d;
u_f = 0.012 .* F;
u_f_2 = .012 .*F2;
A = eval(A);
Af = eval(Af);
u_stress = eval(u_stress);
u_strain = eval(u_strain);

stress = eval(stress);
strain = eval(strain);
U_strain_true = u_strain;
U_Stress_true = eval(U_Stress_true);
Area_reduction = eval(Area_reduction);
U_Area_reduction = eval(U_Area_reduction);

% Area reduction

% t0 = intial thickness
% w0 = intial thickness
% tf = final thickness 
% L = intial length
% wf =final thickness 




% Finding Values

% Finding ultimite stress and strain
ultimite_stress = max(stress); % could also be stress(find(max(stress) == stress))
ultimite_strain = strain(find(max(stress) == stress));

%Rupture stress and strain
Rupture_stress = stress(end);
Rupture_strain = strain(end);

%yield point 0.2 percent offset and Youngs Modulus(E)

p = polyfit(strain(1:50),stress(1:50),1);
E = p(1);
b = stress(1)-E.*(strain(1) .* 1.002);
p(2) = b;
offset = polyval(p,strain(1));
yield_intercept = find(offset >= stress,1);
yield_x = strain(yield_intercept);
yield_y = stress(yield_intercept);
start_index = find(offset >= stress(1),1);

% true values 
true_stress = stress .* (1+strain);
true_strain = log(1+strain);

True_Rupture_stress = true_stress(end);
True_Rupture_strain = true_strain(end);



fprintf('Ultimate Stress: %g +- %g\n',ultimite_stress, u_stress(find(max(stress) == stress)))
fprintf('Ultimate Strain: %g +- %g\n',ultimite_strain, u_strain(find(max(stress) == stress)))
fprintf('Rupture Stress: %g +- %g\n',Rupture_stress, u_stress(end))
fprintf('Rupture Strain: %g +- %g\n',Rupture_strain, u_strain(end))
fprintf('Yield Point: (%g,%g)\n',yield_x,yield_y)



% Resiliance and toughness

i = 1;
area_resilience = 0;
U_resilience = 0;
while i <=50
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain(i+1) - strain(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 
area_resilience = area_resilience + area_trap;

d_u_stress_r = u_stress(i) .* (strain(i+1) - strain(i));
U_resilience = d_u_stress_r + U_resilience;

i = i + 1; 
end

i = 1;
area_toughness = 0;
U_toughness = 0;

while i <= length(stress) - 1
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain(i+1) - strain(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 

area_toughness = area_toughness + area_trap;
d_u_stress_t = u_stress(i) .* (strain(i+1) - strain(i));
U_toughness = d_u_stress_t + U_toughness;

i = i + 1; 
end

% True Uncertainty

rupture_index = length(true_stress);
U_Stress_True_rupture = U_Stress_true(rupture_index);
U_Strain_True_rupture = U_strain_true(rupture_index);

fprintf('Area Resilience: %g +- %g\n',area_resilience, U_resilience)
fprintf('Area Toughness: %g +- %g\n',area_toughness, U_toughness)

fprintf('True Rupture Stress: %g\n',True_Rupture_stress)
fprintf('True Rupture Strain: %g\n',True_Rupture_strain)

fprintf('True Ruptur Stress: %g +- %g\n',True_Rupture_stress, U_Stress_True_rupture)
fprintf('True Ruptur Strain: %g +- %g\n',True_Rupture_strain, U_Strain_True_rupture)





% MMTS

stress_yg_MMTS = stress(1:1:50);
strain_yg_MMTS = strain(1:1:50);


estima = fitlm(strain_yg_MMTS,stress_yg_MMTS);
x1=estima.Coefficients;
youngs_modulus_MMTS=x1{'x1', 'Estimate'};


strain_m_MMTS = strain * 1.002;
strain_m_MMTS = strain_m_MMTS(1:100);
vertical_offset_MMTS = stress(1) - youngs_modulus_MMTS .* strain_m_MMTS(1);

stress_m_MMTS = zeros(100,1);
i = 1;
while i <= 100
stress_m_MMTS(i) = youngs_modulus_MMTS .* strain_m_MMTS(i) + vertical_offset_MMTS;
i = i+1;
end

stress_test = stress (1:100);

index_MMTS = find(stress_m_MMTS > stress_test, 1);
yield_stress_MMTS = stress(index_MMTS);
yield_strain_MMTS = strain(index_MMTS);


oset_plot_MMTS = stress_m_MMTS(1:index_MMTS);
strain_oset_MMTS = strain_m_MMTS(1:index_MMTS);


i = 1;
area_resilience_MMTS = 0;
U_resilience_MMTS = 0;
h1 = 0;
h2 = 0;
b = 0;
area_rect = 0;
area_trap = 0;

while i <= index_MMTS
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain_m_MMTS(i+1) - strain_m_MMTS(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 
area_resilience_MMTS = area_resilience_MMTS + area_trap;

d_u_stress_r_MMTS = u_stress(i) .* (strain_m_MMTS(i+1) - strain_m_MMTS(i));
U_resilience_MMTS = d_u_stress_r_MMTS + U_resilience_MMTS;

i = i + 1; 
end
temp = eval(u_young_MMTS);
temp2 = temp(1:index_MMTS);
u_young_MMTS = mean(temp2)
% Exten
stress_yg_exten = stress(1:1:50);
s_E = info(:,3);
strain_yg_exten = s_E(1:1:50);


estima = fitlm(strain_yg_exten,stress_yg_exten);
x1=estima.Coefficients;
youngs_modulus_extensometer=x1{'x1', 'Estimate'};


strain_m_exten = strain * 1.002;
strain_m_exten = strain_m_exten(1:100);
vertical_offset_exten = stress(1) - youngs_modulus_extensometer .* strain_m_exten(1);

stress_m_exten = zeros(100,1);
i = 1;
while i <= 100
stress_m_exten(i) = youngs_modulus_extensometer .* strain_m_exten(i) + vertical_offset_exten;
i = i+1;
end

stress_test = stress (1:100);

index_exten = find(stress_m_exten > stress_test, 1);
yield_stress_Extensometer = stress(index_exten);
yield_strain_Extensometer = strain(index_exten);


oset_plot_exten = stress_m_exten(1:index_exten);
strain_oset_exten = strain_m_exten(1:index_exten);

i = 1;
area_resilience_EXTEN = 0;
U_resilience_EXTEN = 0;
h1 = 0;
h2 = 0;
b = 0;
area_rect = 0;
area_trap = 0;

while i <= index_exten
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain_m_exten(i+1) - strain_m_exten(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 
area_resilience_EXTEN = area_resilience_EXTEN + area_trap;

d_u_stress_r_EXTEN = u_stress(i) .* (strain_m_exten(i+1) - strain_m_exten(i));
U_resilience_EXTEN = d_u_stress_r_EXTEN + U_resilience_EXTEN;

i = i + 1; 
end
temp3 = eval(u_young_exten);
temp4 = temp3(1:index_exten);
u_young_Exten = mean(temp4)

fprintf('Yield Stress: %g\n',yield_stress_MMTS)
fprintf('Yield Strain: %g\n',yield_strain_MMTS)
fprintf('Area Reduction Percent: %g +- %g\n',Area_reduction, U_Area_reduction)

% poisson Ratio


lateral_strain = ((tf-t0)/t0);
poissons = -(lateral_strain)/(strain(end));
u_poisson = eval(u_poisson);
uid = length(u_poisson);
U_poissons = u_poisson(uid);
fprintf('Poissons Ratio: %g +- %g\n',poissons, U_poissons)

% plotting 

figure(7)
plot(strain,stress,'blue')
hold on
plot(yield_x,yield_y,'diamond','MarkerSize',7,'MarkerFaceColor',[0,0,1])
hold on
plot(strain(start_index:yield_intercept),offset(start_index:yield_intercept),'green')
hold on
plot(ultimite_strain,ultimite_stress,'diamond','MarkerSize',7,'MarkerFaceColor',[0,1,0])
hold on
plot(Rupture_strain,Rupture_stress,'diamond','MarkerSize',7,'MarkerFaceColor',[1,0,0])
title('strain vs Plastic 2')
xlabel('strain')
ylabel('stress(lb/in^2)')
errorbar(strain(1:40:end),stress(1:40:end),u_stress(1:40:end),'vertical','.') % vertical errorbars for stress
errorbar(strain(1:40:end),stress(1:40:end),u_strain(1:40:end),'horizontal','.') % horizontal errorbars for strain
plot(true_strain, true_stress,'red')
hold on
plot(True_Rupture_strain,True_Rupture_stress,'diamond','MarkerSize',7,'MarkerFaceColor',[0.5,0.5,0])
hold on
plot(strain_oset_MMTS,oset_plot_MMTS,'green')
hold on
%plot(yield_strain_Extensometer,yield_stress_Extensometer,'diamond','MarkerSize',7,'MarkerFaceColor',[1,1,1])
hold on
%plot(strain_oset_exten,oset_plot_exten)
hold on
plot(yield_strain_MMTS,yield_stress_MMTS,'diamond','MarkerSize',7,'MarkerFaceColor',[0,1,1])
legend('strain vs stress curve','ultimite stress and strain', 'rupture point', ...
   'uncertianty of stress','uncertianty of strain','true strain vs stress curve','true rupture','0.2% offset line','yield point','Location','southeast')

%fprintf('Youngs Modulus: %g +- %g\n',youngs_modulus_MMTS, u_E)

%% Q2
%% Steel Only
clear
disp('Question 2 Steel Only')
%DATA
syms d F w0 t0 u_c u_tm u_d u_f L tf wf strain_yg_exten u_extens F2 u_f_2
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
stress2 = F2 ./A;
yung = stress2./ strain_yg_exten;
Q = stress .* (1 + strain);
Af = tf.*wf ; % final area
Area_reduction = (1 - (Af./A)).*100 ; % area reduction percent

u_stress = sqrt((diff(stress,F).*u_f).^2 + (diff(stress,t0).*u_c).^2 + (diff(stress,w0).*u_c).^2);
u_strain = sqrt((diff(strain,d).*u_d).^2 + (diff(strain,L).*u_tm).^2);
u_E = sqrt((diff(E,F).*u_f).^2 + (diff(E,t0).*u_c).^2 + (diff(E,w0).*u_c).^2 + ...
    (diff(E,d).*u_d).^2 + (diff(E,L).*u_tm).^2);
u_young_MMTS = sqrt((diff(E,F).*u_f).^2 + (diff(E,t0).*u_c).^2 + (diff(E,w0).*u_c).^2 + ...
    (diff(E,d).*u_d).^2 + (diff(E,L).*u_tm).^2);

u_young_exten = sqrt((diff(yung,F2).*u_f_2).^2 + (diff(yung,t0).* u_c).^2 + (diff(yung,w0).*u_c).^2 + (diff(yung,strain_yg_exten) .* u_extens).^2);

U_Stress_true = sqrt((diff(Q,F).*u_f).^2 + (diff(Q,t0).*u_c).^2 + (diff(Q,w0).*u_c).^2 + ...
    (diff(Q,d).*u_d).^2 + (diff(Q,L).*u_tm).^2);
U_Area_reduction = sqrt((diff(Af,tf).*u_c).^2 + (diff(Af,wf).*u_c).^2 + (diff(A,t0).*u_c).^2 + (diff(A,w0).*u_c).^2);
poissons = -((tf- t0)./t0)./(d./L);
u_poisson = sqrt((diff(poissons,t0) .* u_c).^2 + (diff(poissons,tf) .* u_c).^2 + (diff(poissons,d) .* u_d).^2 + (diff(poissons,L) .* u_tm).^2);

% manual differentiation
%u_stress = sqrt(((1/(t0.*w0)).*u_f).^2 + ((avgF/((t0.^2).*w0).*u_t).^2) + (-avgF/(t0.*w0.^2).*u_w).^2);
%u_strain = sqrt((((1./L).*u_d).^2)+(((avgd./(L.^2)).*u_L).^2));

%Data input
rawdata = importdata('M008_Steel.dat','\t',5);
info = rawdata.data;
d = info(:,2); % displacment data
F = info(:,1); % force data
F2 = F(1:50);
s_E = info(:,3); % strain data for E(not used in this but avaliable)

w0 = 0.523; %in change these values for each sample
t0 = 0.062; %in change these values for each sample
L = 6.5; %in change these values for each sample
tf = 0.056;  % final thickness (change)
Lf = 6.9; % final length 
wf = 0.426; %final width
u_extens = 5 .* 10 .^-8; % Reading on axial strain is x e^-7, half of that

u_c = 0.0005; %in 
u_tm = 1/64; %in
u_d = 0.012 .* d;
u_f = 0.012 .* F;
u_f_2 = .012 .*F2;
A = eval(A);
Af = eval(Af);
u_stress = eval(u_stress);
u_strain = eval(u_strain);
u_E = eval(u_E);
u_E = mean(u_E(1:50)); % to find definite number for youngs modulus
stress = eval(stress);
strain = eval(strain);
U_strain_true = u_strain;
U_Stress_true = eval(U_Stress_true);
Area_reduction = eval(Area_reduction);
U_Area_reduction = eval(U_Area_reduction);

%Area reduction

% t0 = intial thickness
% w0 = intial thickness
% tf = final thickness 
% L = intial length
% wf =final thickness 

%yield point 0.2 percent offset and Youngs Modulus(E)

p = polyfit(strain(1:50),stress(1:50),1);
E = p(1);
b = stress(1)-E.*(strain(1) .* 1.002);
p(2) = b;
offset = polyval(p,strain(1));
yield_intercept = find(offset >= stress,1);
yield_x = strain(yield_intercept);
yield_y = stress(yield_intercept);
start_index = find(offset >= stress(1),1);



% Resiliance and toughness

i = 1;
area_resilience = 0;
U_resilience = 0;
while i <=50
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain(i+1) - strain(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 
area_resilience = area_resilience + area_trap;

d_u_stress_r = u_stress(i) .* (strain(i+1) - strain(i));
U_resilience = d_u_stress_r + U_resilience;

i = i + 1; 
end

i = 1;
area_toughness = 0;
U_toughness = 0;

while i <= length(stress) - 1
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain(i+1) - strain(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 

area_toughness = area_toughness + area_trap;
d_u_stress_t = u_stress(i) .* (strain(i+1) - strain(i));
U_toughness = d_u_stress_t + U_toughness;

i = i + 1; 
end




% MMTS

stress_yg_MMTS = stress(1:1:50);
strain_yg_MMTS = strain(1:1:50);


estima = fitlm(strain_yg_MMTS,stress_yg_MMTS);
x1=estima.Coefficients;
youngs_modulus_MMTS=x1{'x1', 'Estimate'};


strain_m_MMTS = strain * 1.002;
strain_m_MMTS = strain_m_MMTS(1:100);
vertical_offset_MMTS = stress(1) - youngs_modulus_MMTS .* strain_m_MMTS(1);

stress_m_MMTS = zeros(100,1);
i = 1;
while i <= 100
stress_m_MMTS(i) = youngs_modulus_MMTS .* strain_m_MMTS(i) + vertical_offset_MMTS;
i = i+1;
end

stress_test = stress (1:100);

index_MMTS = find(stress_m_MMTS > stress_test, 1);
yield_stress_MMTS = stress(index_MMTS);
yield_strain_MMTS = strain(index_MMTS);


oset_plot_MMTS = stress_m_MMTS(1:index_MMTS);
strain_oset_MMTS = strain_m_MMTS(1:index_MMTS);

i = 1;
area_resilience_MMTS = 0;
U_resilience_MMTS = 0;
h1 = 0;
h2 = 0;
b = 0;
area_rect = 0;
area_trap = 0;

while i <= index_MMTS
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain_m_MMTS(i+1) - strain_m_MMTS(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 
area_resilience_MMTS = area_resilience_MMTS + area_trap;

d_u_stress_r_MMTS = u_stress(i) .* (strain_m_MMTS(i+1) - strain_m_MMTS(i));
U_resilience_MMTS = d_u_stress_r_MMTS + U_resilience_MMTS;

i = i + 1; 
end
temp = eval(u_young_MMTS);
temp2 = temp(1:index_MMTS);
u_young_MMTS = mean(temp2);

Brinell_Hardness_MMTS = youngs_modulus_MMTS ./ 500;

% Exten
stress_yg_exten = stress(1:1:50);
s_E = info(:,3);
strain_yg_exten = s_E(1:1:50);


estima = fitlm(strain_yg_exten,stress_yg_exten);
x1=estima.Coefficients;
youngs_modulus_extensometer=x1{'x1', 'Estimate'};


strain_m_exten = strain * 1.002;
strain_m_exten = strain_m_exten(1:100);
vertical_offset_exten = stress(1) - youngs_modulus_extensometer .* strain_m_exten(1);

stress_m_exten = zeros(100,1);
i = 1;
while i <= 100
stress_m_exten(i) = youngs_modulus_extensometer .* strain_m_exten(i) + vertical_offset_exten;
i = i+1;
end

stress_test = stress (1:100);

index_exten = find(stress_m_exten > stress_test, 1);
yield_stress_Extensometer = stress(index_exten);
yield_strain_Extensometer = strain(index_exten);


oset_plot_exten = stress_m_exten(1:index_exten);
strain_oset_exten = strain_m_exten(1:index_exten);

i = 1;
area_resilience_EXTEN = 0;
U_resilience_EXTEN = 0;
h1 = 0;
h2 = 0;
b = 0;
area_rect = 0;
area_trap = 0;

while i <= index_exten
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain_m_exten(i+1) - strain_m_exten(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 
area_resilience_EXTEN = area_resilience_EXTEN + area_trap;

d_u_stress_r_EXTEN = u_stress(i) .* (strain_m_exten(i+1) - strain_m_exten(i));
U_resilience_EXTEN = d_u_stress_r_EXTEN + U_resilience_EXTEN;

i = i + 1; 
end
temp3 = eval(u_young_exten);
temp4 = temp3(1:index_exten);
u_young_Exten = mean(temp4);

Brinell_Hardness_Exten = youngs_modulus_extensometer ./ 500;

fprintf('Yield Stress Extensometer: %g\n',yield_stress_Extensometer)
fprintf('Yield Strain Extensometer: %g\n',yield_strain_Extensometer)
fprintf('Area Reduction Percent: %g +- %g\n',Area_reduction, U_Area_reduction)


fprintf('Yield Stress MMTS: %g\n',yield_stress_MMTS)
fprintf('Yield Strain MMTS: %g\n',yield_strain_MMTS)

fprintf('Brinell Hardness MMTS: %g\n',Brinell_Hardness_MMTS)
fprintf('Brinell Hardness Extensometer: %g\n',Brinell_Hardness_Exten)

% Plotting 

figure(8)
plot(strain,stress,'blue')
hold on
title('Strain vs Stress of Steel')
xlabel('Strain')
ylabel('Stress(lb/in^2)')
errorbar(strain(1:40:end),stress(1:40:end),u_stress(1:40:end),'vertical','.') % vertical errorbars for stress
errorbar(strain(1:40:end),stress(1:40:end),u_strain(1:40:end),'horizontal','.') % horizontal errorbars for strain
plot(strain_oset_MMTS,oset_plot_MMTS)
hold on
plot(yield_strain_MMTS,yield_stress_MMTS,'diamond','MarkerSize',7,'MarkerFaceColor',[0,1,1])
hold on
plot(strain_oset_exten,oset_plot_exten)
hold on
plot(yield_strain_Extensometer,yield_stress_Extensometer,'diamond','MarkerSize',7,'MarkerFaceColor',[1,1,1])
legend('Strain vs Stress Curve', 'Uncertianty of Stress','Uncertianty of Strain','0.2% Offset Line MMTS','Yield Point MMTS','0.2% Offset Line Extensometer','Yield Point Extensometer','Location','south')

fprintf('Youngs Modulus: %g +- %g\n',youngs_modulus_MMTS, u_E)
fprintf('Modulus Resistance Extensometer: %g\n',area_resilience_EXTEN)
fprintf('Modulus Resistance MMTS: %g\n',area_resilience_MMTS)
%% Q3
%% Both Metals

%Steel---------------------------------------------------------------------------------------------------------------
clear d F w0 t0 u_c u_tm u_d u_f L tf wf strain_yg_exten u_extens F2 u_f_2
disp('Question 3 Both Metals')
disp('Steel')
%DATA
syms d F w0 t0 u_c u_tm u_d u_f L tf wf strain_yg_exten u_extens F2 u_f_2
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
Q = stress .* (1 + strain);
Af = tf.*wf ; % final area
Area_reduction = (1 - (Af./A)).*100 ; % area reduction percent 

u_stress = sqrt((diff(stress,F).*u_f).^2 + (diff(stress,t0).*u_c).^2 + (diff(stress,w0).*u_c).^2);
u_strain = sqrt((diff(strain,d).*u_d).^2 + (diff(strain,L).*u_tm).^2);
u_E = sqrt((diff(E,F).*u_f).^2 + (diff(E,t0).*u_c).^2 + (diff(E,w0).*u_c).^2 + ...
    (diff(E,d).*u_d).^2 + (diff(E,L).*u_tm).^2);
U_Stress_true = sqrt((diff(Q,F).*u_f).^2 + (diff(Q,t0).*u_c).^2 + (diff(Q,w0).*u_c).^2 + ...
    (diff(Q,d).*u_d).^2 + (diff(Q,L).*u_tm).^2);
U_Area_reduction = sqrt((diff(Af,tf).*u_c).^2 + (diff(Af,wf).*u_c).^2 + (diff(A,t0).*u_c).^2 + (diff(A,w0).*u_c).^2);
poissons = -((tf- t0)./t0)./(d./L);
stress2 = F2 ./A;
yung = stress2./ strain_yg_exten;
u_young_MMTS = sqrt((diff(E,F).*u_f).^2 + (diff(E,t0).*u_c).^2 + (diff(E,w0).*u_c).^2 + ...
    (diff(E,d).*u_d).^2 + (diff(E,L).*u_tm).^2);

% manual differentiation
%u_stress = sqrt(((1/(t0.*w0)).*u_f).^2 + ((avgF/((t0.^2).*w0).*u_t).^2) + (-avgF/(t0.*w0.^2).*u_w).^2);
%u_strain = sqrt((((1./L).*u_d).^2)+(((avgd./(L.^2)).*u_L).^2));

%Data input
rawdata = importdata('M008_Steel.dat','\t',5);
info = rawdata.data;
d = info(:,2); % displacment data
F = info(:,1); % force data
s_E = info(:,3); % strain data for E(not used in this but avaliable)

w0 = 0.523; %in change these values for each sample
t0 = 0.062; %in change these values for each sample
L = 6.5; %in change these values for each sample
tf = 0.056;  % final thickness (change)

u_c = 0.0005; %in 
u_tm = 1/64; %in
u_d = 0.012 .* d;
u_f = 0.012 .* F;
A = eval(A);
u_stress = eval(u_stress);
u_strain = eval(u_strain);
u_E = eval(u_E);
u_E = mean(u_E(1:50)); % to find definite number for youngs modulus
stress = eval(stress);
strain = eval(strain);
U_strain_true = u_strain;
U_Stress_true = eval(U_Stress_true);
Af = eval(Af);
Area_reduction = eval(Area_reduction);
U_Area_reduction = eval(U_Area_reduction);

%Finding Values

% Finding ultimite stress and strain
ultimite_stress = max(stress); % could also be stress(find(max(stress) == stress))
ultimite_strain = strain(find(max(stress) == stress));

%Rupture stress and strain
Rupture_stress = stress(end);
Rupture_strain = strain(end);

%yield point 0.2 percent offset and Youngs Modulus(E)

p = polyfit(strain(1:50),stress(1:50),1);
E = p(1);
b = stress(1)-E.*(strain(1) .* 1.002);
p(2) = b;
offset = polyval(p,strain(1));
yield_intercept = find(offset >= stress,1);
yield_x = strain(yield_intercept);
yield_y = stress(yield_intercept);
start_index = find(offset >= stress(1),1);

% true values 
true_stress = stress .* (1+strain);
true_strain = log(1+strain);

True_Rupture_stress = true_stress(end);
True_Rupture_strain = true_strain(end);



fprintf('Ultimate Stress: %g +- %g\n',ultimite_stress, u_stress(find(max(stress) == stress)))
fprintf('Ultimate Strain: %g +- %g\n',ultimite_strain, u_strain(find(max(stress) == stress)))
fprintf('Rupture Stress: %g +- %g\n',Rupture_stress, u_stress(end))
fprintf('Rupture Strain: %g +- %g\n',Rupture_strain, u_strain(end))




% Resiliance and Toughness

i = 1;
area_resilience = 0;
U_resilience = 0;
while i <=50
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain(i+1) - strain(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 
area_resilience = area_resilience + area_trap;

d_u_stress_r = u_stress(i) .* (strain(i+1) - strain(i));
U_resilience = d_u_stress_r + U_resilience;

i = i + 1; 
end

i = 1;
area_toughness = 0;
U_toughness = 0;

while i <= length(stress) - 1
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain(i+1) - strain(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 

area_toughness = area_toughness + area_trap;
d_u_stress_t = u_stress(i) .* (strain(i+1) - strain(i));
U_toughness = d_u_stress_t + U_toughness;

i = i + 1; 
end


fprintf('Area Resilience: %g +- %g\n',area_resilience, U_resilience)
fprintf('Area Toughness: %g +- %g\n',area_toughness, U_toughness)







% MMTS

stress_yg_MMTS = stress(1:1:50);
strain_yg_MMTS = strain(1:1:50);


estima = fitlm(strain_yg_MMTS,stress_yg_MMTS);
x1=estima.Coefficients;
youngs_modulus_MMTS=x1{'x1', 'Estimate'};


strain_m_MMTS = strain * 1.002;
strain_m_MMTS = strain_m_MMTS(1:100);
vertical_offset_MMTS = stress(1) - youngs_modulus_MMTS .* strain_m_MMTS(1);

stress_m_MMTS = zeros(100,1);
i = 1;
while i <= 100
stress_m_MMTS(i) = youngs_modulus_MMTS .* strain_m_MMTS(i) + vertical_offset_MMTS;
i = i+1;
end

stress_test = stress (1:100);

index_MMTS = find(stress_m_MMTS > stress_test, 1);
yield_stress_MMTS = stress(index_MMTS);
yield_strain_MMTS = strain(index_MMTS);


oset_plot_MMTS = stress_m_MMTS(1:index_MMTS);
strain_oset_MMTS = strain_m_MMTS(1:index_MMTS);

i = 1;
area_resilience_MMTS = 0;
U_resilience_MMTS = 0;
h1 = 0;
h2 = 0;
b = 0;
area_rect = 0;
area_trap = 0;

while i <= index_MMTS
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain_m_MMTS(i+1) - strain_m_MMTS(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 
area_resilience_MMTS = area_resilience_MMTS + area_trap;

d_u_stress_r_MMTS = u_stress(i) .* (strain_m_MMTS(i+1) - strain_m_MMTS(i));
U_resilience_MMTS = d_u_stress_r_MMTS + U_resilience_MMTS;

i = i + 1; 
end

temp = eval(u_young_MMTS);
temp2 = temp(1:index_MMTS);
u_young_MMTS = mean(temp2);

% Exten
stress_yg_exten = stress(1:1:50);
s_E = info(:,3);
strain_yg_exten = s_E(1:1:50);


estima = fitlm(strain_yg_exten,stress_yg_exten);
x1=estima.Coefficients;
youngs_modulus_extensometer=x1{'x1', 'Estimate'};


strain_m_exten = strain * 1.002;
strain_m_exten = strain_m_exten(1:100);
vertical_offset_exten = stress(1) - youngs_modulus_extensometer .* strain_m_exten(1);

stress_m_exten = zeros(100,1);
i = 1;
while i <= 100
stress_m_exten(i) = youngs_modulus_extensometer .* strain_m_exten(i) + vertical_offset_exten;
i = i+1;
end

stress_test = stress (1:100);

index_exten = find(stress_m_exten > stress_test, 1);
yield_stress_Extensometer = stress(index_exten);
yield_strain_Extensometer = strain(index_exten);


oset_plot_exten = stress_m_MMTS(1:index_exten);
strain_oset_exten = strain_m_MMTS(1:index_exten);




fprintf('Yield Stress: %g\n',yield_stress_MMTS)
fprintf('Yield Strain: %g\n',yield_strain_MMTS)

% Poisson Ratio


lateral_strain = ((tf-t0)/t0);
poissons = -(lateral_strain)/(strain(end));

fprintf('Poissons Ratio: %g\n',poissons)

% New Variables
AA = strain;
BB = stress;
CC = u_stress;
DD = u_strain;
EE = true_strain;
FF = true_stress;

fprintf('Youngs Modulus: %g +- %g\n',youngs_modulus_MMTS, u_E)

% True Uncertainty

rupture_index = length(true_stress);
U_Stress_True_rupture = U_Stress_true(rupture_index);
U_Strain_True_rupture = U_strain_true(rupture_index);
fprintf('True Rupture Stress: %g +-  %g\n',True_Rupture_stress,U_Stress_True_rupture)
fprintf('True Rupture Strain: %g +-  %g\n',True_Rupture_strain,U_Strain_True_rupture)
%Aluminum------------------------------------------------------------------------------------------------------------
clear d F w0 t0 u_c u_tm u_d u_f L tf wf strain_yg_exten u_extens F2 u_f_2
%DATA
syms d F w0 t0 u_c u_tm u_d u_f L tf wf strain_yg_exten u_extens F2 u_f_2
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
Q = stress .* (1 + strain);
Af = tf.*wf ; % final area
Area_reduction = (1 - (Af./A)).*100 ; % area reduction percent
U_Area_reduction = sqrt((diff(Af,tf).*u_c).^2 + (diff(Af,wf).*u_c).^2 + (diff(A,t0).*u_c).^2 + (diff(A,w0).*u_c).^2);
poissons = -((tf- t0)./t0)./(d./L);

u_stress = sqrt((diff(stress,F).*u_f).^2 + (diff(stress,t0).*u_c).^2 + (diff(stress,w0).*u_c).^2);
u_strain = sqrt((diff(strain,d).*u_d).^2 + (diff(strain,L).*u_tm).^2);
u_E = sqrt((diff(E,F).*u_f).^2 + (diff(E,t0).*u_c).^2 + (diff(E,w0).*u_c).^2 + ...
    (diff(E,d).*u_d).^2 + (diff(E,L).*u_tm).^2);
U_Stress_true = sqrt((diff(Q,F).*u_f).^2 + (diff(Q,t0).*u_c).^2 + (diff(Q,w0).*u_c).^2 + ...
    (diff(Q,d).*u_d).^2 + (diff(Q,L).*u_tm).^2);

% manual differentiation
%u_stress = sqrt(((1/(t0.*w0)).*u_f).^2 + ((avgF/((t0.^2).*w0).*u_t).^2) + (-avgF/(t0.*w0.^2).*u_w).^2);
%u_strain = sqrt((((1./L).*u_d).^2)+(((avgd./(L.^2)).*u_L).^2));
disp('Aluminum')
%Data input
rawdata = importdata('M008_Aluminum.dat','\t',5);
info = rawdata.data;
d = info(:,2); % displacment data
F = info(:,1); % force data
s_E = info(:,3); % strain data for E(not used in this but avaliable)

w0 = 0.518; %in change these values for each sample
t0 = 0.062; %in change these values for each sample
L = 6.5; %in change these values for each sample
tf = 0.0615;  % final thickness (change)
Lf = 6.9; % final length 
wf = 0.426; %final width
u_extens = 5 .* 10 .^-8; % Reading on axial strain is x e^-7, half of that

u_c = 0.0005; %in 
u_tm = 1/64; %in
u_d = 0.012 .* d;
u_f = 0.012 .* F;
A = eval(A);
u_stress = eval(u_stress);
u_strain = eval(u_strain);
u_E = eval(u_E);
u_E = mean(u_E(1:50)); % to find definite number for youngs modulus
stress = eval(stress);
strain = eval(strain);
U_strain_true = u_strain;
U_Stress_true = eval(U_Stress_true);
Af = eval(Af);
Area_reduction = eval(Area_reduction);
U_Area_reduction = eval(U_Area_reduction);

%Finding Values

% Finding ultimite stress and strain
ultimite_stress = max(stress); % could also be stress(find(max(stress) == stress))
ultimite_strain = strain(find(max(stress) == stress));

%Rupture stress and strain
Rupture_stress = stress(end);
Rupture_strain = strain(end);

%yield point 0.2 percent offset and Youngs Modulus(E)

p = polyfit(strain(1:50),stress(1:50),1);
E = p(1);
b = stress(1)-E.*(strain(1) .* 1.002);
p(2) = b;
offset = polyval(p,strain(1));
yield_intercept = find(offset >= stress,1);
yield_x = strain(yield_intercept);
yield_y = stress(yield_intercept);
start_index = find(offset >= stress(1),1);

% true values 
true_stress = stress .* (1+strain);
true_strain = log(1+strain);

True_Rupture_stress = true_stress(end);
True_Rupture_strain = true_strain(end);



fprintf('Ultimate Stress: %g +- %g\n',ultimite_stress, u_stress(find(max(stress) == stress)))
fprintf('Ultimate Strain: %g +- %g\n',ultimite_strain, u_strain(find(max(stress) == stress)))





% Resiliance and toughness

i = 1;
area_resilience = 0;
U_resilience = 0;
while i <=50
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain(i+1) - strain(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 
area_resilience = area_resilience + area_trap;

d_u_stress_r = u_stress(i) .* (strain(i+1) - strain(i));
U_resilience = d_u_stress_r + U_resilience;

i = i + 1; 
end

i = 1;
area_toughness = 0;
U_toughness = 0;

while i <= length(stress) - 1
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain(i+1) - strain(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 

area_toughness = area_toughness + area_trap;
d_u_stress_t = u_stress(i) .* (strain(i+1) - strain(i));
U_toughness = d_u_stress_t + U_toughness;

i = i + 1; 
end



fprintf('True Rupture Stress: %g\n',True_Rupture_stress)
fprintf('True Rupture Strain: %g\n',True_Rupture_strain)





% MMTS

stress_yg_MMTS = stress(1:1:50);
strain_yg_MMTS = strain(1:1:50);


estima = fitlm(strain_yg_MMTS,stress_yg_MMTS);
x1=estima.Coefficients;
youngs_modulus_MMTS=x1{'x1', 'Estimate'};


strain_m_MMTS = strain * 1.002;
strain_m_MMTS = strain_m_MMTS(1:100);
vertical_offset_MMTS = stress(1) - youngs_modulus_MMTS .* strain_m_MMTS(1);

stress_m_MMTS = zeros(100,1);
i = 1;
while i <= 100
stress_m_MMTS(i) = youngs_modulus_MMTS .* strain_m_MMTS(i) + vertical_offset_MMTS;
i = i+1;
end

stress_test = stress (1:100);

index_MMTS = find(stress_m_MMTS > stress_test, 1);
yield_stress_MMTS = stress(index_MMTS);
yield_strain_MMTS = strain(index_MMTS);


oset_plot_MMTS = stress_m_MMTS(1:index_MMTS);
strain_oset_MMTS = strain_m_MMTS(1:index_MMTS);



% Exten
stress_yg_exten = stress(1:1:50);
s_E = info(:,3);
strain_yg_exten = s_E(1:1:50);


estima = fitlm(strain_yg_exten,stress_yg_exten);
x1=estima.Coefficients;
youngs_modulus_extensometer=x1{'x1', 'Estimate'};


strain_m_exten = strain * 1.002;
strain_m_exten = strain_m_exten(1:100);
vertical_offset_exten = stress(1) - youngs_modulus_extensometer .* strain_m_exten(1);

stress_m_exten = zeros(100,1);
i = 1;
while i <= 100
stress_m_exten(i) = youngs_modulus_extensometer .* strain_m_exten(i) + vertical_offset_exten;
i = i+1;
end

stress_test = stress (1:100);

index_exten = find(stress_m_exten > stress_test, 1);
yield_stress_Extensometer = stress(index_exten);
yield_strain_Extensometer = strain(index_exten);


oset_plot_exten = stress_m_MMTS(1:index_exten);
strain_oset_exten = strain_m_MMTS(1:index_exten);




fprintf('Yield Stress: %g\n',yield_stress_MMTS)
fprintf('Yield Strain: %g\n',yield_strain_MMTS)

% Poisson Ratio


lateral_strain = ((tf-t0)/t0);
poissons = -(lateral_strain)/(strain(end));

fprintf('Poissons Ratio: %g\n',poissons)


% Plotting 

figure(9)
plot(strain,stress,'blue')
hold on
plot(AA,BB)
hold on
title('Strain vs Stress of Steel and Aluminum')
xlabel('Strain')
ylabel('Stress(lb/in^2)')
errorbar(strain(1:40:end),stress(1:40:end),u_stress(1:40:end),'vertical','.') % vertical errorbars for stress
errorbar(strain(1:40:end),stress(1:40:end),u_strain(1:40:end),'horizontal','.') % horizontal errorbars for strain
errorbar(AA(1:40:end),BB(1:40:end),CC(1:40:end),'vertical','.') % vertical errorbars for stress
errorbar(AA(1:40:end),BB(1:40:end),DD(1:40:end),'horizontal','.') % horizontal errorbars for strain
plot(true_strain, true_stress)
hold on
plot(EE,FF)
legend('Aluminum Strain vs Stress Curve','Steel Strain vs Stress Curve','Aluminum Uncertianty of Stress','Aluminum Uncertianty of Strain','Steel Uncertianty of Stress','Steel Uncertianty of Strain','Aluminum True Strain vs Stress Curve','Steel True Strain vs Stress Curve','Location','southeast')

fprintf('Youngs Modulus: %g +- %g\n',youngs_modulus_MMTS, u_E)

% True Uncertainty

rupture_index = length(true_stress);
U_Stress_True_rupture = U_Stress_true(rupture_index);
U_Strain_True_rupture = U_strain_true(rupture_index);
fprintf('Rupture Stress: %g +- %g\n',Rupture_stress,U_Stress_True_rupture)
fprintf('Rupture Strain: %g +- %g\n',Rupture_strain,U_Strain_True_rupture)
fprintf('Area Resilience: %g +- %g\n',area_resilience, U_resilience)
fprintf('Area Toughness: %g +- %g\n',area_toughness, U_toughness)
%% Q4
%% Both Carbon Fiber's

%Carbon Fiber 45---------------------------------------------------------------------------------------------------------------
clear d F w0 t0 u_c u_tm u_d u_f L tf wf strain_yg_exten u_extens F2 u_f_2
disp('Carbon Fiber 45')
%DATA
syms d F w0 t0 u_c u_tm u_d u_f L tf wf strain_yg_exten u_extens F2 u_f_2
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
Q = stress .* (1 + strain);

u_stress = sqrt((diff(stress,F).*u_f).^2 + (diff(stress,t0).*u_c).^2 + (diff(stress,w0).*u_c).^2);
u_strain = sqrt((diff(strain,d).*u_d).^2 + (diff(strain,L).*u_tm).^2);
u_E = sqrt((diff(E,F).*u_f).^2 + (diff(E,t0).*u_c).^2 + (diff(E,w0).*u_c).^2 + ...
    (diff(E,d).*u_d).^2 + (diff(E,L).*u_tm).^2);
U_Stress_true = sqrt((diff(Q,F).*u_f).^2 + (diff(Q,t0).*u_c).^2 + (diff(Q,w0).*u_c).^2 + ...
    (diff(Q,d).*u_d).^2 + (diff(Q,L).*u_tm).^2);

% manual differentiation
%u_stress = sqrt(((1/(t0.*w0)).*u_f).^2 + ((avgF/((t0.^2).*w0).*u_t).^2) + (-avgF/(t0.*w0.^2).*u_w).^2);
%u_strain = sqrt((((1./L).*u_d).^2)+(((avgd./(L.^2)).*u_L).^2));

%Data input
rawdata = importdata('M008_Carbon_Fiber_45.dat','\t',5);
info = rawdata.data;
d = info(:,2); % displacment data
F = info(:,1); % force data
s_E = info(:,3); % strain data for E(not used in this but avaliable)

w0 = 0.518; %in change these values for each sample
t0 = 0.060; %in change these values for each sample
L = 6.5; %in change these values for each sample
tf = 0.060;  % final thickness (change)

u_c = 0.0005; %in 
u_tm = 1/64; %in
u_d = 0.012 .* d;
u_f = 0.012 .* F;
A = eval(A);
u_stress = eval(u_stress);
u_strain = eval(u_strain);
u_E = eval(u_E);
u_E = mean(u_E(1:50)); % to find definite number for youngs modulus
stress = eval(stress);
strain = eval(strain);
U_strain_true = u_strain;
U_Stress_true = eval(U_Stress_true);

%Finding Values

% Finding ultimite stress and strain
ultimite_stress = max(stress); % could also be stress(find(max(stress) == stress))
ultimite_strain = strain(find(max(stress) == stress));

%Rupture stress and strain
Rupture_stress = stress(end);
Rupture_strain = strain(end);

%yield point 0.2 percent offset and Youngs Modulus(E)

p = polyfit(strain(1:50),stress(1:50),1);
E = p(1);
b = stress(1)-E.*(strain(1) .* 1.002);
p(2) = b;
offset = polyval(p,strain(1));
yield_intercept = find(offset >= stress,1);
yield_x = strain(yield_intercept);
yield_y = stress(yield_intercept);
start_index = find(offset >= stress(1),1);

% true values 
true_stress = stress .* (1+strain);
true_strain = log(1+strain);

True_Rupture_stress = true_stress(end);
True_Rupture_strain = true_strain(end);



fprintf('Ultimate Stress: %g +- %g\n',ultimite_stress, u_stress(find(max(stress) == stress)))
fprintf('Ultimate Strain: %g +- %g\n',ultimite_strain, u_strain(find(max(stress) == stress)))
fprintf('Rupture Stress: %g +- %g\n',Rupture_stress, u_stress(end))
fprintf('Rupture Strain: %g +- %g\n',Rupture_strain, u_strain(end))


rupture_index = length(true_stress);
U_Stress_True_rupture = U_Stress_true(rupture_index);
U_Strain_True_rupture = U_strain_true(rupture_index);







% Resiliance and toughness

i = 1;
area_resilience = 0;
U_resilience = 0;
while i <=50
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain(i+1) - strain(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 
area_resilience = area_resilience + area_trap;

d_u_stress_r = u_stress(i) .* (strain(i+1) - strain(i));
U_resilience = d_u_stress_r + U_resilience;

i = i + 1; 
end

i = 1;
area_toughness = 0;
U_toughness = 0;

while i <= length(stress) - 1
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain(i+1) - strain(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 

area_toughness = area_toughness + area_trap;
d_u_stress_t = u_stress(i) .* (strain(i+1) - strain(i));
U_toughness = d_u_stress_t + U_toughness;

i = i + 1; 
end



fprintf('Area Toughness: %g +- %g\n',area_toughness, U_toughness)



% MMTS

stress_yg_MMTS = stress(1:1:50);
strain_yg_MMTS = strain(1:1:50);


estima = fitlm(strain_yg_MMTS,stress_yg_MMTS);
x1=estima.Coefficients;
youngs_modulus_MMTS=x1{'x1', 'Estimate'};


strain_m_MMTS = strain * 1.002;
strain_m_MMTS = strain_m_MMTS(1:100);
vertical_offset_MMTS = stress(1) - youngs_modulus_MMTS .* strain_m_MMTS(1);

stress_m_MMTS = zeros(100,1);
i = 1;
while i <= 100
stress_m_MMTS(i) = youngs_modulus_MMTS .* strain_m_MMTS(i) + vertical_offset_MMTS;
i = i+1;
end

stress_test = stress (1:100);

index_MMTS = find(stress_m_MMTS > stress_test, 1);
yield_stress_MMTS = stress(index_MMTS);
yield_strain_MMTS = strain(index_MMTS);


oset_plot_MMTS = stress_m_MMTS(1:index_MMTS);
strain_oset_MMTS = strain_m_MMTS(1:index_MMTS);



% Exten
stress_yg_exten = stress(1:1:50);
s_E = info(:,3);
strain_yg_exten = s_E(1:1:50);


estima = fitlm(strain_yg_exten,stress_yg_exten);
x1=estima.Coefficients;
youngs_modulus_extensometer=x1{'x1', 'Estimate'};


strain_m_exten = strain * 1.002;
strain_m_exten = strain_m_exten(1:100);
vertical_offset_exten = stress(1) - youngs_modulus_extensometer .* strain_m_exten(1);

stress_m_exten = zeros(100,1);
i = 1;
while i <= 100
stress_m_exten(i) = youngs_modulus_extensometer .* strain_m_exten(i) + vertical_offset_exten;
i = i+1;
end

stress_test = stress (1:100);

index_exten = find(stress_m_exten > stress_test, 1);
yield_stress_Extensometer = stress(index_exten);
yield_strain_Extensometer = strain(index_exten);


oset_plot_exten = stress_m_MMTS(1:index_exten);
strain_oset_exten = strain_m_MMTS(1:index_exten);






% Poisson Ratio


lateral_strain = ((tf-t0)/t0);
poissons = -(lateral_strain)/(strain(end));



% New Variables
AA = strain;
BB = stress;
CC = u_stress;
DD = u_strain;



%Carbon Fiber 90------------------------------------------------------------------------------------------------------------
clear d F w0 t0 u_c u_tm u_d u_f L
%DATA
disp('Carbon Fiber 45')
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
Q = stress .* (1 + strain);

u_stress = sqrt((diff(stress,F).*u_f).^2 + (diff(stress,t0).*u_c).^2 + (diff(stress,w0).*u_c).^2);
u_strain = sqrt((diff(strain,d).*u_d).^2 + (diff(strain,L).*u_tm).^2);
u_E = sqrt((diff(E,F).*u_f).^2 + (diff(E,t0).*u_c).^2 + (diff(E,w0).*u_c).^2 + ...
    (diff(E,d).*u_d).^2 + (diff(E,L).*u_tm).^2);
U_Stress_true = sqrt((diff(Q,F).*u_f).^2 + (diff(Q,t0).*u_c).^2 + (diff(Q,w0).*u_c).^2 + ...
    (diff(Q,d).*u_d).^2 + (diff(Q,L).*u_tm).^2);

% manual differentiation
%u_stress = sqrt(((1/(t0.*w0)).*u_f).^2 + ((avgF/((t0.^2).*w0).*u_t).^2) + (-avgF/(t0.*w0.^2).*u_w).^2);
%u_strain = sqrt((((1./L).*u_d).^2)+(((avgd./(L.^2)).*u_L).^2));

%Data input
rawdata = importdata('M008_CF90.dat','\t',5);
info = rawdata.data;
d = info(:,2); % displacment data
F = info(:,1); % force data
s_E = info(:,3); % strain data for E(not used in this but avaliable)

w0 = 0.510; %in change these values for each sample
t0 = 0.072; %in change these values for each sample
L = 6.5; %in change these values for each sample
tf = 0.061;  % final thickness (change)

u_c = 0.0005; %in 
u_tm = 1/64; %in
u_d = 0.012 .* d;
u_f = 0.012 .* F;
A = eval(A);
u_stress = eval(u_stress);
u_strain = eval(u_strain);
u_E = eval(u_E);
u_E = mean(u_E(1:50)); % to find definite number for youngs modulus
stress = eval(stress);
strain = eval(strain);
U_strain_true = u_strain;
U_Stress_true = eval(U_Stress_true);

%Finding Values

% Finding ultimite stress and strain
ultimite_stress = max(stress); % could also be stress(find(max(stress) == stress))
ultimite_strain = strain(find(max(stress) == stress));

%Rupture stress and strain
Rupture_stress = stress(end);
Rupture_strain = strain(end);

%yield point 0.2 percent offset and Youngs Modulus(E)

p = polyfit(strain(1:50),stress(1:50),1);
E = p(1);
b = stress(1)-E.*(strain(1) .* 1.002);
p(2) = b;
offset = polyval(p,strain(1));
yield_intercept = find(offset >= stress,1);
yield_x = strain(yield_intercept);
yield_y = stress(yield_intercept);
start_index = find(offset >= stress(1),1);

% true values 
true_stress = stress .* (1+strain);
true_strain = log(1+strain);

True_Rupture_stress = true_stress(end);
True_Rupture_strain = true_strain(end);



fprintf('Ultimate Stress: %g +- %g\n',ultimite_stress, u_stress(find(max(stress) == stress)))
fprintf('Ultimate Strain: %g +- %g\n',ultimite_strain, u_strain(find(max(stress) == stress)))
fprintf('Rupture Stress: %g +- %g\n',Rupture_stress, u_stress(end))
fprintf('Rupture Strain: %g +- %g\n',Rupture_strain, u_strain(end))




% Resiliance and toughness

i = 1;
area_resilience = 0;
U_resilience = 0;
while i <=50
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain(i+1) - strain(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 
area_resilience = area_resilience + area_trap;

d_u_stress_r = u_stress(i) .* (strain(i+1) - strain(i));
U_resilience = d_u_stress_r + U_resilience;

i = i + 1; 
end

i = 1;
area_toughness = 0;
U_toughness = 0;

while i <= length(stress) - 1
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain(i+1) - strain(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 

area_toughness = area_toughness + area_trap;
d_u_stress_t = u_stress(i) .* (strain(i+1) - strain(i));
U_toughness = d_u_stress_t + U_toughness;

i = i + 1; 
end


fprintf('Area Toughness: %g +- %g\n',area_toughness, U_toughness)







% MMTS

stress_yg_MMTS = stress(1:1:50);
strain_yg_MMTS = strain(1:1:50);


estima = fitlm(strain_yg_MMTS,stress_yg_MMTS);
x1=estima.Coefficients;
youngs_modulus_MMTS=x1{'x1', 'Estimate'};


strain_m_MMTS = strain * 1.002;
strain_m_MMTS = strain_m_MMTS(1:100);
vertical_offset_MMTS = stress(1) - youngs_modulus_MMTS .* strain_m_MMTS(1);

stress_m_MMTS = zeros(100,1);
i = 1;
while i <= 100
stress_m_MMTS(i) = youngs_modulus_MMTS .* strain_m_MMTS(i) + vertical_offset_MMTS;
i = i+1;
end

stress_test = stress (1:100);

index_MMTS = find(stress_m_MMTS > stress_test, 1);
yield_stress_MMTS = stress(index_MMTS);
yield_strain_MMTS = strain(index_MMTS);


oset_plot_MMTS = stress_m_MMTS(1:index_MMTS);
strain_oset_MMTS = strain_m_MMTS(1:index_MMTS);



% Exten
stress_yg_exten = stress(1:1:50);
s_E = info(:,3);
strain_yg_exten = s_E(1:1:50);


estima = fitlm(strain_yg_exten,stress_yg_exten);
x1=estima.Coefficients;
youngs_modulus_extensometer=x1{'x1', 'Estimate'};


strain_m_exten = strain * 1.002;
strain_m_exten = strain_m_exten(1:100);
vertical_offset_exten = stress(1) - youngs_modulus_extensometer .* strain_m_exten(1);

stress_m_exten = zeros(100,1);
i = 1;
while i <= 100
stress_m_exten(i) = youngs_modulus_extensometer .* strain_m_exten(i) + vertical_offset_exten;
i = i+1;
end

stress_test = stress (1:100);

index_exten = find(stress_m_exten > stress_test, 1);
yield_stress_Extensometer = stress(index_exten);
yield_strain_Extensometer = strain(index_exten);


oset_plot_exten = stress_m_MMTS(1:index_exten);
strain_oset_exten = strain_m_MMTS(1:index_exten);






% Poisson Ratio


lateral_strain = ((tf-t0)/t0);
poissons = -(lateral_strain)/(strain(end));



% Plotting 

figure(10)
plot(strain,stress,'blue')
hold on
plot(AA,BB)
hold on
title('Strain vs Stress of Carbon Fiber 45 and Carbon Fiber 90')
xlabel('Strain')
ylabel('Stress(lb/in^2)')
errorbar(strain(1:40:end),stress(1:40:end),u_stress(1:40:end),'vertical','.') % vertical errorbars for stress
errorbar(strain(1:40:end),stress(1:40:end),u_strain(1:40:end),'horizontal','.') % horizontal errorbars for strain
errorbar(AA(1:40:end),BB(1:40:end),CC(1:40:end),'vertical','.') % vertical errorbars for stress
errorbar(AA(1:40:end),BB(1:40:end),DD(1:40:end),'horizontal','.') % horizontal errorbars for strain
legend('CF45 Strain vs Stress Curve','CF90 Strain vs Stress Curve','CF45 Uncertianty of Stress','CF90 Uncertianty of Strain','Location','northeast')

fprintf('Youngs Modulus: %g +- %g\n',youngs_modulus_MMTS, u_E)

% True Uncertainty

rupture_index = length(true_stress);
U_Stress_True_rupture = U_Stress_true(rupture_index);
U_Strain_True_rupture = U_strain_true(rupture_index);

%% Q5
%% Plastic and Aluminum

%Plastic---------------------------------------------------------------------------------------------------------------
clear
disp('Plastic 2')
%DATA
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
Q = stress .* (1 + strain);

u_stress = sqrt((diff(stress,F).*u_f).^2 + (diff(stress,t0).*u_c).^2 + (diff(stress,w0).*u_c).^2);
u_strain = sqrt((diff(strain,d).*u_d).^2 + (diff(strain,L).*u_tm).^2);
u_E = sqrt((diff(E,F).*u_f).^2 + (diff(E,t0).*u_c).^2 + (diff(E,w0).*u_c).^2 + ...
    (diff(E,d).*u_d).^2 + (diff(E,L).*u_tm).^2);
U_Stress_true = sqrt((diff(Q,F).*u_f).^2 + (diff(Q,t0).*u_c).^2 + (diff(Q,w0).*u_c).^2 + ...
    (diff(Q,d).*u_d).^2 + (diff(Q,L).*u_tm).^2);

% manual differentiation
%u_stress = sqrt(((1/(t0.*w0)).*u_f).^2 + ((avgF/((t0.^2).*w0).*u_t).^2) + (-avgF/(t0.*w0.^2).*u_w).^2);
%u_strain = sqrt((((1./L).*u_d).^2)+(((avgd./(L.^2)).*u_L).^2));

%Data input
rawdata = importdata('M008_22_plastic2.dat','\t',5);
info = rawdata.data;
d = info(:,2); % displacment data
F = info(:,1); % force data
s_E = info(:,3); % strain data for E(not used in this but avaliable)

w0 = 0.500; %in change these values for each sample
t0 = 0.137; %in change these values for each sample
L = 6.5; %in change these values for each sample
tf = 0.135;  % final thickness (change)

u_c = 0.0005; %in 
u_tm = 1/64; %in
u_d = 0.012 .* d;
u_f = 0.012 .* F;
A = eval(A);
u_stress = eval(u_stress);
u_strain = eval(u_strain);
u_E = eval(u_E);
u_E = mean(u_E(1:50)); % to find definite number for youngs modulus
stress = eval(stress);
strain = eval(strain);
U_strain_true = u_strain;
U_Stress_true = eval(U_Stress_true);

%Finding Values

% Finding ultimite stress and strain
ultimite_stress = max(stress); % could also be stress(find(max(stress) == stress))
ultimite_strain = strain(find(max(stress) == stress));

%Rupture stress and strain
Rupture_stress = stress(end);
Rupture_strain = strain(end);

%yield point 0.2 percent offset and Youngs Modulus(E)

p = polyfit(strain(1:50),stress(1:50),1);
E = p(1);
b = stress(1)-E.*(strain(1) .* 1.002);
p(2) = b;
offset = polyval(p,strain(1));
yield_intercept = find(offset >= stress,1);
yield_x = strain(yield_intercept);
yield_y = stress(yield_intercept);
start_index = find(offset >= stress(1),1);

% true values 
true_stress = stress .* (1+strain);
true_strain = log(1+strain);

True_Rupture_stress = true_stress(end);
True_Rupture_strain = true_strain(end);



fprintf('Ultimate Stress: %g +- %g\n',ultimite_stress, u_stress(find(max(stress) == stress)))
fprintf('Ultimate Strain: %g +- %g\n',ultimite_strain, u_strain(find(max(stress) == stress)))
fprintf('Rupture Stress: %g +- %g\n',Rupture_stress, u_stress(end))
fprintf('Rupture Strain: %g +- %g\n',Rupture_strain, u_strain(end))




% Resiliance and toughness

i = 1;
area_resilience = 0;
U_resilience = 0;
while i <=50
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain(i+1) - strain(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 
area_resilience = area_resilience + area_trap;

d_u_stress_r = u_stress(i) .* (strain(i+1) - strain(i));
U_resilience = d_u_stress_r + U_resilience;

i = i + 1; 
end

i = 1;
area_toughness = 0;
U_toughness = 0;

while i <= length(stress) - 1
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain(i+1) - strain(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 

area_toughness = area_toughness + area_trap;
d_u_stress_t = u_stress(i) .* (strain(i+1) - strain(i));
U_toughness = d_u_stress_t + U_toughness;

i = i + 1; 
end



fprintf('Area Toughness: %g +- %g\n',area_toughness, U_toughness)







% MMTS

stress_yg_MMTS = stress(1:1:50);
strain_yg_MMTS = strain(1:1:50);


estima = fitlm(strain_yg_MMTS,stress_yg_MMTS);
x1=estima.Coefficients;
youngs_modulus_MMTS=x1{'x1', 'Estimate'};


strain_m_MMTS = strain * 1.002;
strain_m_MMTS = strain_m_MMTS(1:100);
vertical_offset_MMTS = stress(1) - youngs_modulus_MMTS .* strain_m_MMTS(1);

stress_m_MMTS = zeros(100,1);
i = 1;
while i <= 100
stress_m_MMTS(i) = youngs_modulus_MMTS .* strain_m_MMTS(i) + vertical_offset_MMTS;
i = i+1;
end

stress_test = stress (1:100);

index_MMTS = find(stress_m_MMTS > stress_test, 1);
yield_stress_MMTS = stress(index_MMTS);
yield_strain_MMTS = strain(index_MMTS);


oset_plot_MMTS = stress_m_MMTS(1:index_MMTS);
strain_oset_MMTS = strain_m_MMTS(1:index_MMTS);



% Exten
stress_yg_exten = stress(1:1:50);
s_E = info(:,3);
strain_yg_exten = s_E(1:1:50);


estima = fitlm(strain_yg_exten,stress_yg_exten);
x1=estima.Coefficients;
youngs_modulus_extensometer=x1{'x1', 'Estimate'};


strain_m_exten = strain * 1.002;
strain_m_exten = strain_m_exten(1:100);
vertical_offset_exten = stress(1) - youngs_modulus_extensometer .* strain_m_exten(1);

stress_m_exten = zeros(100,1);
i = 1;
while i <= 100
stress_m_exten(i) = youngs_modulus_extensometer .* strain_m_exten(i) + vertical_offset_exten;
i = i+1;
end

stress_test = stress (1:100);

index_exten = find(stress_m_exten > stress_test, 1);
yield_stress_Extensometer = stress(index_exten);
yield_strain_Extensometer = strain(index_exten);


oset_plot_exten = stress_m_MMTS(1:index_exten);
strain_oset_exten = strain_m_MMTS(1:index_exten);






% Poisson Ratio


lateral_strain = ((tf-t0)/t0);
poissons = -(lateral_strain)/(strain(end));



%New Variables
AA = strain;
BB = stress;
CC = u_stress;
DD = u_strain;



% True Uncertainty

rupture_index = length(true_stress);
U_Stress_True_rupture = U_Stress_true(rupture_index);
U_Strain_True_rupture = U_strain_true(rupture_index);

%Aluminum------------------------------------------------------------------------------------------------------------
clear d F w0 t0 u_c u_tm u_d u_f L
%DATA
disp('Aluminum')
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
Q = stress .* (1 + strain);

u_stress = sqrt((diff(stress,F).*u_f).^2 + (diff(stress,t0).*u_c).^2 + (diff(stress,w0).*u_c).^2);
u_strain = sqrt((diff(strain,d).*u_d).^2 + (diff(strain,L).*u_tm).^2);
u_E = sqrt((diff(E,F).*u_f).^2 + (diff(E,t0).*u_c).^2 + (diff(E,w0).*u_c).^2 + ...
    (diff(E,d).*u_d).^2 + (diff(E,L).*u_tm).^2);
U_Stress_true = sqrt((diff(Q,F).*u_f).^2 + (diff(Q,t0).*u_c).^2 + (diff(Q,w0).*u_c).^2 + ...
    (diff(Q,d).*u_d).^2 + (diff(Q,L).*u_tm).^2);

% manual differentiation
%u_stress = sqrt(((1/(t0.*w0)).*u_f).^2 + ((avgF/((t0.^2).*w0).*u_t).^2) + (-avgF/(t0.*w0.^2).*u_w).^2);
%u_strain = sqrt((((1./L).*u_d).^2)+(((avgd./(L.^2)).*u_L).^2));

%Data input
rawdata = importdata('M008_Aluminum.dat','\t',5);
info = rawdata.data;
d = info(:,2); % displacment data
F = info(:,1); % force data
s_E = info(:,3); % strain data for E(not used in this but avaliable)

w0 = 0.518; %in change these values for each sample
t0 = 0.062; %in change these values for each sample
L = 6.5; %in change these values for each sample
tf = 0.0615;  % final thickness (change)

u_c = 0.0005; %in 
u_tm = 1/64; %in
u_d = 0.012 .* d;
u_f = 0.012 .* F;
A = eval(A);
u_stress = eval(u_stress);
u_strain = eval(u_strain);
u_E = eval(u_E);
u_E = mean(u_E(1:50)); % to find definite number for youngs modulus
stress = eval(stress);
strain = eval(strain);
U_strain_true = u_strain;
U_Stress_true = eval(U_Stress_true);

%Finding Values

% Finding ultimite stress and strain
ultimite_stress = max(stress); % could also be stress(find(max(stress) == stress))
ultimite_strain = strain(find(max(stress) == stress));

%Rupture stress and strain
Rupture_stress = stress(end);
Rupture_strain = strain(end);

%yield point 0.2 percent offset and Youngs Modulus(E)

p = polyfit(strain(1:50),stress(1:50),1);
E = p(1);
b = stress(1)-E.*(strain(1) .* 1.002);
p(2) = b;
offset = polyval(p,strain(1));
yield_intercept = find(offset >= stress,1);
yield_x = strain(yield_intercept);
yield_y = stress(yield_intercept);
start_index = find(offset >= stress(1),1);

% true values 
true_stress = stress .* (1+strain);
true_strain = log(1+strain);

True_Rupture_stress = true_stress(end);
True_Rupture_strain = true_strain(end);



fprintf('Ultimate Stress: %g +- %g\n',ultimite_stress, u_stress(find(max(stress) == stress)))
fprintf('Ultimate Strain: %g +- %g\n',ultimite_strain, u_strain(find(max(stress) == stress)))
fprintf('Rupture Stress: %g +- %g\n',Rupture_stress, u_stress(end))
fprintf('Rupture Strain: %g +- %g\n',Rupture_strain, u_strain(end))




% Resiliance and toughness

i = 1;
area_resilience = 0;
U_resilience = 0;
while i <=50
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain(i+1) - strain(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 
area_resilience = area_resilience + area_trap;

d_u_stress_r = u_stress(i) .* (strain(i+1) - strain(i));
U_resilience = d_u_stress_r + U_resilience;

i = i + 1; 
end

i = 1;
area_toughness = 0;
U_toughness = 0;

while i <= length(stress) - 1
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain(i+1) - strain(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 

area_toughness = area_toughness + area_trap;
d_u_stress_t = u_stress(i) .* (strain(i+1) - strain(i));
U_toughness = d_u_stress_t + U_toughness;

i = i + 1; 
end



fprintf('Area Toughness: %g +- %g\n',area_toughness, U_toughness)







% MMTS

stress_yg_MMTS = stress(1:1:50);
strain_yg_MMTS = strain(1:1:50);


estima = fitlm(strain_yg_MMTS,stress_yg_MMTS);
x1=estima.Coefficients;
youngs_modulus_MMTS=x1{'x1', 'Estimate'};


strain_m_MMTS = strain * 1.002;
strain_m_MMTS = strain_m_MMTS(1:100);
vertical_offset_MMTS = stress(1) - youngs_modulus_MMTS .* strain_m_MMTS(1);

stress_m_MMTS = zeros(100,1);
i = 1;
while i <= 100
stress_m_MMTS(i) = youngs_modulus_MMTS .* strain_m_MMTS(i) + vertical_offset_MMTS;
i = i+1;
end

stress_test = stress (1:100);

index_MMTS = find(stress_m_MMTS > stress_test, 1);
yield_stress_MMTS = stress(index_MMTS);
yield_strain_MMTS = strain(index_MMTS);


oset_plot_MMTS = stress_m_MMTS(1:index_MMTS);
strain_oset_MMTS = strain_m_MMTS(1:index_MMTS);



% Exten
stress_yg_exten = stress(1:1:50);
s_E = info(:,3);
strain_yg_exten = s_E(1:1:50);


estima = fitlm(strain_yg_exten,stress_yg_exten);
x1=estima.Coefficients;
youngs_modulus_extensometer=x1{'x1', 'Estimate'};


strain_m_exten = strain * 1.002;
strain_m_exten = strain_m_exten(1:100);
vertical_offset_exten = stress(1) - youngs_modulus_extensometer .* strain_m_exten(1);

stress_m_exten = zeros(100,1);
i = 1;
while i <= 100
stress_m_exten(i) = youngs_modulus_extensometer .* strain_m_exten(i) + vertical_offset_exten;
i = i+1;
end

stress_test = stress (1:100);

index_exten = find(stress_m_exten > stress_test, 1);
yield_stress_Extensometer = stress(index_exten);
yield_strain_Extensometer = strain(index_exten);


oset_plot_exten = stress_m_MMTS(1:index_exten);
strain_oset_exten = strain_m_MMTS(1:index_exten);






% Poisson Ratio


lateral_strain = ((tf-t0)/t0);
poissons = -(lateral_strain)/(strain(end));




% Plotting 

figure(11)
plot(strain,stress,'blue')
hold on
plot(AA,BB)
hold on
title('Strain vs Stress of Plastic and Aluminum')
xlabel('Strain')
ylabel('Stress(lb/in^2)')
errorbar(strain(1:40:end),stress(1:40:end),u_stress(1:40:end),'vertical','.') % vertical errorbars for stress
errorbar(strain(1:40:end),stress(1:40:end),u_strain(1:40:end),'horizontal','.') % horizontal errorbars for strain
errorbar(AA(1:40:end),BB(1:40:end),CC(1:40:end),'vertical','.') % vertical errorbars for stress
errorbar(AA(1:40:end),BB(1:40:end),DD(1:40:end),'horizontal','.') % horizontal errorbars for strain
legend('Aluminum Strain vs Stress Curve','Plastic Strain vs Stress Curve','Aluminum Uncertianty of Stress','Aluminum Uncertainty of Strain','Plastic Uncertainty of Stress','Plastic Uncertianty of Strain','Location','northeast')



% True Uncertainty

rupture_index = length(true_stress);
U_Stress_True_rupture = U_Stress_true(rupture_index);
U_Strain_True_rupture = U_strain_true(rupture_index);

%% Q6
%% All Materials Together

%Steel---------------------------------------------------------------------------------------------------------------
clear
%DATA
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
Q = stress .* (1 + strain);

u_stress = sqrt((diff(stress,F).*u_f).^2 + (diff(stress,t0).*u_c).^2 + (diff(stress,w0).*u_c).^2);
u_strain = sqrt((diff(strain,d).*u_d).^2 + (diff(strain,L).*u_tm).^2);
u_E = sqrt((diff(E,F).*u_f).^2 + (diff(E,t0).*u_c).^2 + (diff(E,w0).*u_c).^2 + ...
    (diff(E,d).*u_d).^2 + (diff(E,L).*u_tm).^2);
U_Stress_true = sqrt((diff(Q,F).*u_f).^2 + (diff(Q,t0).*u_c).^2 + (diff(Q,w0).*u_c).^2 + ...
    (diff(Q,d).*u_d).^2 + (diff(Q,L).*u_tm).^2);

% manual differentiation
%u_stress = sqrt(((1/(t0.*w0)).*u_f).^2 + ((avgF/((t0.^2).*w0).*u_t).^2) + (-avgF/(t0.*w0.^2).*u_w).^2);
%u_strain = sqrt((((1./L).*u_d).^2)+(((avgd./(L.^2)).*u_L).^2));

%Data input
rawdata = importdata('M008_Steel.dat','\t',5);
info = rawdata.data;
d = info(:,2); % displacment data
F = info(:,1); % force data
s_E = info(:,3); % strain data for E(not used in this but avaliable)

w0 = 0.523; %in change these values for each sample
t0 = 0.062; %in change these values for each sample
L = 6.5; %in change these values for each sample
tf = 0.056;  % final thickness (change)

u_c = 0.0005; %in 
u_tm = 1/64; %in
u_d = 0.012 .* d;
u_f = 0.012 .* F;
A = eval(A);
u_stress = eval(u_stress);
u_strain = eval(u_strain);
u_E = eval(u_E);
u_E = mean(u_E(1:50)); % to find definite number for youngs modulus
stress = eval(stress);
strain = eval(strain);
U_strain_true = u_strain;
U_Stress_true = eval(U_Stress_true);

%Finding Values

% Finding ultimite stress and strain
ultimite_stress = max(stress); % could also be stress(find(max(stress) == stress))
ultimite_strain = strain(find(max(stress) == stress));

%Rupture stress and strain
Rupture_stress = stress(end);
Rupture_strain = strain(end);

%yield point 0.2 percent offset and Youngs Modulus(E)

p = polyfit(strain(1:50),stress(1:50),1);
E = p(1);
b = stress(1)-E.*(strain(1) .* 1.002);
p(2) = b;
offset = polyval(p,strain(1));
yield_intercept = find(offset >= stress,1);
yield_x = strain(yield_intercept);
yield_y = stress(yield_intercept);
start_index = find(offset >= stress(1),1);

% true values 
true_stress = stress .* (1+strain);
true_strain = log(1+strain);

True_Rupture_stress = true_stress(end);
True_Rupture_strain = true_strain(end);






% Resiliance and toughness

i = 1;
area_resilience = 0;
U_resilience = 0;
while i <=50
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain(i+1) - strain(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 
area_resilience = area_resilience + area_trap;

d_u_stress_r = u_stress(i) .* (strain(i+1) - strain(i));
U_resilience = d_u_stress_r + U_resilience;

i = i + 1; 
end

i = 1;
area_toughness = 0;
U_toughness = 0;

while i <= length(stress) - 1
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain(i+1) - strain(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 

area_toughness = area_toughness + area_trap;
d_u_stress_t = u_stress(i) .* (strain(i+1) - strain(i));
U_toughness = d_u_stress_t + U_toughness;

i = i + 1; 
end








% MMTS

stress_yg_MMTS = stress(1:1:50);
strain_yg_MMTS = strain(1:1:50);


estima = fitlm(strain_yg_MMTS,stress_yg_MMTS);
x1=estima.Coefficients;
youngs_modulus_MMTS=x1{'x1', 'Estimate'};


strain_m_MMTS = strain * 1.002;
strain_m_MMTS = strain_m_MMTS(1:100);
vertical_offset_MMTS = stress(1) - youngs_modulus_MMTS .* strain_m_MMTS(1);

stress_m_MMTS = zeros(100,1);
i = 1;
while i <= 100
stress_m_MMTS(i) = youngs_modulus_MMTS .* strain_m_MMTS(i) + vertical_offset_MMTS;
i = i+1;
end

stress_test = stress (1:100);

index_MMTS = find(stress_m_MMTS > stress_test, 1);
yield_stress_MMTS = stress(index_MMTS);
yield_strain_MMTS = strain(index_MMTS);


oset_plot_MMTS = stress_m_MMTS(1:index_MMTS);
strain_oset_MMTS = strain_m_MMTS(1:index_MMTS);



% Exten
stress_yg_exten = stress(1:1:50);
s_E = info(:,3);
strain_yg_exten = s_E(1:1:50);


estima = fitlm(strain_yg_exten,stress_yg_exten);
x1=estima.Coefficients;
youngs_modulus_extensometer=x1{'x1', 'Estimate'};


strain_m_exten = strain * 1.002;
strain_m_exten = strain_m_exten(1:100);
vertical_offset_exten = stress(1) - youngs_modulus_extensometer .* strain_m_exten(1);

stress_m_exten = zeros(100,1);
i = 1;
while i <= 100
stress_m_exten(i) = youngs_modulus_extensometer .* strain_m_exten(i) + vertical_offset_exten;
i = i+1;
end

stress_test = stress (1:100);

index_exten = find(stress_m_exten > stress_test, 1);
yield_stress_Extensometer = stress(index_exten);
yield_strain_Extensometer = strain(index_exten);


oset_plot_exten = stress_m_MMTS(1:index_exten);
strain_oset_exten = strain_m_MMTS(1:index_exten);





% Poisson Ratio


lateral_strain = ((tf-t0)/t0);
poissons = -(lateral_strain)/(strain(end));



%New Variables Steel
AA = strain;
BB = stress;
CC = u_stress;
DD = u_strain;



% True Uncertainty

rupture_index = length(true_stress);
U_Stress_True_rupture = U_Stress_true(rupture_index);
U_Strain_True_rupture = U_strain_true(rupture_index);

%Aluminum------------------------------------------------------------------------------------------------------------
clear d F w0 t0 u_c u_tm u_d u_f L
%DATA
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
Q = stress .* (1 + strain);

u_stress = sqrt((diff(stress,F).*u_f).^2 + (diff(stress,t0).*u_c).^2 + (diff(stress,w0).*u_c).^2);
u_strain = sqrt((diff(strain,d).*u_d).^2 + (diff(strain,L).*u_tm).^2);
u_E = sqrt((diff(E,F).*u_f).^2 + (diff(E,t0).*u_c).^2 + (diff(E,w0).*u_c).^2 + ...
    (diff(E,d).*u_d).^2 + (diff(E,L).*u_tm).^2);
U_Stress_true = sqrt((diff(Q,F).*u_f).^2 + (diff(Q,t0).*u_c).^2 + (diff(Q,w0).*u_c).^2 + ...
    (diff(Q,d).*u_d).^2 + (diff(Q,L).*u_tm).^2);

% manual differentiation
%u_stress = sqrt(((1/(t0.*w0)).*u_f).^2 + ((avgF/((t0.^2).*w0).*u_t).^2) + (-avgF/(t0.*w0.^2).*u_w).^2);
%u_strain = sqrt((((1./L).*u_d).^2)+(((avgd./(L.^2)).*u_L).^2));

%Data input
rawdata = importdata('M008_Aluminum.dat','\t',5);
info = rawdata.data;
d = info(:,2); % displacment data
F = info(:,1); % force data
s_E = info(:,3); % strain data for E(not used in this but avaliable)

w0 = 0.518; %in change these values for each sample
t0 = 0.062; %in change these values for each sample
L = 6.5; %in change these values for each sample
tf = 0.0615;  % final thickness (change)

u_c = 0.0005; %in 
u_tm = 1/64; %in
u_d = 0.012 .* d;
u_f = 0.012 .* F;
A = eval(A);
u_stress = eval(u_stress);
u_strain = eval(u_strain);
u_E = eval(u_E);
u_E = mean(u_E(1:50)); % to find definite number for youngs modulus
stress = eval(stress);
strain = eval(strain);
U_strain_true = u_strain;
U_Stress_true = eval(U_Stress_true);

%Finding Values

% Finding ultimite stress and strain
ultimite_stress = max(stress); % could also be stress(find(max(stress) == stress))
ultimite_strain = strain(find(max(stress) == stress));

%Rupture stress and strain
Rupture_stress = stress(end);
Rupture_strain = strain(end);

%yield point 0.2 percent offset and Youngs Modulus(E)

p = polyfit(strain(1:50),stress(1:50),1);
E = p(1);
b = stress(1)-E.*(strain(1) .* 1.002);
p(2) = b;
offset = polyval(p,strain(1));
yield_intercept = find(offset >= stress,1);
yield_x = strain(yield_intercept);
yield_y = stress(yield_intercept);
start_index = find(offset >= stress(1),1);

% true values 
true_stress = stress .* (1+strain);
true_strain = log(1+strain);

True_Rupture_stress = true_stress(end);
True_Rupture_strain = true_strain(end);






% Resiliance and toughness

i = 1;
area_resilience = 0;
U_resilience = 0;
while i <=50
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain(i+1) - strain(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 
area_resilience = area_resilience + area_trap;

d_u_stress_r = u_stress(i) .* (strain(i+1) - strain(i));
U_resilience = d_u_stress_r + U_resilience;

i = i + 1; 
end

i = 1;
area_toughness = 0;
U_toughness = 0;

while i <= length(stress) - 1
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain(i+1) - strain(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 

area_toughness = area_toughness + area_trap;
d_u_stress_t = u_stress(i) .* (strain(i+1) - strain(i));
U_toughness = d_u_stress_t + U_toughness;

i = i + 1; 
end







% MMTS

stress_yg_MMTS = stress(1:1:50);
strain_yg_MMTS = strain(1:1:50);


estima = fitlm(strain_yg_MMTS,stress_yg_MMTS);
x1=estima.Coefficients;
youngs_modulus_MMTS=x1{'x1', 'Estimate'};


strain_m_MMTS = strain * 1.002;
strain_m_MMTS = strain_m_MMTS(1:100);
vertical_offset_MMTS = stress(1) - youngs_modulus_MMTS .* strain_m_MMTS(1);

stress_m_MMTS = zeros(100,1);
i = 1;
while i <= 100
stress_m_MMTS(i) = youngs_modulus_MMTS .* strain_m_MMTS(i) + vertical_offset_MMTS;
i = i+1;
end

stress_test = stress (1:100);

index_MMTS = find(stress_m_MMTS > stress_test, 1);
yield_stress_MMTS = stress(index_MMTS);
yield_strain_MMTS = strain(index_MMTS);


oset_plot_MMTS = stress_m_MMTS(1:index_MMTS);
strain_oset_MMTS = strain_m_MMTS(1:index_MMTS);



% Exten
stress_yg_exten = stress(1:1:50);
s_E = info(:,3);
strain_yg_exten = s_E(1:1:50);


estima = fitlm(strain_yg_exten,stress_yg_exten);
x1=estima.Coefficients;
youngs_modulus_extensometer=x1{'x1', 'Estimate'};


strain_m_exten = strain * 1.002;
strain_m_exten = strain_m_exten(1:100);
vertical_offset_exten = stress(1) - youngs_modulus_extensometer .* strain_m_exten(1);

stress_m_exten = zeros(100,1);
i = 1;
while i <= 100
stress_m_exten(i) = youngs_modulus_extensometer .* strain_m_exten(i) + vertical_offset_exten;
i = i+1;
end

stress_test = stress (1:100);

index_exten = find(stress_m_exten > stress_test, 1);
yield_stress_Extensometer = stress(index_exten);
yield_strain_Extensometer = strain(index_exten);


oset_plot_exten = stress_m_MMTS(1:index_exten);
strain_oset_exten = strain_m_MMTS(1:index_exten);






% Poisson Ratio


lateral_strain = ((tf-t0)/t0);
poissons = -(lateral_strain)/(strain(end));



%New Variables Aluminum
EE = strain;
FF = stress;
GG = u_stress;
HH = u_strain;



% True Uncertainty

rupture_index = length(true_stress);
U_Stress_True_rupture = U_Stress_true(rupture_index);
U_Strain_True_rupture = U_strain_true(rupture_index);

%Carbon Fiber 45---------------------------------------------------------------------------------------------------------------
clear d F w0 t0 u_c u_tm u_d u_f L
%DATA
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
Q = stress .* (1 + strain);

u_stress = sqrt((diff(stress,F).*u_f).^2 + (diff(stress,t0).*u_c).^2 + (diff(stress,w0).*u_c).^2);
u_strain = sqrt((diff(strain,d).*u_d).^2 + (diff(strain,L).*u_tm).^2);
u_E = sqrt((diff(E,F).*u_f).^2 + (diff(E,t0).*u_c).^2 + (diff(E,w0).*u_c).^2 + ...
    (diff(E,d).*u_d).^2 + (diff(E,L).*u_tm).^2);
U_Stress_true = sqrt((diff(Q,F).*u_f).^2 + (diff(Q,t0).*u_c).^2 + (diff(Q,w0).*u_c).^2 + ...
    (diff(Q,d).*u_d).^2 + (diff(Q,L).*u_tm).^2);

% manual differentiation
%u_stress = sqrt(((1/(t0.*w0)).*u_f).^2 + ((avgF/((t0.^2).*w0).*u_t).^2) + (-avgF/(t0.*w0.^2).*u_w).^2);
%u_strain = sqrt((((1./L).*u_d).^2)+(((avgd./(L.^2)).*u_L).^2));

%Data input
rawdata = importdata('M008_Carbon_Fiber_45.dat','\t',5);
info = rawdata.data;
d = info(:,2); % displacment data
F = info(:,1); % force data
s_E = info(:,3); % strain data for E(not used in this but avaliable)

w0 = 0.518; %in change these values for each sample
t0 = 0.060; %in change these values for each sample
L = 6.5; %in change these values for each sample
tf = 0.060;  % final thickness (change)

u_c = 0.0005; %in 
u_tm = 1/64; %in
u_d = 0.012 .* d;
u_f = 0.012 .* F;
A = eval(A);
u_stress = eval(u_stress);
u_strain = eval(u_strain);
u_E = eval(u_E);
u_E = mean(u_E(1:50)); % to find definite number for youngs modulus
stress = eval(stress);
strain = eval(strain);
U_strain_true = u_strain;
U_Stress_true = eval(U_Stress_true);

%Finding Values

% Finding ultimite stress and strain
ultimite_stress = max(stress); % could also be stress(find(max(stress) == stress))
ultimite_strain = strain(find(max(stress) == stress));

%Rupture stress and strain
Rupture_stress = stress(end);
Rupture_strain = strain(end);

%yield point 0.2 percent offset and Youngs Modulus(E)

p = polyfit(strain(1:50),stress(1:50),1);
E = p(1);
b = stress(1)-E.*(strain(1) .* 1.002);
p(2) = b;
offset = polyval(p,strain(1));
yield_intercept = find(offset >= stress,1);
yield_x = strain(yield_intercept);
yield_y = stress(yield_intercept);
start_index = find(offset >= stress(1),1);

% true values 
true_stress = stress .* (1+strain);
true_strain = log(1+strain);

True_Rupture_stress = true_stress(end);
True_Rupture_strain = true_strain(end);






% Resiliance and toughness

i = 1;
area_resilience = 0;
U_resilience = 0;
while i <=50
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain(i+1) - strain(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 
area_resilience = area_resilience + area_trap;

d_u_stress_r = u_stress(i) .* (strain(i+1) - strain(i));
U_resilience = d_u_stress_r + U_resilience;

i = i + 1; 
end

i = 1;
area_toughness = 0;
U_toughness = 0;

while i <= length(stress) - 1
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain(i+1) - strain(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 

area_toughness = area_toughness + area_trap;
d_u_stress_t = u_stress(i) .* (strain(i+1) - strain(i));
U_toughness = d_u_stress_t + U_toughness;

i = i + 1; 
end







% MMTS

stress_yg_MMTS = stress(1:1:50);
strain_yg_MMTS = strain(1:1:50);


estima = fitlm(strain_yg_MMTS,stress_yg_MMTS);
x1=estima.Coefficients;
youngs_modulus_MMTS=x1{'x1', 'Estimate'};


strain_m_MMTS = strain * 1.002;
strain_m_MMTS = strain_m_MMTS(1:100);
vertical_offset_MMTS = stress(1) - youngs_modulus_MMTS .* strain_m_MMTS(1);

stress_m_MMTS = zeros(100,1);
i = 1;
while i <= 100
stress_m_MMTS(i) = youngs_modulus_MMTS .* strain_m_MMTS(i) + vertical_offset_MMTS;
i = i+1;
end

stress_test = stress (1:100);

index_MMTS = find(stress_m_MMTS > stress_test, 1);
yield_stress_MMTS = stress(index_MMTS);
yield_strain_MMTS = strain(index_MMTS);


oset_plot_MMTS = stress_m_MMTS(1:index_MMTS);
strain_oset_MMTS = strain_m_MMTS(1:index_MMTS);



% Exten
stress_yg_exten = stress(1:1:50);
s_E = info(:,3);
strain_yg_exten = s_E(1:1:50);


estima = fitlm(strain_yg_exten,stress_yg_exten);
x1=estima.Coefficients;
youngs_modulus_extensometer=x1{'x1', 'Estimate'};


strain_m_exten = strain * 1.002;
strain_m_exten = strain_m_exten(1:100);
vertical_offset_exten = stress(1) - youngs_modulus_extensometer .* strain_m_exten(1);

stress_m_exten = zeros(100,1);
i = 1;
while i <= 100
stress_m_exten(i) = youngs_modulus_extensometer .* strain_m_exten(i) + vertical_offset_exten;
i = i+1;
end

stress_test = stress (1:100);

index_exten = find(stress_m_exten > stress_test, 1);
yield_stress_Extensometer = stress(index_exten);
yield_strain_Extensometer = strain(index_exten);


oset_plot_exten = stress_m_MMTS(1:index_exten);
strain_oset_exten = strain_m_MMTS(1:index_exten);






% Poisson Ratio


lateral_strain = ((tf-t0)/t0);
poissons = -(lateral_strain)/(strain(end));


%New Variables Carbon Fiber 45
II = strain;
JJ = stress;
KK = u_stress;
LL = u_strain;


% True Uncertainty

rupture_index = length(true_stress);
U_Stress_True_rupture = U_Stress_true(rupture_index);
U_Strain_True_rupture = U_strain_true(rupture_index);

%Carbon Fiber 90------------------------------------------------------------------------------------------------------------
clear d F w0 t0 u_c u_tm u_d u_f L
%DATA
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
Q = stress .* (1 + strain);

u_stress = sqrt((diff(stress,F).*u_f).^2 + (diff(stress,t0).*u_c).^2 + (diff(stress,w0).*u_c).^2);
u_strain = sqrt((diff(strain,d).*u_d).^2 + (diff(strain,L).*u_tm).^2);
u_E = sqrt((diff(E,F).*u_f).^2 + (diff(E,t0).*u_c).^2 + (diff(E,w0).*u_c).^2 + ...
    (diff(E,d).*u_d).^2 + (diff(E,L).*u_tm).^2);
U_Stress_true = sqrt((diff(Q,F).*u_f).^2 + (diff(Q,t0).*u_c).^2 + (diff(Q,w0).*u_c).^2 + ...
    (diff(Q,d).*u_d).^2 + (diff(Q,L).*u_tm).^2);

% manual differentiation
%u_stress = sqrt(((1/(t0.*w0)).*u_f).^2 + ((avgF/((t0.^2).*w0).*u_t).^2) + (-avgF/(t0.*w0.^2).*u_w).^2);
%u_strain = sqrt((((1./L).*u_d).^2)+(((avgd./(L.^2)).*u_L).^2));

%Data input
rawdata = importdata('M008_CF90.dat','\t',5);
info = rawdata.data;
d = info(:,2); % displacment data
F = info(:,1); % force data
s_E = info(:,3); % strain data for E(not used in this but avaliable)

w0 = 0.510; %in change these values for each sample
t0 = 0.072; %in change these values for each sample
L = 6.5; %in change these values for each sample
tf = 0.064;  % final thickness (change)

u_c = 0.0005; %in 
u_tm = 1/64; %in
u_d = 0.012 .* d;
u_f = 0.012 .* F;
A = eval(A);
u_stress = eval(u_stress);
u_strain = eval(u_strain);
u_E = eval(u_E);
u_E = mean(u_E(1:50)); % to find definite number for youngs modulus
stress = eval(stress);
strain = eval(strain);
U_strain_true = u_strain;
U_Stress_true = eval(U_Stress_true);

%Finding Values

% Finding ultimite stress and strain
ultimite_stress = max(stress); % could also be stress(find(max(stress) == stress))
ultimite_strain = strain(find(max(stress) == stress));

%Rupture stress and strain
Rupture_stress = stress(end);
Rupture_strain = strain(end);

%yield point 0.2 percent offset and Youngs Modulus(E)

p = polyfit(strain(1:50),stress(1:50),1);
E = p(1);
b = stress(1)-E.*(strain(1) .* 1.002);
p(2) = b;
offset = polyval(p,strain(1));
yield_intercept = find(offset >= stress,1);
yield_x = strain(yield_intercept);
yield_y = stress(yield_intercept);
start_index = find(offset >= stress(1),1);

% true values 
true_stress = stress .* (1+strain);
true_strain = log(1+strain);

True_Rupture_stress = true_stress(end);
True_Rupture_strain = true_strain(end);






% Resiliance and toughness

i = 1;
area_resilience = 0;
U_resilience = 0;
while i <=50
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain(i+1) - strain(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 
area_resilience = area_resilience + area_trap;

d_u_stress_r = u_stress(i) .* (strain(i+1) - strain(i));
U_resilience = d_u_stress_r + U_resilience;

i = i + 1; 
end

i = 1;
area_toughness = 0;
U_toughness = 0;

while i <= length(stress) - 1
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain(i+1) - strain(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 

area_toughness = area_toughness + area_trap;
d_u_stress_t = u_stress(i) .* (strain(i+1) - strain(i));
U_toughness = d_u_stress_t + U_toughness;

i = i + 1; 
end








% MMTS

stress_yg_MMTS = stress(1:1:50);
strain_yg_MMTS = strain(1:1:50);


estima = fitlm(strain_yg_MMTS,stress_yg_MMTS);
x1=estima.Coefficients;
youngs_modulus_MMTS=x1{'x1', 'Estimate'};


strain_m_MMTS = strain * 1.002;
strain_m_MMTS = strain_m_MMTS(1:100);
vertical_offset_MMTS = stress(1) - youngs_modulus_MMTS .* strain_m_MMTS(1);

stress_m_MMTS = zeros(100,1);
i = 1;
while i <= 100
stress_m_MMTS(i) = youngs_modulus_MMTS .* strain_m_MMTS(i) + vertical_offset_MMTS;
i = i+1;
end

stress_test = stress (1:100);

index_MMTS = find(stress_m_MMTS > stress_test, 1);
yield_stress_MMTS = stress(index_MMTS);
yield_strain_MMTS = strain(index_MMTS);


oset_plot_MMTS = stress_m_MMTS(1:index_MMTS);
strain_oset_MMTS = strain_m_MMTS(1:index_MMTS);



% Exten
stress_yg_exten = stress(1:1:50);
s_E = info(:,3);
strain_yg_exten = s_E(1:1:50);


estima = fitlm(strain_yg_exten,stress_yg_exten);
x1=estima.Coefficients;
youngs_modulus_extensometer=x1{'x1', 'Estimate'};


strain_m_exten = strain * 1.002;
strain_m_exten = strain_m_exten(1:100);
vertical_offset_exten = stress(1) - youngs_modulus_extensometer .* strain_m_exten(1);

stress_m_exten = zeros(100,1);
i = 1;
while i <= 100
stress_m_exten(i) = youngs_modulus_extensometer .* strain_m_exten(i) + vertical_offset_exten;
i = i+1;
end

stress_test = stress (1:100);

index_exten = find(stress_m_exten > stress_test, 1);
yield_stress_Extensometer = stress(index_exten);
yield_strain_Extensometer = strain(index_exten);


oset_plot_exten = stress_m_MMTS(1:index_exten);
strain_oset_exten = strain_m_MMTS(1:index_exten);





% Poisson Ratio


lateral_strain = ((tf-t0)/t0);
poissons = -(lateral_strain)/(strain(end));



%New Variables Carbon Fiber 90
MM = strain;
NN = stress;
OO = u_stress;
PP = u_strain;



% True Uncertainty

rupture_index = length(true_stress);
U_Stress_True_rupture = U_Stress_true(rupture_index);
U_Strain_True_rupture = U_strain_true(rupture_index);

%Plastic---------------------------------------------------------------------------------------------------------------
clear d F w0 t0 u_c u_tm u_d u_f L
%DATA
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
Q = stress .* (1 + strain);

u_stress = sqrt((diff(stress,F).*u_f).^2 + (diff(stress,t0).*u_c).^2 + (diff(stress,w0).*u_c).^2);
u_strain = sqrt((diff(strain,d).*u_d).^2 + (diff(strain,L).*u_tm).^2);
u_E = sqrt((diff(E,F).*u_f).^2 + (diff(E,t0).*u_c).^2 + (diff(E,w0).*u_c).^2 + ...
    (diff(E,d).*u_d).^2 + (diff(E,L).*u_tm).^2);
U_Stress_true = sqrt((diff(Q,F).*u_f).^2 + (diff(Q,t0).*u_c).^2 + (diff(Q,w0).*u_c).^2 + ...
    (diff(Q,d).*u_d).^2 + (diff(Q,L).*u_tm).^2);

% manual differentiation
%u_stress = sqrt(((1/(t0.*w0)).*u_f).^2 + ((avgF/((t0.^2).*w0).*u_t).^2) + (-avgF/(t0.*w0.^2).*u_w).^2);
%u_strain = sqrt((((1./L).*u_d).^2)+(((avgd./(L.^2)).*u_L).^2));

%Data input
rawdata = importdata('M008_22_plastic2.dat','\t',5);
info = rawdata.data;
d = info(:,2); % displacment data
F = info(:,1); % force data
s_E = info(:,3); % strain data for E(not used in this but avaliable)

w0 = 0.506; %in change these values for each sample
t0 = 0.137; %in change these values for each sample
L=6.5; %in change these values for each sample
tf = 0.135;  % final thickness (change)

u_c = 0.0005; %in 
u_tm = 1/64; %in
u_d = 0.012 .* d;
u_f = 0.012 .* F;
A = eval(A);
u_stress = eval(u_stress);
u_strain = eval(u_strain);
u_E = eval(u_E);
u_E = mean(u_E(1:50)); % to find definite number for youngs modulus
stress = eval(stress);
strain = eval(strain);
U_strain_true = u_strain;
U_Stress_true = eval(U_Stress_true);

%Finding Values

% Finding ultimite stress and strain
ultimite_stress = max(stress); % could also be stress(find(max(stress) == stress))
ultimite_strain = strain(find(max(stress) == stress));

%Rupture stress and strain
Rupture_stress = stress(end);
Rupture_strain = strain(end);

%yield point 0.2 percent offset and Youngs Modulus(E)

p = polyfit(strain(1:50),stress(1:50),1);
E = p(1);
b = stress(1)-E.*(strain(1) .* 1.002);
p(2) = b;
offset = polyval(p,strain(1));
yield_intercept = find(offset >= stress,1);
yield_x = strain(yield_intercept);
yield_y = stress(yield_intercept);
start_index = find(offset >= stress(1),1);

% true values 
true_stress = stress .* (1+strain);
true_strain = log(1+strain);

True_Rupture_stress = true_stress(end);
True_Rupture_strain = true_strain(end);






% Resiliance and toughness

i = 1;
area_resilience = 0;
U_resilience = 0;
while i <=50
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain(i+1) - strain(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 
area_resilience = area_resilience + area_trap;

d_u_stress_r = u_stress(i) .* (strain(i+1) - strain(i));
U_resilience = d_u_stress_r + U_resilience;

i = i + 1; 
end

i = 1;
area_toughness = 0;
U_toughness = 0;

while i <= length(stress) - 1
h1 = stress(i);
h2 = stress(i+1);
dh = h2 - h1; 
b = strain(i+1) - strain(i); 
area_rect = h1 * b; 
area_tri = 0.5 * (b * dh); 
area_trap = area_rect + area_tri; 

area_toughness = area_toughness + area_trap;
d_u_stress_t = u_stress(i) .* (strain(i+1) - strain(i));
U_toughness = d_u_stress_t + U_toughness;

i = i + 1; 
end







% MMTS

stress_yg_MMTS = stress(1:1:50);
strain_yg_MMTS = strain(1:1:50);


estima = fitlm(strain_yg_MMTS,stress_yg_MMTS);
x1=estima.Coefficients;
youngs_modulus_MMTS=x1{'x1', 'Estimate'};


strain_m_MMTS = strain * 1.002;
strain_m_MMTS = strain_m_MMTS(1:100);
vertical_offset_MMTS = stress(1) - youngs_modulus_MMTS .* strain_m_MMTS(1);

stress_m_MMTS = zeros(100,1);
i = 1;
while i <= 100
stress_m_MMTS(i) = youngs_modulus_MMTS .* strain_m_MMTS(i) + vertical_offset_MMTS;
i = i+1;
end

stress_test = stress (1:100);

index_MMTS = find(stress_m_MMTS > stress_test, 1);
yield_stress_MMTS = stress(index_MMTS);
yield_strain_MMTS = strain(index_MMTS);


oset_plot_MMTS = stress_m_MMTS(1:index_MMTS);
strain_oset_MMTS = strain_m_MMTS(1:index_MMTS);



% Exten
stress_yg_exten = stress(1:1:50);
s_E = info(:,3);
strain_yg_exten = s_E(1:1:50);


estima = fitlm(strain_yg_exten,stress_yg_exten);
x1=estima.Coefficients;
youngs_modulus_extensometer=x1{'x1', 'Estimate'};


strain_m_exten = strain * 1.002;
strain_m_exten = strain_m_exten(1:100);
vertical_offset_exten = stress(1) - youngs_modulus_extensometer .* strain_m_exten(1);

stress_m_exten = zeros(100,1);
i = 1;
while i <= 100
stress_m_exten(i) = youngs_modulus_extensometer .* strain_m_exten(i) + vertical_offset_exten;
i = i+1;
end

stress_test = stress (1:100);

index_exten = find(stress_m_exten > stress_test, 1);
yield_stress_Extensometer = stress(index_exten);
yield_strain_Extensometer = strain(index_exten);


oset_plot_exten = stress_m_MMTS(1:index_exten);
strain_oset_exten = strain_m_MMTS(1:index_exten);






% Poisson Ratio


lateral_strain = ((tf-t0)/t0);
poissons = -(lateral_strain)/(strain(end));




% Plotting 

figure(12)
plot(AA,BB)%Steel
hold on
plot(EE,FF)%Aluminum
hold on
plot(II,JJ)%CF45
hold on
plot(MM,NN)%CF90
hold on
plot(strain,stress,'blue')%Plastic
hold on
title('Strain vs Stress of All Materials')
xlabel('Strain')
ylabel('Stress(lb/in^2)')
%Steel
errorbar(AA(1:40:end),BB(1:40:end),CC(1:40:end),'vertical','.') % vertical errorbars for stress
errorbar(AA(1:40:end),BB(1:40:end),DD(1:40:end),'horizontal','.') % horizontal errorbars for strain
%Aluminum
errorbar(EE(1:40:end),FF(1:40:end),GG(1:40:end),'vertical','.') % vertical errorbars for stress
errorbar(EE(1:40:end),FF(1:40:end),HH(1:40:end),'horizontal','.') % horizontal errorbars for strain
%CF45
errorbar(II(1:40:end),JJ(1:40:end),KK(1:40:end),'vertical','.') % vertical errorbars for stress
errorbar(II(1:40:end),JJ(1:40:end),LL(1:40:end),'horizontal','.') % horizontal errorbars for strain
%CF90
errorbar(MM(1:40:end),NN(1:40:end),OO(1:40:end),'vertical','.') % vertical errorbars for stress
errorbar(MM(1:40:end),NN(1:40:end),PP(1:40:end),'horizontal','.') % horizontal errorbars for strain
%Plastic
errorbar(strain(1:40:end),stress(1:40:end),u_stress(1:40:end),'vertical','.') % vertical errorbars for stress
errorbar(strain(1:40:end),stress(1:40:end),u_strain(1:40:end),'horizontal','.') % horizontal errorbars for strain
%Legend is commented out because it covers the actual graph. We have a
%screenshot of the legend for the report but uncomment the legend to see it
%in matlab
%legend('Steel Strain vs Stress Curve','Aluminum Strain vs Stress Curve','Carbon Fiber 45 Strain vs Stress Curve','Carbon Fiber 90 Strain vs Stress Curve','Plastic Strain vs Stress Curve','Steel Uncertianty of Stress','Steel Uncertianty of Strain','Aluminum Uncertianty of Stress','Aluminum Uncertianty of Strain','Carbon Fiber 45 Uncertianty of Stress','Carbon Fiber 45 Uncertianty of Strain','Carbon Fiber 90 Uncertianty of Stress','Carbon Fiber 90 Uncertianty of Strain','Plastic Uncertianty of Stress','Plastic Uncertianty of Strain','Location','southeast')

% True Uncertainty

rupture_index = length(true_stress);
U_Stress_True_rupture = U_Stress_true(rupture_index);
U_Strain_True_rupture = U_strain_true(rupture_index);