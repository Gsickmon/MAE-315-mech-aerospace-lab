% Garrett Sickmon
% MAE 315
clear
close all
clc

%% wind tunnel calibration 
clear 
disp('wind tunnel calibration')
uf = 0.5; % Hz
up = 0.0005 .* 133.32; % uncertianty of pressure in pa

syms avg_p_pa p_atm R f temp_K x

rho_eq = (p_atm)./(R.*temp_K); % kg/m^3
velocity_eq = sqrt((2.*avg_p_pa)./rho_eq); % m/s
v_output = 0.003.*(x.^2) + 0.79.*x - 0.5210; % calibration equation

uv_eq = sqrt((diff(velocity_eq,avg_p_pa).*up).^2 + (diff(velocity_eq,f).*uf).^2);

f = [0,5,10,15,20,25,30,35,40]; %frequency in Hz
avg_ptorr =[0,0.0383,0.245,0.707,1.338,2.10,3.14,4.45,6.08];  % average pressure in torr taken from excel spreadsheet theat averaged values internally
avg_p_pa = (133.32.*avg_ptorr); % pa pressure
p_atm = 101325; % pa
% calculated in gage pressure 
temp_F = [73.3,61.5,60.8,61.2,61.7,62.2,62.7,63.9,65.3]; % fahrenheit
temp_C = (temp_F-32).* (5/9); % Celsius
temp_K = temp_C + 273.15; % Kelvin
R = 287; % Ideal gas constant
% Ideal gas law
%rho = (p_atm)./(R.*temp_K); % kg/m^3
% dymanic pressure equation from bernoullis 
%velocity = sqrt((2.*avg_p_pa)./rho); % m/s
%uncertainties
rho = eval(rho_eq);
velocity = eval(velocity_eq);
uv = eval(uv_eq);
%uv = sqrt((up.^2)+(uf.^2));

%plotting
figure(1)
plot(f,velocity,'*')
hold on
errorbar(f(1:1:end),velocity(1:1:end),uf(1:1:end),'horizontal','.')
errorbar(f(1:1:end),velocity(1:1:end),uv(1:1:end),'vertical','.')
grid()
xlabel('Frequency(Hz)')
ylabel('Velocity(m/s)')
title('Frequency Vs. Velocity')
hold on
% regression

[fitted ,s] = polyfit(f,velocity,2); % second order regression
x_calibration = linspace(0,50,100);
y_calibration = fitted(1).*x_calibration.^2+fitted(2).*x_calibration+fitted(3);

plot(x_calibration,y_calibration);

% R^2 value
R_squared = 1 - (s.normr/norm(velocity-mean(velocity)))^2;
hold on
printed_string = sprintf('The second order fit is %.3fx^2 + %.2fx + %.4f',fitted(1),fitted(2),fitted(3));
text(-9,35,printed_string);
printed_string2 = sprintf('Rsquared value is %x',round(R_squared,4));
text(-9,30,printed_string2);
hold off
% Calibration Equation
x = (1:1:50);
%v_output = 0.003.*(x.^2) + 0.79.*x - 0.5210;
v_output = eval(v_output);
fprintf('40Hz velocity: %g m/s +- %g\n',v_output(1,40), uv(1,9))
fprintf('45Hz velocity: %g m/s +- %g\n',v_output(1,45), (uv(1,8)+uv(1,9))./2)

%% vortex shedding 30Hz
clear
disp('vortex shedding 30Hz')
% uncertianty is 1/2 delta t which we caluculate from slides



opts = detectImportOptions('M002_Lab5_30Hz_test1.txt');
opts.Delimiter = {' ',','}; 
Data = readtable('M002_Lab5_30Hz_test1.txt',opts);
t = Data.t;
Amplitude = Data.ExtraVar1;
t = str2double(t); % time
Amplitude = str2double(Amplitude); % amplitude


amp_fft = abs(fft(Amplitude));
sample_rate = 6400;
n_samples = 16384;
deltaf = sample_rate./ n_samples;
%i=1;
f = (0 :(n_samples/2)-1).*deltaf; % mirrors so divide length by two
amp_fft = amp_fft(1:length(f)); % also mirrors so its now the same length as f

%f = zeros(1,length(amp_fft));
%while i < length(amp_fft)
%    f(i) = (deltaf(i).*i);
%    i = i + 1;
%end
%f = f;


% plots
figure(2)
subplot(2,1,1)
plot(t,Amplitude)
hold on 
xlabel('Time(s)')
ylabel('Amplitude(V)')
title('Time vs. Amplitude 30Hz')
hold off

figure(2)
subplot(2,1,2)
plot(f,amp_fft)
xlabel('Frequency(Hz)')
ylabel('Amplitude(V)')
title('Frequency vs. Amplitude 30Hz')
 

% frequency of tunnel is 30Hz

%associted velocity 
x = 30; % 30Hz
velocity = 0.003.*(x.^2) + 0.79.*x - 0.5210;% calibration equation

% shedding frequency
Shedding_f = f(find(max(amp_fft) == amp_fft)); %Hz
uac = 0.5.*(deltaf); % accelerometer uncertianty 

%strouhal number -> fD/U
D = 48.12./1000; %cylinder diameter in meters
Strouhal = (Shedding_f .* D)./velocity;

% Reynold's number 
mu = 1.79E-5; %Dynamic viscosity form fluids text book
rho = 1.2165; % kg/m^3 from first section 30Hz
Re = (rho.*velocity.*D)./mu ;

% Empirical shedding frequency
Empirical_shedding_f = (velocity.*0.198.*(1-(19.7./Re)))./D;

%uncertianties

uac = 0.5.*(deltaf); % accelerometer uncertianty 
uf = 0.5; %hz
uD = 0.005./1000; % caliper uncertianty in meters
syms sD xeq uv Shedding_f1

velocity_eq = 0.003.*(xeq.^2) + 0.79.*xeq - 0.5210;
Strouhal_eq = (Shedding_f1 .* sD)./velocity_eq;
Re_eq = (rho.*velocity_eq.*sD)./mu ;
Empirical_shedding_f_eq = (velocity_eq.*0.198.*(1-(19.7./Re_eq)))./sD;

uv_eq = sqrt((diff(velocity_eq,xeq).*uf).^2);
us_eq = sqrt((diff(Strouhal_eq,sD).*uD).^2 + (diff(Strouhal_eq,xeq).*uf).^2 + (diff(Strouhal_eq,Shedding_f1).*uac).^2);
uRe_eq = sqrt((diff(Re_eq,sD).*uD).^2 + (diff(Re_eq,xeq).*uf).^2); 
ues_eq = sqrt((diff(Empirical_shedding_f_eq,sD).*uD).^2 + (diff(Empirical_shedding_f_eq,xeq).*uf).^2 + (diff(Empirical_shedding_f_eq,Shedding_f1).*uac).^2);
Shedding_f1 = Shedding_f;
xeq = x;
sD = D;

uv = eval(uv_eq);% velocity uncertainty from regrssion equation 
us = eval(us_eq);% strouhal uncertainty
uRe = eval(uRe_eq);% renyolds number uncertainty 
ues = eval(ues_eq);% empircal shedding frequency uncertainty 

fprintf('Frequency of tunnel(30Hz): %gHz\n',x)
fprintf('Associated velocity(30Hz): %g m/s +- %g\n',velocity, uv)
fprintf('Shedding Frequency from Accelerometer(30Hz): %gHz +- %g\n',Shedding_f, uac)
fprintf('Strouhal Number(30Hz): %g +- %g\n',Strouhal, us)
fprintf('Reynold’s Number(30Hz): %g +- %g\n',Re, uRe)
fprintf('Empirical Shedding Frequency(30Hz): %gHz +- %g\n',Empirical_shedding_f, ues)

%% vortex shedding 45Hz
clear
disp('vortex shedding 45Hz')
% uncertianty is 1/2 delta t which we caluculate from slides



opts = detectImportOptions('M002_Lab5_45Hz_test1.txt');
opts.Delimiter = {' ',','}; 
Data = readtable('M002_Lab5_45Hz_test1.txt',opts);
t = Data.t;
Amplitude = Data.ExtraVar1;
t = str2double(t); % time
Amplitude = str2double(Amplitude); % amplitude


amp_fft = abs(fft(Amplitude));
sample_rate = 6400;
n_samples = 16384;
deltaf = sample_rate./ n_samples;
%i=1;
f = (0 :(n_samples/2)-1).*deltaf; % mirrors so divide length by two
amp_fft = amp_fft(1:length(f)); % also mirrors so its now the same length as f

%f = zeros(1,length(amp_fft));
%while i < length(amp_fft)
%    f(i) = (deltaf(i).*i);
%    i = i + 1;
%end
%f = f;


% plots
figure(3)
subplot(2,1,1)
plot(t,Amplitude)
hold on 
xlabel('Time(s)')
ylabel('Amplitude(V)')
title('Time vs. Amplitude 45Hz')
hold off

figure(3)
subplot(2,1,2)
plot(f,amp_fft)
xlabel('Frequency(Hz)')
ylabel('Amplitude(V)')
title('Frequency vs. Amplitude 45Hz')
 

% frequency of tunnel is 30Hz

%associted velocity 
x = 45; % 45Hz
velocity = 0.003.*(x.^2) + 0.79.*x - 0.5210;% calibration equation

% shedding frequency
Shedding_f = f(find(max(amp_fft) == amp_fft)); %Hz
uac = 0.5.*(deltaf); % accelerometer uncertianty 

%strouhal number -> fD/U
D = 48.12./1000; %cylinder diameter in meters
Strouhal = (Shedding_f .* D)./velocity;

% Reynold's number 
mu = 1.79E-5; %Dynamic viscosity from fluids text book
rho = 1.2105; % kg/m^3 from first section 40Hz
Re = (rho.*velocity.*D)./mu ;

% Empirical shedding frequency
Empirical_shedding_f = (velocity.*0.198.*(1-(19.7./Re)))./D;

%uncertianties

uac = 0.5.*(deltaf); % accelerometer uncertianty 
uf = 0.5; %hz
uD = 0.005./1000; % caliper uncertianty in meters
syms sD xeq uv Shedding_f1

velocity_eq = 0.003.*(xeq.^2) + 0.79.*xeq - 0.5210;
Strouhal_eq = (Shedding_f1 .* sD)./velocity_eq;
Re_eq = (rho.*velocity_eq.*sD)./mu ;
Empirical_shedding_f_eq = (velocity_eq.*0.198.*(1-(19.7./Re_eq)))./sD;

uv_eq = sqrt((diff(velocity_eq,xeq).*uf).^2);
us_eq = sqrt((diff(Strouhal_eq,sD).*uD).^2 + (diff(Strouhal_eq,xeq).*uf).^2 + (diff(Strouhal_eq,Shedding_f1).*uac).^2);
uRe_eq = sqrt((diff(Re_eq,sD).*uD).^2 + (diff(Re_eq,xeq).*uf).^2); 
ues_eq = sqrt((diff(Empirical_shedding_f_eq,sD).*uD).^2 + (diff(Empirical_shedding_f_eq,xeq).*uf).^2 + (diff(Empirical_shedding_f_eq,Shedding_f1).*uac).^2);
Shedding_f1 = Shedding_f;
xeq = x;
sD = D;

uv = eval(uv_eq);% velocity uncertainty from regrssion equation 
us = eval(us_eq);% strouhal uncertainty
uRe = eval(uRe_eq);% renyolds number uncertainty 
ues = eval(ues_eq);% empircal shedding frequency uncertainty 

fprintf('Frequency of tunnel(45Hz): %gHz\n',x)
fprintf('Associated velocity(45Hz): %g m/s +- %g\n',velocity, uv)
fprintf('Shedding Frequency from Accelerometer(45Hz): %gHz +- %g\n',Shedding_f, uac)
fprintf('Strouhal Number(45Hz): %g +- %g\n',Strouhal, us)
fprintf('Reynold’s Number(45Hz): %g +- %g\n',Re, uRe)
fprintf('Empirical Shedding Frequency(45Hz): %gHz +- %g\n',Empirical_shedding_f, ues)
%% Cylinder Wake 25Hz
clear
disp('Cylinder Wake 25Hz')
% data in psi convert to pa 1psi = 6894.757 pa
% x axis is pressure tap location in cm
% y axis is velocity in m/s
% avg each pressure rows pressure 

info = importdata('M005_Lab5_Part3_25Hz_without_cylindar','\t',0);

% p1 = info(:,1); p1 is bad data 
%p2 = info(:,2); bad data
p3 = info(:,3);
p4 = info(:,4);
p5 = info(:,5);
p6 = info(:,6);
p7 = info(:,7);
p8 = info(:,8);
p9 = info(:,9);
p10 = info(:,10);
p11 = info(:,11);
p12 = info(:,12);
p13 = info(:,13);
p14 = info(:,14);
p15 = info(:,15);
p16 = info(:,16);
p17 = info(:,17);
p18 = info(:,18);
p19 = info(:,19);
p20 = info(:,20);
p21 = info(:,21);
p22 = info(:,22);
%avg and convert to pa
%p2 = mean(p2).*6894.757; bad data
p3 = mean(p3).*6894.757;
p4 = mean(p4).*6894.757;
p5 = mean(p5).*6894.757;
p6 = mean(p6).*6894.757;
p7 = mean(p7).*6894.757;
p8 = mean(p8).*6894.757;
p9 = mean(p9).*6894.757;
p10 = mean(p10).*6894.757;
p11 = mean(p11).*6894.757;
p12 = mean(p12).*6894.757;
p13 = mean(p13).*6894.757;
p14 = mean(p14).*6894.757;
p15 = mean(p15).*6894.757;
p16 = mean(p16).*6894.757;
p17 = mean(p17).*6894.757;
p18 = mean(p18).*6894.757;
p19 = mean(p19).*6894.757;
p20 = mean(p20).*6894.757;
p21 = mean(p21).*6894.757;
p22 = mean(p22).*6894.757;

p = [p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22];
rho = 1.2177; % at 25Hz
velocity = sqrt((2.*p)./rho); %m/s 

% data with cylinder
infoc = importdata('M005_Lab5_Part3_25Hz_with_cylindar','\t',0);

% p1 = info(:,1); p1 is bad data 
%p2c = infoc(:,2); bad data
p3c = infoc(:,3);
p4c = infoc(:,4);
p5c = infoc(:,5);
p6c = infoc(:,6);
p7c = infoc(:,7);
p8c = infoc(:,8);
p9c = infoc(:,9);
p10c = infoc(:,10);
p11c = infoc(:,11);
p12c = infoc(:,12);
p13c = infoc(:,13);
p14c = infoc(:,14);
p15c = infoc(:,15);
p16c = infoc(:,16);
p17c = infoc(:,17);
p18c = infoc(:,18);
p19c = infoc(:,19);
p20c = infoc(:,20);
p21c = infoc(:,21);
p22c = infoc(:,22);
%avg and convert to pa
%p2c = mean(p2c).*6894.757; bad data
p3c = mean(p3c).*6894.757;
p4c = mean(p4c).*6894.757;
p5c = mean(p5c).*6894.757;
p6c = mean(p6c).*6894.757;
p7c = mean(p7c).*6894.757;
p8c = mean(p8c).*6894.757;
p9c = mean(p9c).*6894.757;
p10c = mean(p10c).*6894.757;
p11c = mean(p11c).*6894.757;
p12c = mean(p12c).*6894.757;
p13c = mean(p13c).*6894.757;
p14c = mean(p14c).*6894.757;
p15c = mean(p15c).*6894.757;
p16c = mean(p16c).*6894.757;
p17c = mean(p17c).*6894.757;
p18c = mean(p18c).*6894.757;
p19c = mean(p19c).*6894.757;
p20c = mean(p20c).*6894.757;
p21c = mean(p21c).*6894.757;
p22c = mean(p22c).*6894.757;

pc = [p3c,p4c,p5c,p6c,p7c,p8c,p9c,p10c,p11c,p12c,p13c,p14c,p15c,p16c,p17c,p18c,p19c,p20c,p21c,p22c];
velocityc = sqrt((2.*pc)./rho); %m/s 



% distance of long gaps = 19.21mm
% distance of short gaps = 9.91
lg = 19.21./1000; %m
sg = 9.91./1000;%m

i = 1;
d = zeros(1,20);
% distance between tubes
while i <= 20
    
    if i <=7
    d(i) = d(i) + lg.*i;
    end
    if i > 7 & i <= 13
        d(i) = d(i) + sg.*i +lg.*7- sg.*7;
    end
    if i > 13
        d(i) = d(i) + lg.*i +d(13) - lg.*13;
    end
    i = i + 1;
end
d=d; % meters (when graphing multiply by 100 to get cm) 

% plotting

figure()
plot(d.*100,velocity)
hold on
plot(d.*100,velocityc)
xlabel('Pressure Rake Spacing (cm)')
ylabel('Velocity(m/s)')
title('Pressure Rake Spacing Vs. Velocity(25Hz)')
legend('Without Cylinder', 'With Cylinder',Location= 'southeast')

%calculating data
x=25; % tunnel frequency
fprintf('Frequency of tunnel(25Hz): %gHz\n',x)

uf = 0.5; %hz
uD = 0.005./1000; % caliper uncertianty in meters
% rho already introduced
mu = 1.79E-5; %Dynamic viscosity
D = 48.12./1000; %cylinder diameter in meters

%trapazoidal integration
Fd_eq = rho.*velocityc.*(velocity - velocityc);
i = 1;
area_resilience = 0;
U_resilience = 0;

% x =d and y = Fd_eq
Area = 0;
while i <= length(d) - 1
Area = Area + ((Fd_eq(i)+Fd_eq(i+1)).*(d(i+1)- d(i)).*0.5); % area of trapazoid

i = i + 1; 
end
Fd = Area;
uf = 0.5; % Hz
up = 0.0005 .* 133.32; % uncertianty of pressure in pa
uL = ((1/64)./12)./ 3.28;

syms sD peq L w A xeq pceq Fd1 uFd

velocity_eq = 0.003.*(xeq.^2) + 0.79.*xeq - 0.5210;
velocity_eq2 = sqrt((2.*peq)./rho); % m/s
velocityc_eq2 = sqrt((2.*pceq)./rho); % m/s
Re_eq = (rho.*velocity_eq.*sD)./mu ;
Fd_eq1 = rho.*velocityc.*(velocity - velocityc);
% w is the width of the waverake 
CD_cyl_eq = Fd./(0.5.*rho.*((mean(velocity).^2).*A));


uRe_eq = sqrt((diff(Re_eq,sD).*uD).^2 + (diff(Re_eq,xeq).*uf).^2 ); 
uFd_eq = sqrt((diff(velocity_eq2,peq).*up).^2 + (diff(velocityc_eq2,pceq).*up).^2);
uCd_eq = sqrt((diff(CD_cyl_eq,peq).*up).^2 + (diff(CD_cyl_eq,sD).*uD).^2 + (diff(CD_cyl_eq,Fd1).*uFd).^2 + (diff(CD_cyl_eq,L).*uL).^2);
sD = D;
peq = p;
pceq = pc;
xeq = x;
Fd1 = Fd;
w = (d(1,end))./100; % width of wave rake in meters
L = (24./12)./ 3.28; % from inchs to meters
A = sD.*L; % frontal area diameter*length

uRe = eval(uRe_eq);% renyolds number uncertainty 
uFd = eval(uFd_eq);
uCd = eval(uCd_eq);
RE = eval(Re_eq);
CD_cyl = eval(CD_cyl_eq);

fprintf('Drag from cylinder(25Hz): %gN +- %g\n',Fd, uFd(1,1))
fprintf('Coefficient of Drag(25Hz): %g +- %g\n',CD_cyl, uCd)
fprintf('Reynolds number(25Hz): %g +- %g\n',RE, uRe)

%% Cylinder Wake 45Hz
clear
disp('Cylinder Wake 45Hz')
% data in psi convert to pa 1psi = 6894.757 pa
% x axis is pressure tap location in cm
% y axis is velocity in m/s
% avg each pressure rows pressure 

info = importdata('M005_Lab5_Part3_45Hz_without_cylindar','\t',0);

% p1 = info(:,1); p1 is bad data 
%p2 = info(:,2); bad data
p3 = info(:,3);
p4 = info(:,4);
p5 = info(:,5);
p6 = info(:,6);
p7 = info(:,7);
p8 = info(:,8);
p9 = info(:,9);
p10 = info(:,10);
p11 = info(:,11);
p12 = info(:,12);
p13 = info(:,13);
p14 = info(:,14);
p15 = info(:,15);
p16 = info(:,16);
p17 = info(:,17);
p18 = info(:,18);
p19 = info(:,19);
p20 = info(:,20);
p21 = info(:,21);
p22 = info(:,22);
%avg and convert to pa
%p2 = mean(p2).*6894.757; bad data
p3 = mean(p3).*6894.757;
p4 = mean(p4).*6894.757;
p5 = mean(p5).*6894.757;
p6 = mean(p6).*6894.757;
p7 = mean(p7).*6894.757;
p8 = mean(p8).*6894.757;
p9 = mean(p9).*6894.757;
p10 = mean(p10).*6894.757;
p11 = mean(p11).*6894.757;
p12 = mean(p12).*6894.757;
p13 = mean(p13).*6894.757;
p14 = mean(p14).*6894.757;
p15 = mean(p15).*6894.757;
p16 = mean(p16).*6894.757;
p17 = mean(p17).*6894.757;
p18 = mean(p18).*6894.757;
p19 = mean(p19).*6894.757;
p20 = mean(p20).*6894.757;
p21 = mean(p21).*6894.757;
p22 = mean(p22).*6894.757;

p = [p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22];
rho = 1.2105; % kg/m^3 from first section 40Hz
velocity = sqrt((2.*p)./rho); %m/s 

% data with cylinder
infoc = importdata('M005_Lab5_Part3_45Hz_with_cylindar','\t',0);

% p1 = info(:,1); p1 is bad data 
%p2c = infoc(:,2); bad data
p3c = infoc(:,3);
p4c = infoc(:,4);
p5c = infoc(:,5);
p6c = infoc(:,6);
p7c = infoc(:,7);
p8c = infoc(:,8);
p9c = infoc(:,9);
p10c = infoc(:,10);
p11c = infoc(:,11);
p12c = infoc(:,12);
p13c = infoc(:,13);
p14c = infoc(:,14);
p15c = infoc(:,15);
p16c = infoc(:,16);
p17c = infoc(:,17);
p18c = infoc(:,18);
p19c = infoc(:,19);
p20c = infoc(:,20);
p21c = infoc(:,21);
p22c = infoc(:,22);
%avg and convert to pa
%p2c = mean(p2c).*6894.757; bad data
p3c = mean(p3c).*6894.757;
p4c = mean(p4c).*6894.757;
p5c = mean(p5c).*6894.757;
p6c = mean(p6c).*6894.757;
p7c = mean(p7c).*6894.757;
p8c = mean(p8c).*6894.757;
p9c = mean(p9c).*6894.757;
p10c = mean(p10c).*6894.757;
p11c = mean(p11c).*6894.757;
p12c = mean(p12c).*6894.757;
p13c = mean(p13c).*6894.757;
p14c = mean(p14c).*6894.757;
p15c = mean(p15c).*6894.757;
p16c = mean(p16c).*6894.757;
p17c = mean(p17c).*6894.757;
p18c = mean(p18c).*6894.757;
p19c = mean(p19c).*6894.757;
p20c = mean(p20c).*6894.757;
p21c = mean(p21c).*6894.757;
p22c = mean(p22c).*6894.757;

pc = [p3c,p4c,p5c,p6c,p7c,p8c,p9c,p10c,p11c,p12c,p13c,p14c,p15c,p16c,p17c,p18c,p19c,p20c,p21c,p22c];
velocityc = sqrt((2.*pc)./rho); %m/s 



% distance of long gaps = 19.21mm
% distance of short gaps = 9.91
lg = 19.21./1000; %m
sg = 9.91./1000;%m

i = 1;
d = zeros(1,20);
% distance between tubes
while i <= 20
    
    if i <=7
    d(i) = d(i) + lg.*i;
    end
    if i > 7 & i <= 13
        d(i) = d(i) + sg.*i +lg.*7- sg.*7;
    end
    if i > 13
        d(i) = d(i) + lg.*i +d(13) - lg.*13;
    end
    i = i + 1;
end
d=d; % meters (when graphing multiply by 100 to get cm) 

% plotting

figure()
plot(d.*100,velocity)
hold on
plot(d.*100,velocityc)
xlabel('Pressure Rake Spacing (cm)')
ylabel('Velocity(m/s)')
title('Pressure Rake Spacing Vs. Velocity(45Hz)')
legend('Without Cylinder', 'With Cylinder',Location= 'southeast')

%calculating data
x=45; % tunnel frequency
fprintf('Frequency of tunnel(45Hz): %gHz\n',x)

uf = 0.5; %hz
uD = 0.005./1000; % caliper uncertianty in meters
% rho already introduced
mu = 1.79E-5; %Dynamic viscosity
D = 48.12./1000; %cylinder diameter in meters

%trapazoidal integration
Fd_eq = rho.*velocityc.*(velocity - velocityc);
i = 1;
area_resilience = 0;
U_resilience = 0;

% x =d and y = Fd_eq
Area = 0;
while i <= length(d) - 1
Area = Area + ((Fd_eq(i)+Fd_eq(i+1)).*(d(i+1)- d(i)).*0.5); % area of trapazoid

i = i + 1; 
end
Fd = Area;
uf = 0.5; % Hz
up = 0.0005 .* 133.32; % uncertianty of pressure in pa
uL = ((1/64)./12)./ 3.28;

syms sD peq L w A xeq pceq Fd1 uFd

velocity_eq = 0.003.*(xeq.^2) + 0.79.*xeq - 0.5210;
velocity_eq2 = sqrt((2.*peq)./rho); % m/s
velocityc_eq2 = sqrt((2.*pceq)./rho); % m/s
Re_eq = (rho.*velocity_eq.*sD)./mu ;
Fd_eq1 = rho.*velocityc.*(velocity - velocityc);
% w is the width of the waverake 
CD_cyl_eq = Fd./(0.5.*rho.*((mean(velocity).^2).*A));


uRe_eq = sqrt((diff(Re_eq,sD).*uD).^2 + (diff(Re_eq,xeq).*uf).^2 ); 
uFd_eq = sqrt((diff(velocity_eq2,peq).*up).^2 + (diff(velocityc_eq2,pceq).*up).^2);
uCd_eq = sqrt((diff(CD_cyl_eq,peq).*up).^2 + (diff(CD_cyl_eq,sD).*uD).^2 + (diff(CD_cyl_eq,Fd1).*uFd).^2 + (diff(CD_cyl_eq,L).*uL).^2);
sD = D;
peq = p;
pceq = pc;
xeq = x;
Fd1 = Fd;
w = (d(1,end))./100; % width of wave rake in meters
L = (24./12)./ 3.28; % from inchs to meters
A = sD.*L; % frontal area diameter*length

uRe = eval(uRe_eq);% renyolds number uncertainty 
uFd = eval(uFd_eq);
uCd = eval(uCd_eq);
RE = eval(Re_eq);
CD_cyl = eval(CD_cyl_eq);

fprintf('Drag from cylinder(45Hz): %gN +- %g\n',Fd, uFd(1,1))
fprintf('Coefficient of Drag(45Hz): %g +- %g\n',CD_cyl, uCd)
fprintf('Reynolds number(45Hz): %g +- %g\n',RE, uRe)

%% Boundary Layer 30Hz 

clear
disp('Boundary Layer 30Hz')
% data in psi convert to pa 1psi = 6894.757 pa
% y axis is pressure tap location in m
% x axis is velocity in m/s
% avg each pressure rows pressure 
% upstream velocity 
info = importdata('m008_upstream_30hz','\t',0);

% p1 = info(:,1); p1 is bad data 
p2 = info(:,2); 
p3 = info(:,3);
p4 = info(:,4);
p5 = info(:,5);
p6 = info(:,6);
p7 = info(:,7);
p8 = info(:,8);
p9 = info(:,9);
p10 = info(:,10);
p11 = info(:,11);
p12 = info(:,12);
p13 = info(:,13);
p14 = info(:,14);
p15 = info(:,15);
p16 = info(:,16);
p17 = info(:,17);
p18 = info(:,18);
p19 = info(:,19);
p20 = info(:,20);
p21 = info(:,21);

%avg and convert to pa
p2 = mean(p2).*6894.757; 
p3 = mean(p3).*6894.757;
p4 = mean(p4).*6894.757;
p5 = mean(p5).*6894.757;
p6 = mean(p6).*6894.757;
p7 = mean(p7).*6894.757;
p8 = mean(p8).*6894.757;
p9 = mean(p9).*6894.757;
p10 = mean(p10).*6894.757;
p11 = mean(p11).*6894.757;
p12 = mean(p12).*6894.757;
p13 = mean(p13).*6894.757;
p14 = mean(p14).*6894.757;
p15 = mean(p15).*6894.757;
p16 = mean(p16).*6894.757;
p17 = mean(p17).*6894.757;
p18 = mean(p18).*6894.757;
p19 = mean(p19).*6894.757;
p20 = mean(p20).*6894.757;
p21 = mean(p21).*6894.757;


p = [p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21];
rho = 1.2165; % kg/m^3 from first section 30Hz
velocity = sqrt((2.*p)./rho); %m/s 

% downstream velocity 
infoc = importdata('m008_downstream','\t',0);

% p1 = info(:,1); p1 is bad data 
p2c = infoc(:,2); 
p3c = infoc(:,3);
p4c = infoc(:,4);
p5c = infoc(:,5);
p6c = infoc(:,6);
p7c = infoc(:,7);
p8c = infoc(:,8);
p9c = infoc(:,9);
p10c = infoc(:,10);
p11c = infoc(:,11);
p12c = infoc(:,12);
p13c = infoc(:,13);
p14c = infoc(:,14);
p15c = infoc(:,15);
p16c = infoc(:,16);
p17c = infoc(:,17);
p18c = infoc(:,18);
p19c = infoc(:,19);
p20c = infoc(:,20);
p21c = infoc(:,21);

%avg and convert to pa
p2c = mean(p2c).*6894.757; 
p3c = mean(p3c).*6894.757;
p4c = mean(p4c).*6894.757;
p5c = mean(p5c).*6894.757;
p6c = mean(p6c).*6894.757;
p7c = mean(p7c).*6894.757;
p8c = mean(p8c).*6894.757;
p9c = mean(p9c).*6894.757;
p10c = mean(p10c).*6894.757;
p11c = mean(p11c).*6894.757;
p12c = mean(p12c).*6894.757;
p13c = mean(p13c).*6894.757;
p14c = mean(p14c).*6894.757;
p15c = mean(p15c).*6894.757;
p16c = mean(p16c).*6894.757;
p17c = mean(p17c).*6894.757;
p18c = mean(p18c).*6894.757;
p19c = mean(p19c).*6894.757;
p20c = mean(p20c).*6894.757;
p21c = mean(p21c).*6894.757;


pc = [p2c,p3c,p4c,p5c,p6c,p7c,p8c,p9c,p10c,p11c,p12c,p13c,p14c,p15c,p16c,p17c,p18c,p19c,p20c,p21c];
velocityc = sqrt((2.*pc)./rho); %m/s 




lg = 0.00852; %m
sg = 0.00364;%m

i = 1;
d = zeros(1,20);
% distance between tubes
while i <= 20 
    
    if i <=10
    d(i) = d(i) + sg.*i;
    end
    if i > 10
        d(i) = d(i) + lg.*i +d(10) - lg.*10;
    end
 
    i = i + 1;
end
d=d; % meters (when graphing multiply by 100 to get cm) 

% plotting

figure()
plot(velocity,d)
hold on
plot(velocityc,d)
grid()
ylabel('Pressure Rake Spacing (m)')
xlabel('Velocity(m/s)')
title('Pressure Rake Spacing Vs. Velocity(30Hz)')
legend('Upstream Velocity', 'Downstream velocity',Location='northwest')

%% Airfoil 
clear 
%Output1 is Angle of Attack (degrees)
%Output2 is the Normal Force (lbf)
%Output3 is the Side Force (lbf) (also, not needed)
%Output4 is the Axial Force (lbf)
%Output5 is the Pitching Moment (lbf-in)

% convert lbs to newtons 1lb = 4.448Newtons
% conver inch lbs to newton meters N m  1 NM = 8.851 inch lbs
%constants 
mu = 1.79E-5; %Dynamic viscosity form fluids text book
x30 = 30;
x45 = 45;
velocity_30 = 0.003.*(x30.^2) + 0.79.*x30 - 0.5210;
velocity_45 = 0.003.*(x45.^2) + 0.79.*x45 - 0.5210;
rho_30 = 1.2165; % kg/m^3 from first section 30Hz
rho_45 = 1.2105; % kg/m^3 
span = (9.875./12)./3.281; % inchs to meters
chord = 89.78./1000; % mm to meters
thickness = 13.78./1000; % mm to meters
A = span.*chord;%platform area
%A_f = thickness.*span; % frontal area

Clean30Hzdata = importdata('M00630hzClean.txt',',',1);
Clean30Hz = Clean30Hzdata.data;
AOA_C30 =  Clean30Hz(:,1); % angle of attack
NF_C30 = (Clean30Hz(:,2)).*4.448; % Normal force(up) newtons
AF_C30 = (Clean30Hz(:,4)).*4.448; % Axial force(front to back) newtons
PM_C30 = (Clean30Hz(:,5))./8.851;% Pitching moment Newton meters

%Lift and drag
Lift_C30 = NF_C30.*cosd(AOA_C30) - AF_C30.*sind(AOA_C30);
Drag_C30 = NF_C30.*sind(AOA_C30) + AF_C30.*cosd(AOA_C30);

%Coefficients of lift, drag, and pithcing moment
Cd_C30 = Drag_C30./(0.5.*rho_30.*((velocity_30).^2).*A);
Cl_C30 = Lift_C30./(0.5.*rho_30.*((velocity_30).^2).*A);
Cm_C30 = PM_C30./(0.5.*rho_30.*((velocity_30).^2).*A);


Dirty30Hzdata = importdata('M00630hzDirty.txt',',',1);
Dirty30Hz = Dirty30Hzdata.data;
AOA_D30 = Dirty30Hz(:,1); % angle of attack
NF_D30 = (Dirty30Hz(:,2)).*4.448; % Normal force(up)
AF_D30 = (Dirty30Hz(:,4)).*4.448; % Axial force(front to back)
PM_D30 = (Dirty30Hz(:,5))./8.851;% Pitching moment

%Lift and drag
Lift_D30 = NF_D30.*cosd(AOA_D30) - AF_D30.*sind(AOA_D30);
Drag_D30 = NF_D30.*sind(AOA_D30) + AF_D30.*cosd(AOA_D30);
LiftoverDrag_C30 = Lift_C30./Drag_C30;
LiftoverDrag_D30 = Lift_D30./Drag_D30;

%Coefficients of lift, drag, and pithcing moment
Cd_D30 = Drag_D30./(0.5.*rho_30.*((velocity_30).^2).*A);
Cl_D30 = Lift_D30./(0.5.*rho_30.*((velocity_30).^2).*A);
Cm_D30 = PM_D30./(0.5.*rho_30.*((velocity_30).^2).*A);

% Plotting clean and dirty 30Hz 
%lift 
figure()
plot(AOA_C30,Lift_C30)
hold on
plot(AOA_D30,Lift_D30)
grid()
ylabel('Lift(N)')
xlabel('Angle of Attack (Degrees)')
title('Angle of Attack vs. Lift for 30Hz')
legend('Clean', 'Dirty',Location='northwest')

%drag
figure()
plot(AOA_C30,Drag_C30)
hold on
plot(AOA_D30,Drag_D30)
grid()
ylabel('Drag(N)')
xlabel('Angle of Attack (Degrees)')
title('Angle of Attack vs. Drag for 30Hz')
legend('Clean', 'Dirty',Location='northwest')

% pithcing moment 
figure()
plot(AOA_C30,PM_C30)
hold on
plot(AOA_D30,PM_D30)
grid()
ylabel('Pitching Moment(N)')
xlabel('Angle of Attack (Degrees)')
title('Angle of Attack vs. Pitching Moment for 30Hz')
legend('Clean', 'Dirty',Location='southwest')

% Coefficent of lift  
figure()
plot(AOA_C30,Cl_C30)
hold on
plot(AOA_D30,Cl_D30)
grid()
ylabel('Coefficent of lift')
xlabel('Angle of Attack (Degrees)')
title('Angle of Attack vs. Coefficent of Lift for 30Hz')
legend('Clean', 'Dirty',Location='northwest')

% Coefficent of Drag  
figure()
plot(AOA_C30,Cd_C30)
hold on
plot(AOA_D30,Cd_D30)
grid()
ylabel('Coefficent of Drag')
xlabel('Angle of Attack (Degrees)')
title('Angle of Attack vs. Coefficent of Drag for 30Hz')
legend('Clean', 'Dirty',Location='northwest')

% Coefficent of Moment  
figure()
plot(AOA_C30,Cm_C30)
hold on
plot(AOA_D30,Cm_D30)
grid()
ylabel('Coefficent of Moment')
xlabel('Angle of Attack (Degrees)')
title('Angle of Attack vs. Coefficent of Moment for 30Hz')
legend('Clean', 'Dirty',Location='southwest')

% Lift/Drag  
figure()
plot(AOA_C30,LiftoverDrag_C30)
hold on
plot(AOA_D30,LiftoverDrag_D30)
grid()
ylabel('Lift/Drag ')
xlabel('Angle of Attack (Degrees)')
title('Angle of Attack vs. Lift/Drag  for 30Hz')
legend('Clean', 'Dirty',Location='northwest')

%Tabulation 30Hz

%renoylds number 
Re_30 = (rho_30.*velocity_30.*chord)./mu;
PeakefficencyClean30 = AOA_C30(find(max(LiftoverDrag_C30) == LiftoverDrag_C30));
PeakefficencyDirty30 = AOA_D30(find(max(LiftoverDrag_D30) == LiftoverDrag_D30));

% Slope of the coefficient of lift vs. angle of attack 

[fitted3 ,s3] = polyfit(AOA_C30(1:7),Cl_C30(1:7),1); % first order regression over linear area before stall
liftcurveslope_30 = fitted3(1);



%45Hz
Clean45Hzdata = importdata('M00645hzClean.txt',',',1);
Clean45Hz = Clean45Hzdata.data;

AOA_C45 =  Clean45Hz(:,1); % angle of attack
NF_C45 = (Clean45Hz(:,2)).*4.448; % Normal force(up)
AF_C45 = (Clean45Hz(:,4)).*4.448; % Axial force(front to back)
PM_C45 = (Clean45Hz(:,5))./8.851;% Pitching moment

%Lift and drag
Lift_C45 = NF_C45.*cosd(AOA_C45) - AF_C45.*sind(AOA_C45);
Drag_C45 = NF_C45.*sind(AOA_C45) + AF_C45.*cosd(AOA_C45);

%Coefficients of lift, drag, and pithcing moment
Cd_C45 = Drag_C45./(0.5.*rho_45.*((velocity_45).^2).*A);
Cl_C45 = Lift_C45./(0.5.*rho_45.*((velocity_45).^2).*A);
Cm_C45 = PM_C45./(0.5.*rho_45.*((velocity_45).^2).*A);

Dirty45Hzdata = importdata('M00645hzDirty.txt',',',1);
Dirty45Hz = Dirty45Hzdata.data;

AOA_D45 = Dirty45Hz(:,1); % angle of attack
NF_D45 = (Dirty45Hz(:,2)).*4.448; % Normal force(up)
AF_D45 = (Dirty45Hz(:,4)).*4.448; % Axial force(front to back)
PM_D45 = (Dirty45Hz(:,5))./8.851;% Pitching moment

%Lift and drag
Lift_D45 = NF_D45.*cosd(AOA_D45) - AF_D45.*sind(AOA_D45);
Drag_D45 = NF_D45.*sind(AOA_D45) + AF_D45.*cosd(AOA_D45);
LiftoverDrag_C45 = Lift_C45./Drag_C45;
LiftoverDrag_D45 = Lift_D45./Drag_D45;

%Coefficients of lift, drag, and pithcing moment
Cd_D45 = Drag_D45./(0.5.*rho_45.*((velocity_45).^2).*A);
Cl_D45 = Lift_D45./(0.5.*rho_45.*((velocity_45).^2).*A);
Cm_D45 = PM_D45./(0.5.*rho_45.*((velocity_45).^2).*A);

% Plotting clean and dirty 45Hz 
%lift 
figure()
plot(AOA_C45,Lift_C45)
hold on
plot(AOA_D45,Lift_D45)
grid()
ylabel('Lift(N)')
xlabel('Angle of Attack (Degrees)')
title('Angle of Attack vs. Lift for 45Hz')
legend('Clean', 'Dirty',Location='northwest')

%drag
figure()
plot(AOA_C45,Drag_C45)
hold on
plot(AOA_D45,Drag_D45)
grid()
ylabel('Drag(N)')
xlabel('Angle of Attack (Degrees)')
title('Angle of Attack vs. Drag for 45Hz')
legend('Clean', 'Dirty',Location='northwest')

% pithcing moment 
figure()
plot(AOA_C45,PM_C45)
hold on
plot(AOA_D45,PM_D45)
grid()
ylabel('Pitching Moment(N)')
xlabel('Angle of Attack (Degrees)')
title('Angle of Attack vs. Pitching Moment for 45Hz')
legend('Clean', 'Dirty',Location='southwest')

% Coefficent of lift  
figure()
plot(AOA_C45,Cl_C45)
hold on
plot(AOA_D45,Cl_D45)
grid()
ylabel('Coefficent of Lift')
xlabel('Angle of Attack (Degrees)')
title('Angle of Attack vs. Coefficent of Lift for 45Hz')
legend('Clean', 'Dirty',Location='northwest')

% Coefficent of Drag  
figure()
plot(AOA_C45,Cd_C45)
hold on
plot(AOA_D45,Cd_D45)
grid()
ylabel('Coefficent of Drag')
xlabel('Angle of Attack (Degrees)')
title('Angle of Attack vs. Coefficent of Drag for 45Hz')
legend('Clean', 'Dirty',Location='northwest')

% Coefficent of Moment  
figure()
plot(AOA_C45,Cm_C45)
hold on
plot(AOA_D45,Cm_D45)
grid()
ylabel('Coefficent of Moment')
xlabel('Angle of Attack (Degrees)')
title('Angle of Attack vs. Coefficent of Moment for 45Hz')
legend('Clean', 'Dirty',Location='southwest')

% Lift/Drag  
figure()
plot(AOA_C45,LiftoverDrag_C45)
hold on
plot(AOA_D45,LiftoverDrag_D45)
grid()
ylabel('Lift/Drag ')
xlabel('Angle of Attack (Degrees)')
title('Angle of Attack vs. Lift/Drag  for 45Hz')
legend('Clean', 'Dirty',Location='northwest')

%Tabulation 30Hz

%renoylds number 
Re_45 = (rho_45.*velocity_45.*chord)./mu;
PeakefficencyClean45 = AOA_C45(find(max(LiftoverDrag_C45) == LiftoverDrag_C45));
PeakefficencyDirty45 = AOA_D45(find(max(LiftoverDrag_D45) == LiftoverDrag_D45));

% Slope of the coefficient of lift vs. angle of attack 

[fitted2 ,s2] = polyfit(AOA_C45(1:7),Cl_C45(1:7),1); % first order regression over linear area before stall
liftcurveslope_45 = fitted2(1);

%All 4 lift drag ratios on one plot

figure()
plot(AOA_C30,LiftoverDrag_C30)
hold on
plot(AOA_D30,LiftoverDrag_D30)
plot(AOA_C45,LiftoverDrag_C45)
plot(AOA_D45,LiftoverDrag_D45)
grid()
ylabel('Lift/Drag ')
xlabel('Angle of Attack (Degrees)')
title('Angle of Attack vs. Lift/Drag  for 30Hz and 45Hz')
legend('Clean30', 'Dirty30','Clean45', 'Dirty45',Location='southeast')

% printing values

disp('Airfoil 30Hz')
fprintf('Fan Frequency(30Hz): %gHz\n',x30)
fprintf('Reynold’s Number(30Hz): %g\n',Re_30)
fprintf('Peak Efficiency clean(30Hz): %g Degrees\n',PeakefficencyClean30)
fprintf('Peak Efficiency dirty(30Hz): %g Degrees\n',PeakefficencyDirty30)
fprintf('Slope of the coefficient of lift vs. angle of attack(30Hz): %gCl/Degrees\n',liftcurveslope_30)

disp('Airfoil 45Hz')
fprintf('Fan Frequency(45Hz): %gHz\n',x45)
fprintf('Reynold’s Number(45Hz): %g\n',Re_45)
fprintf('Peak Efficiency clean(45Hz): %g Degrees\n',PeakefficencyClean45)
fprintf('Peak Efficiency dirty(45Hz): %g Degrees\n',PeakefficencyDirty45)
fprintf('Slope of the coefficient of lift vs. angle of attack(45Hz): %g Cl/Degrees\n',liftcurveslope_45)
