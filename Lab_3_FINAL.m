% Group 1
% MAE 315
clear
close
clc



%% Shaker Tests

%% Shaker 16"
rawdata = importdata('M006_A_shaker_16.txt','\t');
O1s = rawdata(:,4);
F1s = rawdata(:,2); 
I1s = rawdata(:,3);
O1s = O1s(1:8);
I1s = I1s(1:8);
F1s = F1s(1:8);

% Uncertainty
i = 1;
deltf = zeros(1,8);
while i < 7
    deltf(i) = F1s(i+1) - F1s(i);
    i = i + 1;
end
uf = 0.5 .* deltf;

% Theoretical Response data
d = 0.284; %lb/in^3 density
L = 26; % in
w = 25.35 / 25.4; %mm to in
t = 6.55 / 25.4; %mm to in
l = 16; % in
V = L * w * t; % volume
m = (33/140) * d * V * (l/L); % calculated mass lbf
dampning_ratio = 0.04372053143; % from hammer test but calculated is 0.04896 
F0 = mean(I1s).* 8275.4 ; % lbm*s^2

E = 29000 * 1000; %young's modulus ksi to psi
I = (1/12) * w * (t^3); %moment of inertia
k = (3 * E * I) / (l^3); %spring constant
k = k * 32.2 *12; %converted from lbm to lbf
%wn = sqrt(k / m); % theoretical undamped natural frequency
wn = 29; %  expiremental data in hz
%wn = wn ./ (2*pi); % converted from rad/s to Hz
opw = [24:0.05:30]; % operating frequency 
R = (F0./m) ./ sqrt(((((wn*2*pi).^2)-((2*pi*opw).^2)).^2)+(2.*dampning_ratio.*(2*pi*wn).*(opw*2*pi)).^2);

figure(1)
plot(F1s,O1s)
hold on
plot(F1s,I1s)
errorbar(F1s(1:1:end),O1s(1:1:end),uf(1:1:end),'horizontal','.')
errorbar(F1s(1:1:end),I1s(1:1:end),uf(1:1:end),'horizontal','.')
hold on 
plot(opw,R)
axis([23 31 0 210]);
xlabel('Freqency (Hz)')
ylabel('Amplitude (V)')
title('Shaker 16"')
legend('Output','Input','Output Uncertainty','Input Uncertainty','Theoretical Line','Location','northwest')

% Normalized Frequency
NF1 = F1s ./ wn;
i = 1;
ndeltf = zeros(1,8);
while i < 7
    ndeltf(i) = NF1(i+1) - NF1(i);
    i = i + 1;
end
unf = 0.5 .* ndeltf;

figure(2)
plot(NF1,O1s)
hold on
plot(NF1,I1s)
errorbar(NF1(1:1:end),O1s(1:1:end),unf(1:1:end),'horizontal','.')
errorbar(NF1(1:1:end),I1s(1:1:end),unf(1:1:end),'horizontal','.')
hold on 
plot(opw./wn,R)
axis([0.8 1.05 0 210]);
xlabel('Normalized Freqency (Hz)')
ylabel('Amplitude (V)')
title('Normalized Frequency of Shaker 16"')
legend('Output','Input','Output Uncertainty','Input Uncertainty', 'Theoretical Line','Location','northwest')
%% Shaker 20"
rawdata = importdata('M006_A_shaker_20.txt','\t');
O1s = rawdata(:,4);
F1s = rawdata(:,2); 
I1s = rawdata(:,3);
O1s = O1s(1:7);
I1s = I1s(1:7);
F1s = F1s(1:7);

% Uncertainty
i = 1;
deltf = zeros(1,7);
while i < 6
    deltf(i) = F1s(i+1) - F1s(i);
    i = i + 1;
end
uf = 0.5 .* deltf;

% Theoretical Response data
d = 0.284; %lb/in^3 density
L = 26; % in
w = 25.35 / 25.4; %mm to in
t = 6.55 / 25.4; %mm to in
l = 20; % in
V = L * w * t; % volume
m = (33/140) * d * V * (l/L); % calculated mass lbf
dampning_ratio = 0.01961506709; % from hammer test but calculated is 0.04896 
F0 = mean(I1s).* 8275.4 ; % lbm*s^2

E = 29000 * 1000; %young's modulus ksi to psi
I = (1/12) * w * (t^3); %moment of inertia
k = (3 * E * I) / (l^3); %spring constant
k = k * 32.2 *12; %converted from lbm to lbf
%wn = sqrt(k / m); % theoretical undamped natural frequency
wn = 19; %  expiremental data in hz
%wn = wn ./ (2*pi); % converted from rad/s to Hz
opw = [16:0.05:20]; % operating frequency 
R = (F0./m) ./ sqrt(((((wn*2*pi).^2)-((2*pi*opw).^2)).^2)+(2.*dampning_ratio.*(2*pi*wn).*(opw*2*pi)).^2);

figure(3)
plot(F1s,O1s)
hold on
plot(F1s,I1s)
errorbar(F1s(1:1:end),O1s(1:1:end),uf(1:1:end),'horizontal','.')
errorbar(F1s(1:1:end),I1s(1:1:end),uf(1:1:end),'horizontal','.')
hold on 
plot(opw,R)
axis([15 21 0 475]);
xlabel('Freqency (Hz)')
ylabel('Amplitude (V)')
title('Shaker 20"')
legend('Output','Input','Output Uncertainty','Input Uncertainty','Theoretical Line','Location','northwest')

% Normalized Frequency
NF1 = F1s ./ wn;
i = 1;
ndeltf = zeros(1,7);
while i < 6
    ndeltf(i) = NF1(i+1) - NF1(i);
    i = i + 1;
end
unf = 0.5 .* ndeltf;

figure(4)
plot(NF1,O1s)
hold on
plot(NF1,I1s)
errorbar(NF1(1:1:end),O1s(1:1:end),unf(1:1:end),'horizontal','.')
errorbar(NF1(1:1:end),I1s(1:1:end),unf(1:1:end),'horizontal','.')
hold on 
plot(opw./wn,R)
xlabel('Normalized Freqency (Hz)')
ylabel('Amplitude (V)')
title('Normalized Frequency of Shaker 20"')
legend('Output','Input','Output Uncertainty','Input Uncertainty', 'Theoretical Line','Location','northwest')
%% Shaker 24"
rawdata = importdata('M006_A_shaker_24.txt','\t');
O1s = rawdata(:,4);
F1s = rawdata(:,2); 
I1s = rawdata(:,3);
O1s = O1s(1:7);
I1s = I1s(1:7);
F1s = F1s(1:7);

% Uncertainty
i = 1;
deltf = zeros(1,7);
while i < 6
    deltf(i) = F1s(i+1) - F1s(i);
    i = i + 1;
end
uf = 0.5 .* deltf;

% Theoretical Response data
d = 0.284; %lb/in^3 density
L = 26; % in
w = 25.35 / 25.4; %mm to in
t = 6.55 / 25.4; %mm to in
l = 24; % in
V = L * w * t; % volume
m = (33/140) * d * V * (l/L); % calculated mass lbf
dampning_ratio = 0.09149065213; % from hammer test but calculated is 0.04896 
F0 = mean(I1s).* 8275.4 ; % lbm*s^2

E = 29000 * 1000; %young's modulus ksi to psi
I = (1/12) * w * (t^3); %moment of inertia
k = (3 * E * I) / (l^3); %spring constant
k = k * 32.2 *12; %converted from lbm to lbf
%wn = sqrt(k / m); % theoretical undamped natural frequency
wn = 14; %  expiremental data in hz
%wn = wn ./ (2*pi); % converted from rad/s to Hz
opw = [12:0.05:16]; % operating frequency 
R = (F0./m) ./ sqrt(((((wn*2*pi).^2)-((2*pi*opw).^2)).^2)+(2.*dampning_ratio.*(2*pi*wn).*(opw*2*pi)).^2);

figure(5)
plot(F1s,O1s)
hold on
plot(F1s,I1s)
errorbar(F1s(1:1:end),O1s(1:1:end),uf(1:1:end),'horizontal','.')
errorbar(F1s(1:1:end),I1s(1:1:end),uf(1:1:end),'horizontal','.')
hold on 
plot(opw,R)
axis([11 17 0 150]);
xlabel('Freqency (Hz)')
ylabel('Amplitude (V)')
title('Shaker 24"')
legend('Output','Input','Output Uncertainty','Input Uncertainty','Theoretical Line','Location','northwest')

% Normalized Frequency
NF1 = F1s ./ wn;
i = 1;
ndeltf = zeros(1,7);
while i < 6
    ndeltf(i) = NF1(i+1) - NF1(i);
    i = i + 1;
end
unf = 0.5 .* ndeltf;

figure(6)
plot(NF1,O1s)
hold on
plot(NF1,I1s)
errorbar(NF1(1:1:end),O1s(1:1:end),unf(1:1:end),'horizontal','.')
errorbar(NF1(1:1:end),I1s(1:1:end),unf(1:1:end),'horizontal','.')
hold on 
plot(opw./wn,R)
xlabel('Normalized Freqency (Hz)')
ylabel('Amplitude (V)')
title('Normalized Frequency of Shaker 24"')
legend('Output','Input','Output Uncertainty','Input Uncertainty', 'Theoretical Line','Location','northwest')
%% Undamped vs Damped Natural Frequencies at Different Lengths
y1 = 44.5176;
y2 = 44.560;
x1 = 16;

y3 = 13.249;
y4 = 13.2515;
x2 = 20;

y5 = 17.64;
y6 = 17.7143;
x3 = 24;

ux = (1 / 64);
uy = 0.183;

figure(7)
plot(x1,y1,'.',"Color",'red','MarkerSize',20)
hold on
plot(x1,y2,'.',"Color",'blue','MarkerSize',20)
hold on
errorbar(x1(1:1:end),y1(1:1:end),ux(1:1:end),'horizontal','.')
errorbar(x1(1:1:end),y1(1:1:end),uy(1:1:end),'vertical','.')
plot(x2,y3,'.',"Color",'red','MarkerSize',20)
hold on
plot(x2,y4,'.',"Color",'blue','MarkerSize',20)
hold on
plot(x3,y5,'.',"Color",'red','MarkerSize',20)
hold on
plot(x3,y6,'.',"Color",'blue',"MarkerSize",20)
errorbar(x1(1:1:end),y2(1:1:end),ux(1:1:end),'horizontal','.')
errorbar(x1(1:1:end),y2(1:1:end),uy(1:1:end),'vertical','.')
errorbar(x2(1:1:end),y3(1:1:end),ux(1:1:end),'horizontal','.')
errorbar(x2(1:1:end),y3(1:1:end),uy(1:1:end),'vertical','.')
errorbar(x2(1:1:end),y4(1:1:end),ux(1:1:end),'horizontal','.')
errorbar(x2(1:1:end),y4(1:1:end),uy(1:1:end),'vertical','.')
errorbar(x3(1:1:end),y5(1:1:end),ux(1:1:end),'horizontal','.')
errorbar(x3(1:1:end),y5(1:1:end),uy(1:1:end),'vertical','.')
errorbar(x3(1:1:end),y6(1:1:end),ux(1:1:end),'horizontal','.')
errorbar(x3(1:1:end),y6(1:1:end),uy(1:1:end),'vertical','.')
title('Wd vs Wn')
xlabel('Length (in)')
ylabel('Frequency (Hz)')
axis([14 28 10 50])
legend('Wd','Wn','Length Error','Frequency Error','Location','northeast')
%% Hammer Test
%% Hammer 16"

rawdata = importdata('10_28_2024hammer_16_m002adat.txt','\t',1);
info = rawdata.data;
hammer_freq = info(:,1);
hammer_Amp_out = info(:,3);
hammer_Amp_in = info(:,2);
rawdata = importdata('10_28_2024hammer_16_m002afft.txt','\t',1);
info = rawdata.data;
hammer_in = info(:,2);
Time1 = info(:,1); 
hammer_out = info(:,3);

tf_hammer_out = fft(hammer_out);
tf_hammer_in = fft(hammer_in);
Hf = tf_hammer_out./tf_hammer_in;
Mag_response = sqrt(imag(Hf).^2+real(Hf).^2);
Mag_response2 = abs(Hf);
phase = atan(imag(Hf)./real(Hf)).*180./pi;

bipolar = 5;
bitres = 24;
u_V = bipolar ./ (2.^(bitres)-1);

i = 1;
deltf = zeros(1,length(hammer_freq));
while i < length(hammer_freq)
    deltf(i) = hammer_freq(i+1) - hammer_freq(i);
    i = i + 1;
end
uf = 0.5 .* deltf;

i = 1;
deltT = zeros(1,length(Time1));
while i < length(Time1)
    deltT(i) = Time1(i+1) - Time1(i);
    i = i + 1;
end
uT = 0.5 .* deltT;

figure(8)
subplot(2,2,1)
plot(Time1,hammer_Amp_out)
hold on
errorbar(Time1(1:500:end),hammer_Amp_out(1:500:end),u_V,'vertical','.')
errorbar(Time1(1:500:end),hammer_Amp_out(1:500:end),uT(1:500:end),'horizontal','.')
xlabel('Time(s)')
ylabel('Voltage(Vx10^2)')
title('Time vs Voltage Output 16"')
hold off

subplot(2,2,2)
plot(Time1,hammer_Amp_in)
hold on
errorbar(Time1(1:500:end),hammer_Amp_in(1:500:end),u_V,'vertical','.')
errorbar(Time1(1:500:end),hammer_Amp_in(1:500:end),uT(1:500:end),'horizontal','.')
xlabel('Time(s)')
ylabel('Voltage(Vx10^2)')
title('Time vs Voltage Input 16"')
hold off

subplot(2,2,3)
plot(hammer_freq,hammer_out)
hold on
errorbar(hammer_freq(1:400:end),hammer_out(1:400:end),u_V,'vertical','.')
errorbar(hammer_freq(1:400:end),hammer_out(1:400:end),uf(1:400:end),'horizontal','.')
xlabel('Freqency(Hz)')
ylabel('Amplitude(Vx10^2)')
title('Freqency vs Amplitude Output 16"')
axis([0 1500 -75 100]);
hold off

subplot(2,2,4)
plot(hammer_freq,hammer_in)
hold on
errorbar(hammer_freq(1:400:end),hammer_in(1:400:end),u_V,'vertical','.')
errorbar(hammer_freq(1:400:end),hammer_in(1:400:end),uf(1:400:end),'horizontal','.')
xlabel('Freqency(Hz)')
ylabel('Amplitude(Vx10^2)')
title('Freqency vs Amplitude Output 16"')
axis([0 1500 -15 50]);
hold off

figure(9)
subplot (2,1,1)
plot(hammer_freq,Hf)
axis([0 1500 0 200]);
xlabel('Freqency(Hz)')
ylabel('Freqency Response')
title('Freqency vs Freqency Response 16"')

subplot (2,1,2)
plot(hammer_freq,phase)
axis([0 1500 -100 100]);
xlabel('Freqency(Hz)')
ylabel('Phase Lag (Degrees)')
title('Freqency vs Phase Lag 16"')

%% Hammer 20"

rawdata = importdata('10_28_2024hammer_20_m002adat.txt','\t',1);
info = rawdata.data;
hammer_freq_20 = info(:,1);
hammer_Amp_out_20 = info(:,3);
hammer_Amp_in_20 = info(:,2);
rawdata = importdata('10_28_2024hammer_20_m002afft.txt','\t',1);
info = rawdata.data;
hammer_in_20 = info(:,2);
Time1_20 = info(:,1); 
hammer_out_20 = info(:,3);

tf_hammer_out_20 = fft(hammer_out_20);
tf_hammer_in_20 = fft(hammer_in_20);
Hf_20 = tf_hammer_out_20./tf_hammer_in_20;
Mag_response_20 = sqrt(imag(Hf_20).^2+real(Hf_20).^2);
Mag_response2_20 = abs(Hf_20);
phase_20 = atan(imag(Hf_20)./real(Hf_20)).*180./pi;

bipolar = 5;
bitres = 24;
u_V = bipolar ./ (2.^(bitres)-1);

i = 1;
deltf_20 = zeros(1,length(hammer_freq_20));
while i < length(hammer_freq_20)
    deltf_20(i) = hammer_freq_20(i+1) - hammer_freq_20(i);
    i = i + 1;
end
uf_20 = 0.5 .* deltf_20;

i = 1;
deltT_20 = zeros(1,length(Time1_20));
while i < length(Time1_20)
    deltT_20(i) = Time1_20(i+1) - Time1_20(i);
    i = i + 1;
end
uT_20 = 0.5 .* deltT_20;

figure(10)
subplot(2,2,1)
plot(Time1_20,hammer_Amp_out_20)
hold on
errorbar(Time1_20(1:500:end),hammer_Amp_out_20(1:500:end),u_V,'vertical','.')
errorbar(Time1_20(1:500:end),hammer_Amp_out_20(1:500:end),uT_20(1:500:end),'horizontal','.')
xlabel('Time(s)')
ylabel('Voltage(Vx10^2)')
title('Time vs Voltage Output 20"')
axis([-.1 3 0 .6]);
hold off

subplot(2,2,2)
plot(Time1_20,hammer_Amp_in_20)
hold on
errorbar(Time1_20(1:500:end),hammer_Amp_in_20(1:500:end),u_V,'vertical','.')
errorbar(Time1_20(1:500:end),hammer_Amp_in_20(1:500:end),uT_20(1:500:end),'horizontal','.')
xlabel('Time(s)')
ylabel('Voltage(Vx10^2)')
title('Time vs Voltage Input 20"')
axis([0 3 0 .13]);
hold off

subplot(2,2,3)
plot(hammer_freq_20,hammer_out_20)
hold on
errorbar(hammer_freq_20(1:400:end),hammer_out_20(1:400:end),u_V,'vertical','.')
errorbar(hammer_freq_20(1:400:end),hammer_out_20(1:400:end),uf_20(1:400:end),'horizontal','.')
xlabel('Freqency(Hz)')
ylabel('Amplitude(Vx10^2)')
title('Freqency vs Amplitude Output 20"')
axis([0 1500 -25 35]);
hold off

subplot(2,2,4)
plot(hammer_freq_20,hammer_in_20)
hold on
errorbar(hammer_freq_20(1:400:end),hammer_in_20(1:400:end),u_V,'vertical','.')
errorbar(hammer_freq_20(1:400:end),hammer_in_20(1:400:end),uf_20(1:400:end),'horizontal','.')
xlabel('Freqency(Hz)')
ylabel('Amplitude(Vx10^2)')
title('Freqency vs Amplitude Output 20"')
axis([0 1500 -2 15]);
hold off

figure(11)
subplot (2,1,1)
plot(hammer_freq_20,Hf_20)
axis([0 1500 0 170]);
xlabel('Freqency(Hz)')
ylabel('Freqency Response')
title('Freqency vs Freqency Response 20"')

subplot (2,1,2)
plot(hammer_freq_20,phase_20)
axis([0 1500 -100 100]);
xlabel('Freqency(Hz)')
ylabel('Phase Lag (Degrees)')
title('Freqency vs Phase Lag 20"')


%% Hammer 24"

rawdata = importdata('10_28_2024hammer_24_m002adat.txt','\t',1);
info = rawdata.data;
hammer_freq_24 = info(:,1);
hammer_Amp_out_24 = info(:,3);
hammer_Amp_in_24 = info(:,2);
rawdata = importdata('10_28_2024hammer_24_m002afft.txt','\t',1);
info = rawdata.data;
hammer_in_24 = info(:,2);
Time1_24 = info(:,1); 
hammer_out_24 = info(:,3);

tf_hammer_out_24 = fft(hammer_out_24);
tf_hammer_in_24 = fft(hammer_in_24);
Hf_24 = tf_hammer_out_24./tf_hammer_in_24;
Mag_response_24 = sqrt(imag(Hf_24).^2+real(Hf_24).^2);
Mag_response2_24 = abs(Hf_24);
phase_24 = atan(imag(Hf_24)./real(Hf_24)).*180./pi;

bipolar = 5;
bitres = 24;
u_V = bipolar ./ (2.^(bitres)-1);

i = 1;
deltf_24 = zeros(1,length(hammer_freq_24));
while i < length(hammer_freq_24)
    deltf_24(i) = hammer_freq_24(i+1) - hammer_freq_24(i);
    i = i + 1;
end
uf_24 = 0.5 .* deltf_24;

i = 1;
deltT_24 = zeros(1,length(Time1_24));
while i < length(Time1_24)
    deltT_24(i) = Time1_24(i+1) - Time1_24(i);
    i = i + 1;
end
uT_24 = 0.5 .* deltT_24;

figure(12)
subplot(2,2,1)
plot(Time1_24,hammer_Amp_out_24)
hold on
errorbar(Time1_24(1:500:end),hammer_Amp_out_24(1:500:end),u_V,'vertical','.')
errorbar(Time1_24(1:500:end),hammer_Amp_out_24(1:500:end),uT_24(1:500:end),'horizontal','.')
xlabel('Time(s)')
ylabel('Voltage(Vx10^2)')
title('Time vs Voltage Output 24"')
axis([-.1 3 0 .55]);
hold off

subplot(2,2,2)
plot(Time1_24,hammer_Amp_in_24)
hold on
errorbar(Time1_24(1:500:end),hammer_Amp_in_24(1:500:end),u_V,'vertical','.')
errorbar(Time1_24(1:500:end),hammer_Amp_in_24(1:500:end),uT_24(1:500:end),'horizontal','.')
xlabel('Time(s)')
ylabel('Voltage(Vx10^2)')
title('Time vs Voltage Input 24"')
axis([0 3 0 .12]);
hold off

subplot(2,2,3)
plot(hammer_freq_24,hammer_out_24)
hold on
errorbar(hammer_freq_24(1:400:end),hammer_out_24(1:400:end),u_V,'vertical','.')
errorbar(hammer_freq_24(1:400:end),hammer_out_24(1:400:end),uf_24(1:400:end),'horizontal','.')
xlabel('Freqency(Hz)')
ylabel('Amplitude(Vx10^2)')
title('Freqency vs Amplitude Output 24"')
axis([0 1500 -35 50]);
hold off

subplot(2,2,4)
plot(hammer_freq_24,hammer_in_24)
hold on
errorbar(hammer_freq_24(1:400:end),hammer_in_24(1:400:end),u_V,'vertical','.')
errorbar(hammer_freq_24(1:400:end),hammer_in_24(1:400:end),uf_24(1:400:end),'horizontal','.')
xlabel('Freqency(Hz)')
ylabel('Amplitude(Vx10^2)')
title('Freqency vs Amplitude Output 24"')
axis([0 1500 -6 19]);
hold off

figure(13)
subplot (2,1,1)
plot(hammer_freq_24,Hf_24)
axis([0 1500 0 65]);
xlabel('Freqency(Hz)')
ylabel('Freqency Response')
title('Freqency vs Freqency Response 24"')

subplot (2,1,2)
plot(hammer_freq_24,phase_24)
axis([0 1500 -100 100]);
xlabel('Freqency(Hz)')
ylabel('Phase Lag (Degrees)')
title('Freqency vs Phase Lag 24"')

%% Simulation

%% Experiment Validation
rawdata = importdata('10_28_2024caraceni3.1_16_1dat.txt','\t',1);
info = rawdata.data;
Fev16 = info(:,1); 
Vev16 = info(:,2);
rawdata = importdata('10_28_2024caraceni3.1_16_1fft.txt','\t',1);
info = rawdata.data;
Fev16_1 = info(:,1); 
Aev16 = info(:,2);
rawdata = importdata('10_28_2024caraceni3.1_20_1dat.txt','\t',1);
info = rawdata.data;
Fev20 = info(:,1); 
Vev20 = info(:,2);
rawdata = importdata('10_28_2024caraceni3.1_20_1fft.txt','\t',1);
info = rawdata.data;
Fev20_1 = info(:,1); 
Aev20 = info(:,2);
rawdata = importdata('10_28_2024caraceni3.1_24_1dat.txt','\t',1);
info = rawdata.data;
Fev24 = info(:,1); 
Vev24 = info(:,2);
rawdata = importdata('10_28_2024caraceni3.1_24_1fft.txt','\t',1);
info = rawdata.data;
Fev24_1 = info(:,1); 
Aev24 = info(:,2);
figure


subplot(2,2,[1,2])
plot(Fev16_1,Aev16)
hold on
plot(Fev20_1,Aev20)
hold on
plot(Fev24_1,Aev24)
xlabel('Frequency (Hz)')
ylabel('H(f)')
title('Magnitude of Response - Effects of Length')
legend('16 in.','20 in.','24 in.')

subplot(2,2,[3,4])

plot(Fev16,Vev16)
hold on
plot(Fev20,Vev20)
hold on
plot(Fev24,Vev24)
xlabel('Frequency (Hz)')
ylabel('Phase Lag (deg)')
title('Magnitude of Response - Effects of Length')
legend('16 in.','20 in.','24 in.')

%% Effects of Damping
rawdata = importdata('10_28_2024caraceni3.2_20_d1dat.txt','\t',1);
info = rawdata.data;
Fed20_d1 = info(:,1); 
Ved20_d1 = info(:,2);
rawdata = importdata('10_28_2024caraceni3.2_20_d1fft.txt','\t',1);
info = rawdata.data;
Fed20_d1_1 = info(:,1); 
Aed20_d1 = info(:,2);

rawdata = importdata('10_28_2024caraceni3.2_20_d3dat.txt','\t',1);
info = rawdata.data;
Fed20_d3 = info(:,1); 
Ved20_d3 = info(:,2);
rawdata = importdata('10_28_2024caraceni3.2_20_d3fft.txt','\t',1);
info = rawdata.data;
Fed20_d3_1 = info(:,1); 
Aed20_d3 = info(:,2);

rawdata = importdata('10_28_2024caraceni3.2_20_d5dat.txt','\t',1);
info = rawdata.data;
Fed20_d5 = info(:,1); 
Ved20_d5 = info(:,2);
rawdata = importdata('10_28_2024caraceni3.2_20_d5fft.txt','\t',1);
info = rawdata.data;
Fed20_d5_1 = info(:,1); 
Aed20_d5 = info(:,2);

rawdata = importdata('10_28_2024caraceni3.2_20_d7dat.txt','\t',1);
info = rawdata.data;
Fed20_d7 = info(:,1); 
Ved20_d7 = info(:,2);
rawdata = importdata('10_28_2024caraceni3.2_20_d7fft.txt','\t',1);
info = rawdata.data;
Fed20_d7_1 = info(:,1); 
Aed20_d7 = info(:,2);

rawdata = importdata('10_28_2024caraceni3.2_20_d9dat.txt','\t',1);
info = rawdata.data;
Fed20_d9 = info(:,1); 
Ved20_d9 = info(:,2);
rawdata = importdata('10_28_2024caraceni3.2_20_d9fft.txt','\t',1);
info = rawdata.data;
Fed20_d9_1 = info(:,1); 
Aed20_d9 = info(:,2);

figure

subplot(2,2,[1,2])
plot(Fed20_d1_1,Aed20_d1)
hold on
plot(Fed20_d3_1,Aed20_d3)
hold on
plot(Fed20_d5_1,Aed20_d5)
hold on
plot(Fed20_d7_1,Aed20_d7)
hold on
plot(Fed20_d9_1,Aed20_d9)
xlabel('Frequency (Hz)')
ylabel('H(f)')
title('Magnitude of Response - Effects of Damping')
legend('1 kg/s','3 kg/s','5 kg/s','7 kg/s','9 kg/s')

subplot(2,2,[3,4])

plot(Fed20_d1,Ved20_d1)
hold on
plot(Fed20_d3,Ved20_d3)
hold on
plot(Fed20_d5,Ved20_d5)
hold on
plot(Fed20_d7,Ved20_d7)
hold on
plot(Fed20_d9,Ved20_d9)

xlabel('Frequency (Hz)')
ylabel('Phase Lag (deg)')
title('Magnitude of Response - Effects of Damping')
legend('1 kg/s','3 kg/s','5 kg/s','7 kg/s','9 kg/s')
%% Effects of End-Mass

rawdata = importdata('10_28_2024caraceni3.3_20_w0.4dat.txt','\t',1);
info = rawdata.data;
Fem4 = info(:,1); 
Vem4 = info(:,2);
rawdata = importdata('10_28_2024caraceni3.3_20_w0.4fft.txt','\t',1);
info = rawdata.data;
Fem4_1 = info(:,1); 
Aem4 = info(:,2);

rawdata = importdata('10_28_2024caraceni3.3_20_w0.25dat.txt','\t',1);
info = rawdata.data;
Fem25 = info(:,1); 
Vem25 = info(:,2);
rawdata = importdata('10_28_2024caraceni3.3_20_w0.25fft.txt','\t',1);
info = rawdata.data;
Fem25_1 = info(:,1); 
Aem25 = info(:,2);

rawdata = importdata('10_28_2024caraceni3.3_20_w0.65dat.txt','\t',1);
info = rawdata.data;
Fem65 = info(:,1); 
Vem65 = info(:,2);
rawdata = importdata('10_28_2024caraceni3.3_20_w0.65fft.txt','\t',1);
info = rawdata.data;
Fem65_1 = info(:,1); 
Aem65 = info(:,2);

figure

subplot(2,2,[1,2])
plot(Fem25_1,Aem25)
hold on
plot(Fem4_1,Aem4)
hold on
plot(Fem65_1,Aem65)
xlabel('Frequency (Hz)')
ylabel('H(f)')
title('Magnitude of Response - Effects of End Mass')
legend('.25 kg','.4 kg','.65 kg')
subplot(2,2,[3,4])

plot(Fem25,Vem25)
hold on
plot(Fem4,Vem4)
hold on
plot(Fem65,Vem65)
xlabel('Frequency (Hz)')
ylabel('Phase Lag (deg)')
title('Magnitude of Response - Effects of End Mass')
legend('.25 kg','.4 kg','.65 kg')
%% Effects of Material Type

rawdata = importdata('10_28_2024caraceni3.4_20_aldat.txt','\t',1);
info = rawdata.data;
Fet_al = info(:,1); 
Vet_al = info(:,2);
rawdata = importdata('10_28_2024caraceni3.4_20_alfft.txt','\t',1);
info = rawdata.data;
Fet_al_1 = info(:,1); 
Aet_al = info(:,2);

rawdata = importdata('10_28_2024caraceni3.4_20_stainlessdat.txt','\t',1);
info = rawdata.data;
Fet_st = info(:,1); 
Vet_st = info(:,2);
rawdata = importdata('10_28_2024caraceni3.4_20_stainlessfft.txt','\t',1);
info = rawdata.data;
Fet_st_1 = info(:,1); 
Aet_st = info(:,2);

rawdata = importdata('10_28_2024caraceni3.4_20_baselinedat.txt','\t',1);
info = rawdata.data;
Fet_base = info(:,1); 
Vet_base = info(:,2);
rawdata = importdata('10_28_2024caraceni3.4_20_baselinefft.txt','\t',1);
info = rawdata.data;
Fet_base_1 = info(:,1); 
Aet_base = info(:,2);

figure

subplot(2,2,[1,2])
plot(Fet_al_1,Aet_al)
hold on
plot(Fet_st_1,Aet_st)
hold on
plot(Fet_base_1,Aet_base)
xlabel('Frequency (Hz)')
ylabel('H(f)')
title('Magnitude of Response - Effects of Material')
legend('Aluminum','Stainless Steel','Carbon Fiber')
subplot(2,2,[3,4])

plot(Fet_al,Vet_al)
hold on
plot(Fet_st,Vet_st)
hold on
plot(Fet_base,Vet_base)
xlabel('Frequency (Hz)')
ylabel('Phase Lag (deg)')
title('Magnitude of Response - Effects of Material')
legend('Aluminum','Stainless Steel','Carbon Fiber')