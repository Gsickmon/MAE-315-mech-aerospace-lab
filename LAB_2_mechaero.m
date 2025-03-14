% lab 2
%Garrett Sickmon
clc
clear
clear all
close all
%% Plotting graph 1 
clear
%Data input
rawdata = importdata('750Hz_3V_sinewave_0phase_12bit_10BR_25kHzsample_8192samples.txt',' ',9);
info = rawdata.data;
t1 = info(:,1); % time
V1 = info(:,2); % voltage
F1 = info(:,3); % Frequency (Hz)
A1 = info(:,4); % Amplitude

t1 = t1(1:100);
V1 = V1(1:100);
F1 = F1(1:1000);
A1 = A1(1:1000);
% time vs voltage graph
figure(1)
subplot(2,2,1)
plot(t1,V1)
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Simulation data - time')

% Freq graph
figure(1)
subplot(2,2,3)
plot(F1,A1)
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('Simulation data - freq')
% Nyquist graph

F1_sampling = 25; % sample rate is 25000 hz, 25khz
F1_folding = F1_sampling ./ 2;

x = [0 F1_folding 0 F1_folding 0 F1_folding 0 F1_folding]; %100Hz is the Folding frequency in this case
y = [0 0 1 1 2 2 3 3];
figure(1);
subplot(2,2,[2,4])
plot(x,y);
xlabel('Frequency (kHz)')
text(13,0.0,'12.5 kHz') % (point a little past the folding frequnacy, little above line, what the folding frequancy is)
text(-2,0.0,'0kHz')
text(-2,1.0,'25kHz')
text(13,1.0,'37.5 kHz')
text(-2,2.0,'50kHz')
text(13,2.0,'62.5 kHz')
text(-2,3.0,'75kHz')

axis([-3 16 -1 3.5]); % expands the graph
hold on
plot(.750,0,'m*','MarkerSize',18); % shows what the sampling frequancy is
set(text(.750,0.12,'750 Hz'),'Rotation',45)
title('Nyquist Diagram')
hold off

%% Plotting graph 2
clear
%Data input
rawdata = importdata('750Hz_3V_sinewave_0phase_4bit_10BR_25kHzsample_8192samples.txt',' ',9);
info = rawdata.data;
t2 = info(:,1); % time
V2 = info(:,2); % voltage
F2 = info(:,3); % Frequency (Hz)
A2 = info(:,4); % Amplitude

t2 = t2(1:100);
V2 = V2(1:100);
F2 = F2(1:1000);
A2 = A2(1:1000);
% time vs voltage graph
figure(2)
subplot(2,2,1)
plot(t2,V2)
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Simulation data - time')

% Freq graph
figure(2)
subplot(2,2,3)
plot(F2,A2)
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('Simulation data - freq')

% Nyquist graph

F2_sampling = 25; % sample rate is 25000 hz, 25Khz
F2_folding = F2_sampling ./ 2;

x = [0 F2_folding 0 F2_folding 0 F2_folding 0 F2_folding]; %100Hz is the Folding frequency in this case
y = [0 0 1 1 2 2 3 3];
figure(2);
subplot(2,2,[2,4])
plot(x,y);
xlabel('Frequency (kHz)')
text(13,0.0,'12.5 kHz') % (point a little past the folding frequnacy, little above line, what the folding frequancy is)
text(-2,0.0,'0kHz')
text(-2,1.0,'25kHz')
text(13,1.0,'37.5 kHz')
text(-2,2.0,'50kHz')
text(13,2.0,'62.5 kHz')
text(-2,3.0,'75kHz')

axis([-3 16 -1 3.5]); % expands the graph
hold on
plot(.750,0,'m*','MarkerSize',18); % shows what the sampling frequancy is
set(text(.750,0.12,'750 Hz'),'Rotation',45)
title('Nyquist Diagram')
hold off

%% Plotting graph 3
clear
%Data input
rawdata = importdata('750Hz_3V_sinewave_0phase_12bit_1BR_25kHzsample_8192samples.txt',' ',9);
info = rawdata.data;
t3 = info(:,1); % time
V3 = info(:,2); % voltage
F3 = info(:,3); % Frequency (Hz)
A3 = info(:,4); % Amplitude

t3 = t3(1:100);
V3 = V3(1:100);
F3 = F3(1:1000);
A3 = A3(1:1000);
% time vs voltage graph
figure(3)
subplot(2,2,1)
plot(t3,V3)
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Simulation data - time')

% Freq graph
figure(3)
subplot(2,2,3)
plot(F3,A3)
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('Simulation data - freq')

% Nyquist graph

F3_sampling = 25; % sample rate is 25000 hz, 25Khz
F3_folding = F3_sampling ./ 2;

x = [0 F3_folding 0 F3_folding 0 F3_folding 0 F3_folding]; %100Hz is the Folding frequency in this case
y = [0 0 1 1 2 2 3 3];
figure(3);
subplot(2,2,[2,4])
plot(x,y);
xlabel('Frequency (kHz)')
text(13,0.0,'12.5 kHz') % (point a little past the folding frequnacy, little above line, what the folding frequancy is)
text(-2,0.0,'0kHz')
text(-2,1.0,'25kHz')
text(13,1.0,'37.5 kHz')
text(-2,2.0,'50kHz')
text(13,2.0,'62.5 kHz')
text(-2,3.0,'75kHz')

axis([-3 16 -1 3.5]); % expands the graph
hold on
plot(.750,0,'m*','MarkerSize',18); % shows what the sampling frequancy is
set(text(.750,0.12,'750 Hz'),'Rotation',45)
title('Nyquist Diagram')
hold off

%% Plotting graph 4 (3.1.9)
clear
%Data input
rawdata1 = importdata('10_7_2024forschino1dat.txt','\t',1);
rawdata2 = importdata('10_7_2024forschino1fft.txt','\t',1);
info1 = rawdata1.data;
t3 = info1(:,1); % time
V3 = info1(:,2); % voltage
info2 = rawdata2.data;
F3 = info2(:,1); % Frequency (Hz)
A3 = info2(:,2); % Amplitude

t3 = t3(1:100);
V3 = V3(1:100);
F3 = F3(1:1000);
A3 = A3(1:1000);
% time vs voltage graph
figure(4)
subplot(1,2,1)
plot(t3,V3)
axis([-0.0001 0.0041 -3.1 3.1]); % expands the graph
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Simulation data - time')

% Freq graph
figure(4)
subplot(1,2,2)
plot(F3,A3)
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('Simulation data - freq')

%% Plotting graph 5 (3.1.12)
clear
%Data input
rawdata1 = importdata('10_7_2024forschino2dat.txt','\t',1);
rawdata2 = importdata('10_7_2024forschino2fft.txt','\t',1);
info1 = rawdata1.data;
t3 = info1(:,1); % time
V3 = info1(:,2); % voltage
info2 = rawdata2.data;
F3 = info2(:,1); % Frequency (Hz)
A3 = info2(:,2); % Amplitude

t3 = t3(1:100);
V3 = V3(1:100);
F3 = F3(1:1000);
A3 = A3(1:1000);
% time vs voltage graph
figure(5)
subplot(1,2,1)
plot(t3,V3)
axis([-0.0001 0.0041 -3.6 3.6]); % expands the graph
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Simulation data - time')

% Freq graph
figure(5)
subplot(1,2,2)
plot(F3,A3)
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('Simulation data - freq')

%% Plotting graph 6 (3.1.14)
clear
%Data input
rawdata1 = importdata('10_7_2024forschino3dat.txt','\t',1);
rawdata2 = importdata('10_7_2024forschino3fft.txt','\t',1);
info1 = rawdata1.data;
t3 = info1(:,1); % time
V3 = info1(:,2); % voltage
info2 = rawdata2.data;
F3 = info2(:,1); % Frequency (Hz)
A3 = info2(:,2); % Amplitude

t3 = t3(1:100);
V3 = V3(1:100);
F3 = F3(1:1000);
A3 = A3(1:1000);
% time vs voltage graph
figure(6)
subplot(1,2,1)
plot(t3,V3)
axis([-0.0001 0.0041 -3.6 3.6]); % expands the graph
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Simulation data - time')

% Freq graph
figure(6)
subplot(1,2,2)
plot(F3,A3)
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('Simulation data - freq')

%% Plotting graph 7 (3.1.15)
clear
%Data input
rawdata1 = importdata('10_7_2024forschino4dat.txt','\t',1);
rawdata2 = importdata('10_7_2024forschino4fft.txt','\t',1);
info1 = rawdata1.data;
t3 = info1(:,1); % time
V3 = info1(:,2); % voltage
info2 = rawdata2.data;
F3 = info2(:,1); % Frequency (Hz)
A3 = info2(:,2); % Amplitude

t3 = t3(1:100);
V3 = V3(1:100);
F3 = F3(1:1000);
A3 = A3(1:1000);
% time vs voltage graph
figure(7)
subplot(1,2,1)
plot(t3,V3)
axis([-0.0001 0.0041 -3.6 3.6]); % expands the graph
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Simulation data - time')

% Freq graph
figure(7)
subplot(1,2,2)
plot(F3,A3)
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('Simulation data - freq')

%% Plotting graph 8 (3.1.15)
clear
%Data input
rawdata1 = importdata('10_7_2024forschino5dat.txt','\t',1);
rawdata2 = importdata('10_7_2024forschino5fft.txt','\t',1);
info1 = rawdata1.data;
t3 = info1(:,1); % time
V3 = info1(:,2); % voltage
info2 = rawdata2.data;
F3 = info2(:,1); % Frequency (Hz)
A3 = info2(:,2); % Amplitude

t3 = t3(1:100);
V3 = V3(1:100);
F3 = F3(1:1000);
A3 = A3(1:1000);
% time vs voltage graph
figure(8)
subplot(1,2,1)
plot(t3,V3)
axis([-0.0001 0.0041 -3.6 3.6]); % expands the graph
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Simulation data - time')

% Freq graph
figure(8)
subplot(1,2,2)
plot(F3,A3)
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('Simulation data - freq')


%% Plotting graph 9 (3.1.15)
clear
%Data input
rawdata1 = importdata('10_7_2024forschino6dat.txt','\t',1);
rawdata2 = importdata('10_7_2024forschino6fft.txt','\t',1);
info1 = rawdata1.data;
t3 = info1(:,1); % time
V3 = info1(:,2); % voltage
info2 = rawdata2.data;
F3 = info2(:,1); % Frequency (Hz)
A3 = info2(:,2); % Amplitude

t3 = t3(1:100);
V3 = V3(1:100);
F3 = F3(1:1000);
A3 = A3(1:1000);
% time vs voltage graph
figure(9)
subplot(1,2,1)
plot(t3,V3)
axis([-0.0001 0.0041 -3.6 3.6]); % expands the graph
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Simulation data - time')

% Freq graph
figure(9)
subplot(1,2,2)
plot(F3,A3)
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('Simulation data - freq')


%% Plotting graph 10 (4.1.10)
clear
%Data input
rawdata1 = importdata('10_7_2024forschino21dat.txt','\t',1);
rawdata2 = importdata('10_7_2024forschino21fft.txt','\t',1);
info1 = rawdata1.data;
t3 = info1(:,1); % time
V3 = info1(:,2); % voltage
info2 = rawdata2.data;
F3 = info2(:,1); % Frequency (Hz)
A3 = info2(:,2); % Amplitude

t3 = t3(1:100);
V3 = V3(1:100);
F3 = F3(1:500);
A3 = A3(1:500);
% time vs voltage graph
figure(10)
subplot(2,2,1)
plot(t3,V3)
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Simulation data - time')

% Freq graph
figure(10)
subplot(2,2,3)
plot(F3,A3)
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('Simulation data - freq')

% Nyquist graph

F10_sampling = 25; % sample rate is 25000 hz, 25Khz
F10_folding = F10_sampling ./ 2;

x = [0 F10_folding 0 F10_folding 0 F10_folding 0 F10_folding]; %100Hz is the Folding frequency in this case
y = [0 0 1 1 2 2 3 3];
figure(10);
subplot(2,2,[2,4])
plot(x,y);
xlabel('Frequency (kHz)')
text(13,0.0,'12.5 kHz') % (point a little past the folding frequnacy, little above line, what the folding frequancy is)
text(-2,0.0,'0kHz')
text(-2,1.0,'25kHz')
text(13,1.0,'37.5 kHz')
text(-2,2.0,'50kHz')
text(13,2.0,'62.5 kHz')
text(-2,3.0,'75kHz')

axis([-3 16 -1 3.5]); % expands the graph
hold on
plot(.750,0,'m*','MarkerSize',18); % shows what the sampling frequancy is
set(text(.750,0.12,'750 Hz'),'Rotation',45)
title('Nyquist Diagram')
hold off

%% Plotting graph 11 (4.2.3)
clear
%Data input
rawdata1 = importdata('10_7_2024forschino22dat.txt','\t',1);
rawdata2 = importdata('10_7_2024forschino22fft.txt','\t',1);
info1 = rawdata1.data;
t3 = info1(:,1); % time
V3 = info1(:,2); % voltage
info2 = rawdata2.data;
F3 = info2(:,1); % Frequency (Hz)
A3 = info2(:,2); % Amplitude

t3 = t3(1:100);
V3 = V3(1:100);
F3 = F3(1:500);
A3 = A3(1:500);
% time vs voltage graph
figure(11)
subplot(2,2,1)
plot(t3,V3)
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Simulation data - time')

% Freq graph
figure(11)
subplot(2,2,3)
plot(F3,A3)
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('Simulation data - freq')

% Nyquist graph

F11_sampling = 1; % sample rate is 1000 hz, 1Khz
F11_folding = F11_sampling ./ 2;

x = [0 F11_folding 0 F11_folding 0 F11_folding 0 F11_folding]; %100Hz is the Folding frequency in this case
y = [0 0 1 1 2 2 3 3];
figure(11);
subplot(2,2,[2,4])
plot(x,y);
xlabel('Frequency (kHz)')
text(0.5,0.0,'0.5 kHz') % (point a little past the folding frequnacy, little above line, what the folding frequancy is)
text(-0.25,0.0,'0kHz')
text(-0.25,1.0,'1kHz')
text(0.75,1.0,'1.5 kHz')
text(-0.25,2.0,'2kHz')
text(0.75,2.0,'2.5 kHz')
text(-0.25,3.0,'3kHz')
hold on
plot(0.250, 0,'o','MarkerSize',10)
set(text(.250,-0.12,'250 Hz'),'Rotation',-45)
hold on
plot([0.25 0.25],[0 0.5])
axis([-0.35 1 -1 3.5]); % expands the graph
hold on
plot(.250,0.5,'m*','MarkerSize',18); % shows what the sampling frequancy is
set(text(.250,0.62,'750 Hz'),'Rotation',45)
title('Nyquist Diagram')
hold off

%% Plotting graph 12 (4.2.5)
clear
%Data input
rawdata1 = importdata('10_7_2024forschino23dat.txt','\t',1);
rawdata2 = importdata('10_7_2024forschino23fft.txt','\t',1);
info1 = rawdata1.data;
t3 = info1(:,1); % time
V3 = info1(:,2); % voltage
info2 = rawdata2.data;
F3 = info2(:,1); % Frequency (Hz)
A3 = info2(:,2); % Amplitude

t3 = t3(1:100);
V3 = V3(1:100);
F3 = F3(1:500);
A3 = A3(1:500);
% time vs voltage graph
figure(12)
subplot(1,2,1)
plot(t3,V3)
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Simulation data - time')

% Freq graph
figure(12)
subplot(1,2,2)
plot(F3,A3)
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('Simulation data - freq')



%% Plotting graph 13 (4.2.6)
clear
%Data input
rawdata1 = importdata('10_7_2024forschino24dat.txt','\t',1);
rawdata2 = importdata('10_7_2024forschino24fft.txt','\t',1);
info1 = rawdata1.data;
t3 = info1(:,1); % time
V3 = info1(:,2); % voltage
info2 = rawdata2.data;
F3 = info2(:,1); % Frequency (Hz)
A3 = info2(:,2); % Amplitude

t3 = t3(1:100);
V3 = V3(1:100);
F3 = F3(1:500);
A3 = A3(1:500);
% time vs voltage graph
figure(13)
subplot(1,2,1)
plot(t3,V3)
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Simulation data - time')

% Freq graph
figure(13)
subplot(1,2,2)
plot(F3,A3)
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('Simulation data - freq')

%% Plotting graph 14 (4.2.7)
clear
%Data input
rawdata1 = importdata('10_7_2024forschino25dat.txt','\t',1);
rawdata2 = importdata('10_7_2024forschino25fft.txt','\t',1);
info1 = rawdata1.data;
t3 = info1(:,1); % time
V3 = info1(:,2); % voltage
info2 = rawdata2.data;
F3 = info2(:,1); % Frequency (Hz)
A3 = info2(:,2); % Amplitude

t3 = t3(1:100);
V3 = V3(1:100);
F3 = F3(1:500);
A3 = A3(1:500);
% time vs voltage graph
figure(14)
subplot(1,2,1)
plot(t3,V3)
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Simulation data - time')

% Freq graph
figure(14)
subplot(1,2,2)
plot(F3,A3)
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('Simulation data - freq')

%% Plotting graph 15 (4.2.8)
clear
%Data input
rawdata1 = importdata('10_7_2024forschino26dat.txt','\t',1);
rawdata2 = importdata('10_7_2024forschino26fft.txt','\t',1);
info1 = rawdata1.data;
t3 = info1(:,1); % time
V3 = info1(:,2); % voltage
info2 = rawdata2.data;
F3 = info2(:,1); % Frequency (Hz)
A3 = info2(:,2); % Amplitude

t3 = t3(1:100);
V3 = V3(1:100);
F3 = F3(1:500);
A3 = A3(1:500);
% time vs voltage graph
figure(15)
subplot(2,2,1)
plot(t3,V3)
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Simulation data - time')

% Freq graph
figure(15)
subplot(2,2,3)
plot(F3,A3)
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('Simulation data - freq')

% Nyquist graph

F15_sampling = 1; % sample rate is 1000 hz, 1Khz
F15_folding = F15_sampling ./ 2;

x = [0 F15_folding 0 F15_folding 0 F15_folding 0 F15_folding]; %100Hz is the Folding frequency in this case
y = [0 0 1 1 2 2 3 3];
figure(15);
subplot(2,2,[2,4])
plot(x,y);
xlabel('Frequency (kHz)')
text(0.5,0.0,'0.5 kHz') % (point a little past the folding frequnacy, little above line, what the folding frequancy is)
text(-0.25,0.0,'0kHz')
text(-0.25,1.0,'1kHz')
text(0.75,1.0,'1.5 kHz')
text(-0.25,2.0,'2kHz')
text(0.75,2.0,'2.5 kHz')
text(-0.25,3.0,'3kHz')
hold on
plot(0.250, 0,'o','MarkerSize',10)
set(text(.250,-0.12,'250 Hz'),'Rotation',-45)
hold on
plot([0.25 0.25],[0 0.5])
axis([-0.35 1 -1 3.5]); % expands the graph
hold on
plot(.250,0.5,'m*','MarkerSize',18); % shows what the sampling frequancy is
set(text(.250,0.62,'750 Hz'),'Rotation',45)
title('Nyquist Diagram')
hold off

%% Plotting graph 16 (4.2.9)
clear
%Data input
rawdata1 = importdata('10_7_2024forschino27dat.txt','\t',1);
rawdata2 = importdata('10_7_2024forschino27fft.txt','\t',1);
info1 = rawdata1.data;
t3 = info1(:,1); % time
V3 = info1(:,2); % voltage
info2 = rawdata2.data;
F3 = info2(:,1); % Frequency (Hz)
A3 = info2(:,2); % Amplitude

t3 = t3(1:100);
V3 = V3(1:100);
F3 = F3(1:500);
A3 = A3(1:500);
% time vs voltage graph
figure(16)
subplot(2,2,1)
plot(t3,V3)
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Simulation data - time')

% Freq graph
figure(16)
subplot(2,2,3)
plot(F3,A3)
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('Simulation data - freq')

% Nyquist graph

F16_sampling = 1; % sample rate is 1000 hz, 1Khz
F16_folding = F16_sampling ./ 2;

x = [0 F16_folding 0 F16_folding 0 F16_folding 0 F16_folding]; %100Hz is the Folding frequency in this case
y = [0 0 1 1 2 2 3 3];
figure(16);
subplot(2,2,[2,4])
plot(x,y);
xlabel('Frequency (kHz)')
text(0.5,0.0,'0.5 kHz') % (point a little past the folding frequnacy, little above line, what the folding frequancy is)
text(-0.25,0.0,'0kHz')
text(-0.25,1.0,'1kHz')
text(0.75,1.0,'1.5 kHz')
text(-0.25,2.0,'2kHz')
text(0.75,2.0,'2.5 kHz')
text(-0.25,3.0,'3kHz')
hold on
plot(0.250, 0,'o','MarkerSize',10)
set(text(.250,-0.12,'250 Hz'),'Rotation',-45)
hold on
plot([0.25 0.25],[0 0.5])
axis([-0.35 1 -1 3.5]); % expands the graph
hold on
plot(.250,0.5,'m*','MarkerSize',18); % shows what the sampling frequancy is
set(text(.250,0.62,'750 Hz'),'Rotation',45)
title('Nyquist Diagram')
hold off


%% Plotting graph 17 (4.3.5)
clear
%Data input
rawdata1 = importdata('10_7_2024forschino28dat.txt','\t',1);
rawdata2 = importdata('10_7_2024forschino28fft.txt','\t',1);
info1 = rawdata1.data;
t3 = info1(:,1); % time
V3 = info1(:,2); % voltage
info2 = rawdata2.data;
F3 = info2(:,1); % Frequency (Hz)
A3 = info2(:,2); % Amplitude

t3 = t3(1:100);
V3 = V3(1:100);
F3 = F3(1:500);
A3 = A3(1:500);
% time vs voltage graph
figure(17)
subplot(1,2,1)
plot(t3,V3)
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Simulation data - time')

% Freq graph
figure(17)
subplot(1,2,2)
plot(F3,A3)
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('Simulation data - freq')


%% Plotting graph 18 (4.3.6)
clear
%Data input
rawdata1 = importdata('10_7_2024forschino29dat.txt','\t',1);
rawdata2 = importdata('10_7_2024forschino29fft.txt','\t',1);
info1 = rawdata1.data;
t3 = info1(:,1); % time
V3 = info1(:,2); % voltage
info2 = rawdata2.data;
F3 = info2(:,1); % Frequency (Hz)
A3 = info2(:,2); % Amplitude

t3 = t3(1:100);
V3 = V3(1:100);
F3 = F3(1:500);
A3 = A3(1:500);
% time vs voltage graph
figure(18)
subplot(1,2,1)
plot(t3,V3)
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Simulation data - time')

% Freq graph
figure(18)
subplot(1,2,2)
plot(F3,A3)
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('Simulation data - freq')




%% Partial sums sawtooth wave
clear 

syms n t A T w0
bn_f1 = t.*sin(n.*w0.*t); 
ibn = int(bn_f1,t,0,T/2); 

bn = ((8.*A)./(T.^2)).*ibn; 
A = 2.5; %Amplitude of sawtooth wave is 2.5
f = 750; %Frequency of sawtooth wave is 750
T = 1/f; %Calculates period of sawtooth wave
w0 = 2*pi/T; %Calculates
N = 101; %Number Foerier coefficients
n = [1:1:N]; %Find area for each Foerier coefficient
t = [0:0.000001:0.004]; %Time to plot
bn = eval(subs(bn)); %Solve

%Plot the Fourier series for a defined (increasing)
% number of coefficients
figure(19)
tiledlayout (3,2)
for N = [1,3,5,7,9,33]
 X = 0;
 for i = 1:N
 X = X + bn(i).*sin(n(i)*w0.*t);
 i = i+1;
 end
 nexttile
plot(t,X)
 xlim([0,0.004])
 ylim([-4,4])
 title('Sawtooth Partial Sums')
 ylabel('Voltage (V)')
 xlabel('Time(s)')
 leg = legend(num2str(i-1));
 title(leg,'n')
 hold on %Plot them all on the same plot
end



%% Partial sums square wave
clear 

syms n t A T w0
bn_f1 = sin(n.*w0.*t); 
ibn = int(bn_f1,t,0,T/2); 

bn = ((-4.*A)./(T)).*ibn; 
A = 2.5; %Amplitude of square wave is 2.5
f = 750; %Frequency of square wave is 750
T = 1/f; %Calculates period of square wave
w0 = 2*pi/T; %Calculates
N = 101; %Number Foerier coefficients
n = [1:1:N]; %Find area for each Foerier coefficient
t = [0:0.000001:0.004]; %Time to plot
bn = eval(subs(bn)); %Solve

%Plot the Fourier series for a defined (increasing)
% number of coefficients
figure(20)
tiledlayout (3,2)
for N = [1,3,5,7,9,33]
 X = 0;
 for i = 1:N
 X = X + bn(i).*sin(n(i)*w0.*t);
 i = i+1;
 end
 nexttile
plot(t,X)
 xlim([0,0.004])
 ylim([-4,4])
 title('Square Partial Sums')
 ylabel('Voltage (V)')
 xlabel('Time(s)')
 leg = legend(num2str(i-1));
 title(leg,'n')
 hold on %Plot them all on the same plot
end




